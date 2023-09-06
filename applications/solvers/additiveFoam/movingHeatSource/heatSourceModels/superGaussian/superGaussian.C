/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "superGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, superGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::superGaussian::superGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatSourceModels::superGaussian::qDot()
{
    tmp<volScalarField> tqDot
    (
        new volScalarField
        (
            IOobject
            (
                "qDot_",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimPower/dimVolume, 0.0)
        )
    );
    volScalarField& qDot_ = tqDot.ref();

    const scalar power_ = movingBeam_->power();
    
    if (power_ > small)
    {
        const scalar AR = 
            dimensions_.z() / min(dimensions_.x(), dimensions_.y());
        
        const scalar eta_ = absorptionModel_->eta(AR);
        
        //- Calculate volumetric intesity
        const vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

        const dimensionedScalar I0
        (
            "I0",
            dimPower / dimVolume,
            eta_*power_*k_
          / (s.x()*s.y()*s.z()*2.0*pi*Foam::tgamma(3.0/k_))
        );

        //- Calculate heat source weights a cell-faces
        surfaceScalarField factor
        (
            IOobject
            (
                "factor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimless, 0.0)          
        );

        for (label facei=0; facei < mesh_.nInternalFaces(); facei++)
        {
            point d = cmptMag(mesh_.Cf()[facei] - movingBeam_->position());
            
            if (mag(d) < D4sigma()/2.0)
            {
                scalar f = magSqr(cmptDivide(d, s));

                factor[facei] = Foam::exp(-Foam::pow(f, k_/2.0));
            }
        }

        surfaceScalarField::Boundary& bFactor = factor.boundaryFieldRef();
        forAll(bFactor, patchi)
        {
            const pointField& faceCentres = mesh_.Cf().boundaryField()[patchi];
            forAll(faceCentres, facei)
            {
                point d = cmptMag(faceCentres[facei] - movingBeam_->position());
                
                if (mag(d) < D4sigma()/2.0)
                {
                    scalar f = magSqr(cmptDivide(d, s));

                    bFactor[patchi][facei] = Foam::exp(-Foam::pow(f, k_/2.0));
                }
                else
                {
                    bFactor[patchi][facei] = 0.0;
                }
            }
        }

        //- Assemble qDot from normalized power and face weights
        qDot_ = I0 * fvc::average(factor);
        qDot_.correctBoundaryConditions();
        
        //- Perform power correction for coarse meshes
        if (normalize_)
        {
            correctPower(qDot_, eta_);
        }
    }

    return tqDot;
}


scalar Foam::heatSourceModels::superGaussian::D4sigma() const
{
    return 2.0 * max(dimensions_.x(), dimensions_.y());
}


bool Foam::heatSourceModels::superGaussian::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
