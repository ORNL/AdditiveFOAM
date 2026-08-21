/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
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

#include "tabulated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(tabulated, 0);
    addToRunTimeSelectionTable(heatSourceModel, tabulated, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::tabulated::tabulated
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    nSlope_(Zero),
    nIntercept_(Zero),
    k_(1),
    profile_()
{
    nSlope_ =
        heatSourceModelCoeffs_.lookup<scalar>("nSlope");

    nIntercept_ =
        heatSourceModelCoeffs_.lookup<scalar>("nIntercept");

    const scalar minimumDepth =
        heatSourceModelCoeffs_.lookup<scalar>("minimumDepth");

    const fileName fName(heatSourceModelCoeffs_.lookup("file"));

    const fileName tableFile
    (
        mesh_.time().rootPath()
       /mesh_.time().globalCaseName()
       /mesh_.time().constant()
       /fName
    );

    profile_.read(tableFile);

    dimensions_ =
        vector
        (
            0.5*profile_.D4Sigma(),
            0.5*profile_.D4Sigma(),
            minimumDepth
        );

    staticDimensions_ = dimensions_;
    Info<< "Tabulated profile: integral=" << profile_.integral()
        << ", centroid=(" << profile_.centroidX()
        << ' ' << profile_.centroidY() << ") m"
        << ", D4Sigma/major/minor="
        << profile_.D4Sigma() << '/'
        << profile_.D4SigmaMajor() << '/'
        << profile_.D4SigmaMinor() << " m"
        << ", azimuth=" << profile_.azimuth() << " rad" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::tabulated::update()
{
    updateDimensions();

    const scalar a =
        dimensions_.z()
       /min(staticDimensions_.x(), staticDimensions_.y());

    const scalar n =
        min
        (
            max
            (
                nSlope_*std::log2(a) + nIntercept_,
                0.0
            ),
            9.0
        );

    k_ = std::pow(2.0, n);

    const scalar zMax =
        dimensions_.z()
       *std::pow
        (
            invIncGammaRatio_P(1.0/k_, 1.0 - profileTol_)/3.0,
            1.0/k_
        );

    sourceMin_ =
        vector(profile_.x0(), profile_.y0(), -zMax);

    sourceMax_ =
        vector(profile_.x1(), profile_.y1(), 0);

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            profile_.integral()*dimensions_.z()
           *Foam::tgamma(1.0/k_)
           /(k_*std::pow(3.0, 1.0/k_))
        );
}


inline Foam::scalar
Foam::heatSourceModels::tabulated::weight(const vector& r) const
{
    if (r.z() > 0)
    {
        return 0;
    }

    const scalar I = profile_.interpolate(r.x(), r.y());

    const scalar z = -r.z();

    const scalar p =
        std::exp(-3.0*std::pow(z/dimensions_.z(), k_));

    return I*p;
}


bool Foam::heatSourceModels::tabulated::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_.lookup("nSlope") >> nSlope_;
        heatSourceModelCoeffs_.lookup("nIntercept")
            >> nIntercept_;
        const scalar minimumDepth =
            heatSourceModelCoeffs_.lookup<scalar>("minimumDepth");

        const fileName fName(heatSourceModelCoeffs_.lookup("file"));

        const fileName tableFile
        (
            mesh_.time().rootPath()
           /mesh_.time().globalCaseName()
           /mesh_.time().constant()
           /fName
        );

        profile_.read(tableFile);

        staticDimensions_ =
            vector
            (
                0.5*profile_.D4Sigma(),
                0.5*profile_.D4Sigma(),
                minimumDepth
            );

        dimensions_ = staticDimensions_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
