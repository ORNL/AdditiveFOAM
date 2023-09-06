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

#include "movingHeatSourceModel.H"
#include "DynamicList.H"
#include "OFstream.H"
#include "isoPoints.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::movingHeatSourceModel
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (    
        IOobject
        (
            "heatSourceDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sourceNames_(dict_.lookup("sources")),
    qDot_
    (
        IOobject
        (
            "qDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume, 0.0)
    )    
{
    sources_.resize(sourceNames_.size());
        
    //- Create new instance of movingBeam for each beam
    forAll(sources_, i)
    {
        Info << "Adding heatSourceModel for " << sourceNames_[i] << endl;
        sources_.set
        (
            i,
            heatSourceModel::New
            (
                sourceNames_[i],
                dict_,
                mesh_
            ).ptr()
        );
    }
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::~movingHeatSourceModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingHeatSourceModel::adjustDeltaT(scalar& deltaT)
{
    forAll(sources_, i)
    {
        sources_[i].beam().adjustDeltaT(deltaT);
    }
}

bool Foam::movingHeatSourceModel::transientExists()
{
    bool exists_ = false;
    
    forAll(sources_, i)
    {
        if (sources_[i].transient())
        {
            exists_ = true;
            break;
        }
    }

    return exists_;
}

void Foam::movingHeatSourceModel::update()
{
    //- Update the moving heat source shape for each transient case
    if ( transientExists() )
    {
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        forAll(sources_, i)
        {
            if (sources_[i].transient())
            {
                List<point> localPoints = isoPoints(T, sources_[i].isoValue());

                // get the user-defined (static) heat source dimensions
                vector newDimensions = sources_[i].staticDimensions();

                scalar maxDepth = newDimensions.z();

                const scalar searchRadius(cmptMax(newDimensions));

                for (const point& pt : localPoints)
                {
                    vector d = cmptMag(pt - sources_[i].beam().position());

                    scalar r = Foam::sqrt(d.x()*d.x() + d.y()*d.y());
                    
                    if (r <= searchRadius)
                    {
                        maxDepth = max(maxDepth, d.z());
                    }
                }

                reduce(maxDepth, maxOp<scalar>());

                newDimensions[2] = maxDepth;

                // helper function to update the heat source dimensions
                sources_[i].setDimensions(newDimensions);
            }

            Info << sourceNames_[i] << tab << sources_[i].dimensions() << endl;
        }
    }

    //- Subcycle each moving heat source in time and combine into a single field
    qDot_ = dimensionedScalar("Zero", qDot_.dimensions(), 0.0);
    
    forAll(sources_, i)
    {
        if (sources_[i].beam().activePath())
        {
            scalar pathTime = mesh_.time().value();

            const scalar nextTime = pathTime + mesh_.time().deltaTValue();

            const scalar beam_dt = sources_[i].beam().deltaT();

            volScalarField qDoti
            (
                IOobject
                (
                    "qDoti",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("Zero", qDot_.dimensions(), 0.0)
            );
            
            scalar sumWeights = 0.0;

            while ((nextTime - pathTime) > small)
            {
                scalar dt = min(beam_dt, max(0, nextTime - pathTime));

                pathTime += dt;

                sources_[i].beam().move(pathTime);
                                
                qDoti += dt*sources_[i].qDot();

                sumWeights += dt;
            }
            
            qDoti /= sumWeights;
            
            qDot_ += qDoti;
        }
    }
}

// ************************************************************************* //
