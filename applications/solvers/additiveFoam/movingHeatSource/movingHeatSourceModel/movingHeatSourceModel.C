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
    refineDict_(dict_.optionalSubDict("refinementControls")),
    refine_(refineDict_.lookupOrDefault<bool>("refine", false)),
    resolveTail_(false),
    refinementInterval_(0),
    refinementLevel_(0),
    refinementIndex_(0),
    fwdSteps_(0),
    revSteps_(0),
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
    ),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    )
{        
    //- Create new instance of movingBeam for each beam
    sources_.resize(sourceNames_.size());
    
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
    
    if (refine_)
    {
        //- Look up refinementField controls
        resolveTail_ = refineDict_.lookup<bool>("resolveTail");
        refinementInterval_ = refineDict_.lookup<label>("refinementInterval");
        refinementLevel_ = refineDict_.lookup<label>("refinementLevel");
        fwdSteps_ = refineDict_.lookup<label>("fwdSteps");
        revSteps_ = refineDict_.lookup<label>("revSteps");
    }
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::~movingHeatSourceModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::movingHeatSourceModel::refineTime()
{
    //- Return true after specified refinementInterval
    if(mesh_.time().timeIndex() == refinementIndex_)
    {
        Info << "IT'S REFINE TIME" << endl;
    
        refinementIndex_ += refinementInterval_;

        return true;
    }

    //- Also return true for next n time steps, where n = refinementLevel
    else if (mesh_.time().timeIndex()
                 <= refinementIndex_ - refinementInterval_ + refinementLevel_)
    {
        Info << "Performing additional refinement operations." << endl;

        return true;
    }

    else
    {
        return false;
    }
}


void Foam::movingHeatSourceModel::adjustDeltaT(scalar& deltaT)
{
    forAll(sources_, i)
    {
        sources_[i].beam().adjustDeltaT(deltaT);
    }
}


void Foam::movingHeatSourceModel::update()
{
    //- Reset qDot at every time step
    qDot_ = dimensionedScalar("Zero", qDot_.dimensions(), 0.0);

    //- Reset refinement field after desired interval
    bool resetRF = false;
    if ((mesh_.time().timeIndex() == 0)
         ||
        (mesh_.time().timeIndex() == refinementIndex_))
    {
        resetRF = true;
        
        if (resolveTail_)
        {
            const volScalarField& alpha1
                = mesh_.lookupObject<volScalarField>("alpha.solid");
                
            //- Reset to capture melt pool tail
            refinementField_ = pos(1.0 - alpha1);
        }
        
        else
        {
            //- Reset to zero
            refinementField_
                = dimensionedScalar(refinementField_.dimensions(), 0.0);
        }
    }
    
    //- Subcycle each moving heat source in time and combine into a single field
    forAll(sources_, i)
    {
        if (sources_[i].beam().activePath())
        {
            sources_[i].updateDimensions();

            // integrate volumetric heat source over desired time step
            scalar pathTime = mesh_.time().value();

            const scalar nextTime = pathTime + mesh_.time().deltaTValue();

            const scalar beam_dt = sources_[i].beam().deltaT();

            //- Update individual beam heat source contribution
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
            
            //- Update refinementField
            if (refine_ && resetRF)
            {
                //- Create dummy variables for beam position/power during integration
                const scalar currTime = mesh_.time().value();
                const scalar currDt = mesh_.time().deltaTValue();
                vector currPos = Zero;
                scalar currPow = 0.0;
                
                const pointField& points = mesh_.points();
                const vector extend = 1e-10 * vector::one;
                
                //- Buffer number of steps ahead of beam to prevent overrun during
                //  re-refinement of mesh at end of interval
                const int buffer = 5;
                
                //- Integrate forward in time ahead of beam
                for (int j = 0; j < fwdSteps_ + buffer; ++j)
                {
                    //- Update time, position, and power for current step
                    scalar nowTime = currTime + j * currDt;
                    sources_[i].beam().move(currPos, currPow, nowTime);
                    
                    //- Create beam bounding box for current time
                    vector beamDims = sources_[i].dimensions();
                    treeBoundBox beamBb(currPos - beamDims, currPos + beamDims);
                    
                    //- Check if beam bounding box overlaps cells
                    forAll(mesh_.C(), celli)
                    {
                        //- Don't redo checks if cell already marked for refinement
                        if (refinementField_[celli] > 0.0)
                        {
                            continue;
                        }
                        
                        treeBoundBox cellBb(point::max, point::min);
                        
                        const labelList& vertices = mesh_.cellPoints()[celli];
                        
                        //- Create cell bounding box from vertices with furthest
                        //  extent
                        forAll(vertices, k)
                        {
                            cellBb.min()
                                = min(cellBb.min(), points[vertices[k]] - extend);
                            cellBb.max()
                                = max(cellBb.max(), points[vertices[k]] + extend);
                        }

                        //- Mark cell for refinement if its bounding box overlaps
                        //  the beam bounding box
                        if (cellBb.overlaps(beamBb))
                        {
                            refinementField_[celli] = 1.0;
                        }
                    }
                }
                
                //- Integrate backward in time behind beam
                for (int j = 0; j < revSteps_ + buffer; ++j)
                {
                    //- Update time, position, and power for current step
                    scalar nowTime = currTime - j * currDt;
                    sources_[i].beam().move(currPos, currPow, nowTime);
                    
                    //- Create beam bounding box for current time
                    vector beamDims = sources_[i].dimensions();
                    treeBoundBox beamBb(currPos - beamDims, currPos + beamDims);
                    
                    //- Check if beam bounding box overlaps cells
                    forAll(mesh_.C(), celli)
                    {
                        //- Don't redo checks if cell already marked for refinement
                        if (refinementField_[celli] > 0.0)
                        {
                            continue;
                        }
                        
                        treeBoundBox cellBb(point::max, point::min);
                        
                        const labelList& vertices = mesh_.cellPoints()[celli];
                        
                        //- Create cell bounding box from vertices with furthest
                        //  extent
                        forAll(vertices, k)
                        {
                            cellBb.min()
                                = min(cellBb.min(), points[vertices[k]] - extend);
                            cellBb.max()
                                = max(cellBb.max(), points[vertices[k]] + extend);
                        }

                        //- Mark cell for refinement if its bounding box overlaps
                        //  the beam bounding box
                        if (cellBb.overlaps(beamBb))
                        {
                            refinementField_[celli] = 1.0;
                        }
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
