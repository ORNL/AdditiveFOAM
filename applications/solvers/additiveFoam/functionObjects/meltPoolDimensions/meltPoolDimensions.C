/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "meltPoolDimensions.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvc.H"
#include "OSspecific.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meltPoolDimensions, 0);
    
    addToRunTimeSelectionTable
    (
        functionObject,
        meltPoolDimensions,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::meltPoolDimensions::meltPoolDimensions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    T_(mesh_.lookupObject<VolField<scalar>>("T"))
{
    read(dict);
    
    if (Pstream::master())
    {
        const fileName probeDir
        (
            mesh_.time().rootPath()/mesh_.time().globalCaseName()
           /"postProcessing"/"meltPoolDimensions"
        );

        mkDir(probeDir);

        logFilePtrs_.setSize(isoValues_.size());

        forAll(isoValues_, i)
        {
            const fileName logFileName
            (
                probeDir/Foam::name(isoValues_[i]) + ".csv"
            );
            
            Info<< "melt pool dimensions log file name: "
                << logFileName << endl;

            logFilePtrs_.set
            (
                i,
                new OFstream(logFileName)
            );
            
            logFilePtrs_[i] << "Time,Length,Width,Depth" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::meltPoolDimensions::~meltPoolDimensions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::meltPoolDimensions::read(const dictionary& dict)
{
    isoValues_ = dict.lookup<scalarList>("isoValues");
    scanPathAngle_ = dict.lookupOrDefault<scalar>("scanPathAngle", 0.0);
    
    return true;
}

Foam::wordList Foam::functionObjects::meltPoolDimensions::fields() const
{
    return wordList::null();
}

bool Foam::functionObjects::meltPoolDimensions::execute()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    
    const volVectorField& cc = mesh_.C();

    const scalar radians =
        scanPathAngle_ * ( Foam::constant::mathematical::pi / 180.0 );

    const scalar s = sin(radians);
    const scalar c = cos(radians);
  
    List<treeBoundBox> boundBoxes
    (
        isoValues_.size(),
        treeBoundBox(point::max, point::min)
    );

    // check internal faces
    for(label facei=0; facei < mesh_.nInternalFaces(); facei++)
    {        
        const label own = owner[facei];
        const label nei = neighbour[facei];

        scalar minFace = min(T_[own], T_[nei]);
        scalar maxFace = max(T_[own], T_[nei]);

        forAll(isoValues_, i)
        {
            const scalar iso_ = isoValues_[i];

            // update the bounding box
            if ((minFace < iso_) && (maxFace >= iso_))
            {
                vector d = cc[nei] - cc[own];
                
                vector p =
                    cc[own] + d*(iso_ - T_[own])/(T_[nei] - T_[own]);

                vector p_rotated
                (
                    p.x()*c + p.y()*s,
                    p.x()*s + p.y()*c,
                    p.z()
                );
                
                boundBoxes[i].min() =
                    min(p_rotated, boundBoxes[i].min());

                boundBoxes[i].max() =
                    max(p_rotated, boundBoxes[i].max());
            }
        }
    }

    // check boundary faces
    const volScalarField::Boundary& TBf = T_.boundaryField();
    
    forAll(TBf, patchi)
    {   
        const fvPatchScalarField& TPf = TBf[patchi];

        const labelUList& faceCells = TPf.patch().faceCells();

        if (TPf.coupled())
        {
            // processor boundary : interpolate across face
            const vectorField ccn
            (
                cc.boundaryField()[patchi].patchNeighbourField()
            );
            
            const scalarField fn(TPf.patchNeighbourField());

            forAll(faceCells, facei)
            {
                label own = faceCells[facei];

                scalar minFace = min(T_[own], fn[facei]);
                scalar maxFace = max(T_[own], fn[facei]);


                forAll(isoValues_, i)
                {
                    const scalar iso_ = isoValues_[i];
                    
                    // update the bounding box
                    if ((minFace < iso_) && (maxFace >= iso_))
                    {
                        vector d = ccn[facei] -  cc[own];
                        
                        vector p = 
                            cc[own] + d*(iso_ - T_[own])/(fn[facei] - T_[own]);

                        vector p_rotated
                        (
                            p.x()*c + p.y()*s,
                            p.x()*s + p.y()*c,
                            p.z()
                        );
                        
                        boundBoxes[i].min() =
                            min(p_rotated, boundBoxes[i].min());

                        boundBoxes[i].max() =
                            max(p_rotated, boundBoxes[i].max());
                    }
                }
            }
        }
        else
        {
            // physical boundary : take face point if above iso value
            const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];
          
            const scalarField& pif(TPf.patchInternalField());
            
            forAll(faceCells, facei)
            {
                scalar maxFace = max(pif[facei], TPf[facei]);
                
                forAll(isoValues_, i)
                {
                    const scalar iso_ = isoValues_[i];
                    
                    if (maxFace >= iso_)
                    {
                        const vector& p = Cf[facei];
                        
                        vector p_rotated
                        (
                            p.x()*c + p.y()*s,
                            p.x()*s + p.y()*c,
                            p.z()
                        );
                        
                        boundBoxes[i].min() =
                            min(p_rotated, boundBoxes[i].min());

                        boundBoxes[i].max() =
                            max(p_rotated, boundBoxes[i].max());
                    }
                }
            }
        }
    }
    
    forAll(isoValues_, i)
    {
        reduce(boundBoxes[i].min(), minOp<point>());
        reduce(boundBoxes[i].max(), maxOp<point>());
    }
    
    // update the melt pool dimensions log files
    if (Pstream::master())
    {
        forAll(isoValues_, i)
        {
            vector dimensions = max(boundBoxes[i].span(), vector::zero);
            
            OFstream& logFile = logFilePtrs_[i];
                        
            logFile<< mesh_.time().value() << ","
                   << dimensions.x() << ","
                   << dimensions.y() << ","
                   << dimensions.z() << endl;
        }
    }

    return true;
}

bool Foam::functionObjects::meltPoolDimensions::end()
{
    return true;
}


bool Foam::functionObjects::meltPoolDimensions::write()
{
    return true;
}


// ************************************************************************* //
