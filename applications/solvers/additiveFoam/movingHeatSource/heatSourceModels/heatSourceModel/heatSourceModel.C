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

#include "heatSourceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatSourceModel, 0);
    defineRunTimeSelectionTable(heatSourceModel, dictionary);
}

const Foam::word Foam::heatSourceModel::heatSourceDictName
(
    "heatSourceDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::heatSourceModel::createIOobject
(
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dict.name(),
        mesh.time().constant(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModel::heatSourceModel
(
    const word& type,
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sourceName_(sourceName),
    heatSourceDict_(dict),
    sourceDict_(heatSourceDict_.optionalSubDict(sourceName_)),
    heatSourceModelCoeffs_(sourceDict_.optionalSubDict(type + "Coeffs")),
    
    mesh_(mesh),
    scanPatchName_(""),
    scanPatchPoints_(),
    absorptionModel_(nullptr),
    movingBeam_(nullptr)
{
    absorptionModel_ = absorptionModel::New(sourceName_, heatSourceDict_, mesh_);
    movingBeam_ = movingBeam::New(sourceName_, heatSourceDict_, mesh_.time());

    dimensions_ = heatSourceModelCoeffs_.lookup<vector>("dimensions");
    staticDimensions_ = dimensions_;
    
    normalize_ = sourceDict_.lookupOrDefault<Switch>("normalize", true);
    transient_ = heatSourceModelCoeffs_.lookupOrDefault<Switch>("transient", false);

    if (transient_)
    {
        isoValue_ = heatSourceModelCoeffs_.lookup<scalar>("isoValue");
    }
    
    //- Set up scanPatch boundary point list if normalize is enabled
    if (normalize_)
    {
        //- Get name of patch to scan over
        scanPatchName_ = sourceDict_.lookup<word>("scanPatchName");
        
        //- Get instance of polyPatch for scanPatchName
        label scanPatchID = mesh_.boundaryMesh().findPatchID(scanPatchName_);
        const polyPatch& scanPatch = mesh_.boundaryMesh()[scanPatchID];
        
        //- Create list of LOCALLY ordered boundary points
        const labelList& localScanPatchPoints = scanPatch.boundaryPoints();
        labelList procScanPatchPoints;
        
        //- Append each LOCAL boundary point to a list of PROCESSOR boundary points
        forAll (localScanPatchPoints, pointi)
        {
            const labelList& gMeshPoints = scanPatch.meshPoints();
            procScanPatchPoints.append(gMeshPoints[localScanPatchPoints[pointi]]);
        }
        
        //- Create and gather list of lists for storing all processor scan patch points
        labelListList listProcScanPatchPoints(Pstream::nProcs());
        listProcScanPatchPoints[Pstream::myProcNo()] = procScanPatchPoints;
        Pstream::gatherList(listProcScanPatchPoints);
        
        labelList tempList;
        
        //- Loop over each list in the list of lists
        forAll (listProcScanPatchPoints, listi)
        {            
            Info << "Adding points from processor " << listi << endl;
            
            labelList& procListi = listProcScanPatchPoints[listi];
            
            forAll(procListi, pointi)
            {
                tempList.append(procListi[pointi]);
                
                Info << "    Appending procListi entry " << pointi << "with global point value " << procListi[pointi] << "to tempList" << endl;
            }
        }
        
        //- Sort temporary list
        sort(tempList);
        
        Info << "Temp list: " << endl;
        Info << tempList << endl;
        
        //- Remove duplicate entries from tempList
        forAll (tempList, labeli)
        {
            Info << "Checking tempList entry " << labeli << endl;
            //- Check if current entry is equal to previous entry
            if ((labeli > 0) && (tempList[labeli] == tempList[labeli-1]))
            {
                //- Don't do anything
                Info << "Current entry equal to previous entry." << endl;
            }
            //- Check if current entry is equal to next entry
            else if ((labeli < tempList.size()-1) && (tempList[labeli] == tempList[labeli+1]))
            {
                //- Still don't do anything
                Info << "Current entry equal to next entry." << endl;
            }
            //- Append label to real list if not duplicate
            else
            {
                Info << "Appending tempList point " << labeli << " with global point id " << tempList[labeli] << " to scanPatchPoints_" << endl;
                scanPatchPoints_.append(tempList[labeli]);
            }
        }
        
        forAll (scanPatchPoints_, pointi)
        {
            Info << "scanPatchPoint_ index " << pointi << " has global index " << scanPatchPoints_[pointi] << " with coordinates " << mesh_.points()[scanPatchPoints_[pointi]] << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModel::correctPower
(
    volScalarField& qDot,
    const scalar& eta
)
{    
    //- Calculate total resolved power
    scalar resPower = (fvc::domainIntegrate(qDot)).value() / eta;
    
    //- Only normalize if resolved power is > 0
    if (resPower > SMALL)
    {
        //- Initialize boundary test flag
        bool isNearBoundary = false;
        
        //- Get beam position from movingBeam
        vector position = movingBeam_->position();
        
        scalar minWallDist = GREAT;
        
        //- Loop over points to determine if beam center is within 2sigma of
        //  any patch boundary points
        const pointField& points = mesh_.points();
        forAll (scanPatchPoints_, pointi)
        {
            Info << "Checking scanPatchPoints entry " << pointi << " with coordinates " << points[scanPatchPoints_[pointi]] << endl;
            scalar bDist = mag(position - points[scanPatchPoints_[pointi]]);
            
            if (bDist < minWallDist)
            {
                minWallDist = bDist;
            }
            
            //- Set flag and break loop if position is near boundary
            if (bDist < D4sigma()/2.0)
            {
                Info << "Beam center within 2sigma of boundary. Skipping normalization." << endl;
                isNearBoundary = true;
                break;
            }
        }
        
        Info << "Min wall dist = " << minWallDist << endl;
        
        //- Normalize if beam isn't near boundary
        if (!isNearBoundary)
        {
            Info << "Beam center not within 2sigma of boundary. Applying normalization." << endl;
            
            //- Get expected power from movingBeam class
            scalar expPower = movingBeam_->power();
            
            //- Normalize qDot
            qDot *= expPower / resPower;
            qDot.correctBoundaryConditions();
        }
    }
}

bool Foam::heatSourceModel::read()
{
    if (regIOobject::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
