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

    You should have received a copy of the the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::refinementModel

Description
    Base class for models that populate the refinementField used by the
    OpenFOAM dynamic mesh refiner.

    If no refinementModel dictionary is present in heatSourceDict, the
    no-refinement model is selected and no AdditiveFOAM refinement controls
    are required.

    Active refinement models project each heat-source scan path and mark
    cells for refinement using buffers keyed by source name. The refinement
    level is read from topoChanger.maxRefinement in constant/dynamicMeshDict.

SourceFiles
    refinementModel.C

\*---------------------------------------------------------------------------*/

#ifndef refinementModel_H
#define refinementModel_H

#include "labelVector.H"
#include "hexMatcher.H"
#include "dictionary.H"
#include "volFields.H"
#include "movingBeam.H"
#include "treeBoundBox.H"
#include "heatSourceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class refinementModel Declaration
\*---------------------------------------------------------------------------*/

class refinementModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct the base IO object
        IOobject createIOobject
        (
            const dictionary& dict,
            const fvMesh& mesh
        ) const;

        //- Calculate cell axis-aligned bounding boxes
        void calculateCellAABBs(List<treeBoundBox>& cellAABBs) const;

        //- Calculate local scan-path coordinate frame
        void scanPathFrame
        (
            const point& p0,
            const point& p1,
            vector& e0,
            vector& e1,
            vector& e2
        ) const;

        //- Return true if a cell AABB overlaps a yawed path OBB
        bool cellOverlapsOBB
        (
            const treeBoundBox& cellAABB,
            const point& centre,
            const vector& e0,
            const vector& e1,
            const vector& e2,
            const vector& L
        ) const;

        //- Mark/count cells swept by powered path vectors over a time interval
        scalar markScanPath
        (
            const scalar startTime,
            const scalar endTime,
            const List<treeBoundBox>& cellAABBs,
            const bool commit
        );

        //- Return next powered path-vector start or end time after time
        scalar nextPoweredPathEventTime(const scalar time) const;

        //- Read max refinement level from dynamicMeshDict
        label readMaxRefinementLevel() const;

        //- Read source-specific scan-path buffers
        void readSourceBuffers();


protected:

    // Protected Data

        //- PtrList to heatSourceModels
        const PtrList<heatSourceModel>& sources_;

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- Main heat source dictionary
        const dictionary& heatSourceDict_;

        //- Sub-dictionary for refinementModel entries
        dictionary refinementDict_;

        //- Maximum dynamic-mesh refinement level
        label maxRefinementLevel_;

        //- Temperature above which refinement occurs
        scalar refinementTemperature_;

        //- Refinement extent for each source's scan-path OBB slabs
        //  x: endpoint padding along e0, scan direction
        //  y: half-span along e1, transverse direction
        //  z: half-span along e2, depth/vertical direction
        List<vector> sourceBuffers_;

        //- End time of all scan paths
        scalar scanEndTime_;

        //- Marker field for refinement operations
        volScalarField refinementField_;


    // Protected Member Functions

        //- Return summed volume of source refinement buffers
        scalar sourceBufferVolume() const;

        //- Return minimum future scan-path refine volume
        scalar minimumScanPathRefineVolume
        (
            const scalar nBufferVolumes
        ) const;

        //- Mark refinement field by projecting scan path to target volume
        virtual dimensionedScalar markScanPathVolume
        (
            const dimensionedScalar& targetRefineVolume,
            const scalar minScanPathRefineVolume,
            const label maxSearchIter,
            const scalar timeTolerance
        );


public:

    //- Runtime type information
    TypeName("refinementModel");

    //- Default refinement model dictionary name
    static const word refinementModelDictName;

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        refinementModel,
        dictionary,
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (sources, dict, mesh)
    );


    // Constructors

        //- Construct inactive refinement model
        refinementModel
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Construct from components
        refinementModel
        (
            const word& type,
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        refinementModel(const refinementModel&);


    // Selectors

        //- Return a reference to the selected refinement model
        static autoPtr<refinementModel> New
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~refinementModel()
    {}


    // Member Functions

        //- Return access to the refinement marker field
        virtual volScalarField& refinementField()
        {
            return refinementField_;
        }

        //- Return const access to the refinement marker field
        virtual const volScalarField& refinementField() const
        {
            return refinementField_;
        }

        //- Mark refinement field using temperature constraints
        virtual void markTemperature();

        //- Mark refinement field by projecting scan path to specified time
        virtual void markScanPathTime(const scalar& refineTime);

        //- Backwards-compatible wrapper
        virtual void refineUsingTemperature()
        {
            markTemperature();
        }

        //- Backwards-compatible wrapper
        virtual void refineUsingTime(const scalar& refineTime)
        {
            markScanPathTime(refineTime);
        }

        //- Update refinementField
        virtual bool update() = 0;

        //- Read the heat source dictionary
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const refinementModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
