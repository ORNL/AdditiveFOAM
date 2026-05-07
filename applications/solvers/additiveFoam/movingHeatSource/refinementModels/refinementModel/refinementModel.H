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

Class
    Foam::refinementModel::refinementModel

Description
    Abstract base class for refinementModel models

SourceFiles
    refinementModel.C

\*---------------------------------------------------------------------------*/

#ifndef refinementModel_H
#define refinementModel_H

#include "labelVector.H"
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
        bool cellAABBOverlapsOBB
        (
            const treeBoundBox& cellAABB,
            const point& centre,
            const vector& e0,
            const vector& e1,
            const vector& e2,
            const vector& L
        ) const;

        //- Mark/count cells swept by powered path vectors over a time interval
        scalar markScanPathInterval
        (
            const scalar intervalStartTime,
            const scalar intervalEndTime,
            const List<treeBoundBox>& cellAABBs,
            const bool commit
        );

        //- Return next powered path-vector start or end time after time
        scalar nextPoweredPathEventTime(const scalar time) const;


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

        //- Switch for refinementField generation
        const bool refine_;

        //- Number of refinement levels
        label nLevels_;

        //- Temperature above which refinement occurs
        scalar refinementTemperature_;

        //- Refinement extent for scan-path OBB slabs
        //  x: half-span along e0, width/lateral direction
        //  y: endpoint padding along e1, scan direction
        //  z: half-span along e2, depth/vertical direction
        vector buffer_;

        //- End time of last beam
        scalar endTime_;

        //- Marker field for refinement operations
        volScalarField refinementField_;


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

        //- Return the refinement marker field
        virtual volScalarField refinementField()
        {
            return refinementField_;
        }

        //- Set refinementField using temperature constraints
        virtual void refineUsingTemperature();

        //- Refine along scan path until specified time
        virtual void refineUsingTime(const scalar& refinementTime);

        //- Refine along scan path until global refinement volume reaches
        //  provided target. Returns time at which target is reached.
        //  The second argument is kept for backwards-compatible call sites;
        //  use minRefinementVolumeFactor in the refinement dictionary instead.
        virtual dimensionedScalar refineUsingVolume
        (
            const dimensionedScalar& refinementVolume,
            const scalar& unusedMinimumScanPathRefinementVolume = 0.0
        );

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
