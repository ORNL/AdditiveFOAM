/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

Class
    Foam::functionObjects::ExaCA

SourceFiles
    ExaCA.C

\*---------------------------------------------------------------------------*/

#ifndef ExaCA_H
#define ExaCA_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "pointFields.H"
#include "labelVector.H"
#include "DynamicList.H"
#include "Map.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class ExaCA Declaration
\*---------------------------------------------------------------------------*/

class ExaCA
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Temperature field
        const volScalarField& T_;

        //- Cells overlapping the ExaCA box for the current mesh state
        DynamicList<label> overlapCells;

        //- Current mesh cell label to cached cell map id
        Map<label> cellMapIds;

        //- Current mesh data need rebuilding
        bool meshChanged_;

        //- Number of CA points in each direction
        labelVector nPoints_;

        //- Cached CA-point map for one cell
        struct cellMap
        {
            DynamicList<label> pointIds;

            DynamicList<point> positions;

            DynamicList<scalarField> weights;
        };

        //- Thermal event interval in a mapped cell
        struct eventInterval
        {
            label mapId;

            scalar t0;

            scalar t1;

            List<scalar> psi0;

            List<scalar> psi;
        };

        //- Cached cell maps retained by event intervals
        DynamicList<cellMap> cellMaps;

        //- Event intervals collected during the solve
        DynamicList<eventInterval> intervals;

        //- Last melt time for each fixed CA grid-point id
        Map<scalar> meltTimes;

        //- Reduced ExaCA output data: x,y,z,tm,ts,cr
        DynamicList<List<scalar>> data;

        //- ExaCA sampling box
        boundBox box_;

        //- Isotherm value
        scalar isoValue_;

        //- CA grid spacing
        scalar dx_;


    // Private Member Functions

        //- Reset current mesh lookup data
        void resetMeshState();

        //- Update current mesh lookup data
        void updateMeshState();

        //- Return a stable id for a CA grid point
        label pointId
        (
            const label i,
            const label j,
            const label k
        ) const;

        //- Return cached cell map id, constructing the map if needed
        label getCellMap(const label celli);

        //- Construct the CA-point map for one current mesh cell
        cellMap makeCellMap(const label celli) const;

        //- Append an event interval for one current mesh cell
        void appendInterval
        (
            const label celli,
            const pointScalarField& Tp0,
            const pointScalarField& Tp,
            const scalar t0,
            const scalar t
        );

        //- Write reduced ExaCA data
        void writeData() const;


public:

    //- Runtime type information
    TypeName("ExaCA");


    // Constructors

        //- Construct from Time and dictionary
        ExaCA
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ExaCA(const ExaCA&) = delete;


    //- Destructor
    virtual ~ExaCA();


    // Member Functions

        //- Read the ExaCA data
        virtual bool read(const dictionary&);

        //- Initialize overlap-cell list
        virtual void setOverlapCells();

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute
        virtual bool execute();

        //- Interpolate cached event intervals and write ExaCA data
        virtual void interpolate();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write the ExaCA data
        virtual bool write();

        //- Update after mesh points move
        virtual void movePoints(const polyMesh& mesh);

        //- Update after mesh topology changes
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update after mapping from another mesh
        virtual void mapMesh(const polyMeshMap& map);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ExaCA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
