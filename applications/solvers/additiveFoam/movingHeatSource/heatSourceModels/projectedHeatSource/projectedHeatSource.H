/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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
    Foam::heatSourceModels::projectedHeatSource

Description
    Separable heat source constructed from a planar profile and a one-sided
    axial projection.

SourceFiles
    projectedHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef projectedHeatSource_H
#define projectedHeatSource_H

#include "heatSourceModel.H"
#include "heatSourceProfile.H"
#include "heatSourceProjection.H"

namespace Foam
{
namespace heatSourceModels
{

class projectedHeatSource
:
    public heatSourceModel
{
    // Private Data

        //- Planar heat source profile
        autoPtr<heatSourceProfile> profile_;

        //- Axial heat source projection
        autoPtr<heatSourceProjection> projection_;

        //- Maximum power fraction outside each separable component bound
        scalar componentTolerance_;


    // Private Member Functions

        //- Update the projection, source bounds and normalization
        void updateSource();

public:

    //- Runtime type information
    TypeName("projectedHeatSource");


    // Constructors

        projectedHeatSource
        (
            const word& sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~projectedHeatSource()
    {}

    // Member Functions

        virtual void update();

        inline virtual scalar weight(const vector& r) const;

        //- Read the heatSourceDict dictionary
        virtual bool read();
};

} // End namespace heatSourceModels
} // End namespace Foam

#endif

// ************************************************************************* //
