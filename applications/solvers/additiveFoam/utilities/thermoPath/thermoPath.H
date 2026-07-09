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

Class
    Foam::thermoPath

Description
    Phase-change path reader for AdditiveFOAM material properties.

    The path is read from the thermoPath entry in constant/transportProperties.
    If the entry is not present, the legacy constant/thermoPath file is read.

SourceFiles
    thermoPath.C

\*---------------------------------------------------------------------------*/

#ifndef thermoPath_H
#define thermoPath_H

#include "scalarField.H"
#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Time.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermoPath Declaration
\*---------------------------------------------------------------------------*/

class thermoPath
{
    // Private Data

        //- Temperatures from the phase-change path
        scalarField T_;

        //- Solid fractions from the phase-change path
        scalarField alpha_;

        //- Temperature where alpha.solid = 1
        scalar solidus_;

        //- Temperature where alpha.solid = 0
        scalar liquidus_;

        //- Source name used for diagnostics
        fileName sourceName_;


    // Private Member Functions

        //- Read thermoPath from transportProperties or legacy thermoPath file
        void read
        (
            const IOdictionary& transportProperties,
            const Time& runTime
        );

        //- Read legacy constant/thermoPath file
        void readLegacyFile(const Time& runTime);

        //- Set fields from a path list
        void set
        (
            const List<Pair<scalar>>& path,
            const fileName& sourceName
        );

        //- Return temperature at the requested solid fraction
        scalar temperatureAtSolidFraction(const scalar alpha) const;


public:

    // Static Data

        //- Default material dictionary name
        static const word transportPropertiesName;

        //- Thermodynamic path entry name
        static const word thermoPathName;


    // Constructors

        //- Construct from registered transportProperties
        thermoPath(const objectRegistry& db);

        //- Construct from transportProperties
        thermoPath
        (
            const IOdictionary& transportProperties,
            const Time& runTime
        );


    // Member Functions

        //- Return temperatures
        inline const scalarField& x() const
        {
            return T_;
        }

        //- Return solid fractions
        inline const scalarField& y() const
        {
            return alpha_;
        }

        //- Return solidus temperature
        inline scalar solidus() const
        {
            return solidus_;
        }

        //- Return liquidus temperature
        inline scalar liquidus() const
        {
            return liquidus_;
        }

        //- Return source name
        inline const fileName& sourceName() const
        {
            return sourceName_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
