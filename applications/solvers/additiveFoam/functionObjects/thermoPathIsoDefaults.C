/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "thermoPathIsoDefaults.H"
#include "readThermoPath.H"
#include "graph.H"
#include "interpolateXY.H"

namespace Foam
{
namespace functionObjects
{

scalar thermoPathLiquidus(const fvMesh& mesh)
{
    autoPtr<graph> thermoPtr(readThermoPath(mesh));
    const graph& thermo = thermoPtr();

    return interpolateXY(0.0, thermo.y(), thermo.x());
}


scalar thermoPathSolidus(const fvMesh& mesh)
{
    autoPtr<graph> thermoPtr(readThermoPath(mesh));
    const graph& thermo = thermoPtr();

    return interpolateXY(1.0, thermo.y(), thermo.x());
}


scalarList thermoPathTemperatureList(const fvMesh& mesh)
{
    autoPtr<graph> thermoPtr(readThermoPath(mesh));
    const graph& thermo = thermoPtr();

    return scalarList(thermo.x());
}

}
}

// ************************************************************************* //
