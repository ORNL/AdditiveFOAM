/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#ifndef movingHeatSource_H
#define movingHeatSource_H

#include "absorptionModel.H"
#include "heatSourceModel.H"
#include "movingBeam.H"
#include "referenceDimensions.H"
#include "volFields.H"

namespace Foam
{

class movingHeatSource
{
    word name_;
    const fvMesh& mesh_;
    movingBeam beam_;
    autoPtr<absorptionModel> absorption_;
    autoPtr<heatSourceModel> heatSource_;
    referenceDimensions references_;
public:

    movingHeatSource
    (
        const word& name,
        const dictionary& sourceDict,
        const fvMesh& mesh
    );

    movingHeatSource(const movingHeatSource&) = delete;

    inline const word& name() const
    {
        return name_;
    }

    inline movingBeam& beam()
    {
        return beam_;
    }

    inline const movingBeam& beam() const
    {
        return beam_;
    }

    inline const heatSourceModel& heatSource() const
    {
        return heatSource_();
    }

    void update();

    tmp<volScalarField> qDot() const;

    void operator=(const movingHeatSource&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
