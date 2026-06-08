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

\*---------------------------------------------------------------------------*/

#include "movingBeam.H"
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::movingBeam::eps = 1e-10;

namespace Foam
{
    defineTypeNameAndDebug(movingBeam, 0);
    defineRunTimeSelectionTable(movingBeam, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingBeam::movingBeam
(
    const word& sourceName,
    const dictionary& dict,
    const Time& runTime
)
:
    sourceName_(sourceName),
    dict_(dict),
    runTime_(runTime),
    beamDict_(dict_.optionalSubDict(sourceName_)),
    path_(0),
    index_(0),
    position_(Zero),
    power_(Zero),
    endTime_(Zero),
    deltaT_(GREAT),
    hitPathIntervals_(true)
{
    //- Get beam parameters
    deltaT_ = beamDict_.lookupOrDefault<scalar>("deltaT", GREAT);

    hitPathIntervals_ = beamDict_.lookupOrDefault<bool>
    (
        "hitPathIntervals",
        true
    );

    //- Read scan path file
    readPath();

    //- Initialize path index
    if (path_.size())
    {
        index_ = getIndex(runTime_.value());

        move(runTime_.value());
    }

    Info << "Initial path index: " << index_ << endl;

    //- Find path end time
    for (label i = path_.size(); i > 0; --i)
    {
        const label pathi = i - 1;

        if (path_[pathi].power() > eps)
        {
            endTime_ = min(path_[pathi].endTime(), runTime_.endTime().value());
            break;
        }
    }

    Info << "Path end time: " << endTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingBeam::readPath()
{
    const word pName_(beamDict_.lookup("pathName"));

    const fileName pFile_
    (
        runTime_.rootPath()/runTime_.globalCaseName()/runTime_.constant()/pName_
    );

    std::ifstream is(pFile_);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot find file " << pFile_
            << nl
            << exit(FatalError);
    }
    else
    {
        Info << "Reading scan path from: " << pFile_ << endl;
    }

    std::string line;

    // skip the header line
    std::getline(is, line);

    scalar time = Zero;

    point position0 = Zero;

    while (std::getline(is, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream lineStream(line);

        scalar mode = Zero;

        point position = Zero;

        scalar power = Zero;

        scalar parameter = Zero;

        lineStream
            >> mode
            >> position.x()
            >> position.y()
            >> position.z()
            >> power
            >> parameter;

        point startPosition = position0;

        scalar dt = Zero;

        //- Spot mode
        if (mode == 1)
        {
            startPosition = position;

            dt = parameter;
        }
        //- Line mode
        else if (mode == 0)
        {
            dt = mag(position - startPosition)/parameter;
        }

        // add path vectors with non-zero duration to the list
        if (dt > eps)
        {
            path_.append
            (
                pathVector
                (
                    startPosition,  // startPosition
                    position,       // endPosition
                    time,           // startTime
                    time + dt,      // endTime
                    power           // power
                )
            );
        }

        time += dt;
        position0 = position;
    }
}


bool Foam::movingBeam::activePath() const
{
    return ((endTime_ - runTime_.value()) > eps);
}


Foam::label
Foam::movingBeam::getIndex(const scalar time) const
{
    label i = index_;

    const label n = path_.size() - 1;

    // step back path index for safe updating
    for (i = i; i > 0 && (time - path_[i].startTime()) <= eps; --i)
    {}

    // update the path index to the provided time
    for (i = i; i < n && (time - path_[i].endTime()) > eps; ++i)
    {}

    return min(max(i, 0), n);
}


void Foam::movingBeam::adjustDeltaT(scalar& dt) const
{
    if (!(activePath() && hitPathIntervals_ && path_.size()))
    {
        return;
    }

    scalar timeToNextPath = 0;

    label i = index_;

    while (i < path_.size() && timeToNextPath < eps)
    {
        timeToNextPath = max(0, path_[i].endTime() - runTime_.value());

        ++i;
    }

    if (timeToNextPath > eps)
    {
        const scalar nSteps = timeToNextPath/dt;

        if (nSteps < labelMax)
        {
            // allow time step to dilate 1% to hit target path time
            const label nStepsToNextPath = label(max(nSteps, 1) + 0.99);

            // reduce the time step to hit the next path vector end time
            dt = min(timeToNextPath/nStepsToNextPath, dt);
        }
    }
}


void Foam::movingBeam::move(const scalar time)
{
    // update the current index of the path
    index_ = getIndex(time);

    const pathVector& pv = path_[index_];

    //- Update the beam center
    position_ = pv.position(time);

    //- Update the beam power
    power_ = pv.power();
}
// ************************************************************************* //
