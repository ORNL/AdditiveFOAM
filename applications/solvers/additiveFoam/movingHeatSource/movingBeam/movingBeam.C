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

#include "movingBeam.H"

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
    path_(1, segment()),
    index_(0),
    position_(Zero),
    power_(0.0),
    endTime_(0.0),
    deltaT_(GREAT),
    hitPathIntervals_(true)
{
    //- Get beam parameters
    deltaT_ = beamDict_.lookupOrDefault<scalar>("deltaT", GREAT);
    hitPathIntervals_ = beamDict_.lookupOrDefault<bool>("hitPathIntervals", true);
    
    //- Read scan path file
    readPath();
    
    //- Initialize the path index
    index_ = findIndex(runTime_.value());
    
    Info << "Initial path index: " << index_ << endl;
    
    //- Find the beam end time
    for (label i = path_.size() - 1; i > 0; i--)
    {
        if (path_[i].power() > small)
        {
            endTime_ = min(path_[i].time(), runTime_.endTime().value());
            break;
        }
    }

    Info << "Path end time: " << endTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector
Foam::movingBeam::position(const scalar time) const
{
    if (activePath(time))
    {
        const label i = findIndex(time);
        
        if (path_[i].mode() == 1)
        {
            return path_[i].position();
        }
        else
        {
            vector displacement = vector::zero;

            scalar dt = path_[i].time() - path_[i-1].time();

            if (dt > 0)
            {
                const vector dx = path_[i].position() - path_[i-1].position();
                displacement = dx*(time - path_[i-1].time())/dt;
            }

            return path_[i-1].position() + displacement;
        }
    }
    else
    {
        return vector::max;
    }
}


Foam::scalar
Foam::movingBeam::power(const scalar time) const
{
    if (activePath(time))
    {
        const label i = findIndex(time);
        
        if ((time - path_[i-1].time()) > eps)
        {
            return path_[i].power();
        }
        else
        {
            return path_[i-1].power();
        }
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar
Foam::movingBeam::mode(const scalar time) const
{
    if (activePath(time))
    {
        const label i = findIndex(time);
        
        return path_[i].mode();
    }
    else
    {
        return 1.0;
    }
}


Foam::scalar
Foam::movingBeam::parameter(const scalar time) const
{
    if (activePath(time))
    {
        const label i = findIndex(time);
        
        return path_[i].parameter();
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar
Foam::movingBeam::velocity(const scalar time) const
{
    if (activePath(time))
    {
        const label i = findIndex(time);
        
        if (path_[i].mode() == 1)
        {
            return 0.0;
        }
        
        else
        {
            return path_[i].parameter();
        }
    }
    else
    {
        return 0.0;
    }
}


Foam::scalar
Foam::movingBeam::timeToNextPath(const scalar time) const
{
    if (activePath(time))
    {
        label i = findIndex(time);
        
        scalar timeToNextPath = 0;

        while (timeToNextPath < eps)
        {
            timeToNextPath = max(0, path_[i].time() - time);

            i++;

            if (i == path_.size()) 
            {
                break;
            }
        }
        
        return timeToNextPath;
    }
    else
    {
        return GREAT;
    }
}


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

    while (std::getline(is, line))
    {
        if (line.empty())
        {
            continue;
        }

        path_.append(segment(line));
    }

    for (label i=1; i < path_.size(); i++)
    {
        if (path_[i].mode() == 1)
        {
            path_[i].setTime
            (
                path_[i-1].time() + path_[i].parameter()
            );
        }
        else
        {
            scalar d_ = mag(path_[i].position() - path_[i-1].position());
                    
            path_[i].setTime
            (
                path_[i-1].time() + d_/path_[i].parameter()
            );
        }

        Info << i << tab << path_[i].time() << endl;
    }
}


bool Foam::movingBeam::activePath() const
{
    return ((endTime_ - runTime_.value()) > eps);
}


bool Foam::movingBeam::activePath(const scalar time) const
{
    return ((endTime_ - time) > eps);
}


void Foam::movingBeam::move(const scalar time) const
{
    // update the current index of the path
    index_ = findIndex(time);

    const label i = index_;

    // update the beam center
    if (path_[i].mode() == 1)
    {
        position_ = path_[i].position();
    }
    else
    {
        vector displacement = vector::zero;

        scalar dt = path_[i].time() - path_[i-1].time();

        if (dt > 0)
        {
            const vector dx = path_[i].position() - path_[i-1].position();
            displacement = dx*(time - path_[i-1].time())/dt;
        }

        position_ = path_[i-1].position() + displacement;
    }

    // update the beam power
    if ((time - path_[i-1].time()) > eps)
    {
        power_ = path_[i].power();
    }
    else
    {
        power_ = path_[i-1].power();
    }
}


Foam::label
Foam::movingBeam::findIndex(const scalar time) const
{
    label i = index_;

    const label n = path_.size() - 1;

    // step back path index for safe updating
    for (i = i; i > 0 && path_[i].time() > time; --i)
    {}

    // update the path index to the provided time
    for (i = i; i < n && path_[i].time() < time; ++i)
    {}

    // skip any point sources with zero time
    while (i < n)
    {
        if (path_[i].mode() == 1 && path_[i].parameter() == 0)
        {
            ++i;
        }
        else
        {
            break;
        }
    }

    return min(max(i, 0), n);
}


void Foam::movingBeam::adjustDeltaT(scalar& dt) const
{
    if (activePath() && hitPathIntervals_)
    {
        const scalar dtMax = timeToNextPath(runTime_.value());
        const scalar nSteps = dtMax/dt;
                
        if (nSteps < labelMax)
        {
            // allow time step to dilate 1% to hit target path time
            const label nStepsToNextPath = label(max(nSteps, 1) + 0.99);
            dt = min(dtMax/nStepsToNextPath, dt);
        }
    }
}

// ************************************************************************* //
