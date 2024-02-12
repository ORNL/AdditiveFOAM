/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "Timer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Timer::Timer() : totalTime(0.0), name("") {}

Timer::Timer(const std::string& timerName) : totalTime(0.0), name(timerName)
{
    start();
}

void Timer::start()
{
    startTime = std::clock();
}

void Timer::stop()
{
    totalTime += static_cast<double>(std::clock() - startTime) / CLOCKS_PER_SEC;
}

std::string Timer::getName() const
{
    return name;
}

const std::string& Timer::getName()
{
    return name;
}

double Timer::getTotalTime() const
{
    return totalTime;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Timers::Timers(const Time& runTime) : runTime_(runTime) {}

void Timers::createTimer(const std::string& name)
{
    timers.append(Timer(name));
}

void Timers::start(const std::string& name)
{
    bool new_timer = true;

    for (auto& timer : timers)
    {
        if (timer.getName() == name)
        {
            timer.start();
            new_timer = false;
        }
    }

    if (new_timer)
    {
        timers.append(Timer(name));
    }
}

void Timers::stop(const std::string& name)
{
    for (auto& timer : timers)
    {
        if (timer.getName() == name)
        {
            timer.stop();
        }
    }
}

void Timers::write() const
{
    const fileName timerPath(runTime_.rootPath()/runTime_.globalCaseName()/"Profiling");

    mkDir(timerPath);

    OFstream os(timerPath + "/" + "timers_" + Foam::name(Pstream::myProcNo()) + ".csv");

    for (const auto& timer : timers)
    {
        os << word(timer.getName()) << ",";
    }
    os << "elapsedCpuTime" << "\n";

    for (const auto& timer : timers)
    {
        os << timer.getTotalTime() << ",";
    }
    os << runTime_.elapsedCpuTime() << "\n";
}
// ************************************************************************* //
