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
InNamespace
    Foam
Description
    Calculates the cpu time on each processor for an associated task.
SourceFiles
    Timer.C
\*---------------------------------------------------------------------------*/

#ifndef Timer_H
#define Timer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <ctime>
#include <string>
#include "OFstream.H"
#include "Time.H"
#include "DynamicList.H"
#include "OSspecific.H"

using namespace Foam;

class Timer
{
    private:
        std::clock_t startTime;
        double totalTime;
        std::string name;

    public:
        Timer();
        Timer(const std::string& timerName);

        void start();
        void stop();

        std::string getName() const;
        const std::string& getName();
        double getTotalTime() const;
};

class Timers 
{
    private:
        DynamicList<Timer> timers;
        const Time& runTime_;

    public:
        Timers(const Time& runTime);

        void createTimer(const std::string& name);
        void start(const std::string& name);
        void stop(const std::string& name);
        void write() const;
};

#endif // Timer_H
