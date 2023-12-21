#ifndef Timer_H
#define Timer_H

#include <ctime>
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
    Timer() : totalTime(0.0), name("") {}

    Timer(const std::string& timerName)
    :
        totalTime(0.0),
        name(timerName)
    {
        start();
    }

    void start()
    {
        startTime = std::clock();
    }

    void stop()
    {
        totalTime +=
            static_cast<double>(std::clock() - startTime) / CLOCKS_PER_SEC;
    }

    std::string getName() const
    {
        return name;
    }

    const std::string& getName()
    {
        return name;
    }

    double getTotalTime() const
    {
        return totalTime;
    }
};

class Timers 
{
private:
    DynamicList<Timer> timers;
    const Time& runTime_;

public:
    Timers(const Time& runTime) : runTime_(runTime) {};

    void createTimer(const std::string& name)
    {
        timers.append(Timer(name));
    }

    void start(const std::string& name)
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

    void stop(const std::string& name)
    {
        for (auto& timer : timers)
        {
            if (timer.getName() == name)
            {
                timer.stop();
            }
        }
    }

    void write() const
    {
        const fileName timerPath
        (
            runTime_.rootPath()/runTime_.globalCaseName()/"Profiling"
        );

        mkDir(timerPath);

        OFstream os
        (
            timerPath + "/" + "timers_" + Foam::name(Pstream::myProcNo()) + ".csv"
        );

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
};

#endif
