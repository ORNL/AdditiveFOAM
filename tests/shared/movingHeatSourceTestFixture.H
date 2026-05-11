#ifndef movingHeatSourceTestFixture_H
#define movingHeatSourceTestFixture_H

#include <cstdio>
#include <fcntl.h>
#include <memory>
#include <unistd.h>
#include <utility>

#include "Time.H"
#include "IOdictionary.H"
#include "zeroDimensionalFvMesh.H"

namespace additiveFoamTest
{

static const Foam::fileName fixtureRootPath(".");
static const Foam::fileName fixtureCaseName("fixtures/movingHeatSourceCase");

class ScopedStdoutSilencer
{
    int savedFd_;
    int nullFd_;

public:
    ScopedStdoutSilencer()
    :
        savedFd_(-1),
        nullFd_(-1)
    {
        std::fflush(stdout);

        savedFd_ = dup(STDOUT_FILENO);
        nullFd_ = open("/dev/null", O_WRONLY);

        if (savedFd_ >= 0 && nullFd_ >= 0)
        {
            dup2(nullFd_, STDOUT_FILENO);
        }
    }

    ~ScopedStdoutSilencer()
    {
        std::fflush(stdout);

        if (savedFd_ >= 0)
        {
            dup2(savedFd_, STDOUT_FILENO);
            close(savedFd_);
        }

        if (nullFd_ >= 0)
        {
            close(nullFd_);
        }
    }

    ScopedStdoutSilencer(const ScopedStdoutSilencer&) = delete;
    ScopedStdoutSilencer& operator=(const ScopedStdoutSilencer&) = delete;
};

template<class Fn>
inline auto suppressStdout(Fn&& fn) -> decltype(fn())
{
    ScopedStdoutSilencer silencer;
    return fn();
}

inline std::unique_ptr<Foam::Time> makeTime()
{
    return suppressStdout
    (
        []()
        {
            return std::unique_ptr<Foam::Time>
            (
                new Foam::Time
                (
                    Foam::Time::controlDictName,
                    fixtureRootPath,
                    fixtureCaseName,
                    false
                )
            );
        }
    );
}

inline Foam::IOdictionary makeHeatSourceDict(Foam::Time& runTime)
{
    return Foam::IOdictionary
    (
        Foam::IOobject
        (
            "heatSourceDict",
            runTime.constant(),
            runTime,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        )
    );
}

inline Foam::fvMesh makeFixtureMesh(Foam::Time& runTime)
{
    return Foam::zeroDimensionalFvMesh(runTime);
}

} // namespace additiveFoamTest

#endif
