#include "doctest.h"

#include "movingBeam.H"
#include "movingHeatSourceTestFixture.H"

namespace
{

bool scalarClose
(
    const Foam::scalar lhs,
    const Foam::scalar rhs,
    const Foam::scalar tol = 1e-9
)
{
    return Foam::mag(lhs - rhs) <= tol;
}

Foam::movingBeam makeBeam(Foam::Time& runTime, const char* sourceName)
{
    Foam::IOdictionary heatSourceDict(additiveFoamTest::makeHeatSourceDict(runTime));
    return additiveFoamTest::suppressStdout
    (
        [&]()
        {
            return Foam::movingBeam(sourceName, heatSourceDict, runTime);
        }
    );
}

} // namespace

TEST_CASE("Foam::movingBeam computes path times and activity from the fixture scan path")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    CHECK_EQ(beam.findIndex(0.0), 1);
    CHECK_EQ(beam.findIndex(1.5), 1);
    CHECK_EQ(beam.findIndex(1.500001), 2);
    CHECK_EQ(beam.findIndex(3.000001), 3);
    CHECK(beam.activePath());

    runTime->setTime(3.6, 0);
    CHECK(!beam.activePath());
}

TEST_CASE("Foam::movingBeam skips zero-duration point sources when locating the active segment")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "skipBeam"));

    CHECK_EQ(beam.findIndex(0.0), 2);
    CHECK_EQ(beam.findIndex(0.5), 2);
}

TEST_CASE("Foam::movingBeam move interpolates travel segments and switches power on boundaries")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    beam.move(0.0);
    CHECK_EQ(beam.position().x(), Foam::scalar(1.0));
    CHECK_EQ(beam.power(), Foam::scalar(0.0));

    beam.move(2.25);
    CHECK(scalarClose(beam.position().x(), 2.5));
    CHECK(scalarClose(beam.position().y(), 0.0));
    CHECK(scalarClose(beam.power(), 200.0));

    beam.move(3.0);
    CHECK(scalarClose(beam.position().x(), 4.0));
    CHECK(scalarClose(beam.power(), 200.0));

    beam.move(3.000001);
    CHECK(scalarClose(beam.position().x(), 4.0));
    CHECK(scalarClose(beam.power(), 50.0));
}

TEST_CASE("Foam::movingBeam adjustDeltaT lands on the next path interval when enabled")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    Foam::scalar dt = 1.0;
    beam.adjustDeltaT(dt);

    CHECK(scalarClose(dt, 0.75));
}

TEST_CASE("Foam::movingBeam adjustDeltaT leaves the timestep unchanged when path hits are disabled")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "noHitBeam"));

    Foam::scalar dt = 1.0;
    beam.adjustDeltaT(dt);

    CHECK(scalarClose(dt, 1.0));
}
