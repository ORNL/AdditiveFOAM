#include <gtest/gtest.h>

#include "movingBeam.H"
#include "movingHeatSourceTestFixture.H"

namespace
{

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

TEST(movingBeamTests, computesPathTimesAndActivityFromTheFixtureScanPath)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    EXPECT_EQ(beam.findIndex(0.0), 1);
    EXPECT_EQ(beam.findIndex(1.5), 1);
    EXPECT_EQ(beam.findIndex(1.500001), 2);
    EXPECT_EQ(beam.findIndex(3.000001), 3);
    EXPECT_TRUE(beam.activePath());

    runTime->setTime(3.6, 0);
    EXPECT_FALSE(beam.activePath());
}

TEST(movingBeamTests, skipsZeroDurationPointSourcesWhenLocatingTheActiveSegment)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "skipBeam"));

    EXPECT_EQ(beam.findIndex(0.0), 2);
    EXPECT_EQ(beam.findIndex(0.5), 2);
}

TEST(movingBeamTests, moveInterpolatesTravelSegmentsAndSwitchesPowerOnBoundaries)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    beam.move(0.0);
    EXPECT_DOUBLE_EQ(beam.position().x(), Foam::scalar(1.0));
    EXPECT_DOUBLE_EQ(beam.power(), Foam::scalar(0.0));

    beam.move(2.25);
    EXPECT_NEAR(beam.position().x(), 2.5, 1e-9);
    EXPECT_NEAR(beam.position().y(), 0.0, 1e-9);
    EXPECT_NEAR(beam.power(), 200.0, 1e-9);

    beam.move(3.0);
    EXPECT_NEAR(beam.position().x(), 4.0, 1e-9);
    EXPECT_NEAR(beam.power(), 200.0, 1e-9);

    beam.move(3.000001);
    EXPECT_NEAR(beam.position().x(), 4.0, 1e-9);
    EXPECT_NEAR(beam.power(), 50.0, 1e-9);
}

TEST(movingBeamTests, adjustDeltaTLandsOnTheNextPathIntervalWhenEnabled)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "testBeam"));

    Foam::scalar dt = 1.0;
    beam.adjustDeltaT(dt);

    EXPECT_NEAR(dt, 0.75, 1e-9);
}

TEST(movingBeamTests, adjustDeltaTLeavesTheTimestepUnchangedWhenPathHitsAreDisabled)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::movingBeam beam(makeBeam(*runTime, "noHitBeam"));

    Foam::scalar dt = 1.0;
    beam.adjustDeltaT(dt);

    EXPECT_NEAR(dt, 1.0, 1e-9);
}
