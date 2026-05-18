#include <string>

#include <gtest/gtest.h>
#include "segment.H"

TEST(segmentTests, defaultConstructionYieldsAZeroedPointSource)
{
    Foam::segment seg;

    EXPECT_DOUBLE_EQ(seg.mode(), Foam::scalar(1));
    EXPECT_DOUBLE_EQ(seg.position().x(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.position().y(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.position().z(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.power(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.parameter(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.time(), Foam::scalar(0));
}

TEST(segmentTests, parsesASpaceDelimitedSegmentDefinition)
{
    Foam::segment seg(std::string("0 1 2 3 400 5"));

    EXPECT_DOUBLE_EQ(seg.mode(), Foam::scalar(0));
    EXPECT_DOUBLE_EQ(seg.position().x(), Foam::scalar(1));
    EXPECT_DOUBLE_EQ(seg.position().y(), Foam::scalar(2));
    EXPECT_DOUBLE_EQ(seg.position().z(), Foam::scalar(3));
    EXPECT_DOUBLE_EQ(seg.power(), Foam::scalar(400));
    EXPECT_DOUBLE_EQ(seg.parameter(), Foam::scalar(5));
}
