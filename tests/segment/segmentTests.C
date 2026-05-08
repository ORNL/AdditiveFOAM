#include <string>

#include "doctest.h"
#include "segment.H"

TEST_CASE("Foam::segment default construction yields a zeroed point source")
{
    Foam::segment seg;

    CHECK_EQ(seg.mode(), Foam::scalar(1));
    CHECK_EQ(seg.position().x(), Foam::scalar(0));
    CHECK_EQ(seg.position().y(), Foam::scalar(0));
    CHECK_EQ(seg.position().z(), Foam::scalar(0));
    CHECK_EQ(seg.power(), Foam::scalar(0));
    CHECK_EQ(seg.parameter(), Foam::scalar(0));
    CHECK_EQ(seg.time(), Foam::scalar(0));
}

TEST_CASE("Foam::segment parses a space-delimited segment definition")
{
    Foam::segment seg(std::string("0 1 2 3 400 5"));

    CHECK_EQ(seg.mode(), Foam::scalar(0));
    CHECK_EQ(seg.position().x(), Foam::scalar(1));
    CHECK_EQ(seg.position().y(), Foam::scalar(2));
    CHECK_EQ(seg.position().z(), Foam::scalar(3));
    CHECK_EQ(seg.power(), Foam::scalar(400));
    CHECK_EQ(seg.parameter(), Foam::scalar(5));
}
