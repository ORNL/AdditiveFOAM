#include <gtest/gtest.h>

#include "OStringStream.H"
#include "graph.H"
#include "interpolateXY.H"

namespace
{

} // namespace

TEST(utilityTests, interpolateXYHandlesExactHitsInterpolationClampingAndUnsortedXData)
{
    Foam::scalarField xOld(3);
    xOld[0] = 3.0;
    xOld[1] = 1.0;
    xOld[2] = 2.0;

    Foam::scalarField yOld(3);
    yOld[0] = 30.0;
    yOld[1] = 10.0;
    yOld[2] = 20.0;

    EXPECT_NEAR(Foam::interpolateXY(2.0, xOld, yOld), 20.0, 1e-9);
    EXPECT_NEAR(Foam::interpolateXY(1.5, xOld, yOld), 15.0, 1e-9);
    EXPECT_NEAR(Foam::interpolateXY(0.0, xOld, yOld), 10.0, 1e-9);
    EXPECT_NEAR(Foam::interpolateXY(4.0, xOld, yOld), 30.0, 1e-9);

    const Foam::labelPair labels = Foam::interpolateXYLabels(1.5, xOld, yOld);
    EXPECT_EQ(labels.first(), 1);
    EXPECT_EQ(labels.second(), 2);
}

TEST(utilityTests, graphWordifyNormalizesLabelsAndYReturnsTheOnlyCurve)
{
    Foam::scalarField x(2);
    x[0] = 0.0;
    x[1] = 1.0;

    Foam::scalarField y(2);
    y[0] = 2.0;
    y[1] = 3.0;

    EXPECT_EQ(Foam::graph::wordify("Melt Pool (mm)"), Foam::word("Melt_Pool__mm"));

    Foam::graph g("title", "x", "Melt Pool (mm)", x, y);

    EXPECT_DOUBLE_EQ(g.y()[0], Foam::scalar(2.0));
    EXPECT_DOUBLE_EQ(g.y()[1], Foam::scalar(3.0));
}

TEST(utilityTests, graphWriteTableEmitsTheStoredXYPairs)
{
    Foam::scalarField x(2);
    x[0] = 0.0;
    x[1] = 1.0;

    Foam::scalarField y(2);
    y[0] = 2.0;
    y[1] = 3.0;

    Foam::graph g("title", "x", "y", x, y);
    Foam::OStringStream os;
    g.writeTable(os);

    const std::string output = os.str();

    EXPECT_NE(output.find("0"), std::string::npos);
    EXPECT_NE(output.find("1"), std::string::npos);
    EXPECT_NE(output.find("2"), std::string::npos);
    EXPECT_NE(output.find("3"), std::string::npos);
}
