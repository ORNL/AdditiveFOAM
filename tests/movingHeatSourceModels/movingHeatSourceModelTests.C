#include <cmath>

#include <gtest/gtest.h>

#include "KellyAbsorption.H"
#include "constantAbsorption.H"
#include "modifiedSuperGaussian.H"
#include "movingHeatSourceTestFixture.H"
#include "projectedGaussian.H"
#include "superGaussian.H"

namespace
{

Foam::IOdictionary makeHeatSourceDict(Foam::Time& runTime)
{
    return additiveFoamTest::makeHeatSourceDict(runTime);
}

Foam::scalar expectedKellyEta
(
    const Foam::word& geometry,
    const Foam::scalar aspectRatio,
    const Foam::scalar eta0,
    const Foam::scalar etaMin
)
{
    if (aspectRatio <= 1.0)
    {
        return etaMin;
    }

    const Foam::scalar theta = Foam::atan(1.0 / aspectRatio);

    Foam::scalar F = 0.0;
    Foam::scalar G = 0.0;

    if (geometry == "cone")
    {
        F = 0.25 * (3.0 * Foam::sin(theta) - Foam::sin(3.0 * theta));
        G = 1.0 / (1.0 + Foam::sqrt(1.0 + Foam::pow(aspectRatio, 2)));
    }
    else
    {
        F = 0.5 * (1.0 - Foam::cos(2.0 * theta));
        G = 0.5 / (1.0 + aspectRatio);
    }

    return eta0 * (1.0 + (1.0 - eta0) * (G - F))
         / (1.0 - (1.0 - eta0) * (1.0 - G));
}

Foam::scalar projectedK
(
    const Foam::scalar aspectRatio,
    const Foam::scalar A,
    const Foam::scalar B
)
{
    const Foam::scalar n = Foam::min(Foam::max(A * std::log2(aspectRatio) + B, 0.0), 9.0);
    return std::pow(2.0, n);
}

} // namespace

TEST(movingHeatSourceModelTests, constantAbsorptionReturnsTheConfiguredEtaForAnyAspectRatio)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::absorptionModels::constant model("testBeam", heatSourceDict, mesh);

    EXPECT_NEAR(model.eta(0.5), 0.35, 1e-9);
    EXPECT_NEAR(model.eta(7.5), 0.35, 1e-9);
}

TEST(movingHeatSourceModelTests, KellyAbsorptionMatchesTheConeAndCylinderFormulasAndRespectsEtaMin)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::absorptionModels::Kelly cone("kellyConeBeam", heatSourceDict, mesh);
    Foam::absorptionModels::Kelly cylinder("kellyCylinderBeam", heatSourceDict, mesh);

    const Foam::scalar aspectRatio = 2.0;

    EXPECT_NEAR(cone.eta(aspectRatio), expectedKellyEta("cone", aspectRatio, 0.45, 0.15), 1e-9);
    EXPECT_NEAR(cylinder.eta(aspectRatio), expectedKellyEta("cylinder", aspectRatio, 0.45, 0.15), 1e-9);
    EXPECT_NEAR(cone.eta(1.0), 0.15, 1e-9);
    EXPECT_NEAR(cylinder.eta(0.75), 0.15, 1e-9);
}

TEST(movingHeatSourceModelTests, superGaussianWeightIsCenteredAndSymmetricAndV0MatchesTheNormalizationFormula)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::heatSourceModels::superGaussian model
    (
        additiveFoamTest::suppressStdout
        (
            [&]()
            {
                return Foam::heatSourceModels::superGaussian("testBeam", heatSourceDict, mesh);
            }
        )
    );

    EXPECT_NEAR(model.weight(Foam::vector::zero), 1.0, 1e-9);
    EXPECT_NEAR(model.weight(Foam::vector(1.0, 0.0, 0.0)), model.weight(Foam::vector(-1.0, 0.0, 0.0)), 1e-9);
    EXPECT_LT(model.weight(Foam::vector(1.0, 0.0, 0.0)), 1.0);

    const Foam::scalar k = 2.0;
    const Foam::scalar a = Foam::pow(2.0, 1.0 / k);
    const Foam::vector s = Foam::vector(2.0, 2.0, 4.0) / a;
    const Foam::scalar expectedV0 =
        (2.0 / 3.0) * s.x() * s.y() * s.z()
      * Foam::constant::mathematical::pi * Foam::tgamma(1.0 + 3.0 / k);

    EXPECT_NEAR(model.V0().value(), expectedV0, 1e-9);
}

TEST(movingHeatSourceModelTests, modifiedSuperGaussianTruncatesBeyondTheBeamDepthAndRemainsSymmetricInPlane)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::heatSourceModels::modifiedSuperGaussian model
    (
        additiveFoamTest::suppressStdout
        (
            [&]()
            {
                return Foam::heatSourceModels::modifiedSuperGaussian("modifiedBeam", heatSourceDict, mesh);
            }
        )
    );

    EXPECT_NEAR(model.weight(Foam::vector::zero), 1.0, 1e-9);
    EXPECT_NEAR(model.weight(Foam::vector(0.5, 0.0, 1.0)), model.weight(Foam::vector(-0.5, 0.0, 1.0)), 1e-9);
    EXPECT_NEAR(model.weight(Foam::vector(0.0, 0.0, 4.0)), 0.0, 1e-9);
    EXPECT_GT(model.V0().value(), 0.0);
}

TEST(movingHeatSourceModelTests, projectedGaussianClampsTheDerivedExponentAndDecaysAwayFromTheCenter)
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::heatSourceModels::projectedGaussian model
    (
        additiveFoamTest::suppressStdout
        (
            [&]()
            {
                return Foam::heatSourceModels::projectedGaussian("projectedBeam", heatSourceDict, mesh);
            }
        )
    );

    EXPECT_NEAR(model.weight(Foam::vector::zero), 1.0, 1e-9);
    EXPECT_LT(model.weight(Foam::vector(0.25, 0.0, 0.0)), 1.0);
    EXPECT_LT(model.weight(Foam::vector(0.0, 0.0, 16.0)), model.weight(Foam::vector(0.0, 0.0, 1.0)));

    const Foam::scalar aspectRatio = 16.0 / Foam::min(1.0, 2.0);
    const Foam::scalar k = projectedK(aspectRatio, 3.0, 0.0);
    const Foam::scalar expectedV0 =
        0.5 * Foam::constant::mathematical::pi * 1.0 * 2.0 * 16.0
      * Foam::tgamma(1.0 / k) / (k * std::pow(3.0, 1.0 / k));

    EXPECT_NEAR(k, 512.0, 1e-9);
    EXPECT_NEAR(model.V0().value(), expectedV0, 1e-9);
}
