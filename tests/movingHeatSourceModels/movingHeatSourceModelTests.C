#include <cmath>

#include "doctest.h"

#include "KellyAbsorption.H"
#include "constantAbsorption.H"
#include "modifiedSuperGaussian.H"
#include "movingHeatSourceTestFixture.H"
#include "projectedGaussian.H"
#include "superGaussian.H"

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

TEST_CASE("constant absorption returns the configured eta for any aspect ratio")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::absorptionModels::constant model("testBeam", heatSourceDict, mesh);

    CHECK(scalarClose(model.eta(0.5), 0.35));
    CHECK(scalarClose(model.eta(7.5), 0.35));
}

TEST_CASE("Kelly absorption matches the cone and cylinder formulas and respects etaMin")
{
    auto runTime = additiveFoamTest::makeTime();
    Foam::fvMesh mesh(additiveFoamTest::makeFixtureMesh(*runTime));
    Foam::IOdictionary heatSourceDict(makeHeatSourceDict(*runTime));

    Foam::absorptionModels::Kelly cone("kellyConeBeam", heatSourceDict, mesh);
    Foam::absorptionModels::Kelly cylinder("kellyCylinderBeam", heatSourceDict, mesh);

    const Foam::scalar aspectRatio = 2.0;

    CHECK(scalarClose(cone.eta(aspectRatio), expectedKellyEta("cone", aspectRatio, 0.45, 0.15)));
    CHECK(scalarClose(cylinder.eta(aspectRatio), expectedKellyEta("cylinder", aspectRatio, 0.45, 0.15)));
    CHECK(scalarClose(cone.eta(1.0), 0.15));
    CHECK(scalarClose(cylinder.eta(0.75), 0.15));
}

TEST_CASE("superGaussian weight is centered and symmetric and V0 matches the normalization formula")
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

    CHECK(scalarClose(model.weight(Foam::vector::zero), 1.0));
    CHECK(scalarClose(model.weight(Foam::vector(1.0, 0.0, 0.0)), model.weight(Foam::vector(-1.0, 0.0, 0.0))));
    CHECK(model.weight(Foam::vector(1.0, 0.0, 0.0)) < 1.0);

    const Foam::scalar k = 2.0;
    const Foam::scalar a = Foam::pow(2.0, 1.0 / k);
    const Foam::vector s = Foam::vector(2.0, 2.0, 4.0) / a;
    const Foam::scalar expectedV0 =
        (2.0 / 3.0) * s.x() * s.y() * s.z()
      * Foam::constant::mathematical::pi * Foam::tgamma(1.0 + 3.0 / k);

    CHECK(scalarClose(model.V0().value(), expectedV0));
}

TEST_CASE("modifiedSuperGaussian truncates beyond the beam depth and remains symmetric in-plane")
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

    CHECK(scalarClose(model.weight(Foam::vector::zero), 1.0));
    CHECK(scalarClose(model.weight(Foam::vector(0.5, 0.0, 1.0)), model.weight(Foam::vector(-0.5, 0.0, 1.0))));
    CHECK(scalarClose(model.weight(Foam::vector(0.0, 0.0, 4.0)), 0.0));
    CHECK(model.V0().value() > 0.0);
}

TEST_CASE("projectedGaussian clamps the derived exponent and decays away from the center")
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

    CHECK(scalarClose(model.weight(Foam::vector::zero), 1.0));
    CHECK(model.weight(Foam::vector(0.25, 0.0, 0.0)) < 1.0);
    CHECK(model.weight(Foam::vector(0.0, 0.0, 16.0)) < model.weight(Foam::vector(0.0, 0.0, 1.0)));

    const Foam::scalar aspectRatio = 16.0 / Foam::min(1.0, 2.0);
    const Foam::scalar k = projectedK(aspectRatio, 3.0, 0.0);
    const Foam::scalar expectedV0 =
        0.5 * Foam::constant::mathematical::pi * 1.0 * 2.0 * 16.0
      * Foam::tgamma(1.0 / k) / (k * std::pow(3.0, 1.0 / k));

    CHECK(scalarClose(k, 512.0));
    CHECK(scalarClose(model.V0().value(), expectedV0));
}
