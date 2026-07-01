/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
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

\*---------------------------------------------------------------------------*/

#include "nLightAFX.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(nLightAFX, 0);
    addToRunTimeSelectionTable(heatSourceModel, nLightAFX, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::nLightAFX::nLightAFX
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    const dictionary innerDict
    (
        heatSourceModelCoeffs_.optionalSubDict("inner")
    );

    const dictionary outerDict
    (
        heatSourceModelCoeffs_.optionalSubDict("outer")
    );

    alpha_ = heatSourceModelCoeffs_.lookup<scalar>("alpha");

    r0_ = innerDict.lookup<scalar>("radius");
    sigma0_ = innerDict.lookup<scalar>("sigma");
    A0_ = innerDict.lookup<scalar>("A");
    B0_ = innerDict.lookup<scalar>("B");

    r1_ = outerDict.lookup<scalar>("radius");
    sigma1_ = outerDict.lookup<scalar>("sigma");
    A1_ = outerDict.lookup<scalar>("A");
    B1_ = outerDict.lookup<scalar>("B");

    const scalar d = dimensions_.z();

    const scalar x =
        max
        (
            d/min(staticDimensions_.x(), staticDimensions_.y()),
            0.001
        );

    const scalar p0 = min(max(A0_*std::log2(x) + B0_, 0.0), 9.0);
    const scalar p1 = min(max(A1_*std::log2(x) + B1_, 0.0), 9.0);

    n0_ = Foam::pow(2.0, p0);
    n1_ = Foam::pow(2.0, p1);

    a0_ = ai(r0_, sigma0_, d, n0_);
    a1_ = ai(r1_, sigma1_, d, n1_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::nLightAFX::weight(const vector& x)
{
    const scalar z = Foam::mag(x.z());

    const scalar r = Foam::sqrt(x.x()*x.x() + x.y()*x.y());

    const scalar d = dimensions_.z();

    const scalar x0 =
        std::exp(-0.5*Foam::sqr((r - r0_)/sigma0_))
      + std::exp(-0.5*Foam::sqr((r + r0_)/sigma0_));

    const scalar x1 =
        std::exp(-0.5*Foam::sqr((r - r1_)/sigma1_))
      + std::exp(-0.5*Foam::sqr((r + r1_)/sigma1_));

    const scalar s0 =
        std::exp(-3.0*Foam::pow(mag(z/d), n0_));

    const scalar s1 =
        std::exp(-3.0*Foam::pow(mag(z/d), n1_));

    return ((1.0 - alpha_)*x0*s0*a1_) + (alpha_*x1*s1*a0_);
}


inline Foam::scalar
Foam::heatSourceModels::nLightAFX::ai
(
    scalar x,
    scalar s,
    scalar d,
    scalar n
)
{
    const scalar t1 =
        2.0*pi*s*d*Foam::tgamma(1.0/n)
      / (n*Foam::pow(3.0, 1.0/n));

    const scalar t2 =
        2.0*s*std::exp(-0.5*Foam::sqr(x/s))
      + Foam::sqrt(2.0*pi)*x*Foam::erf(x/(Foam::sqrt(2.0)*s));

    return t1*t2;
}


inline Foam::dimensionedScalar
Foam::heatSourceModels::nLightAFX::V0()
{
    const scalar d = dimensions_.z();

    const scalar x =
        max
        (
            d/min(staticDimensions_.x(), staticDimensions_.y()),
            0.001
        );

    const scalar p0 = min(max(A0_*std::log2(x) + B0_, 0.0), 9.0);
    const scalar p1 = min(max(A1_*std::log2(x) + B1_, 0.0), 9.0);

    n0_ = Foam::pow(2.0, p0);
    n1_ = Foam::pow(2.0, p1);

    a0_ = ai(r0_, sigma0_, d, n0_);
    a1_ = ai(r1_, sigma1_, d, n1_);

    return dimensionedScalar("V0", dimVolume, a0_*a1_);
}


bool Foam::heatSourceModels::nLightAFX::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        const dictionary innerDict
        (
            heatSourceModelCoeffs_.optionalSubDict("inner")
        );

        const dictionary outerDict
        (
            heatSourceModelCoeffs_.optionalSubDict("outer")
        );

        heatSourceModelCoeffs_.lookup("alpha") >> alpha_;

        innerDict.lookup("radius") >> r0_;
        innerDict.lookup("sigma") >> sigma0_;
        innerDict.lookup("A") >> A0_;
        innerDict.lookup("B") >> B0_;

        outerDict.lookup("radius") >> r1_;
        outerDict.lookup("sigma") >> sigma1_;
        outerDict.lookup("A") >> A1_;
        outerDict.lookup("B") >> B1_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
