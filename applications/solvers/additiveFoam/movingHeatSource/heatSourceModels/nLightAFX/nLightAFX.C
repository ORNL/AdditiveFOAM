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
    alpha_ = heatSourceModelCoeffs_.lookup<scalar>("alpha");

    r0_ = heatSourceModelCoeffs_.lookup<scalar>("r0");
    sigma0_ = heatSourceModelCoeffs_.lookup<scalar>("sigma0");

    r1_ = heatSourceModelCoeffs_.lookup<scalar>("r1");
    sigma1_ = heatSourceModelCoeffs_.lookup<scalar>("sigma1");

    nSlope_ = heatSourceModelCoeffs_.lookup<scalar>("nSlope");
    nIntercept_ = heatSourceModelCoeffs_.lookup<scalar>("nIntercept");

    const scalar d = dimensions_.z();

    const scalar x =
        max
        (
            d/min(staticDimensions_.x(), staticDimensions_.y()),
            0.001
        );

    const scalar n =
        min
        (
            max
            (
                nSlope_*std::log2(x)
              + nIntercept_,
                0.0
            ),
            9.0
        );

    k_ = Foam::pow(2.0, n);

    a0_ = ai(r0_, sigma0_, d, k_);
    a1_ = ai(r1_, sigma1_, d, k_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::nLightAFX::weight(const vector& x)
{
    if (x.z() > 0)
    {
        return 0;
    }

    const scalar z = -x.z();

    const scalar r = Foam::sqrt(x.x()*x.x() + x.y()*x.y());

    const scalar d = dimensions_.z();

    const scalar x0 =
        std::exp(-0.5*Foam::sqr((r - r0_)/sigma0_))
      + std::exp(-0.5*Foam::sqr((r + r0_)/sigma0_));

    const scalar x1 =
        std::exp(-0.5*Foam::sqr((r - r1_)/sigma1_))
      + std::exp(-0.5*Foam::sqr((r + r1_)/sigma1_));

    const scalar s =
        std::exp(-3.0*Foam::pow(z/d, k_));

    return s*((1.0 - alpha_)*x0 + alpha_*x1*a0_/a1_);
}


inline Foam::scalar
Foam::heatSourceModels::nLightAFX::ai
(
    scalar x,
    scalar s,
    scalar d,
    scalar k
)
{
    const scalar t1 =
        2.0*pi*s*d*Foam::tgamma(1.0/k)
      / (k*Foam::pow(3.0, 1.0/k));

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

    const scalar n =
        min
        (
            max
            (
                nSlope_*std::log2(x)
              + nIntercept_,
                0.0
            ),
            9.0
        );

    k_ = Foam::pow(2.0, n);

    a0_ = ai(r0_, sigma0_, d, k_);
    a1_ = ai(r1_, sigma1_, d, k_);

    sourceMin_.z() = -1.5*d;
    sourceMax_.z() = 0;

    return dimensionedScalar("V0", dimVolume, a0_);
}


bool Foam::heatSourceModels::nLightAFX::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        heatSourceModelCoeffs_.lookup("alpha") >> alpha_;

        heatSourceModelCoeffs_.lookup("r0") >> r0_;
        heatSourceModelCoeffs_.lookup("sigma0") >> sigma0_;

        heatSourceModelCoeffs_.lookup("r1") >> r1_;
        heatSourceModelCoeffs_.lookup("sigma1") >> sigma1_;

        heatSourceModelCoeffs_.lookup("nSlope") >> nSlope_;
        heatSourceModelCoeffs_.lookup("nIntercept") >> nIntercept_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
