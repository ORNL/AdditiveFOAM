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
    k_(1),
    J0_(Zero),
    J1_(Zero),
    A0_(Zero),
    A1_(Zero)
{
    alpha_ = heatSourceModelCoeffs_.lookup<scalar>("alpha");

    r0_ = heatSourceModelCoeffs_.lookup<scalar>("r0");
    sigma0_ = heatSourceModelCoeffs_.lookup<scalar>("sigma0");

    r1_ = heatSourceModelCoeffs_.lookup<scalar>("r1");
    sigma1_ = heatSourceModelCoeffs_.lookup<scalar>("sigma1");

    nSlope_ = heatSourceModelCoeffs_.lookup<scalar>("nSlope");
    nIntercept_ = heatSourceModelCoeffs_.lookup<scalar>("nIntercept");

    J0_ = J(0, r0_, sigma0_);
    J1_ = J(0, r1_, sigma1_);

    updateBounds();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar Foam::heatSourceModels::nLightAFX::J
(
    const scalar R,
    const scalar r,
    const scalar sigma
) const
{
    return
        Foam::sqr(sigma)
       *(
            std::exp(-0.5 * Foam::sqr((R - r) / sigma))
          + std::exp(-0.5 * Foam::sqr((R + r) / sigma))
        )
      + r*sigma*Foam::sqrt(0.5*pi)
       *(
            Foam::erfc((R - r)/(Foam::sqrt(2.0)*sigma))
          - Foam::erfc((R + r)/(Foam::sqrt(2.0)*sigma))
        );
}


inline Foam::scalar Foam::heatSourceModels::nLightAFX::cutoff
(
    const scalar r,
    const scalar sigma,
    const scalar J0,
    const scalar epsilon
) const
{
    // J(r + sigma*u) <= C*exp(-u^2/2)
    const scalar C =
        2.0*Foam::sqr(sigma)
      + r*sigma*Foam::sqrt(0.5*pi);

    const scalar u2 = 2.0*max(std::log(C/(epsilon*J0)), 0.0);

    return r + sigma*Foam::sqrt(u2);
}


void Foam::heatSourceModels::nLightAFX::updateBounds()
{
    // (1 - epsilon)^2 = (1 - profileTol_)
    const scalar epsilon =
        profileTol_/(1.0 + Foam::sqrt(1.0 - profileTol_));

    scalar rMax =
        max
        (
            cutoff(r0_, sigma0_, J0_, epsilon),
            cutoff(r1_, sigma1_, J1_, epsilon)
        );

    scalar rMin = 0;

    // Bisect to radial tolerance
    const label nIter = ceil(std::log2(1/small));

    for (label i=0; i<nIter && rMax-rMin>1e-8; ++i)
    {
        const scalar R = 0.5*(rMin + rMax);

        const scalar epsilonR =
            (1.0 - alpha_)
           *J(R, r0_, sigma0_)/J0_
          + alpha_*J(R, r1_, sigma1_)/J1_;

        (epsilonR > epsilon ? rMin : rMax) = R;
    }

    sourceMin_.x() = -rMax;
    sourceMin_.y() = -rMax;
    sourceMax_.x() = rMax;
    sourceMax_.y() = rMax;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::nLightAFX::update()
{
    updateDimensions();

    const scalar d = dimensions_.z();

    const scalar a =
        d/min(staticDimensions_.x(), staticDimensions_.y());

    const scalar n =
        min
        (
            max
            (
                nSlope_*std::log2(a)
              + nIntercept_,
                0.0
            ),
            9.0
        );

    k_ = Foam::pow(2.0, n);

    const scalar Ap =
        d*Foam::tgamma(1.0/k_)
       /(k_*Foam::pow(3.0, 1.0/k_));

    A0_ = 2.0*pi*J0_*Ap;
    A1_ = 2.0*pi*J1_*Ap;

    // (1 - epsilon)^2 = (1 - profileTol_)
    const scalar epsilon =
        profileTol_/(1.0 + Foam::sqrt(1.0 - profileTol_));

    const scalar zMax =
        d
       *Foam::pow
        (
            invIncGammaRatio_P(1.0/k_, 1.0 - epsilon)/3.0,
            1.0/k_
        );

    sourceMin_.z() = -zMax;
    sourceMax_.z() = 0;

    V0_ = dimensionedScalar("V0", dimVolume, A0_);
}

inline Foam::scalar
Foam::heatSourceModels::nLightAFX::weight(const vector& x) const
{
    if (x.z() > 0)
    {
        return 0;
    }

    const scalar z = -x.z();

    const scalar r = Foam::sqrt(x.x()*x.x() + x.y()*x.y());

    const scalar d = dimensions_.z();

    const scalar I0 =
        std::exp(-0.5*Foam::sqr((r - r0_)/sigma0_))
      + std::exp(-0.5*Foam::sqr((r + r0_)/sigma0_));

    const scalar I1 =
        std::exp(-0.5*Foam::sqr((r - r1_)/sigma1_))
      + std::exp(-0.5*Foam::sqr((r + r1_)/sigma1_));

    const scalar p = std::exp(-3.0*Foam::pow(z/d, k_));

    return p*((1.0 - alpha_)*I0 + alpha_*I1*A0_/A1_);
}


bool Foam::heatSourceModels::nLightAFX::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_.lookup("alpha") >> alpha_;

        heatSourceModelCoeffs_.lookup("r0") >> r0_;
        heatSourceModelCoeffs_.lookup("sigma0") >> sigma0_;

        heatSourceModelCoeffs_.lookup("r1") >> r1_;
        heatSourceModelCoeffs_.lookup("sigma1") >> sigma1_;

        heatSourceModelCoeffs_.lookup("nSlope") >> nSlope_;
        heatSourceModelCoeffs_.lookup("nIntercept") >> nIntercept_;

        J0_ = J(0, r0_, sigma0_);
        J1_ = J(0, r1_, sigma1_);

        updateBounds();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
