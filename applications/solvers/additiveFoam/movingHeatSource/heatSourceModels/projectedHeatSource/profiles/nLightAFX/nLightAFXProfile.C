/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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

#include "nLightAFXProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace heatSourceProfiles
{
    defineTypeNameAndDebug(nLightAFXProfile, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProfile,
        nLightAFXProfile,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;


Foam::heatSourceProfiles::nLightAFXProfile::nLightAFXProfile
(
    const dictionary& dict,
    const fvMesh&,
    const scalar epsilon
)
:
    alpha_(0),
    r0_(0),
    sigma0_(0),
    r1_(0),
    sigma1_(0),
    J0_(0),
    J1_(0),
    integral_(0),
    D4Sigma_(0),
    profileBb_(point::zero, point::zero)
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");

    alpha_ = coeffs.lookup<scalar>("alpha");

    r0_ = coeffs.lookup<scalar>("r0");
    sigma0_ = coeffs.lookup<scalar>("sigma0");

    r1_ = coeffs.lookup<scalar>("r1");
    sigma1_ = coeffs.lookup<scalar>("sigma1");

    J0_ = J(0, r0_, sigma0_);
    J1_ = J(0, r1_, sigma1_);

    integral_ = 2.0*pi*J0_;

    const scalar r2 =
        (1.0 - alpha_)*K(r0_, sigma0_)/J0_
      + alpha_*K(r1_, sigma1_)/J1_;

    D4Sigma_ = 4.0*Foam::sqrt(0.5*r2);

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
            (1.0 - alpha_)*J(R, r0_, sigma0_)/J0_
          + alpha_*J(R, r1_, sigma1_)/J1_;

        (epsilonR > epsilon ? rMin : rMax) = R;
    }

    profileBb_ =
        boundBox
        (
            point(-rMax, -rMax, 0),
            point(rMax, rMax, 0)
        );
}


inline Foam::scalar Foam::heatSourceProfiles::nLightAFXProfile::J
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


inline Foam::scalar Foam::heatSourceProfiles::nLightAFXProfile::K
(
    const scalar r,
    const scalar sigma
) const
{
    const scalar f = std::exp(-0.5*Foam::sqr(r/sigma));
    const scalar e = Foam::erf(r/(Foam::sqrt(2.0)*sigma));

    return
        2.0*Foam::sqr(sigma)
       *(Foam::sqr(r) + 2.0*Foam::sqr(sigma))*f
      + r*sigma*Foam::sqrt(2.0*pi)
       *(Foam::sqr(r) + 3.0*Foam::sqr(sigma))*e;
}


inline Foam::scalar Foam::heatSourceProfiles::nLightAFXProfile::cutoff
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


inline Foam::scalar Foam::heatSourceProfiles::nLightAFXProfile::weight
(
    const scalar x,
    const scalar y
) const
{
    const scalar r = Foam::sqrt(x*x + y*y);

    const scalar I0 =
        std::exp(-0.5*Foam::sqr((r - r0_)/sigma0_))
      + std::exp(-0.5*Foam::sqr((r + r0_)/sigma0_));

    const scalar I1 =
        std::exp(-0.5*Foam::sqr((r - r1_)/sigma1_))
      + std::exp(-0.5*Foam::sqr((r + r1_)/sigma1_));

    return (1.0 - alpha_)*I0 + alpha_*I1*J0_/J1_;
}

// ************************************************************************* //
