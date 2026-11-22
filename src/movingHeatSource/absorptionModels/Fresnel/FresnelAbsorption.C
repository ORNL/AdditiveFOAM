/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#include "FresnelAbsorption.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include <complex>

namespace Foam
{
namespace absorptionModels
{
    defineTypeNameAndDebug(Fresnel, 0);
    addToRunTimeSelectionTable(absorptionModel, Fresnel, dictionary);
}
}


Foam::absorptionModels::Fresnel::Fresnel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    n_(0),
    k_(0)
{
    if (dict.found("n") && dict.found("k"))
    {
        n_ = dict.lookup<scalar>("n");
        k_ = dict.lookup<scalar>("k");
    }
    else
    {
        const IOdictionary& transportProperties =
            mesh.lookupObject<IOdictionary>("transportProperties");

        if (!transportProperties.found("optical"))
        {
            FatalIOErrorInFunction(dict)
                << "The Fresnel absorption model needs n and k, in the "
                << "absorption dictionary or in the optical sub-dictionary "
                << "of transportProperties" << exit(FatalIOError);
        }

        const dictionary& optical = transportProperties.subDict("optical");
        n_ = optical.lookup<scalar>("n");
        k_ = optical.lookup<scalar>("k");
    }

    Info<< "    Fresnel absorption: n = " << n_ << ", k = " << k_
        << ", normal-incidence absorptivity " << eta(0, 1.0) << endl;
}


Foam::scalar Foam::absorptionModels::Fresnel::eta
(
    const scalar,
    const scalar cosTheta
) const
{
    // Fresnel coefficients of an absorbing medium, with the complex
    // refraction angle from Snell's law; unpolarized light
    const scalar cosI = min(max(cosTheta, scalar(0)), scalar(1));
    const std::complex<scalar> n2(n_, k_);
    const std::complex<scalar> sin2I(1.0 - sqr(cosI), 0.0);
    const std::complex<scalar> cosT = std::sqrt(1.0 - sin2I/(n2*n2));

    const std::complex<scalar> rs = (cosI - n2*cosT)/(cosI + n2*cosT);
    const std::complex<scalar> rp = (n2*cosI - cosT)/(n2*cosI + cosT);

    return 1.0 - 0.5*(std::norm(rs) + std::norm(rp));
}


// ************************************************************************* //
