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

Class
    Foam::absorptionModels::Fresnel

Description
    Absorptivity of a metal surface from the Fresnel equations of an
    absorbing medium with complex refractive index n + ik, for unpolarized
    light at the local incidence angle. The constants are those of the metal
    at the wavelength of the beam: given as n and k in the absorption
    dictionary, or taken from the optical sub-dictionary of
    transportProperties when the absorption dictionary has none.

    A surface model: the beam is absorbed at the first metal surface it
    reaches, with the normal of the metal fraction there.

SourceFiles
    FresnelAbsorption.C

\*---------------------------------------------------------------------------*/

#ifndef FresnelAbsorption_H
#define FresnelAbsorption_H

#include "absorptionModel.H"

namespace Foam
{
namespace absorptionModels
{

class Fresnel
:
    public absorptionModel
{
    //- Real and imaginary parts of the refractive index
    scalar n_;
    scalar k_;

public:

    TypeName("Fresnel");

    Fresnel(const dictionary& dict, const fvMesh& mesh);

    virtual ~Fresnel()
    {}

    //- Normal-incidence absorptivity, for a volumetric call
    virtual scalar eta(const scalar aspectRatio) const
    {
        return eta(aspectRatio, 1.0);
    }

    //- Absorptivity at the incidence cosine cosTheta
    virtual scalar eta(const scalar aspectRatio, const scalar cosTheta) const;

    virtual word scope() const
    {
        return "surface";
    }
};

}
}

#endif
