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

#ifndef absorptionModel_H
#define absorptionModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class absorptionModel
{
public:

    TypeName("absorptionModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        absorptionModel,
        dictionary,
        (const dictionary& dict, const fvMesh& mesh),
        (dict, mesh)
    );

    absorptionModel()
    {}

    absorptionModel(const absorptionModel&) = delete;

    static autoPtr<absorptionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~absorptionModel()
    {}

    //- Absorptivity of a volumetric source, from the melt pool aspect ratio
    virtual scalar eta(const scalar aspectRatio) const = 0;

    //- Absorptivity of a surface hit at the incidence cosine cosTheta; a
    //  model without an angle dependence returns eta(aspectRatio)
    virtual scalar eta(const scalar aspectRatio, const scalar cosTheta) const
    {
        return eta(aspectRatio);
    }

    //- The deposition the model is written for: "surface" (a flux on the
    //  first metal surface), "volumetric" (a source in the metal), or "any"
    virtual word scope() const
    {
        return "any";
    }

    void operator=(const absorptionModel&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
