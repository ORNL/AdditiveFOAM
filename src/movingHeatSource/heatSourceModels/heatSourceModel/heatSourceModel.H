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

Class
    Foam::heatSourceModel

Description
    Common interface for volumetric heat-source distributions. Profile metrics
    always describe the intensity footprint in the beam plane at z=0.

SourceFiles
    heatSourceModel.C
    heatSourceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatSourceModel_H
#define heatSourceModel_H

#include "autoPtr.H"
#include "profileMetrics.H"
#include "boundBox.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "labelVector.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class heatSourceModel
{
protected:

    labelVector nPoints_;

    scalar tolerance_;

    scalar depth_;

    scalar sourceDepth_;

    boundBox bounds_;

    dimensionedScalar V0_;

    heatSourceModel(const dictionary& dict);

public:

    TypeName("heatSourceModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        heatSourceModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );

    heatSourceModel(const heatSourceModel&) = delete;

    static autoPtr<heatSourceModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~heatSourceModel()
    {}

    virtual void update
    (
        const scalar depth,
        const scalar aspectRatio
    ) = 0;

    virtual scalar weight(const vector& r) const = 0;

    virtual const profileMetrics& metrics() const = 0;

    inline scalar depth() const
    {
        return depth_;
    }

    inline const labelVector& nPoints() const
    {
        return nPoints_;
    }

    inline const boundBox& bounds() const
    {
        return bounds_;
    }

    inline const dimensionedScalar& V0() const
    {
        return V0_;
    }

    void operator=(const heatSourceModel&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
