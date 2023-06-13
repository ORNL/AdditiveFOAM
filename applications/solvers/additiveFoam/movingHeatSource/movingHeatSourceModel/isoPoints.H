/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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

Function
    Foam::isoPoints

Description
    Function for fast determination of isocontour points given a scalar field
    and value.

\*---------------------------------------------------------------------------*/

#ifndef isoPoints_H
#define isoPoints_H

List<point> isoPoints(const volScalarField& field, const scalar& iso)
{
    // set local reference to mesh data
    const fvMesh& mesh = field.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& cc = mesh.C();

    DynamicList<point> localPoints;

    // isocontour across internal faces
    for(label facei=0; facei < mesh.nInternalFaces(); facei++)
    {        
        const label own = owner[facei];
        const label nei = neighbour[facei];

        scalar minFace = min(field[own], field[nei]);
        scalar maxFace = max(field[own], field[nei]);

        if ((minFace < iso) && (maxFace >= iso))
        {
            vector d = cc[nei] - cc[own];
            vector p = cc[own] + d*(iso - field[own])/(field[nei] - field[own]);
            
            localPoints.append(p);
        }
    }

    // isocontour across processor and boundary faces
    const volScalarField::Boundary& fieldBf = field.boundaryField();

    forAll(fieldBf, patchi)
    {   
        const fvPatchScalarField& fieldPf = fieldBf[patchi];

        const labelUList& faceCells = fieldPf.patch().faceCells();

        if (fieldPf.coupled())
        {
            // processor boundary : interpolate across face
            const vectorField ccn(cc.boundaryField()[patchi].patchNeighbourField());
            const scalarField fn(fieldPf.patchNeighbourField());

            forAll(faceCells, facei)
            {
                label own = faceCells[facei];

                scalar minFace = min(field[own], fn[facei]);
                scalar maxFace = max(field[own], fn[facei]);

                if ((minFace < iso) && (maxFace >= iso))
                {
                    vector d = ccn[facei] -  cc[own];
                    vector p = cc[own] + d*(iso - field[own])/(fn[facei] - field[own]);
                    
                    localPoints.append(p);
                }
            }
        }
    }

    localPoints.shrink();

    return localPoints;
}

#endif
