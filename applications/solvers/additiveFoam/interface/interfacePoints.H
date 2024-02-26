vector interfacePoints(const volScalarField& field, const scalar& iso)
{
    const fvMesh& mesh = field.mesh();

    // set local reference to mesh data
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& cc = mesh.C();

    DynamicList<vector> positions(mesh.nFaces());

    // check internal faces
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
            positions.append(p);
        }
    }

    // check boundary faces
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
                    positions.append(p);
                }
            }
        }
        else
        {
            // physical boundary : take face centre in liquid
            const vectorField fc(mesh.boundaryMesh()[patchi].faceCentres());
            const scalarField pif(fieldPf.patchInternalField());

            forAll(fc, i)
            {
                if (pif[i] >= iso)
                {
                    positions.append(fc[i]);
                }
            }
        }
    }

    positions.shrink();

    boundBox isoBb(positions);

    return returnReduce(isoBb.span(), maxOp<vector>());
}
