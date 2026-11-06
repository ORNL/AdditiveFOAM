/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

Application
    smoothSetFields

Description
    Set volScalarField values from a setFieldsDict with volume fractions
    instead of the 0/1 cell selection of setFields.

    Reads the same dictionary as setFields (default system/setFieldsDict,
    -dict to select another) and supports the sphereToCell (centre, radius)
    and boxToCell (box) regions with volScalarFieldValue entries. Regions are
    applied in order and a later region overrides an earlier one, as in
    setFields. A cell whose vertices all lie inside a region takes its value;
    a cell cut by the region surface is sub-sampled nSamples^3 times over its
    bounding box and takes the mean of the sample values. The sub-sampling is
    exact for axis-aligned hexahedral cells (blockMesh).

    nSamples is read from the dictionary (default 8) or the -nSamples option.

Usage
    smoothSetFields [-dict <file>] [-nSamples <n>] [-parallel]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "systemDict.H"
#include "treeBoundBox.H"
#include "meshSearch.H"
#include "Map.H"
#include "Tuple2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A region shape: inside test and bounding box
class shape
{
    bool sphere_;
    point centre_;
    scalar radiusSqr_;
    treeBoundBox bb_;

public:

    shape(const word& type, const dictionary& dict)
    :
        sphere_(type == "sphereToCell"),
        centre_(Zero),
        radiusSqr_(0)
    {
        if (sphere_)
        {
            centre_ = dict.lookup<point>("centre");
            const scalar r = dict.lookup<scalar>("radius");
            radiusSqr_ = sqr(r);
            bb_ = treeBoundBox(centre_ - vector::one*r, centre_ + vector::one*r);
        }
        else if (type == "boxToCell")
        {
            bb_ = treeBoundBox(dict.lookup("box"));
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Region type " << type << " is not supported; "
                << "use sphereToCell or boxToCell" << exit(FatalIOError);
        }
    }

    const treeBoundBox& bb() const
    {
        return bb_;
    }

    bool inside(const point& p) const
    {
        return sphere_ ? magSqr(p - centre_) <= radiusSqr_ : bb_.contains(p);
    }
};


//- Read the volScalarFieldValue entries of a fieldValues list
List<Tuple2<word, scalar>> readFieldValues(const dictionary& dict)
{
    ITstream& is = dict.lookup("fieldValues");
    DynamicList<Tuple2<word, scalar>> values;

    is.readBegin("fieldValues");
    token t(is);
    while (!(t.isPunctuation() && t.pToken() == token::END_LIST))
    {
        is.putBack(t);
        const word type(is);
        const word name(is);
        const scalar value(readScalar(is));

        if (type != "volScalarFieldValue")
        {
            FatalIOErrorInFunction(dict)
                << "Field value type " << type << " is not supported; "
                << "use volScalarFieldValue" << exit(FatalIOError);
        }
        values.append(Tuple2<word, scalar>(name, value));
        is >> t;
    }

    return List<Tuple2<word, scalar>>(values);
}


int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "addRegionOption.H"
    argList::addOption
    (
        "nSamples",
        "label",
        "sub-samples per direction in a cut cell (default 8)"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createRegionMeshNoChangers.H"

    const dictionary setFieldsDict(systemDict("setFieldsDict", args, mesh));

    const label n =
        args.optionLookupOrDefault<label>
        (
            "nSamples",
            setFieldsDict.lookupOrDefault<label>("nSamples", 8)
        );
    const label nSub = n*n*n;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    // Fields set by any region
    wordHashSet fieldNames;
    forAll(regions, regioni)
    {
        const List<Tuple2<word, scalar>> values
        (
            readFieldValues(regions[regioni].dict())
        );
        forAll(values, i)
        {
            fieldNames.insert(values[i].first());
        }
    }

    const pointField& points = mesh.points();
    const labelListList& cellPoints = mesh.cellPoints();
    const indexedOctree<treeDataCell>& tree = meshSearch::New(mesh).cellTree();

    forAllConstIter(wordHashSet, fieldNames, iter)
    {
        const word& fieldName = iter.key();

        volScalarField field
        (
            IOobject
            (
                fieldName,
                runTime.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // Value of a cell wholly inside the last region that covered it,
        // and the sample values of cells cut by a region surface
        scalarField uniform(field.primitiveField());

        if (setFieldsDict.found("defaultFieldValues"))
        {
            ITstream& is = setFieldsDict.lookup("defaultFieldValues");
            is.readBegin("defaultFieldValues");
            token t(is);
            while (!(t.isPunctuation() && t.pToken() == token::END_LIST))
            {
                is.putBack(t);
                const word type(is);
                const word name(is);
                const scalar value(readScalar(is));
                if (name == fieldName)
                {
                    uniform = value;
                }
                is >> t;
            }
        }

        Map<scalarField> samples;

        label nFull = 0, nCut = 0;

        forAll(regions, regioni)
        {
            const entry& region = regions[regioni];
            const List<Tuple2<word, scalar>> values
            (
                readFieldValues(region.dict())
            );

            bool found = false;
            scalar value = 0;
            forAll(values, i)
            {
                if (values[i].first() == fieldName)
                {
                    found = true;
                    value = values[i].second();
                }
            }
            if (!found)
            {
                continue;
            }

            const shape s(region.keyword(), region.dict());

            const labelList cells(tree.findBox(s.bb()));

            forAll(cells, i)
            {
                const label celli = cells[i];
                const labelList& cp = cellPoints[celli];

                bool allInside = true;
                boundBox cellBb(points, cp, false);
                forAll(cp, pi)
                {
                    if (!s.inside(points[cp[pi]]))
                    {
                        allInside = false;
                        break;
                    }
                }

                if (allInside)
                {
                    uniform[celli] = value;
                    samples.erase(celli);
                    nFull++;
                    continue;
                }

                const vector span(cellBb.span());
                scalarField cellSamples;
                bool cut = false;

                label k = 0;
                for (label i0 = 0; i0 < n; i0++)
                {
                    for (label i1 = 0; i1 < n; i1++)
                    {
                        for (label i2 = 0; i2 < n; i2++)
                        {
                            const point p
                            (
                                cellBb.min()
                              + cmptMultiply
                                (
                                    span,
                                    vector(i0 + 0.5, i1 + 0.5, i2 + 0.5)/n
                                )
                            );

                            if (s.inside(p))
                            {
                                if (!cut)
                                {
                                    cut = true;
                                    if (!samples.found(celli))
                                    {
                                        samples.insert
                                        (
                                            celli,
                                            scalarField(nSub, uniform[celli])
                                        );
                                    }
                                }
                                samples[celli][k] = value;
                            }
                            k++;
                        }
                    }
                }

                if (cut)
                {
                    nCut++;
                }
            }
        }

        scalarField& f = field.primitiveFieldRef();
        f = uniform;
        forAllConstIter(Map<scalarField>, samples, sIter)
        {
            f[sIter.key()] = average(sIter());
        }
        field.correctBoundaryConditions();

        Info<< "Field " << fieldName << ": "
            << returnReduce(nFull, sumOp<label>()) << " whole-cell and "
            << returnReduce(nCut, sumOp<label>()) << " cut-cell region hits, "
            << returnReduce(samples.size(), sumOp<label>())
            << " cells with a fractional value, integral "
            << gSum(f*mesh.V().primitiveField()) << nl << endl;

        field.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
