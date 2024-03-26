// write out interface dimensions
Info<< "interface span " << "(" << Tliq.value() << "): "
    << interfacePoints(T, Tliq.value()) << endl;

Info<< "interface span " << "(" << Tsol.value() << "): "
    << interfacePoints(T, Tsol.value()) << endl;

// write out the maximum temperature
Info<< "maximum temperature: " << gMax(T.internalField()) << endl;

volScalarField G("G", mag(fvc::grad(T)));

forAll(mesh.cells(), celli)
{
    if ( (T[celli] <= Tliq.value()) && (T.oldTime()[celli] > Tliq.value()) )
    {
        const point pt = mesh.C()[celli];

        interfaceData.append(
            {
                pt[0],
                pt[1],
                pt[2],
                R[celli] / G[celli],
                G[celli]
            });
    }
}
