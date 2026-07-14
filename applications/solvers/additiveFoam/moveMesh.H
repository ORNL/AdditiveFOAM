if (pimple.firstPimpleIter() || pimple.moveMeshOuterCorrectors())
{
    // Move the mesh
    mesh.move();

    if (mesh.changing())
    {
        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;

        if (correctPhi || mesh.poly().topoChanged())
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            correctUphiBCs(U, phi, true);

            fv::correctPhi
            (
                phi,
                U,
                p_rgh,
                autoPtr<volScalarField>(),
                autoPtr<volScalarField>(),
                pressureReference,
                pimple
            );

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        if (checkMeshCourantNo)
        {
            const scalarField sumPhi
            (
                fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
            );

            const scalar meshCoNum
            (
                0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue()
            );

            const scalar meanMeshCoNum
            (
                0.5
               *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
               *runTime.deltaTValue()
            );

            Info<< "Mesh Courant Number mean: " << meanMeshCoNum
                << " max: " << meshCoNum << endl;
        }
    }
}
