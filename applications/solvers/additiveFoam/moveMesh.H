// Taken from: Foam::solvers::incompressibleFluid::moveMesh()
//if (pimple.firstPimpleIter() || pimple.moveMeshOuterCorrectors())
{
    // Move the mesh
    Info << "about to move mesh" << endl;
    mesh.move();
    Info << "just moved mesh" << endl;

    if (mesh.changing())
    {
        MRF.update();

        if (correctPhi || mesh.topoChanged())
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            correctUphiBCs(U, phi, true);

            fv::correctPhi
            (
                phi,
                U,
                p,
                autoPtr<volScalarField>(),
                autoPtr<volScalarField>(),
                pressureReference,
                pimple
            );

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        // Taken from: Foam::solvers::fluidSolver::meshCourantNo
        if (checkMeshCourantNo)
        {
            const scalarField sumPhi
            (
                fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
            );

            const scalar meshCoNum
            (
                0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue()
            );

            const scalar meanMeshCoNum
            (
                0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue()
            );

            Info<< "Mesh Courant Number mean: " << meanMeshCoNum
                << " max: " << meshCoNum << endl;
        }
    }
}
