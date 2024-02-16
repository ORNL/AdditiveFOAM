{
    fvScalarMatrix TEqn(T, dimPower);

    if (explicitSolve)
    {
        TEqn =
        (
            rho*Cp*(fvm::ddt(T) + fvc::div(phi, T))
          - fvc::laplacian(fvc::interpolate(kappa), T)
          - sources.qDot()
        );
    }
    else
    {
        TEqn =
        (
            rho*Cp*(fvm::ddt(T) + fvm::div(phi, T))
          - fvm::laplacian(fvc::interpolate(kappa), T)
          - sources.qDot()
        );
    }

    //- solve thermodynamic coupling using root-finding algorithm
    volScalarField dFdT
    (
        IOobject
        (
            "dFdT",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dFdT", dimless/dimTemperature, 0.0)
    );

    dimensionedScalar dt(runTime.deltaT());

    volScalarField T0("T0", T);
    volScalarField alpha10("alpha10", alpha1);

    for (int tCorr=0; tCorr < nThermoCorr; tCorr++)
    {
        #include "thermo/thermoSource.H"

        //- solve the energy equation
        solve
        (
            TEqn
         ==
            fvm::Sp(rho*Lf*dFdT/dt, T)
          + rho*Lf*(alpha10 - alpha1.oldTime() - dFdT*T0)/dt
        );

        T.correctBoundaryConditions();

        //- update solid fraction via Taylor's series expansion
        alpha1 = min(max(alpha10 + dFdT*(T - T0), 0.0), 1.0);
        alpha1.correctBoundaryConditions();

        //- check convergence of iterative procedure
        scalar residual =
            gMax(mag(alpha1.internalField() - alpha10.internalField())());

        Info<< "Thermo: iteration " << tCorr
            << " residual: " << residual << endl;

        if (residual < thermoTol)
        {
            break;
        }
    }

    // set solid fraction consistent with temperature
    forAll(mesh.cells(), cellI)
    {
        scalar alpha1_ = interpolateXY(T[cellI], thermo.x(), thermo.y());
        alpha1[cellI] = min(max(alpha1_, 0.0), 1.0);
    }

    alpha1.correctBoundaryConditions();
}
