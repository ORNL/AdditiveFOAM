{
    #include "thermo/thermoScheme.H"

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
    
    volScalarField T0("T0", T);
    
    for (int tCorr=0; tCorr < nThermoCorr; tCorr++)
    {
        #include "thermo/thermoSource.H"

        //- solve the energy equation
        solve
        (
            TEqn
         ==
            rDeltaT*rho*Lf*(fvm::Sp(dFdT, T) - dFdT*T0)
          + rho*Lf*fvc::ddt(alpha1)
        );

        T.correctBoundaryConditions();

        //- update solid fraction via Taylor's series expansion
        volScalarField alpha10("alpha10", alpha1);
        alpha1 = min(max(alpha10 + dFdT*(T - T0), 0.0), 1.0);
        alpha1.correctBoundaryConditions();

        //- check convergence of iterative procedure
        scalar residual =
            gMax(mag(alpha1.internalField() - alpha10.internalField())());

        Info<< "Thermo: iteration " << tCorr
            << " residual: " << residual << endl;

        if ((residual < thermoTol) && (tCorr > 0))
        {
            break;
        }
    }
}
