{
    fvScalarMatrix TEqn(T, dimPower);

    if (explicitSolve)
    {
        TEqn =
        (
            rho*Cp*(fvm::ddt(T) + fvc::div(phi, T))
          - fvc::laplacian(kappa, T)
          - sources.qDot()
        );
    }
    else
    {
        TEqn =
        (
            rho*Cp*(fvm::ddt(T) + fvm::div(phi, T))
          - fvm::laplacian(kappa, T)
          - sources.qDot()
        );
    }
           
    //- set source term linearization weights for time discretization scheme  
    scalar coefft = 1;
    const dimensionedScalar deltaT  = runTime.deltaT();
    
    const word ddtScheme = mesh.schemes().ddt(T.name());
    if ((ddtScheme != "Euler") && (ddtScheme != "backward"))
    {
        FatalErrorInFunction
            << "Unsupported ddtScheme " << ddtScheme
            << ". Choose either Euler or backwards."
            << abort(FatalError);
    }
    else if (ddtScheme == "backward")
    {
        const dimensionedScalar deltaT0 = runTime.deltaT0();
        coefft += (deltaT/(deltaT + deltaT0)).value();    
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
    
    volScalarField T0("T0", T);
    
    for (int tCorr=0; tCorr < nThermoCorr; tCorr++)
    {
        #include "thermo/thermoSource.H"

        //- solve the energy equation
        solve
        (
            TEqn
         ==
            coefft*rho*Lf/deltaT*(fvm::Sp(dFdT, T) - dFdT*T0)
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
