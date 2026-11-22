{
    //- Energy equation in conservative enthalpy form, linearized about the
    //  state of the last projection: the implicit terms in T are the Newton
    //  corrections of the explicit enthalpy terms, and vanish at convergence.
    //  The conductivity of the step is that of its old state. The explicit
    //  solution is forward Euler with the fluxes of the step at the old time;
    //  otherwise the ddt scheme of T applies to the sensible and the latent
    //  heat alike
    const dimensionedScalar rDeltaT = 1.0/runTime.deltaT();

    tmp<volScalarField> tFluxExplicit;

    if (explicitSolve)
    {
        tFluxExplicit =
            fvc::div(fvc::flux(phi, rhoH, "div(phi,T)"))
          - fvc::laplacian(kappa, T)
          - sources.qDot();
    }

    //- Terms of the step that do not depend on the enthalpy state, assembled
    //  once: convection of the sensible heat with the heat capacity of the
    //  old state, conduction with its conductivity, and the heat source. The
    //  convective flux of the enthalpy is corrected to the state of each
    //  iteration below
    tmp<surfaceScalarField> tRhoCpPhi;
    tmp<fvScalarMatrix> tTEqnStep;

    if (!explicitSolve)
    {
        tRhoCpPhi = phi*linearInterpolate(rhoCp);

        tTEqnStep =
            fvm::div(tRhoCpPhi(), T, "div(phi,T)")
          - fvm::laplacian(kappa, T)
          - sources.qDot();
    }

    //- the implicit temperature limit is assembled only when one is set
    const bool limitT = Tmax.value() < vGreat;

    phaseChange.begin(nThermoCorr, thermoTol);

    while (phaseChange.loop())
    {
        fvScalarMatrix TEqn(T, dimPower);

        if (explicitSolve)
        {
            fv::EulerDdtScheme<scalar> euler(mesh);

            TEqn =
            (
                euler.fvmDdt(rhoCp, T) - euler.fvcDdt(rhoCp, T)
              + euler.fvcDdt(rhoH)
              + tFluxExplicit()
             ==
                phaseChange.latentSource(true)
            );
        }
        else
        {
            TEqn =
            (
                fvm::ddt(rhoCp, T) - fvc::ddt(rhoCp, T) + fvc::ddt(rhoH)
              + tTEqnStep()
              - fvc::div(tRhoCpPhi(), T, "div(phi,T)")
              + fvc::div(fvc::flux(phi, rhoH, "div(phi,T)"))
             ==
                phaseChange.latentSource()
            );
        }

        //- optional implicit limiting of temperature field
        if (limitT)
        {
            const volScalarField A
            (
                "A",
                dimensionedScalar(rhoCp.dimensions(), 1e15)*pos(T - Tmax)
            );

            TEqn += rDeltaT*(fvm::Sp(A, T) - A*Tmax);
        }

        TEqn.solve();

        phaseChange.project(sensible);

        #include "thermo/updateEnthalpy.H"
    }

    #include "thermo/updateProperties.H"
}
