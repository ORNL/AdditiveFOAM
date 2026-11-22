{
    const surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

    //- The solid fraction follows the metal advected into the cell: the
    //  solid fraction of the metal is kept, which on an isothermal step of
    //  the path is a state and not a function of the temperature
    alphaSolid.oldTime();

    forAll(alphaSolid, celli)
    {
        const scalar alpha1cOld =
            min(max(alpha1.oldTime()[celli], scalar(0)), scalar(1));

        const scalar fs =
            alpha1cOld > small
          ? min(alphaSolid[celli]/alpha1cOld, scalar(1))
          : thermo.fs(T[celli]);

        alphaSolid[celli] = alpha1c[celli]*fs;
    }
    alphaSolid.correctBoundaryConditions();

    #include "thermo/updateProperties.H"

    //- Energy equation in conservative enthalpy form, linearized about the
    //  state of the last property update: the implicit terms in T are the
    //  Newton corrections of the explicit enthalpy terms, and vanish at
    //  convergence
    phaseChange.begin(nThermoCorr, thermoTol);

    const labelList surfaceCells(sources.surfaceCells(alpha1));

    while (phaseChange.loop())
    {
        #include "thermo/surfaceCooling.H"

        const surfaceScalarField rhoCpPhi
        (
            "rhoCpPhi",
            rho1*alphaPhi1*linearInterpolate(Cp1)
          + rho2*alphaPhi2*linearInterpolate(Cp2)
        );

        // Upwind in both convection terms: their high-order parts cancel
        const fv::gaussConvectionScheme<scalar> upwindCp
        (
            mesh,
            rhoCpPhi,
            tmp<surfaceInterpolationScheme<scalar>>
            (
                new upwind<scalar>(mesh, rhoCpPhi)
            )
        );

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
          - fvc::ddt(rhoCp, T) + fvc::ddt(rhoH)
          + upwindCp.fvmDiv(rhoCpPhi, T)
          - upwindCp.fvcDiv(rhoCpPhi, T)
          + fvc::div(rhoHPhi)
          - fvm::laplacian(kappa, T)
         ==
            sources.qDot()
          + phaseChange.latentSource()
          - surfaceCoolingSu
          - fvm::Sp(surfaceCoolingSp, T)
        );

        TEqn.relax();
        TEqn.solve();

        phaseChange.project(sensible);

        #include "thermo/updateProperties.H"
    }

    if (evaporationCoeff > 0 || emissivity > 0)
    {
        Info<< "Surface cooling: evaporation " << evaporationPower
            << " W, radiation " << radiationPower
            << " W, area " << surfaceArea << " m2" << endl;
    }

    // Solid metal does not move, so the melted metal of a cell is its largest
    // liquid fraction, limited by the metal that remains
    alphaMelted = min(alpha1c, max(alphaMelted, alpha1c - alphaSolid));
    alphaMelted.correctBoundaryConditions();
}
