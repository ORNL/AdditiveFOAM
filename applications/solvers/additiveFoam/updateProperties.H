// store old time values
kappa.oldTime();

forAll(mesh.cells(), cellI)
{
    scalar alpha1_ = alpha1[cellI];
    scalar alpha2_ = 1.0 - alpha1[cellI];

    //- Discrete powder-liquid transition model
    alpha3[cellI] = (alpha2_ > 0.0) ? 0.0 : alpha3[cellI];
    scalar alpha3_ = alpha3[cellI];

    scalar T1 = min(max(T[cellI], 300.0), Tsol.value());
    scalar T2 = min(max(T[cellI], Tliq.value()), 2.0*Tliq.value());

    //- Update solid-liquid mixture properties
    scalar Cp12 = alpha1_*Cp1.value(T1) + alpha2_*Cp2.value(T2);
    scalar kappa12 = alpha1_*kappa1.value(T1) + alpha2_*kappa2.value(T2);

    //- Update solid-liquid-powder mixture properties
    kappa[cellI] = (1.0 - alpha3_)*kappa12 + alpha3_*kappa3.value(T1);
    Cp[cellI] = (1.0 - alpha3_)*Cp12 + alpha3_*Cp3.value(T1);
}

alpha3.correctBoundaryConditions();
Cp.correctBoundaryConditions();
kappa.correctBoundaryConditions();

// smooth discontinuities across continuum powder interface
kappa = 0.5*fvc::average(kappa + kappa.oldTime());
