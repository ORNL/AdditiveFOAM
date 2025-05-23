//- Construct the temperature equation for the chosen time integration scheme
word ddtScheme = mesh.schemes().ddt(T.name());

if
(
    (ddtScheme != "Euler")
 && (ddtScheme != "backward")
 && (ddtScheme != "CrankNicolson")
)
{
    FatalErrorInFunction
        << "Unsupported ddtScheme " << ddtScheme
        << ". Choose Euler, CrankNicolson, or backward."
        << abort(FatalError);
}

fvScalarMatrix TEqn(T, dimPower);
dimensionedScalar rDeltaT = 1.0 / runTime.deltaT();
   
if (explicitSolve)
{
    // Forward Euler is the only explicit method supported.
    ddtScheme = "Euler";
    
    TEqn =
    (
        rho*Cp*(fv::EulerDdtScheme<scalar>(mesh).fvmDdt(T) + fvc::div(phi, T))
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

    scalar coefft = 1.0;
    
    if (ddtScheme == "backward")
    {
        const dimensionedScalar deltaT = runTime.deltaT();
        const dimensionedScalar deltaT0 = runTime.deltaT0();
        
        coefft = 1.0 + (deltaT/(deltaT + deltaT0)).value();    
    }
    
    if (ddtScheme == "CrankNicolson")
    {
        scalar ocCoeff =
            refCast<const fv::CrankNicolsonDdtScheme<scalar>>
            (
                fv::ddtScheme<scalar>::New
                (
                    mesh,
                    mesh.schemes().ddt(ddtScheme)
                )()
            ).ocCoeff();

        coefft = 1.0 + ocCoeff;
    }
    
    Info << "thero time coefficient (coefft): " << coefft << endl;
    
    rDeltaT *= coefft;
}
