"""Case generator and reader for the thermocapillary slot tests.

freeSurface: additiveVoFFoam, liquid layer under a gas layer in a closed slot,
             the Sen and Davis surface deflection.
flatSurface: additiveFoam, the liquid alone under the Marangoni slip
             condition, the flat-surface return flow.

Default parameters are those of the Sen and Davis benchmark (Wirth et al.
2022, table 2): L = 1 mm, liquid depth h = 0.2 mm, gas height 0.3 mm,
tau = |dsigma/dT| |dT/dx| = 0.4 Pa, sigma = 0.01 N/m at the mean temperature,
zero gravity. The thermal diffusivity is raised so that the Marangoni number
is 2e-3 and the surface temperature stays linear, as the reference assumes
(a 100 times smaller diffusivity changes the deflection by under 1%).

The default free-surface mesh is the 100 x 240 mesh of the original
marangoni-test tutorial, with the default interfaceCompression (MULES)
transport, run to 0.1 s.
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

P = dict(L=1e-3, h=2e-4, H=3e-4, rho1=1562.5, mu1=5e-3, rho2=0.05,
         mu2=2.83e-5, sigma=0.01, dSigmadT=-1e-4, Thot=302.0, Tcold=298.0,
         alpha_thermal=1.6e-3, endTime=0.1, writeInterval=0.01)


def tau(p=P):
    return p['dSigmadT']*(p['Tcold'] - p['Thot'])/p['L']


def sigma_of_x(x, p=P):
    T = p['Thot'] + (p['Tcold'] - p['Thot'])*x/p['L']
    return p['sigma'] + p['dSigmadT']*(T - 0.5*(p['Thot'] + p['Tcold']))


HEADER = '''FoamFile
{{
    format      ascii;
    class       {cls};
    object      {obj};
}}
'''


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def field(case, name, cls, dims, internal, bcs):
    body = '\n'.join(f'    {k}\n    {{\n' + ''.join(f'        {e};\n' for e in v) + '    }'
                     for k, v in bcs.items())
    write(case/'0'/name, HEADER.format(cls=cls, obj=name)
          + f'dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n{body}\n}}\n')


def mesh(case, L, Ly, nx, ny, dz):
    write(case/'system/blockMeshDict', HEADER.format(cls='dictionary', obj='blockMeshDict') + f'''
convertToMeters 1;
vertices
(
    (0 0 0) ({L} 0 0) ({L} {Ly} 0) (0 {Ly} 0)
    (0 0 {dz}) ({L} 0 {dz}) ({L} {Ly} {dz}) (0 {Ly} {dz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1) );
edges ();
boundary
(
    bottom {{ type wall; faces ((0 1 5 4)); }}
    top    {{ type wall; faces ((3 7 6 2)); }}
    left   {{ type wall; faces ((0 4 7 3)); }}
    right  {{ type wall; faces ((1 2 6 5)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
''')


def temperature(case, nx, ny, L, p):
    x = (np.arange(nx) + 0.5)*L/nx
    T = np.tile(p['Thot'] + (p['Tcold'] - p['Thot'])*x/L, ny)
    values = 'nonuniform List<scalar>\n%d\n(\n%s\n)' % (T.size, '\n'.join(f'{t:.9f}' for t in T))
    field(case, 'T', 'volScalarField', '[0 0 0 1 0 0 0]', values,
          {'left': ['type fixedValue', f'value uniform {p["Thot"]}'],
           'right': ['type fixedValue', f'value uniform {p["Tcold"]}'],
           '"(top|bottom)"': ['type zeroGradient'], 'frontAndBack': ['type empty']})


def controls(case, app, p, extra=''):
    write(case/'system/controlDict', HEADER.format(cls='dictionary', obj='controlDict') + f'''
application     {app};
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {p["endTime"]};
deltaT          1e-7;
writeControl    adjustableRunTime;
writeInterval   {p["writeInterval"]};
purgeWrite      0;
writeFormat     binary;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           0.5;
maxAlphaCo      1;
maxDeltaT       1e-4;
{extra}
''')
    write(case/'system/decomposeParDict', HEADER.format(cls='dictionary', obj='decomposeParDict')
          + 'numberOfSubdomains 1;\nmethod simple;\nsimpleCoeffs { n (1 1 1); }\n')


def zero_beam(case):
    write(case/'constant/heatSourceDict', HEADER.format(cls='dictionary', obj='heatSourceDict') + '''
sources
{
    beam
    {
        path            scanPath;
        widthReference  D4Sigma;
        absorption { model constant; eta 0; }
        heatSource
        {
            model       superGaussian;
            radius      (1e-4 1e-4);
            depth       1e-5;
            k           2;
            nPoints     (10 10 10);
        }
    }
}
refinement { model none; }
''')
    write(case/'constant/scanPath', 'Mode X Y Z Power Param\n1 0 0 0 0 0\n1 0 0 0 0 1\n')
    write(case/'constant/g', HEADER.format(cls='uniformDimensionedVectorField', obj='g')
          + 'dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n')


def material(p, k1):
    liquidus = 250.0
    sigma0 = p['sigma'] - p['dSigmadT']*(0.5*(p['Thot'] + p['Tcold']) - liquidus)
    cp = 500.0
    return f'''
solid  {{ kappa ({k1} 0 0); Cp ({cp} 0 0); }}
liquid {{ kappa ({k1} 0 0); Cp ({cp} 0 0); }}
powder {{ kappa ({k1} 0 0); Cp ({cp} 0 0); }}
rho     [1 -3 0 0 0 0 0]    {p["rho1"]};
mu      [1 -1 -1 0 0 0 0]   {p["mu1"]};
beta    [0 0 0 -1 0 0 0]    0;
Lf      [0 2 -2 0 0 0 0]    1e3;
emissivity  0;
dSigmadT    {p["dSigmadT"]};
sigma       {sigma0};
thermoPath
(
    (200 1)
    ({liquidus} 0)
);
'''


def vof_case(case, nx=100, p=P, slip_walls=False, k_liquid=None, k_gas=None, schemes=None, ny=240,
             maxCapillaryCo=10):
    """additiveVoFFoam slot, interface on a cell face. Square cells by default;
    ny gives anisotropic cells (the original tutorial uses 100 x 240)."""
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    L, h, H = p['L'], p['h'], p['H']
    dx = L/nx
    ny = ny or int(round((h + H)/dx))
    dy = (h + H)/ny
    assert abs(round(h/dy)*dy - h) < 1e-6*h, 'the interface must lie on a cell face'
    mesh(case, L, h + H, nx, ny, dx)
    cp1, cp2 = 500.0, 520.0
    k1 = k_liquid if k_liquid is not None else p['alpha_thermal']*p['rho1']*cp1
    k2 = k_gas if k_gas is not None else p['alpha_thermal']*p['rho2']*cp2
    write(case/'constant/transportProperties', HEADER.format(cls='dictionary', obj='transportProperties')
          + material(p, k1) + f'''
gas
{{
    rho     {p["rho2"]};
    mu      {p["mu2"]};
    Cp      ({cp2} 0 0);
    kappa   ({k2} 0 0);
}}
''')
    zero_beam(case)
    controls(case, 'additiveVoFFoam', p, f'maxCapillaryCo  {maxCapillaryCo};')
    src = Path(__file__).resolve().parents[2]/'tutorials/additiveVoFFoam/tabulated/system'
    shutil.copy(src/'fvSchemes', case/'system/fvSchemes')
    shutil.copy(src/'fvSolution', case/'system/fvSolution')
    s = (case/'system/fvSolution').read_text()
    s = re.sub(r'thermoTolerance\s+[^;]+;', 'thermoTolerance     1e-8;', s)
    s = re.sub(r'nOuterCorrectors\s+\d+;', 'nOuterCorrectors    1;', s)
    s = re.sub(r'nCorrectors\s+\d+;', 'nCorrectors         3;', s)
    (case/'system/fvSolution').write_text(s)
    if schemes:
        t = (case/'system/fvSchemes').read_text()
        for a, b in schemes.items():
            t = t.replace(a, b)
        (case/'system/fvSchemes').write_text(t)
    y = (np.arange(ny) + 0.5)*dy
    alpha = np.repeat((y < h).astype(float), nx)
    field(case, 'alpha.metal', 'volScalarField', '[0 0 0 0 0 0 0]',
          'nonuniform List<scalar>\n%d\n(\n%s\n)' % (alpha.size, '\n'.join(f'{a:g}' for a in alpha)),
          {'"(left|right|top|bottom)"': ['type zeroGradient'], 'frontAndBack': ['type empty']})
    temperature(case, nx, ny, L, p)
    side = ['type slip'] if slip_walls else ['type noSlip']
    field(case, 'U', 'volVectorField', '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)',
          {'"(top|bottom)"': ['type noSlip'], '"(left|right)"': side, 'frontAndBack': ['type empty']})
    field(case, 'p_rgh', 'volScalarField', '[1 -1 -2 0 0 0 0]', 'uniform 0',
          {'"(left|right|top|bottom)"': ['type fixedFluxPressure', 'value uniform 0'],
           'frontAndBack': ['type empty']})
    return dict(case=str(case), solver='additiveVoFFoam', nx=nx, ny=ny, dx=dx, dy=dy, slip_walls=slip_walls,
                k_liquid=k1, k_gas=k2, **p)


def run(case, solver, env=None, log='log.solver'):
    """blockMesh and the solver, serial, in the environment given as a
    shell snippet (default: the current one)."""
    pre = (env + '; ') if env else ''
    cmd = f'{pre}blockMesh > log.blockMesh 2>&1 && {solver} > {log} 2>&1'
    return subprocess.run(['bash', '-c', cmd], cwd=case).returncode


def read(path, comps=1):
    t = Path(path).read_bytes()
    m = re.search(rb'internalField\s+nonuniform\s+List<(?:scalar|vector)>\s*\n?(\d+)\s*\n?\(', t)
    if not m:
        u = re.search(rb'internalField\s+uniform\s+\(?([^;)]+)\)?;', t).group(1).decode().split()
        return np.array([float(v) for v in u])
    n = int(m.group(1))
    s = m.end()
    if re.search(rb'format\s+binary', t):
        a = np.frombuffer(t[s:s + 8*comps*n], dtype='<f8').copy()
    else:
        a = np.array(t[s:t.index(b'\n)', s)].decode().replace('(', ' ').replace(')', ' ').split(), float)
    return a.reshape(n, comps) if comps > 1 else a


def times(case):
    return sorted((d for d in Path(case).iterdir()
                   if re.fullmatch(r'[0-9.e+-]+', d.name) and float(d.name) > 0),
                  key=lambda d: float(d.name))


def vof_profile(case, meta, tdir):
    """Surface height minus the initial depth per column, from the metal
    volume in the column: eta = sum(alpha) dy - h."""
    nx, ny, dx = meta['nx'], meta['ny'], meta['dx']
    dy = meta.get('dy', dx)
    alpha = read(tdir/'alpha.metal').reshape(ny, nx)
    U = read(tdir/'U', 3)
    return dict(t=float(tdir.name), x=(np.arange(nx) + 0.5)*dx,
                eta=alpha.sum(axis=0)*dy - meta['h'],
                umax=float(np.abs(U).max()) if U.ndim == 2 else 0.0,
                U=U.reshape(ny, nx, 3) if U.ndim == 2 else None,
                alpha=alpha, T=read(tdir/'T').reshape(ny, nx) if (tdir/'T').exists() else None)


def af_case(case, nx=100, p=dict(P, alpha_thermal=1.6e-5), endTime=0.1):
    """additiveFoam slot: the liquid alone under the Marangoni slip condition."""
    case = Path(case)
    if case.exists():
        shutil.rmtree(case)
    L, h = p['L'], p['h']
    dx = L/nx
    ny = int(round(h/dx))
    mesh(case, L, h, nx, ny, dx)
    k1 = p['alpha_thermal']*p['rho1']*500.0
    write(case/'constant/transportProperties', HEADER.format(cls='dictionary', obj='transportProperties')
          + material(p, k1))
    zero_beam(case)
    q = dict(p, endTime=endTime, writeInterval=endTime/10)
    controls(case, 'additiveFoam', q, 'maxDi           10;\nmaxPhaseChangeCo 1;')
    write(case/'system/fvSchemes', HEADER.format(cls='dictionary', obj='fvSchemes') + '''
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss upwind; div(phi,U) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; laplacian(kappa,T) Gauss harmonic corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''')
    write(case/'system/fvSolution', HEADER.format(cls='dictionary', obj='fvSolution') + '''
solvers
{
    p_rgh { solver GAMG; tolerance 1e-10; relTol 0; smoother DIC; }
    p_rghFinal { $p_rgh; }
    "U.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; }
    "T.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-15; relTol 0; minIter 1; maxIter 50; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors     3;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
    nThermoCorrectors   20;
    thermoTolerance     1e-8;
}
relaxationFactors { equations { ".*" 1; } }
''')
    temperature(case, nx, ny, L, p)
    field(case, 'U', 'volVectorField', '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)',
          {'"(left|right|bottom)"': ['type noSlip'],
           'top': ['type marangoni', f'dSigmadT {p["dSigmadT"]}', 'value uniform (0 0 0)'],
           'frontAndBack': ['type empty']})
    field(case, 'p_rgh', 'volScalarField', '[0 2 -2 0 0 0 0]', 'uniform 0',
          {'"(left|right|top|bottom)"': ['type fixedFluxPressure', 'rho rhok', 'value uniform 0'],
           'frontAndBack': ['type empty']})
    field(case, 'alpha.solid', 'volScalarField', '[0 0 0 0 0 0 0]', 'uniform 0',
          {'"(left|right|top|bottom)"': ['type zeroGradient'], 'frontAndBack': ['type empty']})
    return dict(case=str(case), solver='additiveFoam', nx=nx, ny=ny, dx=dx, **p)


def af_fields(case, meta, tdir):
    nx, ny, dx = meta['nx'], meta['ny'], meta['dx']
    U = read(tdir/'U', 3).reshape(ny, nx, 3)
    return dict(t=float(tdir.name), x=(np.arange(nx) + 0.5)*dx, y=(np.arange(ny) + 0.5)*dx,
                u=U[..., 0], v=U[..., 1])


def run_parallel(case, solver, n, env=None, log='log.solver'):
    """blockMesh, decomposition into n subdomains (scotch) and the solver in parallel."""
    write(Path(case)/'system/decomposeParDict', HEADER.format(cls='dictionary', obj='decomposeParDict')
          + f'numberOfSubdomains {n};\nmethod scotch;\n')
    pre = (env + '; ') if env else ''
    cmd = (f'{pre}blockMesh > log.blockMesh 2>&1 && decomposePar > log.decomposePar 2>&1 && '
           f'mpirun -np {n} --bind-to none {solver} -parallel > {log} 2>&1 && reconstructPar > log.reconstructPar 2>&1')
    return subprocess.run(['bash', '-c', cmd], cwd=case).returncode
