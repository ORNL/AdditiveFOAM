"""Reference solutions for the thermocapillary slot of Sen and Davis (1982).

A liquid layer of depth h fills a slot of length L between rigid side walls on
a no-slip bottom. Its upper surface carries the Marangoni stress
tau = dsigma/dT dT/dx of a linear temperature field. Gravity is zero and the
capillary number is small, so the surface stays close to flat and its
deflection eta(x) follows from the normal-stress balance of the flow under a
flat surface:

    sigma(x) eta'' = -(N(x) - c),   eta'(0) = eta'(L) = 0 (90 degree contact),
    integral of eta = 0 (mass),

with N the jump of normal stress across the surface. An optional passive gas
layer of height H sits under a no-slip lid.

Two references:

thin_slot(...)   Sen and Davis leading order in the aspect ratio A = h/L:
                 the core return flow u = (tau/mu)(3y^2/(4h) - y/2) and the
                 cubic deflection eta = tau L^3/(32 h sigma) (3 xi - xi^3),
                 xi = (2x - L)/L, with the exact two-layer factor for a gas
                 layer.

stokes(...)      The same small-capillary-number problem at finite aspect
                 ratio, solved without the thin-slot approximation: Stokes
                 flow in the liquid (and gas) on a staggered grid with the
                 shear jump on the flat interface, then the normal-stress
                 balance for eta. It contains the corner regions that the
                 thin-slot expansion leaves out and is second order in the
                 grid spacing.

flat_surface(...) The liquid alone under a flat Marangoni-slip surface: the
                 problem the Marangoni boundary condition of additiveFoam
                 solves (stokes with no gas layer).
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def two_layer_factor(h, H, mu1, mu2):
    """Ratio of the thin-slot pressure-jump gradient with a gas layer of
    height H under a no-slip lid to that of the liquid alone, 3 tau/(2h).
    Both layers carry zero net flux; derived in the README."""
    if H is None or H <= 0 or mu2 == 0:
        return 1.0
    return (mu1 - mu2*h**2/H**2)/(mu1 + mu2*h/H)


def thin_slot(x, L, h, tau, sigma, mu1, H=None, mu2=0.0):
    """Leading-order deflection of Sen and Davis with uniform sigma."""
    xi = (2*np.asarray(x) - L)/L
    f = two_layer_factor(h, H, mu1, mu2)
    return f*tau*L**3/(32*h*sigma)*(3*xi - xi**3)


def core_velocity(y, h, tau, mu):
    """Thin-slot return flow in the liquid, y from the bottom."""
    y = np.asarray(y)
    return tau/mu*(3*y**2/(4*h) - y/2)


def deflection(x, N, sigma):
    """Solve sigma(x) eta'' = -(N - c) on cell centres x with eta' = 0 at
    both walls and zero mean; c makes the problem solvable."""
    n = len(x)
    dx = x[1] - x[0]
    sigma = np.broadcast_to(sigma, (n,)).astype(float)
    c = np.sum(N/sigma)/np.sum(1/sigma)
    rhs = -(N - c)/sigma
    A = np.zeros((n + 1, n))
    for i in range(n):
        if i > 0:
            A[i, i - 1] += 1/dx**2
            A[i, i] -= 1/dx**2
        if i < n - 1:
            A[i, i + 1] += 1/dx**2
            A[i, i] -= 1/dx**2
    A[n, :] = 1.0
    b = np.append(rhs, 0.0)
    eta = np.linalg.lstsq(A, b, rcond=None)[0]
    return eta


def stokes(L, h, tau, mu1, nx, ny1, H=0.0, mu2=0.0, sigma=None, ny2=None):
    """Stokes flow in a slot with a flat interface at y = h carrying the
    shear jump tau, liquid below (viscosity mu1), optional gas above up to a
    lid at y = h + H (viscosity mu2). No-slip on the bottom, the side walls
    and the lid. Returns the fields, the interface normal-stress jump and,
    when sigma is given (scalar or a function of x), the deflection.

    Without gas the interface is the top boundary with the imposed shear,
    which is the flat-surface problem of the Marangoni boundary condition."""
    gas = H > 0
    dx = L/nx
    dy = h/ny1
    if gas:
        ny2 = ny2 if ny2 is not None else int(round(H/dy))
        if abs(ny2*dy - H) > 1e-9*H:
            raise ValueError('H must be a multiple of h/ny1')
    else:
        ny2 = 0
    ny = ny1 + ny2
    js = ny1                                    # interface face index
    mu = np.where(np.arange(ny) < ny1, mu1, mu2)

    # unknown numbering: u at faces i = 1..nx-1, rows j; v at faces
    # j = 1..ny-1 except js, columns i; p in every cell
    uid = -np.ones((ny, nx + 1), int)
    vid = -np.ones((ny + 1, nx), int)
    k = 0
    for j in range(ny):
        for i in range(1, nx):
            uid[j, i] = k
            k += 1
    for j in range(1, ny):
        if gas and j == js:
            continue
        for i in range(nx):
            vid[j, i] = k
            k += 1
    pid = (np.arange(ny*nx) + k).reshape(ny, nx)
    size = k + ny*nx

    rows, cols, vals = [], [], []
    rhs = np.zeros(size)

    def add(r, c, v):
        rows.append(r)
        cols.append(c)
        vals.append(v)

    # x-momentum: -mu u_xx - d(shear)/dy + dp/dx = 0
    for j in range(ny):
        m = mu[j]
        for i in range(1, nx):
            r = uid[j, i]
            diag = 2*m/dx**2
            for ii in (i - 1, i + 1):
                if 0 < ii < nx:
                    add(r, uid[j, ii], -m/dx**2)
            # bottom face of this control volume
            if j == 0:
                diag += 2*m/dy**2                      # wall: u = 0 at y = 0
            elif gas and j == js:
                # interface below a gas row: shear from the gas side
                g = 2*mu1*mu2/((mu1 + mu2)*dy)
                diag += g/dy
                add(r, uid[j - 1, i], -g/dy)
                rhs[r] += -mu2/(mu1 + mu2)*tau/dy
            else:
                diag += m/dy**2
                add(r, uid[j - 1, i], -m/dy**2)
            # top face of this control volume
            if j == ny - 1:
                if gas:
                    diag += 2*m/dy**2                  # lid
                else:
                    rhs[r] += tau/dy                   # imposed shear
            elif gas and j == js - 1:
                g = 2*mu1*mu2/((mu1 + mu2)*dy)
                diag += g/dy
                add(r, uid[j + 1, i], -g/dy)
                rhs[r] += mu1/(mu1 + mu2)*tau/dy
            else:
                diag += m/dy**2
                add(r, uid[j + 1, i], -m/dy**2)
            add(r, r, diag)
            add(r, pid[j, i], 1/dx)
            add(r, pid[j, i - 1], -1/dx)

    # y-momentum at interior horizontal faces of each layer
    for j in range(1, ny):
        if gas and j == js:
            continue
        m = mu[j - 1] if j <= ny1 else mu[j]
        for i in range(nx):
            r = vid[j, i]
            diag = 2*m/dy**2
            for jj in (j - 1, j + 1):
                if vid[jj, i] >= 0:
                    add(r, vid[jj, i], -m/dy**2)
            for ii in (i - 1, i + 1):
                if 0 <= ii < nx:
                    add(r, vid[j, ii], -m/dx**2)
                    diag += m/dx**2
                else:
                    diag += 2*m/dx**2                  # side wall: v = 0
            add(r, r, diag)
            add(r, pid[j, i], 1/dy)
            add(r, pid[j - 1, i], -1/dy)

    # continuity, one gauge per closed layer
    gauges = {(0, 0)}
    if gas:
        gauges.add((js, 0))
    for j in range(ny):
        for i in range(nx):
            r = pid[j, i]
            if (j, i) in gauges:
                add(r, r, 1.0)
                continue
            if uid[j, i + 1] >= 0:
                add(r, uid[j, i + 1], 1/dx)
            if uid[j, i] >= 0:
                add(r, uid[j, i], -1/dx)
            if vid[j + 1, i] >= 0:
                add(r, vid[j + 1, i], 1/dy)
            if vid[j, i] >= 0:
                add(r, vid[j, i], -1/dy)

    A = coo_matrix((vals, (rows, cols)), shape=(size, size)).tocsr()
    sol = spsolve(A, rhs)

    u = np.zeros((ny, nx + 1))
    v = np.zeros((ny + 1, nx))
    for j in range(ny):
        u[j, 1:nx] = sol[uid[j, 1:nx]]
    for j in range(1, ny):
        if vid[j, 0] >= 0:
            v[j] = sol[vid[j]]
    p = sol[pid]
    div = np.diff(u, axis=1)/dx + np.diff(v, axis=0)/dy

    # normal-stress jump N = (p - 2 mu v_y)_liquid - (p - 2 mu v_y)_gas
    pl = 1.5*p[js - 1] - 0.5*p[js - 2]
    vyl = (4*v[js - 1] - v[js - 2])/(2*dy)            # v = 0 at js
    N = pl + 2*mu1*vyl
    N = pl - 2*mu1*(-vyl)
    N = pl - 2*mu1*((3*0 - 4*v[js - 1] + v[js - 2])/(2*dy))
    if gas:
        pg = 1.5*p[js] - 0.5*p[js + 1]
        vyg = (-3*0 + 4*v[js + 1] - v[js + 2])/(2*dy)
        N = N - (pg - 2*mu2*vyg)

    x = (np.arange(nx) + 0.5)*dx
    yc = (np.arange(ny) + 0.5)*dy
    out = dict(x=x, y=yc, xf=np.arange(nx + 1)*dx, yf=np.arange(ny + 1)*dy,
               u=u, v=v, p=p, N=N, divergence=float(np.abs(div).max()),
               residual=float(np.abs(A@sol - rhs).max()),
               surface_u=0.5*(u[js - 1, :-1] + u[js - 1, 1:]) if gas else None)
    if not gas:
        # surface velocity of the flat slip surface: u_top = u + tau dy/(2 mu)
        uc = 0.5*(u[ny - 1, :-1] + u[ny - 1, 1:])
        out['surface_u'] = uc + tau*dy/(2*mu1)
    else:
        ul = 0.5*(u[js - 1, :-1] + u[js - 1, 1:])
        ug = 0.5*(u[js, :-1] + u[js, 1:])
        out['surface_u'] = (tau*dy/2 + mu1*ul + mu2*ug)/(mu1 + mu2)
    if sigma is not None:
        s = sigma(x) if callable(sigma) else np.full(nx, float(sigma))
        out['eta'] = deflection(x, N, s)
        e = out['eta']
        out['endpoints'] = np.array([(9*e[0] - e[1])/8, (9*e[-1] - e[-2])/8])
    return out


def flat_surface(L, h, tau, mu, nx, ny):
    """The liquid alone under a flat surface carrying the Marangoni shear:
    the problem of the Marangoni slip boundary condition."""
    return stokes(L, h, tau, mu, nx, ny)


if __name__ == '__main__':
    # self-checks: thin-slot limit, grid convergence, two-layer factor
    tau, sigma, mu1, mu2 = 0.4, 0.01, 0.005, 2.83e-5
    h = 2e-4
    print('A = h/L   endpoint / thin-slot endpoint (liquid only, ny1 = 16 per h)')
    for L in (1e-3, 2e-3, 4e-3):
        nx = int(round(16*L/h))
        r = stokes(L, h, tau, mu1, nx, 16, sigma=sigma)
        ref = thin_slot(np.array([L]), L, h, tau, sigma, mu1)[0]
        print(f'{h/L:6.3f}   {r["endpoints"][1]/ref:.5f}   div {r["divergence"]:.1e}')
    print('grid convergence, A = 0.2, liquid and gas (H = 0.3 mm)')
    for n1 in (8, 16, 32):
        r = stokes(1e-3, h, tau, mu1, 5*n1, n1, H=3e-4, mu2=mu2, sigma=sigma)
        print(f'ny1 {n1:3d}: endpoint {r["endpoints"][1]*1e6:.4f} um')
    print('two-layer factor', two_layer_factor(h, 3e-4, mu1, mu2))
