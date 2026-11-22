"""AdditiveFOAM against the thermocapillary slot references.

    python3 plot.py <runs directory> [<output html>]

The runs directory holds case folders written by slot.vof_case (names
starting with "vof-") and slot.af_case ("af-"), each with its meta.json.
Writes one self-contained HTML page with the figures and the numbers.
"""

import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference as ref
import slot

plt.rcParams.update({'font.size': 9})
COLORS = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']


def png(fig):
    b = io.BytesIO()
    fig.savefig(b, format='png', dpi=110)
    plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()


def load(runs, prefix):
    out = {}
    for d in sorted(Path(runs).glob(prefix + '*')):
        if (d/'meta.json').exists() and slot.times(d):
            out[d.name] = (d, json.loads((d/'meta.json').read_text()))
    return out


def reference_deflection(p, n1=64):
    """Finite-aspect reference on two grids and its Richardson value."""
    tau = slot.tau(p)
    rs = [ref.stokes(p['L'], p['h'], tau, p['mu1'], int(n1*p['L']/p['h']*f), int(n1*f), H=p['H'],
                     mu2=p['mu2'], sigma=lambda x: slot.sigma_of_x(x, p)) for f in (0.5, 1.0)]
    x = rs[1]['x']
    e0 = np.interp(x, rs[0]['x'], rs[0]['eta'])
    order = 1.15                            # measured, corner-limited (README)
    eta = rs[1]['eta'] + (rs[1]['eta'] - e0)/(2**order - 1)
    ends = rs[1]['endpoints'] + (rs[1]['endpoints'] - rs[0]['endpoints'])/(2**order - 1)
    return x, eta, ends, rs[1]


def main(runs, html, extras=(), estimators=None, extra_captions={}):
    p = slot.P
    tau = slot.tau(p)
    figs, rows = [], []

    # ---- 1. the reference itself -------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True)
    A, ratio = [], []
    for L in (1e-3, 2e-3, 4e-3, 8e-3):
        r = ref.stokes(L, p['h'], tau, p['mu1'], int(16*L/p['h']), 16, sigma=p['sigma'])
        A.append(p['h']/L)
        ratio.append(r['endpoints'][1]/ref.thin_slot(np.array([L]), L, p['h'], tau, p['sigma'], p['mu1'])[0])
    ax[0].plot(A, ratio, 'o-')
    ax[0].axhline(1, color='k', lw=0.6)
    ax[0].set_xlabel('aspect ratio A = h/L')
    ax[0].set_ylabel('endpoint / Sen-Davis thin-slot endpoint')
    ax[0].set_title('Finite-aspect reference tends to Sen and Davis as A -> 0')
    ns, ends = [8, 16, 32, 64], []
    for n in ns:
        r = ref.stokes(p['L'], p['h'], tau, p['mu1'], 5*n, n, H=p['H'], mu2=p['mu2'], sigma=p['sigma'])
        ends.append(r['endpoints'][1])
    ends = np.array(ends)
    ax[1].loglog([p['h']/n*1e6 for n in ns[:-1]], np.abs(ends[:-1] - ends[-1])*1e6, 'o-', label='endpoint')
    ax[1].set_xlabel('reference grid spacing (um)')
    ax[1].set_ylabel('|endpoint - finest| (um)')
    ax[1].set_title('Reference grid convergence (corner-limited, order ~1.2)')
    figs.append(('The reference', png(fig),
                 'Left: the exact small-capillary-number Stokes solution of the slot at finite aspect ratio '
                 'approaches the Sen and Davis thin-slot deflection linearly in A; at the benchmark A = 0.2 the '
                 'endpoint is 1.5% below the thin-slot value. Right: grid convergence of the reference itself. '
                 'The order is limited by the stress singularity where the Marangoni shear meets the no-slip wall; '
                 'the Richardson value used below is good to about 0.2%.'))

    # ---- 2. flat surface, additiveFoam Marangoni condition -----------------
    af = load(runs, 'af-')
    if af:
        r = ref.flat_surface(p['L'], p['h'], tau, p['mu1'], 400, 80)
        uc = np.vstack([np.zeros(r['u'].shape[1]), r['u'], np.zeros(r['u'].shape[1])])
        yc = np.concatenate([[0], r['y'], [p['h']]])
        # top boundary value from the imposed shear: u_top = u_row + tau dy/(2 mu)
        uc[-1] = r['u'][-1] + tau*(r['y'][1] - r['y'][0])/(2*p['mu1'])
        interp = RegularGridInterpolator((yc, r['xf']), uc)
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
        yy = np.linspace(0, p['h'], 200)
        ax[0].plot(ref.core_velocity(yy, p['h'], tau, p['mu1'])*1e3, yy*1e6, 'k:', lw=1, label='Sen-Davis core flow')
        ax[0].plot(interp(np.c_[yy, np.full_like(yy, p['L']/2)])*1e3, yy*1e6, 'k-', lw=1, label='finite-aspect Stokes')
        errs = []
        for i, (name, (d, meta)) in enumerate(sorted(af.items(), key=lambda kv: kv[1][1]['nx'])):
            f = slot.af_fields(d, meta, slot.times(d)[-1])
            nx = meta['nx']
            um = 0.5*(f['u'][:, nx//2 - 1] + f['u'][:, nx//2])
            ax[0].plot(um*1e3, f['y']*1e6, 'o', ms=3, color=COLORS[i], label=f'additiveFoam {nx} x {meta["ny"]}')
            X, Y = np.meshgrid(f['x'], f['y'])
            uref = interp(np.c_[Y.ravel(), X.ravel()]).reshape(Y.shape)
            e = np.sqrt(np.mean((f['u'] - uref)**2))/np.sqrt(np.mean(uref**2))
            emid = np.abs(um - uref[:, nx//2]).max()/np.abs(uref).max()
            errs.append((meta['dx'], e, emid))
            rows.append(('flat surface (additiveFoam)', f'{nx} x {meta["ny"]}', f'{f["t"]*1e3:.0f} ms',
                         f'u field L2 {e:.4f}', f'mid-slot profile max {emid:.4f}'))
            ax[1].plot(f['x']*1e3, f['u'][-1]*1e3, '-', color=COLORS[i], label=f'{nx} x {meta["ny"]}, top cell row')
            ax[1].plot(f['x']*1e3, uref[-1]*1e3, ':', color=COLORS[i])
        ax[0].set_xlabel('u (mm/s) at mid-slot')
        ax[0].set_ylabel('height above the bottom (um)')
        ax[0].legend(fontsize=7)
        ax[0].set_title('Return flow under the Marangoni condition')
        ax[1].set_xlabel('x (mm)')
        ax[1].set_ylabel('u (mm/s)')
        ax[1].set_title('Top cell row; dotted: reference at the same points')
        ax[1].legend(fontsize=7)
        e = np.array(errs)
        ax[2].loglog(e[:, 0]*1e6, e[:, 1], 'o-', label='u field, relative L2')
        ax[2].loglog(e[:, 0]*1e6, e[:, 2], 's-', label='mid-slot profile, max')
        ax[2].loglog(e[:, 0]*1e6, e[0, 1]*(e[:, 0]/e[0, 0])**2, 'k--', lw=0.7, label='second order')
        ax[2].set_xlabel('cell size (um)')
        ax[2].set_ylabel('error')
        ax[2].legend(fontsize=7)
        ax[2].set_title('Convergence')
        figs.append(('Flat surface: the additiveFoam Marangoni condition', png(fig),
                     'additiveFoam with the marangoni slip condition on the top of a closed slot, L = 1 mm, '
                     'h = 0.2 mm, tau = 0.4 Pa, zero gravity, steady state. The reference is the Stokes flow '
                     'of the same slot with the imposed surface shear, solved without the thin-slot '
                     'approximation; the Sen and Davis core profile u = (tau/mu)(3y^2/4h - y/2) is shown for '
                     'comparison and holds away from the walls.'))

    # ---- 3. free surface, additiveVoFFoam deflection -----------------------
    vof = load(runs, 'vof-')
    if vof:
        xr, er, ends_r, r1 = reference_deflection(p)
        thin = ref.thin_slot(xr, p['L'], p['h'], tau, p['sigma'], p['mu1'], p['H'], p['mu2'])
        main_runs = {k: v for k, v in vof.items() if '-' not in k[4:]}
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
        ax[0].plot(xr*1e3, er*1e6, 'k-', lw=1.2, label='finite-aspect reference')
        ax[0].plot(xr*1e3, thin*1e6, 'k:', lw=1, label='Sen-Davis thin slot')
        conv = []
        for i, (name, (d, meta)) in enumerate(sorted(main_runs.items(), key=lambda kv: kv[1][1]['nx'])):
            ts = slot.times(d)
            pr = slot.vof_profile(d, meta, ts[-1])
            ax[0].plot(pr['x']*1e3, pr['eta']*1e6, 'o-', ms=2.5, lw=0.8, color=COLORS[i],
                       label=f'{meta.get("scheme", "MULES")} {meta["nx"]} x {meta["ny"]}, {pr["t"]*1e3:.0f} ms')
            hist = [slot.vof_profile(d, meta, t) for t in ts]
            ax[1].plot([h['t']*1e3 for h in hist], [(h['eta'][-1] - h['eta'][0])/(ends_r[1] - ends_r[0]) for h in hist],
                       'o-', ms=3, color=COLORS[i], label=f'{meta["nx"]} x {meta["ny"]}')
            e_ref = np.interp(pr['x'], xr, er)
            l2 = np.sqrt(np.mean((pr['eta'] - e_ref)**2))/np.sqrt(np.mean(e_ref**2))
            conv.append((meta['dx'], l2))
            e = pr['eta']
            yl, yr = (p['h'] + (9*e[0] - e[1])/8)/p['L'], (p['h'] + (9*e[-1] - e[-2])/8)/p['L']
            rows.append(('free surface (additiveVoFFoam)', f'{meta["nx"]} x {meta["ny"]}', f'{pr["t"]*1e3:.0f} ms',
                         f'profile L2 {l2:.3f}', f'wall heights y = Y/w {yl:.4f} / {yr:.4f} (reference {(p["h"] + ends_r[0])/p["L"]:.4f} / {(p["h"] + ends_r[1])/p["L"]:.4f}; Sen and Davis in Wirth 2022 table 6: 0.188 / 0.213)'))
        ax[0].set_xlabel('x (mm)')
        ax[0].set_ylabel('surface deflection (um)')
        ax[0].set_title('Steady surface deflection')
        ax[0].legend(fontsize=7)
        ax[1].axhline(1, color='k', lw=0.6)
        ax[1].set_xlabel('time (ms)')
        ax[1].set_ylabel('end-to-end rise / reference')
        ax[1].set_title('Approach to steady state')
        ax[1].legend(fontsize=7)
        c = np.array(conv)
        ax[2].loglog(c[:, 0]*1e6, c[:, 1], 'o-')
        ax[2].set_xlabel('cell size (um)')
        ax[2].set_ylabel('profile error, relative L2')
        ax[2].set_title('Error at the last write')
        figs.append(('Free surface: the Sen and Davis deflection', png(fig),
                     'additiveVoFFoam, liquid 0.2 mm under 0.3 mm of gas in a closed 1 mm slot, zeroGradient alpha (90 degree) '
                     'contact, zero gravity, tau = 0.4 Pa, sigma = 0.01 N/m. The deflection is measured from '
                     'the metal volume in each column.'))


    # ---- 4. geometric transport: steady state -------------------------------
    for extra in extras:
        ex = load(extra, 'vof-')
        if not ex:
            continue
        xr, er, ends_r, _ = reference_deflection(p)
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
        ax[0].plot(xr*1e3, er*1e6, 'k-', lw=1.2, label='finite-aspect reference')
        for i, (name, (d, meta)) in enumerate(sorted(ex.items(), key=lambda kv: kv[1][1]['nx'])):
            ts = slot.times(d)
            hist = [slot.vof_profile(d, meta, t) for t in ts]
            pr = hist[-1]
            ax[0].plot(pr['x']*1e3, pr['eta']*1e6, 'o-', ms=2.5, lw=0.8, color=COLORS[i],
                       label=f'{meta.get("scheme", "MULES")} {meta["nx"]} x {meta["ny"]}, {pr["t"]*1e3:.0f} ms')
            ax[1].plot([q['t']*1e3 for q in hist], [(q['eta'][-1] - q['eta'][0])/(ends_r[1] - ends_r[0]) for q in hist],
                       'o-', ms=3, color=COLORS[i], label=f'{meta.get("scheme", "MULES")} {meta["nx"]}')
            band = [((q['alpha'] > 0.01) & (q['alpha'] < 0.99)).sum(axis=0).mean() for q in hist]
            ax[2].plot([q['t']*1e3 for q in hist], band, 'o-', ms=3, color=COLORS[i], label=f'{meta.get("scheme", "MULES")} {meta["nx"]}')
            e_ref = np.interp(pr['x'], xr, er)
            l2 = np.sqrt(np.mean((pr['eta'] - e_ref)**2))/np.sqrt(np.mean(e_ref**2))
            rows.append((f'free surface, {meta.get("scheme", "MULES")} ({Path(extra).name})', f'{meta["nx"]} x {meta["ny"]}', f'{pr["t"]*1e3:.0f} ms',
                         f'profile L2 {l2:.3f}', f'end-to-end rise {(pr["eta"][-1] - pr["eta"][0])/(ends_r[1] - ends_r[0]):.3f} of the reference'))
        ax[0].set_xlabel('x (mm)'); ax[0].set_ylabel('surface deflection (um)'); ax[0].legend(fontsize=7); ax[0].set_title('Deflection')
        ax[1].axhline(1, color='k', lw=0.6); ax[1].set_xlabel('time (ms)'); ax[1].set_ylabel('end-to-end rise / reference'); ax[1].legend(fontsize=7)
        ax[1].set_title('Approach to steady state')
        ax[2].set_xlabel('time (ms)'); ax[2].set_ylabel('interface cells per column'); ax[2].legend(fontsize=7); ax[2].set_title('Interface band width')
        figs.append((f'Long runs: {Path(extra).name}', png(fig), extra_captions.get(Path(extra).name, '')))

    # ---- 5. the curvature estimator on the exact interface -------------------
    if estimators and Path(estimators).exists():
        est = json.loads(Path(estimators).read_text())
        labels = ['additiveVoFFoam SSF estimator', 'without curvature smoothing', 'without alpha smoothing', 'without either', 'reconstructed distance function']
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), constrained_layout=True)
        nx_show = '80'
        for k, lab in enumerate(labels):
            if k == 2:
                continue
            q = est[nx_show][k]
            ax[0].plot(np.array(q['x'])*1e3, q['eff'], '.-', ms=3, lw=0.7, label=lab)
        ax[0].plot(np.array(est[nx_show][0]['x'])*1e3, est[nx_show][0]['exact'], 'k-', lw=1.2, label='exact')
        ax[0].set_ylim(-400, 400); ax[0].set_xlabel('x (mm)'); ax[0].set_ylabel('curvature seen by the surface force (1/m)')
        ax[0].set_title(f'Exact interface, {nx_show} x {int(nx_show)//2}'); ax[0].legend(fontsize=6.5)
        ns = sorted(est, key=int)
        for k, lab in enumerate(labels):
            if k == 2:
                continue
            slope, err = [], []
            for n in ns:
                q = est[n][k]; x = np.array(q['x']); e = np.array(q['eff']); ke = np.array(q['exact'])
                m = (x > 0.2*p['L']) & (x < 0.8*p['L'])
                slope.append(np.polyfit(ke[m], e[m], 1)[0]); err.append(np.sqrt(np.mean((e[m] - ke[m])**2))/np.sqrt(np.mean(ke[m]**2)))
            dxs = [p['L']/int(n)*1e6 for n in ns]
            ax[1].plot(dxs, slope, 'o-', label=lab); ax[2].loglog(dxs, err, 'o-', label=lab)
        ax[1].axhline(1, color='k', lw=0.6); ax[1].set_xscale('log'); ax[1].set_xlabel('cell size (um)'); ax[1].set_ylabel('mean curvature / exact, middle 60%')
        ax[1].set_title('Bias'); ax[1].legend(fontsize=6.5)
        ax[2].set_xlabel('cell size (um)'); ax[2].set_ylabel('pointwise error, relative L2'); ax[2].set_title('Scatter'); ax[2].legend(fontsize=6.5)
        figs.append(('The curvature estimator, isolated', png(fig),
                     'The surface force of additiveVoFFoam reproduced outside the solver and applied to the exact '
                     'reference interface at rest, so no flow or transport enters. The true curvature is about '
                     '150 1/m (kappa dx = 0.001 to 0.007), which makes this a far harder curvature test than a '
                     'droplet. The solver estimator is biased by 5 to 50% and its bias does not converge; '
                     'removing the curvature smoothing removes most of the bias; a reconstructed distance '
                     'function (a first, unoptimised implementation) has the right mean but large pointwise '
                     'scatter. The scatter of every estimator that differentiates volume fractions twice grows '
                     'with refinement here.'))

    # ---- page ------------------------------------------------------------
    table = ''.join(f'<tr>{"".join(f"<td>{c}</td>" for c in row)}</tr>' for row in rows)
    body = ''.join(f'<section><h2>{t}</h2><figure><img src="data:image/png;base64,{b}" alt="{t}">'
                   f'<figcaption>{c}</figcaption></figure></section>' for t, b, c in figs)
    page = f'''<title>Thermocapillary Slot</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>
:root {{ --bg:#f7f6f2; --fg:#1c1e22; --muted:#5c6170; --line:#d8d6ce; --panel:#fff; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#16181d; --fg:#e8e6e0; --muted:#a3a7b2; --line:#343842; --panel:#1f2229; color-scheme:dark }} }}
:root[data-theme="dark"] {{ --bg:#16181d; --fg:#e8e6e0; --muted:#a3a7b2; --line:#343842; --panel:#1f2229; color-scheme:dark }}
body {{ background:var(--bg); color:var(--fg); font-family:system-ui,-apple-system,"Segoe UI",sans-serif; line-height:1.45; margin:0;
  padding-left:max(16px,env(safe-area-inset-left,0px)); padding-right:max(16px,env(safe-area-inset-right,0px)); padding-block:24px 64px; }}
main {{ max-width:1050px; margin:0 auto; display:grid; gap:28px; }}
h1 {{ font-size:1.7rem; margin:0 0 4px }} h2 {{ font-size:1.15rem; margin:0 0 8px; border-bottom:1px solid var(--line); padding-bottom:4px }}
p, figcaption {{ max-width:85ch }} .lede, figcaption {{ color:var(--muted) }} figcaption {{ font-size:.88rem }}
figure {{ margin:0; display:grid; gap:6px }} img {{ width:100%; height:auto; background:#fff; border:1px solid var(--line); border-radius:4px }}
.tbl {{ overflow-x:auto; border:1px solid var(--line); border-radius:6px; background:var(--panel) }}
table {{ border-collapse:collapse; width:100%; font-size:.86rem }} td, th {{ padding:6px 10px; border-bottom:1px solid var(--line); text-align:left; white-space:nowrap }}
</style>
<main>
<section><h1>Thermocapillary Slot</h1>
<p class="lede">AdditiveFOAM against the Sen and Davis (1982) thermocapillary slot: the flat-surface return flow of the
additiveFoam Marangoni condition and the surface deflection of additiveVoFFoam, each against an exact finite-aspect
reference. Generated by test/thermocapillarySlot/plot.py.</p></section>
{body}
<section><h2>Numbers</h2><div class="tbl"><table><tr><th>case</th><th>mesh</th><th>time</th><th>error</th><th></th></tr>{table}</table></div></section>
</main>
'''
    Path(html).write_text(page)
    print('wrote', html, len(page)//1024, 'kB')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('runs'); ap.add_argument('html', nargs='?', default='thermocapillarySlot.html')
    ap.add_argument('--extra', action='append', default=[], help='further run directories (e.g. long or MPLIC runs)')
    ap.add_argument('--estimators', help='JSON written by the curvature-estimator study')
    a = ap.parse_args()
    main(a.runs, a.html, a.extra, a.estimators)
