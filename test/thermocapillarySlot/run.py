"""Run the thermocapillary slot tests and write the comparison page.

    python3 run.py <output directory> [--np N] [--endTime T]

Builds and runs, in the current OpenFOAM-14 and AdditiveFOAM environment:

  vof-100x240   additiveVoFFoam, the Sen and Davis surface deflection on the
                100 x 240 mesh of the original marangoni-test tutorial,
                interfaceCompression (MULES) transport;
  af-100x20     additiveFoam, the flat-surface return flow under the
                marangoni slip condition;

then writes <output directory>/thermocapillarySlot.html with plot.py.
"""

import argparse
import json
from pathlib import Path

import slot
import plot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('output')
    ap.add_argument('--np', type=int, default=1, help='MPI ranks per case (default serial)')
    ap.add_argument('--endTime', type=float, default=None, help='override the end time (s), e.g. for a smoke test')
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    p = dict(slot.P)
    if a.endTime:
        p.update(endTime=a.endTime, writeInterval=a.endTime/2)
    cases = [('vof-100x240', 'additiveVoFFoam', lambda c: slot.vof_case(c, p=p)),
             ('af-100x20', 'additiveFoam', lambda c: slot.af_case(c, p=dict(p, alpha_thermal=1.6e-5),
                                                                  endTime=p['endTime']))]
    for name, solver, make in cases:
        case = out/name
        meta = make(case)
        (case/'meta.json').write_text(json.dumps(meta))
        rc = (slot.run_parallel(case, solver, a.np) if a.np > 1 else slot.run(case, solver))
        print(name, 'exit', rc, flush=True)
        if rc:
            raise SystemExit(f'{name} failed: see {case}/log.solver')
    plot.main(str(out), str(out/'thermocapillarySlot.html'))


if __name__ == '__main__':
    main()
