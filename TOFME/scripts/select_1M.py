#!/usr/bin/env python3
"""Select ~N events of P24iy 19.6 GeV picoDst spread over the run range: nRuns runs evenly spaced
in run order, nPerRun largest files per run. Writes filelists/me_1M/job_XXX.list (one file per job)
and filelists/me_1M_files.txt (run file events)."""
import collections, os, sys

target = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
nPerRun = int(sys.argv[2]) if len(sys.argv) > 2 else 3
base = "/gpfs01/star/pwg/dongx/tof/TOFME/filelists"
xrd = "root://xrdstar.rcf.bnl.gov:1095/"
runs = collections.defaultdict(list)
for line in open(f"{base}/local_P24iy_all.txt"):
    run, path, fn, ev, node = line.split()
    runs[int(run)].append((int(ev), f"{path}/{fn}"))
good = sorted(r for r, fl in runs.items() if sum(1 for e, _ in fl if e > 10000) >= nPerRun)
# events per run from nPerRun largest files
perRun = {r: sorted(runs[r], reverse=True)[:nPerRun] for r in good}
avg = sum(sum(e for e, _ in perRun[r]) for r in good) / len(good)
nRuns = min(len(good), int(round(target / avg)) + 1)
step = len(good) / nRuns
chosen = [good[int(i * step)] for i in range(nRuns)]
os.makedirs(f"{base}/me_1M", exist_ok=True)
tot = 0; j = 0
with open(f"{base}/me_1M_files.txt", "w") as out:
    for r in chosen:
        for e, f in perRun[r]:
            out.write(f"{r} {xrd}{f} {e}\n")
            with open(f"{base}/me_1M/job_{j:03d}.list", "w") as jl:
                jl.write(f"{xrd}{f}\n")
            j += 1; tot += e
print(f"runs with >= {nPerRun} files >10k events: {len(good)} (run {good[0]}..{good[-1]})")
print(f"chosen runs: {nRuns}, files/jobs: {j}, events: {tot}")
print("runs:", " ".join(str(r) for r in chosen))
