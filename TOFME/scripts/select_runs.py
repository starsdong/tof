#!/usr/bin/env python3
"""Select ~target events of P24iy 19.6 GeV picoDst as one job per run with >= minEv events per job.
Runs in the StRefMultCorr bad-run list are skipped. The catalog list has one line per replica: files are
deduplicated by name, and files with a replica on a node in avoidNodes (stalled XRootD server) are dropped. For each run the largest files are taken until
the job has >= minEv events (at most maxFiles files). Runs are chosen evenly spaced in run order.
usage: select_runs.py <target> <minEv> <tag> [avoidNode,...]
writes filelists/<tag>/job_<run>.list, filelists/<tag>_files.txt (run url events), filelists/<tag>_jobids.txt"""
import collections, os, re, sys

target = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
minEv = int(sys.argv[2]) if len(sys.argv) > 2 else 60_000
tag = sys.argv[3] if len(sys.argv) > 3 else "me_runs"
avoidNodes = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else set()
maxFiles = 8
base = "/gpfs01/star/pwg/dongx/tof/TOFME/filelists"
xrd = "root://xrdstar.rcf.bnl.gov:1095/"
badh = open("/star/nfs4/AFS/star/packages/SL24y/StRoot/StRefMultCorr/BadRun.h").read()
blk = badh[badh.index("badrun_refmult_2019["):]
blk = blk[blk.index("{"):blk.index("};")]
bad = set(int(x) for x in re.findall(r"\b(20\d{6})\b", blk))
files = {}          # unique file name -> (run, events, path, set of nodes)
for line in open(f"{base}/local_P24iy_all.txt"):
    run, path, fn, ev, node = line.split()
    if fn not in files:
        files[fn] = [int(run), int(ev), f"{path}/{fn}", set()]
    files[fn][3].add(node.split(".")[0])
runs = collections.defaultdict(list)
nAvoid = 0
for fn, (run, ev, full, nodes) in files.items():
    if nodes & avoidNodes:
        nAvoid += 1
        continue
    runs[run].append((ev, full))
sel = {}
for r, fl in runs.items():
    if r in bad or not (20056032 <= r <= 20093036):
        continue
    tot, take = 0, []
    for e, f in sorted(fl, reverse=True)[:maxFiles]:
        take.append((e, f)); tot += e
        if tot >= minEv:
            break
    if tot >= minEv:
        sel[r] = take
good = sorted(sel)
avg = sum(sum(e for e, _ in sel[r]) for r in good) / len(good)
nRuns = min(len(good), int(round(target / avg)) + (1 if target % avg else 0))
step = len(good) / nRuns
chosen = [good[int(i * step)] for i in range(nRuns)]
os.makedirs(f"{base}/{tag}", exist_ok=True)
tot = 0; nf = 0
with open(f"{base}/{tag}_files.txt", "w") as out, open(f"{base}/{tag}_jobids.txt", "w") as ids:
    for r in chosen:
        with open(f"{base}/{tag}/job_{r}.list", "w") as jl:
            for e, f in sel[r]:
                out.write(f"{r} {xrd}{f} {e}\n"); jl.write(f"{xrd}{f}\n"); tot += e; nf += 1
        ids.write(f"{r}\n")
print(f"catalog unique files {len(files)}, dropped (replica on {sorted(avoidNodes)}) {nAvoid}")
print(f"bad runs (StRefMultCorr 2019 list): {len(bad)}; eligible runs with >= {minEv} events in <= {maxFiles} files: {len(good)}")
print(f"chosen runs/jobs: {nRuns}, files: {nf}, events: {tot}, events/job: {tot/nRuns:.0f}")
print("runs:", " ".join(str(r) for r in chosen))
