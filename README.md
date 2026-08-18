# NanoAOD Calib_ChargeTagger studies

Date: 2026-07-27

This repo studies `.root` outputs from [Calib_ChargeTagger_JH](https://github.com/jiashu-huang/Calib_ChargeTagger_JH).


## Set up 

Make sure you initiate `cmsenv` before you run!
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/jhuan166/Vcb/CMSSW_15_1_0_patch4/src
cmsenv
```

To recompile the C++:
```bash
make-clean
make
```

## Design

Provide one executable command-line entry point:
```bash
./calib-study <input-address-txt> <production-tag> [tree-name] [--Vcb]
```
Arguments:
- `input-address-txt`: a `.txt` containing one input `.root` file address per line.
- `production-tag`: a label identifying the NanoAOD production under study. 
    This value is used as the output subdirectory name, so it should be treated 
    as a filesystem path component rather than free-form prose.
- `tree-name`: optional ROOT `TTree` name. Default: `Events`.
- `--Vcb`: optional flag that additionally runs studies 9–11. It can appear
  before or after the optional `tree-name`.

Runtime behavior:
1. Read every non-comment `.root` path from `input-address-txt` and chain the requested `TTree` across them.
2. Open the chained input and retrieve the requested `TTree`.
3. In this repo, create `./outputs/<production_tag>/` if it does not already exist.
4. Run the inspection studies listed below.
5. Write all generated text reports into `./outputs/<production_tag>/`.

Every study is filled in a **single pass** over the tree, and only the branches
the studies actually read are enabled, so the cost of a run is one read of the
input rather than one read per plot.

### Splitting a run across machines

The studies are all histogram fills, so they split cleanly: each input file can
be histogrammed on its own, and the resulting histograms added. Two optional
mode flags expose that split.

```bash
# map: fill the histograms for one input list and write them to a partial file.
#      No reports and no plots are produced.
./calib-study <input-address-txt> <production-tag> [tree-name] [--Vcb] --partial <file.root>

# reduce: add every .root partial in a directory, then write the reports and plots.
./calib-study <input-address-txt> <production-tag> [tree-name] [--Vcb] --merge <partials-dir>
```

A partial file holds one named `TH1D` per study histogram, plus `cutflow_events`
and `cutflow_weight` (one labelled bin per cutflow step) and a
`meta_has_final_weight` flag. `--merge` adds them bin by bin, so the merged
result is identical to running the whole input list in one go.

Pass the same `--Vcb` setting to the map and merge steps. If the map jobs ran
without it, they never booked the studies 9–11 histograms and a `--Vcb` merge
stops with a missing-histogram error; in the other direction the merge warns and
plots only studies 3–8. `./calib-study-condor` always passes the same flag to
both, so this only comes up when running the two steps by hand.

`--merge` still needs `<input-address-txt>`, because studies 1 and 2 describe the
branch schema and entry 0. Those are shared by every file in a production, so
only the first file in the list is opened for them.

### Run through HTCondor

```bash
./calib-study-condor <input-address-txt> <production-tag> [tree-name] [--Vcb]
```

This submits a Condor **DAG** that splits the work across the input files:

- one **map** job per `.root` file in `<input-address-txt>`, each writing
  `outputs/<production-tag>/condor/partials/partial-NNN.root`;
- one **merge** job, released by DAGMan only once every map job has succeeded,
  which adds the partials and writes the reports and plots.

The map jobs run concurrently, so the wall time is roughly the time for the
single largest input file plus a few seconds to merge, instead of the sum over
all files.

Every worker initializes the CMSSW release named in
[`cmsenv-source.txt`](./cmsenv-source.txt), changes to this repository, and runs
the equivalent `./calib-study` command. Workers read the input list and ROOT
files directly; no files are transferred to the Condor sandbox.

Results are written to `./outputs/<production-tag>/`. The generated worker
scripts, submit descriptions, per-file input lists, partial histograms, and
Condor log/stdout/stderr files are kept in
`./outputs/<production-tag>/condor/`. Re-submitting the same production tag
clears `inputs/` and `partials/` first, so a stale partial from an earlier
submission is never summed into the merge.

Watch the DAG with `condor_q`, and check
`outputs/<production-tag>/condor/calib-study.dag.dagman.out` if a node fails.
If a map job fails, DAGMan writes a rescue DAG next to the `.dag` file;
re-running `./calib-study-condor` with the same tag starts over from scratch
instead of resuming, since it clears the partials.

To fall back to the old behavior of one job that processes the whole list
sequentially, add `--single-job`:

```bash
./calib-study-condor <input-address-txt> <production-tag> [tree-name] [--Vcb] --single-job
```

`<input-address-txt>` may contain multiple `.root` files. The files are treated
as one dataset in either mode, so each study produces one combined set of plots.

You can edit the CPU, memory and disk usages [at the top of the
script](calib-study-condor); the map and merge jobs are requested separately.

### Purity at the 0.5 cut

Every `PosvsNeg` plot (studies 6 and 9) cuts the score at 0.5 and quotes, for the
region above the cut:

```
purity        = M / (M + N)
contamination = N / (M + N)
```

where **M is the red `b~ l- v~` yield** and **N is the blue `b l+ v` yield**, both
weighted by `finalWeight`. The line is drawn under the x-axis, and the cut itself
is marked with a dashed vertical line.

0.5 is exactly a bin edge of the 100-bin `[0, 1]` score axis, so the cut splits
the histogram cleanly rather than slicing a bin. The overflow bin is included so
that a score of exactly 1 is not dropped.

This is computed from the merged histograms at plotting time, so re-quoting it at
a different working point costs one `--merge` over the existing partials — the
input `.root` files are never reopened.

### Event counts on the PNetB plots

Every `PNetB` plot (studies 7, 8, 10 and 11) is normalized, so the absolute scale
is printed under the axis instead:

```
entries = 157077,   weighted yield = 317.59   (16.07 outside the axis range, not shown)
```

- `entries` is the raw fill count — every event that reached the plot, unweighted.
  Use it to judge the statistical uncertainty of the shape.
- `weighted yield` is the `finalWeight`-weighted sum over the drawn bins. This is
  the quantity the normalization divides by.
- The parenthetical appears only when some events land outside `[0, 1]`. In the
  Summer24 sample these are all **underflow** — jets carrying a sentinel PNetB
  score rather than a real one. They are filled, so `entries` counts them, but
  they are neither drawn nor included in the normalization, which is why they are
  reported separately rather than folded into the yield.

Note that `entries` therefore counts slightly more events than the drawn
distribution contains. The partial files store weighted bin contents, so an exact
unweighted count of only the in-range events is not recoverable from them; it
would need a re-run of the map jobs.

### Studies in `./calib-study`

One should create a C++ implementation of these studies and compile them so they
can be used through `./calib-study`.

Notes:
* If ROOT software is needed, use the CMS environment provided in [cmsenv-source.txt](./cmsenv-source.txt).
* Title each plot using the default format `<branch_name> in <production_tag>, post processor`, unless specified below.

1. Write the branch schema for the requested `TTree` to `fields.txt`. 
    The report should list every branch name and its ROOT/C++ data type.

2. Write a dump of entry index `0` from the requested `TTree` to
   `0th-event.txt`. The report should include branch names and the stored values
   for that entry.

3. Perform a cutflow, and print out the efficiencies of each step both with 
    respect to the last step, and with respect to the initial input.
    The efficiencies should be calculated using the branch `finalWeight` 
    if available; if not, weight each MC event as one.

    - Event must have at least four jets, each with `PT >= 20 GeV` and `abs(eta) < 2.5`
    - Each gen-level quark from the semileptonic top decay must have a matched jet,
        i.e. `ak4MatchedHadB_*`, `ak4MatchedLepB_*`, `ak4MatchedHadQ1_*`, `ak4MatchedHadQ2_*`
        should each have a `1` ("matched") entry.
    - The four jets matched to gen-level quarks, obtained in the step above, should be distinct.

    (For the physicist: Calib_ChargeTagger_JH already selects a trigger lepton 
    per event that passes through the activation threshold and abs(eta) reqs.)

    Output the results of the cutflow in a .txt file.

The following studies should be performed on events that passed the cutflow above.

4. Plot out the invariant mass of the two jets matched to `ak4MatchedHadQ1` 
    and `ak4MatchedHadQ2` in a histogram binned from 0 to 160 GeV, bin width 
    10 GeV. Weight each event by `finalWeight` if available.

5. Plot out the invariant mass of the three jets matched to `ak4MatchedHadQ1`, 
    `ak4MatchedHadQ2`, and `ak4MatchedHadB` in a histogram binned from 80 to 260 GeV, 
    bin width 10 GeV. Weight each event by `finalWeight` if available.

6. Plot out the `ak4JetParTPosvsNeg` score of the jet matched to `ak4MatchedLepB`
    in six plots. Each plot contains a PT region in 
    `[20, 30), [30, 50), [50, 80), [80, 110), [110, 150), [>150]`.

    Inside each plot, separate the events into two sets according to the trigger
    lepton's sign. If the trigger lepton is positive, call the set "b l+ v" and 
    fill it semitransparent blue, and if the trigger lepton is negative, call 
    the set "b~ l- v~" and fill it semitransparent red. Overlay the two 
    histograms in the same plot without stacking. Bin from 0 to 1 with width 0.01.

    Title each plot as "ParTPosvsNeg score of LepB jet, <PT bin>". Output the 
    plots to a folder named "LepB_jet_PosvsNeg" under the output subdirectory.

    Mark the 0.5 cut with a dashed line, and print the purity and contamination
    of the region above it under the axis. See
    [Purity at the 0.5 cut](#purity-at-the-05-cut) below.
    
7. Plot out the `ak4JetbtagPNetB` score of the jet matched to `ak4MatchedLepB`
    in six plots. Each plot contains a PT region in 
    `[20, 30), [30, 50), [50, 80), [80, 110), [110, 150), [>150]`.

    Inside each plot, bin according to the following bins 
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the 
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of LepB jet, <PT bin>". Output the 
    plots to a folder named "LepB_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield under the axis. See
    [Event counts on the PNetB plots](#event-counts-on-the-pnetb-plots).

8. Plot out the `ak4JetbtagPNetB` score of the jet matched to `ak4MatchedHadB`
    in six plots. Each plot contains a PT region in 
    `[20, 30), [30, 50), [50, 80), [80, 110), [110, 150), [>150]`.

    Inside each plot, bin according to the following bins 
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the 
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of HadB jet, <PT bin>". Output the 
    plots to a folder named "HadB_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield as in study 7.

The following studies should be performed on events that passed the cutflow
above, along with a `--Vcb` tag in the `./calib-study` input.

9. Plot out the `ak4JetParTPosvsNeg` score of the jet matched to the charm jet
    (check if `GenQ1PdgId`, `GenQ2PdgId` is 4 or -4) in `ak4MatchedHadQ1` or `ak4MatchedHadQ2`.

    Plot in six plots. Each plot contains a PT region in
    `[20, 30), [30, 50), [50, 80), [80, 110), [110, 150), [>150]`.

    Inside each plot, separate the events into two sets according to the trigger
    lepton's sign. If the trigger lepton is positive, call the set "b l+ v" and
    fill it semitransparent blue, and if the trigger lepton is negative, call
    the set "b~ l- v~" and fill it semitransparent red. Overlay the two
    histograms in the same plot without stacking. Bin from 0 to 1 with width 0.01.

    Title each plot as "ParTPosvsNeg score of C jet, <PT bin>". Output the
    plots to a folder named "C_jet_PosvsNeg" under the output subdirectory.

    Mark the 0.5 cut and print the purity and contamination as in study 6.

10. Plot out the `ak4JetbtagPNetB` score of the jet matched to the charm jet
    (check if `GenQ1PdgId`, `GenQ2PdgId` is 4 or -4) in `ak4MatchedHadQ1` or `ak4MatchedHadQ2`.

    Inside each plot, bin according to the following bins
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of C jet, <PT bin>". Output the
    plots to a folder named "C_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield as in study 7.

11. Plot out the `ak4JetbtagPNetB` score of the jet matched to the bottom jet
    (check if `GenQ1PdgId`, `GenQ2PdgId` is 5 or -5) in `ak4MatchedHadQ1` or `ak4MatchedHadQ2`.

    Inside each plot, bin according to the following bins
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of WB jet, <PT bin>". Output the
    plots to a folder named "WB_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield as in study 7.

The following study should be performed with the `--Vcb` tag but on events of 
at least four jets, each with `PT >= 20 GeV` and `abs(eta) < 2.5`. No extra 
cutflow is needed.

12. Select events with only 4 jets, and three of them should have a `PNetB` 
    score >= 0.7862 (XT). Plot out the DeltaR between the fourth, unselected jet
    and the gen-level charm quark (PID = +-4).

    Bin from 0 to 3.2 with bin width 0.1. Title 
    `DeltaR between c quark and c-jet-candidate in 4 jet events`.

    In the legends, show `Selected events: <histogram-count>`. Weight histogram
    by `finalWeight`.

