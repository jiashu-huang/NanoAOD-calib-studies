# NanoAOD Calib_ChargeTagger studies

Updated: 2026-09-16

This repo studies `.root` outputs from [Calib_ChargeTagger_JH](https://github.com/jiashu-huang/Calib_ChargeTagger_JH).

It produces a cutflow and nominal physics plots. It does not validate the full
processor or evaluate systematic variations. See the [audit](docs/audit-20260909.md)
for the checks, selection details, and remaining limits.

## Set up 

Make sure you initiate `cmsenv` before you run!
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/jhuan166/Vcb/CMSSW_15_1_0_patch4/src
cmsenv
```

To recompile the C++:
```bash
make clean
make
# Small ROOT fixtures; no production data needed:
python3 tests/test_calib_study.py
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

`histograms.root` retains the unscaled histograms, including underflow, overflow,
and squared weight sums, before the plotting code normalizes the PNetB shapes.
Every file within one input chain must have the same branch names and types,
including consistent `finalWeight` availability. Duplicate inputs, unreadable
files/events, non-finite values read by a study, and invalid or duplicate jet
assignments stop the run with a nonzero exit code.

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

Pass the same `--Vcb` setting to the map and merge steps. A mismatch is an error.
Mixing weighted and unweighted partials or different histogram axes also fails.
Partials made before the September 2026 audit must be regenerated: the score=1
handling changed, and new partials carry a `meta_format_v2` marker.
For `--Vcb`, partials must also contain the UpType histograms split by W daughter
flag. Regenerate older, unsplit partials with the updated executable before
merging; their combined histograms cannot recover the decay components.

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

The reported cut is `score >= 0.5`; 0.5 is exactly a bin edge. Scores equal to 1
are placed in the last visible bin in both ParT and PNetB plots. Sentinel and
out-of-range scores stay in underflow/overflow and are excluded from the purity.
The quoted fractions use signed MC weights, not unweighted counts.

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

### The up-type W prong

Studies 9 and 10 score **the up-type quark of the hadronic W** — the `u` or `c`
of the `W -> q qbar'` pair. That prong is named by a single branch: the skimmer
stores the two W daughters by *type*, not by generator child order, so `GenQ1`
is always the down-type quark (`|pdgId|` 1, 3, 5 — odd) and `GenQ2` always the
up-type one (2, 4 — even). Charge conservation gives a W exactly one of each, so
the split is total. A W that somehow yielded two same-type quarks would leave
the missing side at `PAD_VAL` rather than slot the wrong quark into `GenQ2`.
This is enforced in `_split_w_quarks_by_type` in `Calib_ChargeTagger_JH`.

So the jet is simply the one at `GenHadQ2JetIdx`. No per-event `|pdgId| == 4`
check is needed, and there is no "both quarks match / neither matches" case to
fall back on.

**This is not a change of selection for `TTtoLNuCB`.** In `W -> cb` the up-type
prong *is* the charm, in every event — verified over all 5 648 574 events of
`prod_20260816`: `|GenQ1PdgId|` is 5 and `|GenQ2PdgId|` is 4, 100% of the time.
The plots that used to be labelled "C jet" are filled from exactly the same jets
as before; only the name changed.

What it buys is `TTtoLNu2Q`, where the hadronic W goes to `ud`, `us`, `cd` or
`cs` and only **half** the events contain a charm at all. Selecting on
`|pdgId| == 4` there would silently drop that half; selecting the up-type prong
keeps every event and plots the `u` and the `c` prong together, which is the
quantity a charge tagger actually sees.

Because the histogram names changed with the labels (`c_part_posvsneg_*` ->
`uptype_part_posvsneg_*`, `pnetb_C_*` -> `pnetb_UpType_*`), partial files
written by an older build cannot be `--merge`d by this one. Re-run the map jobs
rather than mixing the two.

### Studies in `./calib-study`

Studies 1–11 are implemented in `CalibStudy.cc`. Study 12 remains a proposal.

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
    - Each gen-level quark from the semileptonic top decay must have a matched
      saved jet: `GenHadBJetIdx`, `GenLepBJetIdx`, `GenHadQ1JetIdx`, and
      `GenHadQ2JetIdx` must each be nonnegative and name an existing jet slot.
    - The four assignments must be distinct. This is an input consistency check
      that aborts on a violation, rather than a fourth cutflow bin.

    The 20 GeV requirement counts jets in the event; it is not reapplied to each
    truth-matched jet. A matched jet with 15 < pT < 20 GeV can enter the mass
    plots, but cannot enter any of the tagger plots' pT bins. Only saved jet slots
    are available here (currently at most ten); this is not a count of all input jets.

    (For the physicist: Calib_ChargeTagger_JH already selects a trigger lepton 
    per event that passes through the activation threshold and abs(eta) reqs.)

    Output the results of the cutflow in a .txt file.

The following studies should be performed on events that passed the cutflow above.

4. Plot out the invariant mass of the jets at `GenHadQ1JetIdx`
    and `GenHadQ2JetIdx` in a histogram binned from 0 to 160 GeV, bin width
    10 GeV. Weight each event by `finalWeight` if available.

5. Plot out the invariant mass of the jets at `GenHadQ1JetIdx`,
    `GenHadQ2JetIdx`, and `GenHadBJetIdx` in a histogram binned from 80 to 260 GeV,
    bin width 10 GeV. Weight each event by `finalWeight` if available.

6. Plot out the `ak4JetParTPosvsNeg` score of the jet at `GenLepBJetIdx`
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
    
7. Plot out the `ak4JetbtagPNetB` score of the jet at `GenLepBJetIdx`
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

8. Plot out the `ak4JetbtagPNetB` score of the jet at `GenHadBJetIdx`
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

9. Plot out the `ak4JetParTPosvsNeg` score of the jet matched to the **up-type
    quark of the hadronic W**, i.e. the jet at `GenHadQ2JetIdx`. See
    [The up-type W prong](#the-up-type-w-prong) for why that branch alone
    identifies it, and why this is the old charm-jet selection on `TTtoLNuCB`.

    Plot in six plots. Each plot contains a PT region in
    `[20, 30), [30, 50), [50, 80), [80, 110), [110, 150), [>150]`.

    Inside each plot, separate the events into two sets according to the trigger
    lepton's sign. If the trigger lepton is positive, call the set "b l+ v" and
    fill it semitransparent blue, and if the trigger lepton is negative, call
    the set "b~ l- v~" and fill it semitransparent red, both with alpha 0.35.
    Split each color by the saved W daughter flags and stack those components
    within that color. Overlay the blue and red stacks, each starting at zero.
    Use the same fill pattern for a decay mode in both colors: diagonal hatching
    for `WtoUD` and `WtoCD` in opposite directions, vertical hatching for `WtoUS`,
    horizontal hatching for `WtoCS`, and solid fill for `WtoCB`.
    Bin from 0 to 1 with width 0.01.

    The saved branches are `GenWtoUD`, `GenWtoCD`, `GenWtoUS`, `GenWtoCS`, and
    **`GenWtoBC`**, which is labelled **`WtoCB`** in these plots. `GenWtoUB` is
    also retained if present in the selected events. Modes with no entries
    across all six pT regions are omitted from every legend, so `TTtoLNuCB`
    shows only `WtoCB`, while `TTtoLNu2Q` shows `WtoUD`, `WtoCD`, `WtoUS`, and
    `WtoCS`. Each plotted jet must have exactly one saved W daughter flag set.
    The histograms retain signed `finalWeight` values, and each legend lists
    the weighted yield by color and decay mode. Purity uses the combined red
    and blue yields as before.

    Title each plot as "ParTPosvsNeg score of UpType jet, <PT bin>". Output the
    plots to a folder named "UpType_jet_PosvsNeg" under the output subdirectory.

    Mark the 0.5 cut and print the purity and contamination as in study 6.

10. Plot out the `ak4JetbtagPNetB` score of the jet matched to the **up-type
    quark of the hadronic W**, i.e. the jet at `GenHadQ2JetIdx`, as in study 9.

    Inside each plot, bin according to the following bins
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of UpType jet, <PT bin>". Output the
    plots to a folder named "UpType_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield as in study 7.

11. Plot out the `ak4JetbtagPNetB` score of the jet matched to the bottom jet
    (check if `GenQ1PdgId`, `GenQ2PdgId` is 5 or -5) in `GenHadQ1JetIdx` or
    `GenHadQ2JetIdx`. That is the down-type prong of the W restricted to the b,
    so its population depends on whether the sample contains a W decay to a b.
    Inspect the resulting yield rather than inferring it from the dataset label.

    Inside each plot, bin according to the following bins
    `[0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0]`.
    For `0.047, 0.245, 0.6734, 0.7862, 0.961`, show `L, M, T, XT, XXT` on the
    bin boundaries. Weight by `finalWeight` if available, and normalize.

    Title each plot as "PNetB score of WB jet, <PT bin>". Output the
    plots to a folder named "WB_jet_PNetB" under the output subdirectory.

    Print the entry count and weighted yield as in study 7.

**Not implemented:** the following proposed study would use the `--Vcb` tag on events of
at least four jets, each with `PT >= 20 GeV` and `abs(eta) < 2.5`. No extra 
cutflow is needed.

12. Select events with only 4 jets, and three of them should have a `PNetB` 
    score >= 0.7862 (XT). Plot out the DeltaR between the fourth, unselected jet
    and the gen-level charm quark (PID = +-4).

    Bin from 0 to 3.2 with bin width 0.1. Title 
    `DeltaR between c quark and c-jet-candidate in 4 jet events`.

    In the legends, show `Selected events: <histogram-count>`. Weight histogram
    by `finalWeight`.
