# NanoAOD Calib_ChargeTagger studies

Date: 2026-07-24

This repo studies `.root` outputs from [Calib_ChargeTagger_JH](https://github.com/jiashu-huang/Calib_ChargeTagger_JH).

## Design

Provide one executable command-line entry point:
```bash
./calib-study <input-address-txt> <production-tag> <OPTIONAL:tree-name>
```
Arguments:
- `input-address-txt`: a `.txt` containing the input `.root` file's address.
- `production-tag`: a label identifying the NanoAOD production under study. 
    This value is used as the output subdirectory name, so it should be treated 
    as a filesystem path component rather than free-form prose.
- `tree-name`: optional ROOT `TTree` name. Default: `Events`.

Runtime behavior:
1. Locate one `.root` file at the first line of `input_path` (omit `#` lines)
2. Open the ROOT file and retrieve the requested `TTree`.
3. In this repo, create `./outputs/<production_tag>/` if it does not already exist.
4. Run the inspection studies listed below.
5. Write all generated text reports into `./outputs/<production_tag>/`.

### Studies in `./calib-study`

One should create a C++ implementation of these studies and compile them so they
can be used through `./calib-study`.

Notes:
* If ROOT software is needed, use the CMS environment provided in [cmsenv-source.txt](./cmsenv-source.txt).
* Title each plot using the default format `<branch_name> in <production_tag>, post CMSSW_CHARGE`, unless specified below.

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

4. Plot out the invariant mass of the two jets matched to `ak4MatchedHadQ1` 
    and `ak4MatchedHadQ2` in a histogram binned from 0 to 160 GeV, bin width 
    10 GeV. Weight each event by `finalWeight` if available.

5. Plot out the invariant mass of the three jets matched to `ak4MatchedHadQ1`, 
    `ak4MatchedHadQ2`, and `ak4MatchedHadB` in a histogram binned from 80 to 260 GeV, bin width 10 GeV. Weight each event by `finalWeight` if available.
