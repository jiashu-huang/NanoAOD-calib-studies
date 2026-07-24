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
