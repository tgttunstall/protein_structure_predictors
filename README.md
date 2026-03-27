# Protein Structure Predictor Wrappers

Python wrappers for submitting protein mutation jobs to mCSM, mCSM-LIG, and DynaMut, then fetching and formatting results.

## Setup

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Project directories (created already): `input/`, `input/pdb/`, `output/`, `results/`.

## Input schema

Base columns (all tools): `mutation_id`, `gene`, `pdb_id`, `chain`, `residue_number`, `wt_aa`, `mut_aa`.

Additional columns per tool:
- mCSM-LIG: `ligand_id`, `wt_affinity` (nM)

## Usage

Run commands with `PYTHONPATH=src` (or install the package).

### mCSM
- Submit grouped by PDB: `PYTHONPATH=src python -m psp.mcsm submit --input input/mcsm_muts.csv --output output/mcsm_jobs.csv --pdb-dir input/pdb`
- Fetch results: `PYTHONPATH=src python -m psp.mcsm fetch --jobs output/mcsm_jobs.csv --results-dir results/mcsm`
- Format to CSV: `PYTHONPATH=src python -m psp.mcsm format --results-dir results/mcsm --output output/mcsm_results.csv`

### mCSM-PPI
- Submit grouped by PDB: `PYTHONPATH=src python -m psp.mcsm_ppi submit --input input/mcsm_muts.csv --output output/mcsm_ppi_jobs.csv --pdb-dir input/pdb`
- Fetch results: `PYTHONPATH=src python -m psp.mcsm_ppi fetch --jobs output/mcsm_ppi_jobs.csv --results-dir results/mcsm_ppi`
- Format to CSV: `PYTHONPATH=src python -m psp.mcsm_ppi format --results-dir results/mcsm_ppi --output output/mcsm_ppi_results.csv`

### mCSM-LIG
- Submit (one job per row): `PYTHONPATH=src python -m psp.mcsm_lig submit --input input/mcsm_lig_muts.csv --output output/mcsm_lig_jobs.csv --pdb-dir input/pdb`
- Fetch: `PYTHONPATH=src python -m psp.mcsm_lig fetch --jobs output/mcsm_lig_jobs.csv --results-dir results/mcsm_lig`
- Format: `PYTHONPATH=src python -m psp.mcsm_lig format --results-dir results/mcsm_lig --output output/mcsm_lig_results.csv`

### mCSM-NA
- Submit grouped by PDB and nucleic acid type: `PYTHONPATH=src python -m psp.mcsm_na submit --input input/mcsm_na_muts.csv --output output/mcsm_na_jobs.csv --pdb-dir input/pdb`
- Fetch: `PYTHONPATH=src python -m psp.mcsm_na fetch --jobs output/mcsm_na_jobs.csv --results-dir results/mcsm_na`
- Format: `PYTHONPATH=src python -m psp.mcsm_na format --results-dir results/mcsm_na --output output/mcsm_na_results.csv`

### DynaMut
- Submit grouped by PDB and chain: `PYTHONPATH=src python -m psp.dynamut submit --input input/dynamut_2xb7.csv --output output/dynamut_jobs.csv --pdb-dir input/pdb`
- Fetch: `PYTHONPATH=src python -m psp.dynamut fetch --jobs output/dynamut_jobs.csv --results-dir results/dynamut`
- Format: `PYTHONPATH=src python -m psp.dynamut format --results-dir results/dynamut --output output/dynamut_results.csv`

### DynaMut2
- Submit grouped by PDB: `PYTHONPATH=src python -m psp.dynamut2 submit --input input/mcsm_muts.csv --output output/dynamut2_jobs.csv`
- Fetch: `PYTHONPATH=src python -m psp.dynamut2 fetch --jobs output/dynamut2_jobs.csv --results-dir results/dynamut2`
- Format: `PYTHONPATH=src python -m psp.dynamut2 format --results-dir results/dynamut2 --output output/dynamut2_results.csv`

## Notes
- PDB files are cached in `input/pdb/`; they are fetched from RCSB if missing.
- Fetch steps skip existing result files unless `--force` is set.
- Format steps pull the first table from each HTML result; if the upstream layout changes, adjust parsers accordingly.
- For DynaMut, `mutation_id` is just the mutation (e.g., `F1174C`); the submitted mutation-list file has no header and one mutation per line.
- For mCSM-NA, include `nucleic_acid_type` (dsDNA/ssDNA/RNA) in the input CSV; jobs are grouped by PDB ID and nucleic acid type. The mutation-list sent to the server has no header and one mutation per line.
- For DynaMut2, mutation-list sent to the server has no header and one mutation per line in the form `Chain Mutation` (e.g., `A R282W`); PDB is provided as accession.
