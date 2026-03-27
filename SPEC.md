# Bionformatics tool wrappers

Write a set of wrapper scripts for the following online tools:

* mCSM
* mCSM-LIG
* Dynamut

## Submission stage

Each script should read from a source CSV in the `input/` directory containing the following fields:


mutation_id: Mutation ID
gene: Gene
pdb_id: PDB ID
chain: PDB Chain
residue_number: PDB Residue Number
wt_aa: Wild Type Amino Acid in PDB
mut_aa: Mutant Amino Acid (Not in PDB)

An example of this CSV is:

```
mutation_id,gene,pdb_id,chain,residue_number,wt_aa,mut_aa
A R282W,TP53,2OCJ,A,282,R,W
A Y236F,TP53,2OCJ,A,236,Y,F
A N239Y,TP53,2OCJ,A,239,N,Y
A C242S,TP53,2OCJ,A,242,C,S
```

The URLs and submission guidelines for each of these tools are:

* mCSM
URL: `https://biosig.lab.uq.edu.au/mcsm/stability`
Notes: Accepts multiple mutations. Use the "Mutation list" submission method. Read the submission page thoroughly in order to determine the correct format for submission.

* mCSM-LIG
URL: `https://biosig.lab.uq.edu.au/mcsm_lig/prediction`
Notes: Only accepts a single mutation at a time, with Chain, 3-Letter ligand ID, and Wild-Type affinity. Read the submission page thoroughly in order to determine the correct format for submission.

* Dynamut
URL: `https://biosig.lab.uq.edu.au/dynamut/prediction`
Notes: Multiple mutations. Use the "Mutation List" submission method. Read the submission page thoroughly in order to determine the correct format for submission.

After submitting the jobs to each tool, retrieve the Job ID and Job URL for each entry, and store these in a tool-specific CSV in the `output/` directory.

## Retrieval Stage

Using the set of tool-specific CSVs in the `output/` directory, retrieve (via web scraping where necessary) the content from each result URL and store in the `results/` directory, taking care not to overwrite anything.

## Formatting Stage

Using the results from the the `results/` directory, convert each of the downloaded files to CSV.
