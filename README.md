# fruitfly
Utilities and CLI for molbloom filters.

THIS IS A WORK IN PROGRESS. Several changes are expected.

**Note:** To date this is just a convenient script with a CLI that wraps
around `molbloom` to help finding matches of compounds in `.bloom` files,
given an input SMILES (`.smi`) file.

## Installing dependencies
The recommended way is to use conda/mamba (or similar) to create an environment
with the required tools. 

After you clone this repository, you can create the conda environment using 
the following command (note that we recommend `mamba` but you can change the
command to use `conda` instead of `mamba`, if desired):

```bash
mamba env create -f environmentl.yml -n fruitfly-env
```

## How to use

To use it you would just activate the created environment with 
`conda activate fruitfly-env`, and call the script like the following in the
same directory where you cloned this repository:

```bash
python fruitfly.py find-matches --glob-exp "../*.bloom" --smiles-file matching_examples.sm
i --out-file results.json
```

Please note that the quotes in the `glob-exp` are required for the shell to not interpret
wildcards (`*`).

Additional arguments can be specified, for more information please use 
`python fruitfly find-matches --help`.

## Results
The results are stored in a generated `json` file, which has as keys the base names
of the filters you used, and the matching smiles as values for each of the filters.

This is convenient in order to know from which bloom filter the matches were obtained.

Additionally, if you specified the `--include-unmatched True` argument, a key `unmatched`
will be generated with the unmatched SMILES as the values.


## Future work
In the near future, commands to easily create your own filters from `.smi` files will
be added. Improving the CLI and overall modularity of the project is also planned.
