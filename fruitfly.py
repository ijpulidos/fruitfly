import glob
import bz2
import pickle
import random
from functools import partial
from itertools import chain
import multiprocessing as mp

import click

# Helper functions
import re
def parse_humanized_string(humanized_string):
    """Parse humanized strings like '128M' to '16K' to a machine-readable numeric value."""
    humanized_string = humanized_string.strip().upper()
    multipliers = {'K': 10**3, 'M': 10**6, 'B': 10**9}

    for suffix, multiplier in multipliers.items():
        if humanized_string.upper().endswith(suffix):
            numeric_value = float(humanized_string[:-1]) * multiplier
            return int(numeric_value)

    # If no valid suffix is found, assume it's already a numeric value
    try:
        return int(humanized_string)
    except ValueError:
        raise ValueError("Invalid humanized string format")

def detect_humanized_number(some_string):
    """Detects numeric humanized expression in some string"""
    regex = r"[-+]?(?:\d*\.*\d+)[kKmMbB]"
    results = re.search(regex, some_string)
    try:
        first_result = results.group(0)
    except AttributeError:
        raise ValueError("No numeric humanized expression found in string.")
    return first_result

@click.command()
@click.option("--file", help="Path to bz2 file to be processed.")
@click.option("--num-lines", help="Number of random lines to select from file.")
def sample_bz2(file_path, num_lines):
    """
    Reads a file in bz2 compressed format, selects random lines, and returns them.

    Useful for extracting smiles from Enamine REAL database files.

    Parameters
    ----------
    file_path : str
        The path to the bz2 compressed file to be processed.
    num_lines : int
        The number of random lines to select from the file.

    Returns
    -------
    List[str]
        A list containing randomly selected lines from the file.
    """
    total_lines = parse_humanized_string(detect_humanized_number(file_path))
    print(f"Processing file: {file_path}")
    with bz2.BZ2File(file_path, "r") as afile:
        # lines = afile.readlines()

        chosen_lines = []
        # chosen_lines.append(next(afile))  # Always append first line
        num = int(total_lines/num_lines)
        for aline in afile:
            if random.randrange(num):
                continue
            chosen_lines.append(aline)

    print(f"Finished processing file: {file_path}")
    return chosen_lines


@click.group()
def cli():
    pass


def bloom_match(bloom_file: str, smiles_iterable):
    """
    Filters SMILES strings using a Bloom filter.

    Parameters
    ----------
    bloom_file : str
        The path to the Bloom filter file.
    smiles_iterable : Iterable[str]
        An iterable containing SMILES strings to be filtered.

    Returns
    -------
    Dict[str, List[str]]
        A dictionary containing filtered SMILES strings.
        The key is the base name of the Bloom filter file, and the value is a list of SMILES strings that passed
        through the filter.
    """
    from pathlib import Path
    from molbloom import BloomFilter
    bf = BloomFilter(bloom_file)
    results = []
    for smile in smiles_iterable:
        if smile in bf:
            results.append(smile)

    base_name = Path(bloom_file).stem
    data = {base_name: results}

    return data

@click.command()
@click.option("--glob-exp", default="*.bloom", help="glob expression to match bloom filter files. Ex: 'filters/*.bloom' .")
@click.option("--num-procs", default=-1, help="Number of subprocesses to use. Default is to choose automatically based on number of bloom filters.")
@click.option("--smiles-file", help="SMILES file in .smi format to use for matches.")
@click.option("--out-file", default="results.json", help="Output JSON file to store the results. Default: results.json.")
@click.option("--include-unmatched", default=False, help="Include unmatched SMILES in the results. Default: False.")
def find_matches(glob_exp: str, num_procs: int, smiles_file: str, out_file="results.json", include_unmatched: bool=False):
    import json
    filter_files = glob.glob(glob_exp)
    if not filter_files:
        raise ValueError(f"No bloom filters found using the glob expression: {glob_exp}")
    if num_procs < 1:
        # Match number of processes with number of filters
        n_processes = len(filter_files)
    else:
        n_processes = num_procs

    with open(smiles_file, "r") as in_file:
        smiles = []
        for line in in_file:
            smiles.append(line.strip())

    with mp.Pool(n_processes) as pool:
        partial_f = partial(bloom_match, smiles_iterable=smiles)
        results = pool.map(partial_f, filter_files, chunksize=1)

    # Merge dictionaries
    results_dict = {}
    for result in results:
        results_dict.update(result)

    # Add unmatched
    if include_unmatched:
        results_dict["unmatched"] = []
        matched = list(chain(*results_dict.values()))
        # print(matched)
        for smile in smiles:
            if smile not in matched:
                # print(f"SMILES {smile} not found.")
                # break
                results_dict["unmatched"].append(smile)

    with open(out_file, "w") as out_file:
        json.dump(results_dict, out_file, indent=4)

if __name__ == "__main__":
    cli.add_command(find_matches)
    cli()
    """
    n_processes = 9
    file_paths = glob.glob("../*cxsmiles.bz2")
    # file_paths = ["../Enamine_REAL_HAC_11_21_666M_CXSMILES.cxsmiles.bz2"]
    n_compounds_from_each = 10000  # Number of mols to extract from each db file
    with mp.Pool(n_processes) as pool:
        partial_f = partial(sample_bz2, num_lines=n_compounds_from_each)
        results = pool.map(partial_f, file_paths, chunksize=1)

    # Pickle/serialize the results
    with open("results.pickle", "wb") as out_file:
        pickle.dump(results, out_file, pickle.HIGHEST_PROTOCOL)
    """
