#!/usr/bin/env python3

"""
Author: Rebecca Locke
Date: 03/05/24

Script to parse the output from Sambamba and return genes with any exons
which are below 30x coverage.

Usage: python find_low_coverage_genes.py -i <sambamba_output.txt>
"""

import argparse
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments given

    Returns
    -------
    argparse.Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        description='Information required to find low coverage genes'
    )

    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Input which is .txt output of Sambamba'
    )

    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help=(
            'Custom name for output report (otherwise will name according to '
            'input file)'
        )
    )

    args = parser.parse_args()

    return args


def read_in_sambamba_file(input_file) -> pd.DataFrame:
    """
    Read in the Sambamba output file to a pandas df and split
    the GeneSymbol;Accession column into two separate columns

    Parameters
    ----------
    input_file : str
        Sambamba per-exon coverage file

    Returns
    -------
    coverage_df : pd.DataFrame
        _description_
    """
    coverage_df = pd.read_csv(input_file, delim_whitespace=True)

    coverage_df[['GeneSymbol', 'Accession']] = coverage_df[
        'GeneSymbol;Accession'
    ].str.split(';', expand=True)

    return coverage_df


def main():
    args = parse_args()
    sambamba_df = read_in_sambamba_file(args.input)


if __name__ == '__main__':
    main()
