#!/usr/bin/env python3

"""
Author: Rebecca Locke
Date: 03/05/24

Script to parse the output from Sambamba and return genes which are below 30x
coverage.

Usage:
python find_low_coverage_genes.py [-h]
-i <sambamba_file.txt> \
[-o <name_of_output_file.csv>] \
[-t <100>]
"""

import argparse
import numpy as np
import pandas as pd
import sys

from pathlib import Path


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
        help='.txt file input which is output of Sambamba'
    )

    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='Name for output report of genes below 30x threshold (optional)'
    )

    parser.add_argument(
        '-t',
        '--threshold',
        type=int,
        help='Get genes < this percentage threshold at 30x (default 100pc)',
        default=100
    )

    args = parser.parse_args()

    return args


def read_in_sambamba_file(input_file) -> pd.DataFrame:
    """
    Read in the Sambamba output file to a pandas df and split the
    GeneSymbol;Accession column into two separate columns so we can do
    checks and group by gene / transcript

    Parameters
    ----------
    input_file : str
        Sambamba per-exon coverage .txt file

    Returns
    -------
    sambamba_subset : pd.DataFrame
        pandas df with subset of Sambamba info, with columns split out

    Example output format:
    +-------------+---------------+-------------+------------+-------------+--------------+
    | #chromosome | StartPosition | EndPosition | GeneSymbol | Accession   | percentage30 |
    +-------------+---------------+-------------+------------+-------------+--------------+
    | 1           | 26126711      | 26126914    | SELENON    | NM_020451.2 | 49.2611      |
    | 1           | 26127523      | 26127661    | SELENON    | NM_020451.2 | 100.0000     |
    +-------------+---------------+-------------+------------+-------------+--------------+
    """
    # Read in .txt file, accounting for whitespace
    sambamba_df = pd.read_csv(input_file, delim_whitespace=True)

    # Remove any full duplicate rows
    sambamba_df.drop_duplicates(keep='first', inplace=True)

    # Subset to useful columns
    try:
        sambamba_subset = sambamba_df[[
            '#chromosome', 'StartPosition', 'EndPosition', 'GeneSymbol;Accession', 'percentage30'
        ]]
    except KeyError:
        print(
            "Expecting columns #chromosome, StartPosition, EndPosition, "
            "GeneSymbol;Accession, percentage30. Please check input file "
            "has these columns"
        )
        sys.exit(1)

    # Split out GeneSymbol;Accession column into two columns
    sambamba_subset[['GeneSymbol', 'Accession']] = sambamba_subset[
        'GeneSymbol;Accession'
    ].str.split(';', expand=True)

    return sambamba_subset


def calculate_gene_30x_coverage(sambamba_df) -> pd.DataFrame:
    """
    Calculate the coverage of each gene at 30x based on exon 30x coverage

    Parameters
    ----------
    sambamba_df : pd.DataFrame
        pandas df with subset of Sambamba info, with columns split out

    Returns
    -------
    gene_cov_30x_df: pd.DataFrame
        pandas df with each gene + transcript and the coverage at 30x

    Example output format:
    +------------+----------------+--------------------+
    | GeneSymbol | Accession      | AvgGeneCoverage30x |
    +------------+----------------+--------------------+
    | ACTA1      | NM_001100.3    | 100.000000         |
    | B3GALNT2   | NM_152490.4    | 100.000000         |
    | MSTO1      | NM_018116.3    | 98.6452609633718   |
    | NEB        | NM_001271208.1 | 90.2228444732621   |
    +------------+----------------+--------------------+
    """
    # Add column representing length of each exon
    sambamba_df['ExonLength'] = (
        sambamba_df['EndPosition'] - sambamba_df['StartPosition']
    )

    # Calculate fraction of gene each exon covers wrt all exon lengths summed
    # multiply each exon's gene fraction by the percentage30 column
    # then average this over the gene
    gene_cov_30x_df = sambamba_df.groupby(['GeneSymbol', 'Accession']).apply(
        lambda row: np.average(row['percentage30'], weights=row['ExonLength'])
    ).to_frame('AvgGeneCoverage30x').reset_index()

    return gene_cov_30x_df


def main():
    """
    Main function to generate gene 30x coverage report output at threshold
    given
    """
    args = parse_args()
    sambamba_df = read_in_sambamba_file(args.input)
    gene_coverage_30x = calculate_gene_30x_coverage(sambamba_df)

if __name__ == '__main__':
    main()
