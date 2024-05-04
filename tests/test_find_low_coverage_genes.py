import os
import pandas as pd
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
import find_low_coverage_genes as flcg


class TestReadInFile():
    """
    Functions for testing reading in the Sambamba file
    """
    bad_format = os.path.join(
        TEST_DATA_DIR,
        "NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output_bad_format.txt"
    )

    def test_read_in_sambamba_file(self):
        """
        Check that error raised as expected when Sambamba file doesn't have
        columns we expect
        """
        expected_msg = (
            "Expecting columns #chromosome, StartPosition, EndPosition, "
            "GeneSymbol;Accession, percentage30. Please check input file"
        )

        with pytest.raises(KeyError, match=expected_msg):
            flcg.read_in_sambamba_file(self.bad_format)


class TestFilterLowCoverage():
    """
    Functions for testing getting low coverage genes and thresholds used
    """
    # Make example df to test with
    example_df = pd.DataFrame.from_dict({
        'GeneSymbol': ['GENE1', 'GENE2', 'GENE3', 'GENE4'],
        'Accession': ['NM_12345.1', 'NM12246.1', 'NM17596.1', 'NM2345.1'],
        'AvgGeneCoverage30x': [100.00, 98.64526096337, 92.1447575, 96.4475757]
    })

    def test_get_genes_below_30x_threshold_at_100pc(self):
        """
        Check the correct 3 genes are returned which are under 100% coverage
        """
        below_100pc_at_30x = flcg.get_genes_below_30x_threshold(
            self.example_df, 100
        )
        below_100pc_at_30x.equals(
            pd.DataFrame.from_dict({
                'GeneSymbol': ['GENE2', 'GENE3', 'GENE4'],
                'Accession': ['NM12246.1', 'NM17596.1', 'NM2345.1'],
                'AvgGeneCoverage30x': [98.65, 92.14, 96.44]
            })
        )

    def test_get_genes_below_30x_threshold_at_98pc(self):
        """
        Check the correct 2 genes are returned which are under 98% coverage
        """
        below_98pc_at_30x = flcg.get_genes_below_30x_threshold(
            self.example_df, 98
        )
        below_98pc_at_30x.equals(
            pd.DataFrame.from_dict({
                'GeneSymbol': ['GENE3', 'GENE4'],
                'Accession': ['NM17596.1', 'NM2345.1'],
                'AvgGeneCoverage30x': [92.14, 96.44]
            })
        )

    def test_get_genes_below_30x_threshold_at_96pc(self):
        """
        Check only one gene left which is under 96% at 30x coverage
        """
        below_96pc_at_30x = flcg.get_genes_below_30x_threshold(
            self.example_df, 96
        )
        below_96pc_at_30x.equals(
            pd.DataFrame.from_dict({
                'GeneSymbol': ['GENE3'],
                'Accession': ['NM17596.1'],
                'AvgGeneCoverage30x': [92.14]
            })
        )


    def test_no_genes_under_30x_threshold(self):
        """
        Check df empty when no genes are under the threshold given are found
        """
        below_92pc_at_30x = flcg.get_genes_below_30x_threshold(
            self.example_df, 92
        )
        assert below_92pc_at_30x.empty
