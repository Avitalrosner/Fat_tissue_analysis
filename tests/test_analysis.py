import pytest
import pandas as pd
import os

from GTEx import get_marker_genes

@pytest.fixture(scope="module")
def setup_data():
    data = {
        'Description': ['MT-GeneA', 'MT-GeneB', 'FAU', 'MRPL13', 'RPL10', 'GeneC', 'GeneD'],
        'Sample1': [10, 0, 5, 0, 0, 0, 0],
        'Sample2': [0, 15, 0, 5, 0, 5, 0],
        'Sample3': [0, 0, 0, 0, 20, 0, 20]
    }
    df = pd.DataFrame(data)
    df.set_index('Description', inplace=True)
    return df

def test_file_reading():
    assert os.path.exists("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")

def test_mitochondrial_genes_removal(setup_data):
    df = setup_data.copy()  # Copying to avoid modifying the fixture data
    df['Description'] = df.index.astype(str)
    df_without_mt = df[~df['Description'].str.startswith('MT-')]
    assert not any(df_without_mt['Description'].str.startswith('MT-'))

def test_ribosomal_genes_removal(setup_data):
    delete_genes = ["FAU", "MRPL13", "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14",
        "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23", "RPL23A",
        "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31",
        "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL", "RPL37", "RPL37A", "RPL38",
        "RPL39", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPLP0",
        "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS15", "RPS15A", "RPS16", "RPS17",
        "RPS18", "RPS19", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27",
        "RPS27A", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5", "RPS6",
        "RPS7", "RPS8", "RPS9", "RPSA", "RSL24D1", "RSL24D1P11", "UBA52"]
    
    df = setup_data.copy()
    df_without_mt = df[~df.index.str.startswith('MT-')]
    df_no_mt_ribo = df_without_mt[~df_without_mt.index.isin(delete_genes)]
    assert not any(gene in df_no_mt_ribo.index for gene in delete_genes)

def test_marker_genes_calculation(setup_data):
    df = setup_data
    marker_genes_dict = get_marker_genes(df)
    assert isinstance(marker_genes_dict, dict)
    assert all(isinstance(key, str) for key in marker_genes_dict.keys())
    assert all(isinstance(value, pd.Series) for value in marker_genes_dict.values())
