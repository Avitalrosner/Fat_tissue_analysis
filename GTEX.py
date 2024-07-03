#!/usr/bin/env python
# coding: utf-8

# In[35]:


## Finding MARKER GENES from GTEx


# In[36]:


# Import modules

import pandas as pd
import os


# In[37]:


# Getting the dataset to analyze

GTEx_file = os.path.join(os.path.dirname(__file__), '..', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')


# In[38]:


# Converting to DataFrame

df = pd.read_csv(GTEx_file, sep='\t', skiprows=2)
df


# In[39]:


# Deleting the mitochondrial genes from the dataset

df['Description'] = df['Description'].astype(str)
df_without_mt = df[~df['Description'].str.startswith('MT-')]
df_without_mt


# In[40]:


#Delete all genes that are considered ribo 

delete_genes = ["FAU", "MRPL13", "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14",
    "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23", "RPL23A",
    "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL3", "RPL30", "RPL31",
    "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL", "RPL37", "RPL37A", "RPL38",
    "RPL39", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPLP0",
    "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS15", "RPS15A", "RPS16", "RPS17",
    "RPS18", "RPS19", "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27",
    "RPS27A", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5", "RPS6",
    "RPS7", "RPS8", "RPS9", "RPSA", "RSL24D1", "RSL24D1P11", "UBA52"]

df_no_mt_ribo = df_without_mt[~df_without_mt.iloc[:, 1].isin(delete_genes)]
df_no_mt_ribo


# In[41]:


# Checking if the change worked - expecting *not* to find the gene

gene_name = "RPL22"
is_in_description = gene_name in df_no_mt_ribo.iloc[:, 1].values
is_in_description

# Indeed got "False"


# In[42]:


# Organizing the dataframe to the wanted index

df_no_mt_ribo.set_index("Description", inplace=True)
df_no_mt_ribo


# In[43]:


# Dropping unnecessary coulombs

df_filtered = df_no_mt_ribo.drop("Name", axis = 1)
df_filtered


# In[44]:


# Making 3 sub-dataframes for each fat tissue.
# Each sub-dataframe doesn't have the other 2 fat tissues. 

df_only_adipose_subcutaneous = df_filtered.drop(['Adipose - Visceral (Omentum)', 'Breast - Mammary Tissue'], axis = 1)
df_only_adipose_visceral = df_filtered.drop(['Adipose - Subcutaneous', 'Breast - Mammary Tissue'], axis = 1)
df_only_breast = df_filtered.drop(["Adipose - Visceral (Omentum)", 'Adipose - Subcutaneous'], axis = 1)


# In[45]:


# Creating a function to get marker genes for each sample in a dictionary

def get_marker_genes(df, FC_THRESH=2, MAX_EXP_THRESH=1e-4, PN=1e-9):

    df = df / df.sum(axis=0)
    marker_genes_dict = dict()

    # Loop through each tissue/sample
    for ii, sample in enumerate(df.columns):
        sample_expression = df[sample]
        other_samples_max_expression = df[df.columns[df.columns != sample]].max(axis=1)
        all_samples_max_exp = df.max(axis=1)

        # Calculate the ratio of expression in the current sample vs other samples
        tissue_ratio = (sample_expression[all_samples_max_exp > MAX_EXP_THRESH] + PN) / (other_samples_max_expression[all_samples_max_exp > MAX_EXP_THRESH] + PN)
        marker_genes_dict[sample] = tissue_ratio[tissue_ratio > FC_THRESH]

    return marker_genes_dict


# In[46]:


# To see all marker genes for all tissues:

#get_marker_genes(df_filtered)


# In[47]:


# Applying the marker gene function to the 3 sub-dataframes

only_adipose_visceral = get_marker_genes(df_only_adipose_visceral)
only_adipose_subcutaneous = get_marker_genes(df_only_adipose_subcutaneous)
only_breast = get_marker_genes(df_only_breast)


# In[48]:


#Making a list for each fat tissue, with gene_name and exp_value

vis = list(only_adipose_visceral.values())
sub = list(only_adipose_subcutaneous.values())
bre = list(only_breast.values())


# In[49]:


# Converting to lists to have only gene names 

lst_vis = []
df_split_vis = pd.DataFrame(vis)
for gene in df_split_vis:
    lst_vis.append(gene)

lst_sub = []
df_split_sub = pd.DataFrame(sub)
for gene in df_split_sub:
    lst_sub.append(gene)

lst_bre = []
df_split_bre = pd.DataFrame(bre)
for gene in df_split_bre:
    lst_bre.append(gene)


# In[50]:


#Combine visceral (vis) and subcutaneous(sub), and delete genes from breast(bre)

union_vis_sub = lst_sub + lst_vis         
union_vis_sub = list(set(union_vis_sub)) 
sorted_lst = sorted(union_vis_sub)
intersect_bre_union = list(set(union_vis_sub).intersection(set(lst_bre)))
vis_genes_not_in_intersection = [gene for gene in lst_vis if gene not in intersect_bre_union]
sub_genes_not_in_intersection = [gene for gene in lst_sub if gene not in intersect_bre_union]


# In[51]:


# Saving the final vis and sub lists, that don't have breast

subdirectory = "C:/python/Fat_tissue_analysis"
if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)
file_path_vis = os.path.join(subdirectory, "vis.csv")
df_vis = pd.DataFrame(vis_genes_not_in_intersection, columns=["Gene"])
df_vis.to_csv(file_path_vis, index=False)

file_path_sub = os.path.join(subdirectory, "sub.csv")
df_sub = pd.DataFrame(sub_genes_not_in_intersection, columns=["Gene"])
df_sub.to_csv(file_path_sub, index=False)

# The output files are in the repository.

