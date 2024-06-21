# Fat tissue analysis for cfRNA data

#
## Background:

Although human beings seem very different from the outside, when it comes to our biological info, aka DNA - we are very similar to each other, and mostly we express the same genes respectivally to each cell type.
With that being said, external factor like our life-style, the food we eat and such, can have an impact of the expression levels of some genes. 

There is a cool noninvasive method called [cfRNA](https://www.illumina.com/science/genomics-research/articles/circulating-rna-sequencing-enables-noninvasive-monitoring-dynamic-changes-in-human-health.html), with which we can research the RNA content in the blood.
Meaning that after analysing the results of the RNA seq, we might be able to find interesting differences between blood samples of different patients.

The processed in simple terms goes like this:

1. Cells die due to [Apoptosis / necrosis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1087413/) and spill their content to the blood.
2. Another option - cells secrete RNA via [extra-cellular vesicle](https://www.nature.com/articles/s41577-022-00763-8).
3. Either way - there are RNA molecules circulating freely in the blood stream.
4. Blood sample can be harvest, processed, and RNA material is sent to sequencing.

![Workflow](https://ars.els-cdn.com/content/image/1-s2.0-S1471491421000022-gr1_lrg.jpg)

[More information can be found here](https://www.sciencedirect.com/science/article/pii/S1471491421000022?via%3Dihub#f0005)

#
## Project goal:

This project is a part of a bigger project that looks for differences in the gene expression of human blood samples, that are processed using the cfRNA method.

Spesifically here, I want to create an analysis based on the [GTEx database](https://gtexportal.org/home/) - that contains information of RNA seq data. The RNA was harvest from biopsies of most tissues in the body.

Here I'm interested in finding the marker genes of the fat tissues - both Adipose Subcutaneous and Adipose Visceral.
The markers I find will be used for analysing the RNA seq data that would be received from our patients blood samples that went through the cfRNA processing.

This way, I might be able to find interesting phenomenas in regards to the fat tissues - perhaps the fat marker genes are represented more in people with obisity, or represented less in those who suffer from anorexia.

#
## Technical info:

### The analysis is planned to work this way:

1. Downloading the GTEx database.
2. Importing relevant python libraries and modules.
3. Uploading the GTEx database and reading it as pandas table.
4. Filtering the genes - deleting mitochondrial genes and ribo genes.
5. Adjusting the table - setting gene names as indexes, dropping the ENSG names.
6. Creating 3 sub tables, each contains 2 out of 3 fat tissues - visceral, subcutaneous and brest.
7. Defining a function for finding marker genes.
8. Applying the function for the 3 sub tables, to receive the marker genes.
9. Combining the marker genes from visceral and subcutaneous adipose tissues.
10. Finding the intersection between the union list, and the marker genes in breast tissue.
10. Deleting from both the sub and vis lists the marker genes that appear in the breast tissue.
11. Extracting the marker genes for vis and sub to a csv file.

   
### Installing the dependencies:
Ensure you have the following installed:
- Python 3.10 or later
- Jupyter Notebook
- pandas
- pytest
- nbimporter

The needed installations can be found in `requirements.txt`
```
pip install -r requirements.txt
```


### Running the program:
Clone this repository (https://github.com/Avitalrosner/Fat_tissue_analysis.git), open the  jupyter notebook provided - GTEx.ipynb,
and press on Run all as shown here:
 ![Run all](./run_code.JPG)


### Testing the program:
Run-
```
pytest
```
The output should look like this:
![Output of test](./output_of_test.JPG)


The project was originally implemented as part of the [Python programming course](https://github.com/szabgab/wis-python-course-2024-04) at the [Weizmann Institute of Science](https://www.weizmann.ac.il/) taught by [Gabor Szabo](https://szabgab.com/).


