## **ICE**
### interoperable, integrated CTD<sup>2</sup> cancer cell-line computational environment exemplar

#### *Paul Clemons, Vlado Dancik, Joshua Dempster, Eugene Douglass, Kwat Medetgul-Ernar, Becky Kuzma, Philip Montgomery, Sandrine Muller, Brittany Petros*

##### ICE is an environment to contain structured datasets in a modeling framework for computational biologists that was developed as part of a [CTD<sup>2</sup>](https://ocg.cancer.gov/programs/ctd2) Trans-Network Collaboration Project. ICE uses the [DepMap quarterly-release schedule](https://depmap.org/portal/) as a foundation to harmonize entity names and pre-align matrices of binary and numeric data to a common space of cell lines.

#### 1. `MATLAB` reference implementation of ICE build

##### We provide a `MATLAB` reference implementation of ICE onboarding scripts, but we note that the data files resulting from running these scripts are platform-independent. This project started as a way to make it easy to consume and align data from large-scale profiling projects initially performed or seeded at the Broad Institute:

 - [Cancer Cell-Line Encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle)
 - [Project Achilles](https://depmap.org/portal/achilles/)
 - [Cancer Therapeutics Response Portal (CTRP)](https://portals.broadinstitute.org/ctrp/)
 - [Cancer Dependency Map (DepMap)](https://depmap.org/portal/)

##### The reference implementation was built for use at the Broad Institute, but hardened to provide a facile starting point for other computational biologists working with the same data. The reference implementation uses three setup functions (`iceopts.m`, `icenv.m`, `chkfile.m`) and a helper file (`tree.csv`) to structure a filesystem in which to run onboarding scripts (`code/build/m`), of which there are several (see below).

##### To run an ICE build on a local copy of this repository, the `data` folder should be populated with downloaded data files (***see Section 3***), after which `build.m` *should* run to completion within about a kilosecond depending on the local system. After a build, the resulting `MAT` files should load within a few seconds.

###### *At the moment, `MATLAB 2018b` is the only version of `MATLAB` explicitly tested with this code, and we note that `MATLAB 2020a` will* **not** *work due to changes in the built-in function `readtable.m` that we will address in a later revision of ICE.*

##### Each script saves output both as `MAT` files and as platform-independent `CSV` files, each representing both data matrices (as sparse triples in the `CSV` files) and rich metadata tables describing row or column entities. Filenaming is controlled to reflect both semantic content and the "shape" of data: for matrices, the numbers of rows and columns; for metadata tables, the number of records.

#### 2. Cell-line names and annotations

##### As a primary source of cell-line names and annotations, we extracted metadata from an internal Broad Institute database named Cancer Dependencies Database (CDDB). CDDB reflects synchronization of cell-line information between the above projects, among others (see below), and communicates with an internal sample-level database based on `ArxSpan` software. CDDB stores names, including synyonyms, and multiple types of annotation such as primary disease, histology, and demographic information. Annotations and synonyms have types that reflect their origins, and synonyms are assigned relative priorities.

##### To ensure that we expose only public annotations, we whitelisted term mappings imported to CDDB from ArxSpan output files to ensure that any annotations from this internal source are also available in other public DepMap files, regardless of their annotation type. We used the most recent public DepMap sample information file (`sample_info.csv`) to ensure that each referenced cell line was present in the namespace constructed from CDDB. For transparency of this process, we include the `Pipeline Pilot` script we used to perform this process (`code/build/pp`), even though CDDB itself is not public.

##### The output of this process is a set of four files that represent a value-add on top of several sources of public information, including the four Broad Institute projects listed above, plus additional public cell-line annotation data from `COSMIC`, `Cell Line Data Base`, and the CTD<sup>2</sup> `Dashboard`. These four files are included directly in this repository (`data\cddb`), which represent (as of **DepMap 20Q2**):

 - `cell-line-best-names.csv`: numeric IDs for 3,302 distinct cell lines mapped to their highest-priority name
 - `cell-line-synon-names.csv`: 18,011 distinct validated names and synonyms for 3,302 cell lines
 - `cell-line-anno-names.csv`: 706 distinct term-type + term-name pairs used for annotating cell lines
 - `map-cell-anno-names.csv`: 18,558 distinct mappings of term-type + term-name pairs mapped to cell lines

##### This dataset uses the Project Achilles name assignment (of the form `ACH-999999`) to index cell lines into a common namespace. For clarity of representation across datasets, we enforce that column indices correspond to the numeric portion of the Project Achilles name for all downstream data matrices with a cell-line dimension, including 82 "missing" Project Achilles names due to deprecated, merged, or redirected name assignments.

#### 3. Datasets addressed by ICE build

##### Not included in this repository (but necessary to run an ICE build) are the source data files themselves, representing different types of data collected on various subsets of the cell lines. Here, we provide links to the data files, whenever possible using a link provided by DepMap to minimize the number of different locations. Regardless of point-of-origin, users downloading these files are responsible to abide by the license terms provided by the original source. The license for this repository covers only the code in this repository and the four data files described explicitly in this readme file.

 - [`data\sample_info.csv`](https://ndownloader.figshare.com/files/22629137)
 - [`data\ctrp\v10.M1.informer_set.txt`](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv1.0_2013_pub_Cell_154_1151/CTRPv1.0_2013_pub_Cell_154_1151.zip)
 - [`data\ctrp\v10.D3.area_under_conc_curve.txt`](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv1.0_2013_pub_Cell_154_1151/CTRPv1.0_2013_pub_Cell_154_1151.zip)
 - [`data\ctrp\v20.meta.per_compound.txt`](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip)
 - [`data\ctrp\v20.meta.per_cell_line.txt`](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip)
 - [`data\ctrp\v20.data.curves_post_qc.txt`](https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip)
 - [`data\ctrp\new-abs-auc-with-qc.txt`](https://github.com/remontoire-pac/ctrp-reference/blob/master/auc/new-abs-auc-with-qc.txt)
 - [`data\depmap\sanger-dose-response.csv`](https://depmap.org/portal/download/api/download/external?file_name=processed_portal_downloads%2Fgdsc-drug-set-export-658c.5%2Fsanger-dose-response.csv)
 - [`data\depmap\primary-screen-replicate-collapsed-logfold-change.csv`](https://ndownloader.figshare.com/files/20237709)
 - [`data\depmap\Achilles_gene_effect.csv`](https://ndownloader.figshare.com/files/22629068)
 - [`data\depmap\Achilles_gene_dependency.csv`](https://ndownloader.figshare.com/files/22629071)
 - [`data\depmap\D2_combined_gene_dep_scores.csv`](https://ndownloader.figshare.com/files/13515395)
 - [`data\depmap\CCLE_expression.csv`](https://ndownloader.figshare.com/files/22897976)
 - [`data\depmap\total-proteome-_v2-normalized-protein-abundance.csv`](https://cds.team/taiga/dataset/total-proteome--5c50/2#)
 - [`data\depmap\CCLE_RPPA_20181003.csv`](https://depmap.org/portal/download/api/download/external?file_name=ccle%2Fccle_2019%2FCCLE_RPPA_20181003.csv)
 - [`data\depmap\CCLE_metabolomics_20190502.csv`](https://depmap.org/portal/download/api/download?file_name=ccle%2Fccle_2019%2FCCLE_metabolomics_20190502.csv&bucket=depmap-external-downloads)
 - [`data\depmap\CCLE_gene_cn.csv`](https://ndownloader.figshare.com/files/22629107)

#### 4. Output of ICE build (links to Figshare; ZIP files there are now versioned by date)

- [`CSV file data structure (platform-independent)`](https://figshare.com/s/37699d7af138e2b643c3)
- [`MAT file data structure (MATLAB)`](https://figshare.com/s/ee840d9eb16d10ae7ab8)
