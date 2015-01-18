# Taxonomic Diversity and Distinctness Measures Applied to William's Lake Soil Metagenomes
Niels Hanson, Aria Hahn  
January 16, 2015  



# Summary

This document details the implementation of the Clarke's and Warwick's Diversity and Distinctness measures using the functional annotation to reference genomes [[@Warwick:1995vf, @Pienkowski:1998um, @Clarke:1999tx]]. However, the current implementation in the Vegan package does not take into account incomplete taxonomies [[@Oksanen:2007vs]], a common situation when using functional-genes to infer taxonomy via the Lowest Common Ancestor method in meta'omic samples [[@Huson:2007jl]]. Moreover, we would like to take advantage of properties and intuition built into the Weighted Taxonomic Distance (WTD) [[@Hanson:2014bz]], which creates a hard separation between distances across sub-trees. To accomplish this we implement the calculation of distances from partial taxonomies, allowing the efficient calculation of WTD matrices. Moreover, we implement a number of analytical utility scripts to format taxonomic lineages the MetaPathways pipeline to create compatible distance and sample abundance matrix. Using the William's Lake (WL) soil metagenomes as a test case, we using that both taxonomic diversity measures, DeltaStar and DeltaPlus, show a positive correlation with increasing depth.

This document shows the development of the method step-by-step. First it describes the original behavior of the implemented taxonomic diversity measures as implemented in the Vegan package, demonstrating their incompatibility with partial taxonomies. Next, we describe the implementation of a matrix form calculation of partial taxonomies, and through weighting demonstrate the implementation of the WTD. Next, we describe the development and use of utility scripts for use with MetaPathways output for creating taxonomies with respect to samples, pathways, and reactions, using presence/absence, ORF counts, or a normalized read-mapping measure, Reads per Kilobase per Million Mapped reads. Finally, we demonstrate the use of this implementation using 25 metagenomic samples from the William's Lake long term soil productivity (LTSP) site, and show that both diversity metrics reflect the increase in taxonomic diversity with increasing depth in the natural, undisturbed treatment, while a flat or negative trend is observed in the clear-cut treatment. 

## Contents

This analysis consists of the following files and data products:

* [scripts/](): Directory containing perl and python scripts to handle MetaPathways output and ePGDBs
     * [extract_pathway_table_from_pgdb.pl](): Slightly-modified version of the MetaPathways script to extract pathways in the 'long' format with pathways and reactions as variables
     * [link_refseq_taxa_to_pathway.py](): Script to link pathway tables with RefSeq Best-HIT, Lowest Common Ancestor (LCA), LCA* taxonomies, as well as RPKM values. This script requires the parsed B/LAST outputs, RPKM files, produced by MetaPathways, and the 'long' pathway tables produced by `extract_pathway_table_from_pgdb.pl` from the ePGDBs. Large tables can be gziped (.gz) for convenience.
     * [extract_LTSP_WL_pwys.sh](): Shell script uses `extract_pathway_table_from_pgdb.pl` to extract WL pathways from Pathway Tools
     * [run_link_refseq_taxa_to_pathway.sh](): Shell script runs `link_refseq_taxa_to_pathway.py` on the WL samples
     * [python_resources/](): Folder containing classes and methods for running `link_refseq_taxa_to_pathway.py`
          * [ncbi_taxonomy_tree.txt](): Slightly-modified NCBI Taxonomy Database Hierarchy.
          * [ncbi.map](): Preferred taxonomy names from NCBI IDs.
          * [LCAStar.py](): Class and functions for calculating LCA, MEGAN Taxonomy, and LCA*.

* [annotations/](): Directory containing RefSeq annotations (`.LASTout.parsed.txt.gz`), RPKM values (`.orf_rpkm.txt.gz`), and long-style pathway tables (`.long.pwy.txt.gz`) for the WL LTSP samples.
* [r_scripts/](): Directory of r functions.
     * [my_taxa2dist.R](): Modification of the `taxa2dist()` to handle partial taxonomies and the WTD
* [data/](): Processed data ready for analysis in R
     * [simple_test.txt](): Very basic taxonomy matrix with incomplete taxonomies for testing our implementation
     * [Phylodisttest.txt](): Matrix of partial taxonomies from the WL LTSP samples 
     * [LTSP_WL_MappingandInfo_clean.csv](): Matrix of metadata, chemistry and experimental conditions for the WL samples
          * [pwy_taxa/](): long-style pathway tables with LCA taxonomies and RPKM values from the LTSP samples
     * [temp/](): Directory to store previous `.Rdata` files for fast processing. Remove contents of this directory to reboot the analysis and reprocess taxonomic distance matrix
* [figures/](): Directory for pdf figures in the analysis.

## R Libraries

* Load R libraries


```r
library(vegan)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
figure_dir <- paste(wd, "figures", sep="/")
```

# Implementation

The initial implementation in VEGAN does allow for incomplete taxonomies. In this section we demonstrate this behavior using the Dutch Dune dataset and a subset of the WL LTSP partial taxonomies. The Dutch Dune Meadow dataset has full taxonomies and is not particularly large compared to the number of unique taxonomies found in environmental datasets. We observe that the current implementation of `taxa2dist`, while improved to some degree class-based scaling, produces distances that do not respect partial taxonomies.

* Load the Dutch Dune Meadow data.


```r
data(dune)
data(dune.taxon)
```

* Run through the example data


```r
# calculate the taxonomic distance
dune.taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
# plot the tree
dune.taxontree <- hclust(dune.taxdis)
plot(dune.taxontree)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-4-1.png) 

* Plotting the matrix as a heatmap will show the structure in more detail


```r
pheatmap(dune.taxdis, treeheight_row=150, treeheight_col=150)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-5-1.png) 

* Vegan's implementation of the taxonomic distance


```r
mod <- taxondive(dune, dune.taxdis)
summary(mod)
```

```
##            Delta  Delta*  Delta+ sd(Delta+) z(Delta+) Pr(>|z|)   
## 1        25.0089 32.1543 51.5455     9.5637   -2.8405 0.004504 **
## 2        60.4931 66.3497 66.6869     4.6664   -2.5769 0.009970 **
## 3        46.5985 51.7024 70.7475     4.6664   -1.7067 0.087881 . 
## 4        58.1988 63.1763 73.4033     3.5073   -1.5135 0.130150   
## 5        72.0452 76.9903 77.6024     3.2187   -0.3446 0.730386   
## 6        76.4148 83.1205 80.8430     4.2170    0.5054 0.613259   
## 7        70.2500 75.4752 78.3974     3.5073   -0.0896 0.928626   
## 8        59.1259 63.4363 74.7865     3.8361   -1.0232 0.306220   
## 9        56.9481 60.9854 71.1597     3.5073   -2.1532 0.031303 * 
## 10       68.9021 74.5133 77.7617     3.8361   -0.2476 0.804432   
## 11       76.2014 85.1259 82.5379     5.2089    0.7346 0.462605   
## 12       69.5554 77.7922 82.6136     5.2089    0.7491 0.453792   
## 13       55.2961 62.9232 76.6566     4.6664   -0.4404 0.659657   
## 14       77.6087 89.2500 88.1818     6.7459    1.4039 0.160362   
## 15       74.9245 84.2485 86.5584     5.8818    1.3341 0.182178   
## 16       57.8435 66.5389 73.3604     5.8818   -0.9098 0.362936   
## 17       64.8225 72.4081 70.9740     6.7459   -1.1470 0.251377   
## 18       76.7314 85.7730 81.4015     5.2089    0.5164 0.605570   
## 19       73.4487 81.3182 83.3712     5.2089    0.8945 0.371029   
## 20       78.0762 87.0634 87.9221     5.8818    1.5659 0.117368   
## Expected 73.2888 70.1816 78.7116                                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

* Try it out with some partial taxonomies


```r
taxon.test = read.table(paste(wd, "data/Phylodisttest.txt", sep="/"), header=TRUE, row.names="row.names", sep ="\t")
```

* Build a tree based on the distances as calculated via the original implementation of taxa2dist in Vegan. We will vary the step lengths between successive levels relative to proportional loss of the number of distinct classes.


```r
taxon.test.matrix <- taxa2dist(taxon.test, varstep=TRUE)
taxontree <- hclust(taxon.test.matrix)
plot(taxontree)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-8-1.png) 

* Structure with number of distinct class step lengths.


```r
pheatmap(taxon.test.matrix, treeheight_row=150, treeheight_col=150)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-9-1.png) 

```r
pdf(file = paste(figure_dir, "WL_matrix_class_step.pdf", sep="/"), width=12, height=12)
pheatmap(as.matrix(taxon.test.matrix), treeheight_row=150, treeheight_col=150)
dev.off()
```

```
## pdf 
##   2
```

* Uniform step lengths


```r
taxon.test.matrix2 <- taxa2dist(taxon.test)
taxontree2 <- hclust(taxon.test.matrix2)
plot(taxontree2)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-10-1.png) 

* Structure with uniform step lengths.


```r
pheatmap(taxon.test.matrix2, treeheight_row=150, treeheight_col=150)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-11-1.png) 

## Implementing WTD and Partial Taxonomy

This demonstrates our implementation of the weighted taxonomic distance and partial taxonomies in `my_taxa2dist.R`, a slightly modified version of the original. 

*TODO: Add some details of the matrix form exclusive OR and the outer product form of the taxonomic distance matrix.*

* source the implementation of `my_taxa2dist`:


```r
source(paste(wd, "r_scripts/my_taxa2dist.R", sep="/"))
```

* we'll do a simple test to do a sanity check of the implementation


```r
simple_test1 <- read.table(paste(wd, "data/simple_test.txt", sep="/"), header=TRUE, row.names="row.names", sep ="\t")
simple_test1_dist <- my_taxa2dist(simple_test1, check=FALSE, wtd=TRUE)
simple_test1_dist
```

```
##          A;B;C;D  A;B;E;-
## A;B;E;- 16.12903         
## A;F;-;- 35.48387 32.25806
```

* plot the simple tree


```r
plot(hclust(simple_test1_dist))
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-14-1.png) 

* We'll do another quick test with the LTSP data


```r
taxon.test.wtd <- my_taxa2dist(taxon.test, check=FALSE, wtd=TRUE)
pheatmap(as.matrix(taxon.test.wtd), treeheight_row=150, treeheight_col=150)
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-15-1.png) 

```r
# create pdf
pdf(file = paste(figure_dir, "WL_matrix_wtd.pdf", sep="/"), width=9.9603, height=9.614)
pheatmap(as.matrix(taxon.test.wtd), treeheight_row=150, treeheight_col=150)
dev.off()
```

```
## pdf 
##   2
```

* lets take a look at the Dune data again, looks about the same


```r
# calculate the taxonomic distance
dune.taxdis.wtd <- my_taxa2dist(dune.taxon, check=FALSE, varstep=FALSE, wtd=TRUE)
```

```
## Warning in my_taxa2dist(dune.taxon, check = FALSE, varstep = FALSE, wtd =
## TRUE): you used 'check=FALSE' and some distances are zero -- was this
## intended?
```

```r
# plot the tree
dune.taxontree.wtd <- hclust(dune.taxdis.wtd)
par(mfrow=c(1,2))
plot(dune.taxontree, main="Dune (Original)")
plot(dune.taxontree.wtd, main="Dune (WTD)" )
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-16-1.png) 

```r
par(mfrow=c(1,1))

# create pdf
pdf(file = paste(figure_dir, "dune_wtd_compare.pdf", sep="/"), width=16, height=10)
par(mfrow=c(1,2))
plot(dune.taxontree, main="Dune (Original)")
plot(dune.taxontree.wtd, main="Dune (WTD)" )
par(mfrow=c(1,1))
dev.off()
```

```
## pdf 
##   2
```

```r
mod <- taxondive(dune, dune.taxdis.wtd)
summary(mod)
```

```
##            Delta  Delta*  Delta+ sd(Delta+) z(Delta+) Pr(>|z|)  
## 1        12.8852 16.5666 40.9524    12.2475   -2.1681  0.03015 *
## 2        51.6454 56.6454 56.8607     5.7431   -1.8536  0.06380 .
## 3        35.9178 39.8519 62.2222     5.7431   -0.9200  0.35757  
## 4        48.7606 52.9309 64.5096     4.2466   -0.7056  0.48044  
## 5        62.8979 67.2152 68.1319     3.8791    0.1614  0.87181  
## 6        68.3046 74.2986 71.5152     5.1594    0.7771  0.43712  
## 7        61.3187 65.8796 68.8645     4.2466    0.3199  0.74903  
## 8        47.2487 50.6932 64.5022     4.6680   -0.6435  0.51991  
## 9        45.2704 48.4798 60.7245     4.2466   -1.5969  0.11028  
## 10       60.3175 65.2295 68.5426     4.6680    0.2221  0.82426  
## 11       68.2284 76.2191 70.7231     6.4524    0.4986  0.61806  
## 12       61.3579 68.6239 73.5450     6.4524    0.9359  0.34930  
## 13       44.9495 51.1494 65.6085     5.7431   -0.3304  0.74111  
## 14       70.4854 81.0582 78.4580     8.4818    1.2912  0.19662  
## 15       67.1937 75.5556 77.0975     7.3378    1.3071  0.19116  
## 16       46.3985 53.3734 62.8118     7.3378   -0.6397  0.52236  
## 17       55.8125 62.3438 60.3175     8.4818   -0.8475  0.39671  
## 18       67.9600 75.9681 69.3122     6.4524    0.2799  0.77953  
## 19       62.7001 69.4180 71.2522     6.4524    0.5806  0.56151  
## 20       71.1043 79.2889 80.2721     7.3378    1.7398  0.08190 .
## Expected 63.5233 60.8301 67.5059                                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Taxonomy and Pathways from MetaPathways

* Slightly modified MetaPathways utility script `extract_pathway_table_from_pgdb.pl` that produces long-table output labeling pathways and reactions
* wrote `scripts/link_refseq_taxa_to_pwy.py` to combine output of pathways in metapathways with parsedblast and RPKM values

```
# Run Pathway Tools in API mode
./pathway-tools/pathway-tools -api

# navigate to the scripts directory 
cd .../scripts/

# Extract pathways
./extract_LTSP_WL_pwys.sh

# Calculate LCA and prepare pathway taxonomies
./run_link_refseq_taxa_to_pathway.sh
```

# Analysis of William's Lake

Here we utilize our implementation to analyze the William's Lake metagenomes by sample and pathway.

## Diversity by Sample

In order to calculate the diversity measures, we need to create a distance matrix and a sample by taxa abundance matrix. To do this we will parse all the extract taxonomies and create a non-redundant list of taxonomies, and then use these to create the distance matrix for our samples using the WTD. From here we construct two sample by taxa abundance matrices; the first based on the raw ORF counts, and the second based on the RPKM normalized abundance measure. In both cases we show that the taxonomic diversity and distinctness measures, delta* and delta+, show a strong positive correlation with depth in the Natural samples, while the treated samples show no clear or negative correlations.

* Various functions for creating properly formatted taxa matrix from partial taxonomies for `my_taxa2dist()`


```r
# fills in the rest of the NA
create_rows <- function(x, max) {
  c(x, rep(NA, max - length(x)))
}

# creates dataframe delimiting the taxonomy hierarchy
create_taxa_df <- function(x) {
  split_taxa <- sapply(x, strsplit, ";")
  lengths <- sapply(split_taxa, length)
  max_length <- max(lengths)
  rows <- lapply(split_taxa, create_rows, max_length)
  rows <- lapply(rows, rev) # reverse rows
  taxa_df <- t(data.frame(rows))
  row.names(taxa_df) <- x
  taxa_df[is.na(taxa_df)] <- "unclassified"
  return(taxa_df)
}
```

* Create a non-redundant list of unique taxonomies for each of the annotations


```r
pwy_taxa_dir <- paste(wd, "data/pwy_taxa/", sep="/")
temp_dir <- paste(wd, "data", "temp", sep="/")

# select files from the pwy_taxa_dir
pwy_taxa_files <- list.files(pwy_taxa_dir, pattern="*taxa.txt.gz")
rdata_files <- list.files(temp_dir, pattern="*.Rdata")

# create non-redundant list of taxonomies for all samples
nr_taxa_list <- c()
if(!("nr_taxa_list.Rdata" %in% rdata_files)) { 
  for (pt_file in pwy_taxa_files) {
    print(pt_file)
    df <- read.table(paste(pwy_taxa_dir, pt_file, sep="/"), header=TRUE, sep="\t",
                   colClasses= c('character', 'character', 'character', 'character', 'character',
                                 'numeric', 'numeric', 'numeric', 'character', 'character', 'numeric'), quote = "")
    my_taxa <- unique(df$TAXONOMY)
    nr_taxa_list <- c(nr_taxa_list, my_taxa)
  }
  # save to disk
  nr_taxa_list <- unique(nr_taxa_list)
  save(nr_taxa_list, file=paste(temp_dir, "nr_taxa_list.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "nr_taxa_list.Rdata", sep="/"))
}
# remove NA class
nr_taxa_list <- nr_taxa_list[!is.na(nr_taxa_list)]
```

* Construct the WTD distance matrix from these annotations


```r
if(!("taxdis.Rdata" %in% rdata_files)) { 
  nr_taxa_list <- unique(nr_taxa_list)
  split_taxa <- sapply(nr_taxa_list, strsplit, ";")
  lengths <- sapply(split_taxa, length)
  max_length <- max(lengths)
  test_df <- create_taxa_df(nr_taxa_list)
  taxdis <- my_taxa2dist(test_df, check=FALSE, varstep=FALSE, wtd=TRUE)

  # save to disk
  save(taxdis, file=paste(temp_dir, "taxdis.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "taxdis.Rdata", sep="/"))
}
```

* Construct the sample by taxa matrix (ORF counts)


```r
sample_taxa_matrix_orfs <- NULL
if(!("sample_taxa_matrix_orfs.Rdata" %in% rdata_files)) { 
  for (pt_file in pwy_taxa_files) {
    print(pt_file)
    df <- read.table(paste(pwy_taxa_dir, pt_file, sep="/"), header=TRUE, sep="\t",
                    colClasses= c('character', 'character', 'character', 'character', 'character',
                                  'numeric', 'numeric', 'numeric', 'character', 'character', 'numeric'), quote = "")
    my_row <- table(df$TAXONOMY)[nr_taxa_list]
    my_row[is.na(my_row)] <- 0
    sample_taxa_matrix_orfs <- rbind(sample_taxa_matrix_orfs,my_row)
  }
  # set sample as row_names
  colnames(sample_taxa_matrix_orfs) <- nr_taxa_list
  row.names(sample_taxa_matrix_orfs) <- sub("-scaffolds.long.pwy.taxa.txt.gz", "", pwy_taxa_files)
  
  save(sample_taxa_matrix_orfs, file=paste(temp_dir, "sample_taxa_matrix_orfs.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "sample_taxa_matrix_orfs.Rdata", sep="/"))
}
```

* Calculate the diversity scores (ORF Counts)


```r
if(!("diversity_result_orfs.Rdata" %in% rdata_files)) { 
  diversity_result_orfs <- taxondive(sample_taxa_matrix_orfs, taxdis)
  save(diversity_result_orfs, file=paste(temp_dir, "diversity_result_orfs.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "diversity_result_orfs.Rdata", sep="/"))
}
div_res_orfs <- summary(diversity_result_orfs)
div_res_orfs
```

```
##             Delta   Delta*   Delta+ sd(Delta+) z(Delta+)  Pr(>|z|)    
## a26980   2.713658 3.300335 5.569843   0.133225   -3.7608 0.0001694 ***
## a26981   2.832185 3.350476 5.456805   0.076916   -7.9836 1.421e-15 ***
## a26982   2.213846 3.182478 5.776865   0.097220   -3.0242 0.0024932 ** 
## a26983   2.264232 3.155152 5.821844   0.074784   -3.3300 0.0008685 ***
## a26984   2.722351 3.328560 5.525936   0.102338   -5.3249 1.010e-07 ***
## a26985   2.406815 3.232813 5.740493   0.130431   -2.5330 0.0113092 *  
## a26986   2.328577 3.174704 5.740546   0.080343   -4.1115 3.932e-05 ***
## a26987   2.650290 3.373605 5.730377   0.065297   -5.2146 1.842e-07 ***
## a26988   2.752728 3.376664 5.606916   0.086292   -5.3766 7.590e-08 ***
## a26989   2.546945 3.268091 5.731587   0.068330   -4.9655 6.854e-07 ***
## a26990   2.365828 3.123147 5.609957   0.078196   -5.8944 3.761e-09 ***
## a26991   2.582217 3.153291 5.416265   0.110195   -5.9405 2.842e-09 ***
## a26992   2.518166 3.078194 5.284151   0.081861   -9.6105 < 2.2e-16 ***
## a26993   2.250203 2.884950 5.431156   0.076353   -8.3784 < 2.2e-16 ***
## a26994   2.265440 2.968183 5.516938   0.091396   -6.0608 1.354e-09 ***
## a26995   2.334446 3.105256 5.676712   0.089167   -4.4205 9.848e-06 ***
## a26996   2.573047 3.288550 5.710418   0.063524   -5.6743 1.392e-08 ***
## a26997   2.622502 3.162706 5.253316   0.082588   -9.8992 < 2.2e-16 ***
## a26998   2.426347 2.925116 5.252682   0.107952   -7.5792 3.476e-14 ***
## a26999   2.317200 3.001527 5.585082   0.103890   -4.6760 2.925e-06 ***
## a27000   2.314597 3.112346 5.721840   0.093060   -3.7506 0.0001764 ***
## a27001   2.620753 3.321069 5.748672   0.074039   -4.3518 1.350e-05 ***
## a27002   2.632780 3.173803 5.473587   0.089111   -6.7027 2.046e-11 ***
## a27003   2.474919 3.114082 5.639768   0.079971   -5.3908 7.015e-08 ***
## a27004   2.272844 2.985303 5.588673   0.102924   -4.6851 2.799e-06 ***
## a27005   2.601804 3.321572 5.712673   0.075330   -4.7551 1.984e-06 ***
## Expected 2.540804 2.539189 6.070875                                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

* Load the William's Lake metadata


```r
wl_meta_df <- read.table(paste(wd, "data", "LTSP_WL_MappingandInfo_clean.csv", sep="/"), sep=",", header=TRUE)
```

* Create summary Data Frame (ORF Counts)


```r
div_res_orfs_matrix <- as.matrix(div_res_orfs)
div_res_orfs_matrix <- cbind(div_res_orfs_matrix[-nrow(div_res_orfs_matrix),], as.matrix(wl_meta_df)) # combine with df
master_df_orfs <- data.frame(delta=as.numeric(div_res_orfs_matrix[,1]),
           delta_star=as.numeric(div_res_orfs_matrix[,2]),
           delta_plus=as.numeric(div_res_orfs_matrix[,3]),
           delta_plus_sd=as.numeric(div_res_orfs_matrix[,4]),
           delta_plus_z=as.numeric(div_res_orfs_matrix[,5]),
           delta_plus_pr=as.numeric(div_res_orfs_matrix[,6]),
           horizon=as.character(div_res_orfs_matrix[,8]),
           treatment=as.character(div_res_orfs_matrix[,9]),
           depth=as.numeric(div_res_orfs_matrix[,10]),
           pH=as.numeric(div_res_orfs_matrix[,11]),
           totalC_percent=as.numeric(div_res_orfs_matrix[,12]),
           totalN_percent=as.numeric(div_res_orfs_matrix[,13]),
           totalorganicN_ppm=as.numeric(div_res_orfs_matrix[,14]),
           nh4_ppm=as.numeric(div_res_orfs_matrix[,15]),
           no3_ppm=as.numeric(div_res_orfs_matrix[,16]))
```

* Depth (cm) by taxonomic diversity (ORF Counts)


```r
p <- ggplot(master_df_orfs, aes(x=depth, y=delta_star)) 
p <- p + geom_smooth(method="glm", size=1.5) 
p <- p + geom_point(size=3.5, aes(col=horizon)) 
p <- p + theme_bw() 
p <- p + facet_wrap(~treatment)
p <- p + ggtitle("Williams Lake (ORF Counts)")
p
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-24-1.png) 

```r
pdf(file = paste(figure_dir, "wl_delta_star_orfs.pdf", sep="/"), width=9.832, height=4.399)
p
dev.off()
```

```
## pdf 
##   2
```


```r
# delta plus (ORF count)
# quartz(width=9.832677, height=4.399606)
p <- ggplot(master_df_orfs, aes(x=depth, y=delta_plus)) 
p <- p + geom_smooth(formula=y ~ x, method="glm", size=1.5) 
p <- p + geom_point(size=3.5, aes(col=horizon)) 
p <- p + facet_wrap(~treatment) 
p <- p + theme_bw()
p <- p + ggtitle("Williams Lake (ORF Counts)")
p
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-25-1.png) 

```r
pdf(file = paste(figure_dir, "wl_delta_plus_orfs.pdf", sep="/"), width=9.832, height=4.399)
p
dev.off()
```

```
## pdf 
##   2
```

* Construct the sample by taxa matrix (RPKM)


```r
div_res_rpkm_matrix <- NULL
if(!("div_res_rpkm_matrix.Rdata" %in% rdata_files)) { 
  for (pt_file in pwy_taxa_files) {
    print(pt_file)
    df <- read.table(paste(pwy_taxa_dir, pt_file, sep="/"), header=TRUE, sep="\t",
                   colClasses= c('character', 'character', 'character', 'character', 'character',
                                 'numeric', 'numeric', 'numeric', 'character', 'character', 'numeric'), quote = "")
    
    res <- df %>% 
    select(TAXONOMY, RPKM) %>%
    group_by(TAXONOMY) %>% 
    summarize( rpkm=sum(RPKM) )
    
    my_row <- res$rpkm
    names(my_row) <- res$TAXONOMY
    my_row <- my_row[nr_taxa_list]
    my_row[is.na(my_row)] <- 0
    div_res_rpkm_matrix <- rbind(div_res_rpkm_matrix,my_row)
  }
  
  # set sample as row_names
  colnames(div_res_rpkm_matrix) <- nr_taxa_list
  row.names(div_res_rpkm_matrix) <- sub("-scaffolds.long.pwy.taxa.txt.gz", "", pwy_taxa_files)
  
  save(div_res_rpkm_matrix, file=paste(temp_dir, "div_res_rpkm_matrix.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "div_res_rpkm_matrix.Rdata", sep="/"))
}
```

* RPKM diversity results


```r
if(!("diversity_result_rpkm.Rdata" %in% rdata_files)) { 
  diversity_result_rpkm <- taxondive(div_res_rpkm_matrix, taxdis)
  save(diversity_result_rpkm, file=paste(temp_dir, "diversity_result_rpkm.Rdata", sep="/"))
} else {
  load(paste(temp_dir, "diversity_result_rpkm.Rdata", sep="/"))
}
div_res_rpkm <- summary(diversity_result_rpkm)
div_res_rpkm
```

```
##             Delta   Delta*   Delta+ sd(Delta+) z(Delta+)  Pr(>|z|)    
## a26980   2.758167 3.302781 5.569843   0.133225   -3.7608 0.0001694 ***
## a26981   2.810832 3.302222 5.456805   0.076916   -7.9836 1.421e-15 ***
## a26982   2.207982 3.153793 5.776865   0.097220   -3.0242 0.0024932 ** 
## a26983   2.241932 3.125682 5.821844   0.074784   -3.3300 0.0008685 ***
## a26984   2.735179 3.305028 5.525936   0.102338   -5.3249 1.010e-07 ***
## a26985   2.392132 3.107883 5.740493   0.130431   -2.5330 0.0113092 *  
## a26986   2.329099 3.134762 5.740546   0.080343   -4.1115 3.932e-05 ***
## a26987   2.648708 3.371166 5.730377   0.065297   -5.2146 1.842e-07 ***
## a26988   2.792178 3.356256 5.606916   0.086292   -5.3766 7.590e-08 ***
## a26989   2.424724 3.158521 5.731587   0.068330   -4.9655 6.854e-07 ***
## a26990   2.596629 3.277595 5.609957   0.078196   -5.8944 3.761e-09 ***
## a26991   2.466538 2.956193 5.416265   0.110195   -5.9405 2.842e-09 ***
## a26992   2.438620 2.927395 5.284151   0.081861   -9.6105 < 2.2e-16 ***
## a26993   2.164103 2.744631 5.431156   0.076353   -8.3784 < 2.2e-16 ***
## a26994   2.197456 2.820860 5.516938   0.091396   -6.0608 1.354e-09 ***
## a26995   2.328675 3.069486 5.676712   0.089167   -4.4205 9.848e-06 ***
## a26996   2.568293 3.297345 5.710418   0.063524   -5.6743 1.392e-08 ***
## a26997   2.610637 3.106320 5.253316   0.082588   -9.8992 < 2.2e-16 ***
## a26998   2.257050 2.702555 5.252682   0.107952   -7.5792 3.476e-14 ***
## a26999   2.269697 2.921286 5.585082   0.103890   -4.6760 2.925e-06 ***
## a27000   2.287615 3.076943 5.721840   0.093060   -3.7506 0.0001764 ***
## a27001   2.593318 3.282045 5.748672   0.074039   -4.3518 1.350e-05 ***
## a27002   2.571294 3.081967 5.473587   0.089111   -6.7027 2.046e-11 ***
## a27003   2.388313 2.982336 5.639768   0.079971   -5.3908 7.015e-08 ***
## a27004   2.231190 2.892473 5.588673   0.102924   -4.6851 2.799e-06 ***
## a27005   2.552190 3.331558 5.712673   0.075330   -4.7551 1.984e-06 ***
## Expected 2.518899 2.517737 6.070875                                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

* RPKM result matrix


```r
res_matrix_rpkm <- as.matrix(div_res_rpkm)
res_matrix_rpkm <- cbind(res_matrix_rpkm[-nrow(res_matrix_rpkm),], as.matrix(wl_meta_df)) # combine with df
master_df_rpkm <- data.frame(delta=as.numeric(res_matrix_rpkm[,1]),
           delta_star=as.numeric(res_matrix_rpkm[,2]),
           delta_plus=as.numeric(res_matrix_rpkm[,3]),
           delta_plus_sd=as.numeric(res_matrix_rpkm[,4]),
           delta_plus_z=as.numeric(res_matrix_rpkm[,5]),
           delta_plus_pr=as.numeric(res_matrix_rpkm[,6]),
           horizon=as.character(res_matrix_rpkm[,8]),
           treatment=as.character(res_matrix_rpkm[,9]),
           depth=as.numeric(res_matrix_rpkm[,10]),
           pH=as.numeric(res_matrix_rpkm[,11]),
           totalC_percent=as.numeric(res_matrix_rpkm[,12]),
           totalN_percent=as.numeric(res_matrix_rpkm[,13]),
           totalorganicN_ppm=as.numeric(res_matrix_rpkm[,14]),
           nh4_ppm=as.numeric(res_matrix_rpkm[,15]),
           no3_ppm=as.numeric(res_matrix_rpkm[,16]))
```

* Depth by taxonomic diversity (RPKM)


```r
# delta_star (RPKM count)
p <- ggplot(master_df_rpkm, aes(x=depth, y=delta_star)) 
p <- p + geom_smooth(method="glm", size=1.5) 
p <- p + geom_point(size=3.5, aes(col=horizon)) 
p <- p + theme_bw() 
p <- p + facet_wrap(~treatment)
p <- p + ggtitle("William's Lake (RPKM)")
p
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-29-1.png) 

```r
pdf(file = paste(figure_dir, "wl_delta_star_rpkm.pdf", sep="/"), width=9.832, height=4.399)
p
dev.off()
```

```
## pdf 
##   2
```


```r
# delta plus (RPKM count)
p <- ggplot(master_df_rpkm, aes(x=depth, y=delta_plus)) 
p <- p + geom_smooth(formula=y ~ x, method="glm", size=1.5) 
p <- p + geom_point(size=3.5, aes(col=horizon)) 
p <- p + facet_wrap(~treatment) 
p <- p + theme_bw()
p <- p + ggtitle("William's Lake (RPKM)")
p
```

![](./taxonomic_variance_notes_files/figure-html/unnamed-chunk-30-1.png) 

```r
pdf(file = paste(figure_dir, "wl_delta_plus_rpkm.pdf", sep="/"), width=9.832, height=4.399)
p
dev.off()
```

```
## pdf 
##   2
```

## Diversity by Pathway

* Construct the Sample, Pathway, Measure, and Values table 
     * Calculate diversity per pathway


```r
# cleans up summary results from taxondive
construct_summary_matrix <- function(taxon_dive_summary, clean=TRUE) {
  x <- data.frame(delta=taxon_dive_summary[,1],
                 delta_star=taxon_dive_summary[,2],
                 delta_plus=taxon_dive_summary[,3],
                 delta_plus_sd=taxon_dive_summary[,4],
                 delta_plus_z=taxon_dive_summary[,5],
                 delta_plus_pr=taxon_dive_summary[,6])
  # save expected values
  delta_e <- x[nrow(x),1]
  delta_star_e <- x[nrow(x),2]
  delta_plus_e <- x[nrow(x),3]
  
  x$PWY_NAME = rownames(x)
  if (clean) {
    # removes NAs, NaN, Inf
    x <- x[complete.cases(x),]
    # must be positive 
    positive_set <- (x[,1] >= 0 & x[,2] >= 0 & x[,3] >= 0)
    x <- x[positive_set,]
  }
  
  # bind results together
  res <- list(df=x, delta_e=delta_e, delta_star_e=delta_star_e, delta_plus_e=delta_plus_e)
  return(res)
}
```


```r
load(paste(temp_dir, "pathway_master_table.Rdata", sep="/"))
# clean up names
pathway_master_table$SAMPLE <- sub("-scaffolds", "", pathway_master_table$SAMPLE)
pathway_master_table$SAMPLE <- sub("a", "A", pathway_master_table$SAMPLE)
```
```


```r
if("pathway_master_table.Rdata" %in% rdata_files) {
  load(paste(temp_dir, "pathway_master_table.Rdata", sep="/"))
} else {
  pathway_master_table = NULL
  for (pt_file in pwy_taxa_files) {
    print(pt_file)
    # load dataframe
    df <- read.table(paste(pwy_taxa_dir, pt_file, sep="/"), header=TRUE, sep="\t",
                    colClasses= c('character', 'character', 'character', 'character', 'character',
                                 'numeric', 'numeric', 'numeric', 'character', 'character', 'numeric'), quote = "")
    ## calculate diversity per pathway
    
    # summarize pathways by taxonomy abundance
    res <- df %>% 
        select(SAMPLE, PWY_NAME, PWY_COMMON_NAME, TAXONOMY, NUM_REACTIONS, NUM_COVERED_REACTIONS, ORF_COUNT, RPKM) %>%
        group_by(SAMPLE, PWY_NAME, PWY_COMMON_NAME, TAXONOMY) %>%
        filter( is.na(TAXONOMY) != TRUE) %>%
        summarize( ORF_COUNT=max(ORF_COUNT), RPKM=sum(RPKM)  )
    
    # convert into dataframe for taxonomic calculation
    res2 <- dcast(res, PWY_NAME ~ TAXONOMY, value.var='RPKM', sum)
    rownames(res2) <- res2$PWY_NAME
    res2 <- res2[-1]

    # prepare taxonomic distance matrix
    my_taxa_list <- unique(colnames(res2))
    split_taxa <- sapply(my_taxa_list, strsplit, ";")
    lengths <- sapply(split_taxa, length)
    max_length <- max(lengths)
    test_df <- create_taxa_df(my_taxa_list)
    my_taxdis <- my_taxa2dist(test_df, check=FALSE, varstep=FALSE, wtd=TRUE)
    
    # calculate statistics
    test <- taxondive(res2, my_taxdis)
    summary_res <- summary(test)
    
    # clean up summary
    test_list <- construct_summary_matrix(summary_res)

    # create a data frame of sample, pathway, measure, value
    res3 <- df %>% 
        select(SAMPLE, PWY_NAME, PWY_COMMON_NAME, TAXONOMY, NUM_REACTIONS, NUM_COVERED_REACTIONS, ORF_COUNT, RPKM) %>%
        group_by(SAMPLE, PWY_NAME, PWY_COMMON_NAME) %>%
        filter( is.na(TAXONOMY) != TRUE) %>%
        summarize( ORF_COUNT=max(ORF_COUNT), RPKM=sum(RPKM)  )

    # bind new statistics to pathway
    res4 <- merge(res3, test_list$df, by="PWY_NAME")
    
    # add to master table
    pathway_master_table <- rbind(pathway_master_table, res4)
  }
  
  # save when complete
  save(pathway_master_table, file=paste(temp_dir, "pathway_master_table.Rdata", sep="/"))
}

# clean up names
pathway_master_table$SAMPLE <- sub("-scaffolds", "", pathway_master_table$SAMPLE)
pathway_master_table$SAMPLE <- sub("a", "A", pathway_master_table$SAMPLE)
```

* Now with the pathway master table constructed we can compute a few things






```r
# join by metadata
names(wl_meta_df)[names(wl_meta_df) == 'Library'] <- "SAMPLE"
pathway_master_table <- merge(pathway_master_table, wl_meta_df, "SAMPLE")

res<-pathway_master_table %>%
  select(SAMPLE, PWY_NAME, PWY_COMMON_NAME, ORF_COUNT, delta_star) %>%
  group_by(SAMPLE, PWY_NAME) %>%
  summarize(PWY_COMMON_NAME, min_delta_star = min(delta_star, na.rm = TRUE), ORF_COUNT) %>%
  filter(ORF_COUNT > 10) %>%
  filter(min_rank(min_delta_star) < 5)
```

# Future Ideas

Ideas for future work.

* Bootstrap p-values for delta_star (could be time-consuming)
* Investigate measures in the context of the pathway

# References

Section for work cited.
