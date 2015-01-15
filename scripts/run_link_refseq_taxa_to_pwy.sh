#!/bin/bash

# script to run link_refseq_taxa_to_pwy.py with pathway long tables
samples=(
a26980-scaffolds 
a26981-scaffolds 
a26982-scaffolds 
a26983-scaffolds 
a26984-scaffolds 
a26985-scaffolds 
a26986-scaffolds 
a26987-scaffolds 
a26988-scaffolds 
a26989-scaffolds 
a26990-scaffolds 
a26991-scaffolds 
a26992-scaffolds 
a26993-scaffolds 
a26994-scaffolds 
a26995-scaffolds 
a26996-scaffolds 
a26997-scaffolds 
a26998-scaffolds 
a26999-scaffolds 
a27000-scaffolds 
a27001-scaffolds 
a27002-scaffolds 
a27003-scaffolds 
a27004-scaffolds 
a27005-scaffolds 
)

# location of important directories
annotation_dir=/Users/nielsh/Dropbox/projects/TaxonomicVariance/annotations
python_resources=/Users/nielsh/Dropbox/projects/TaxonomicVariance/scripts/python_resources
output_dir=/Users/nielsh/Dropbox/projects/TaxonomicVariance/data/pwy_taxa

# create table for each sample
for sample in ${samples[@]}
do
  python2.7 link_refseq_taxa_to_pwy.py --parsed_blast ${annotation_dir}/${sample}.refseq-nr-2014-01-18.LASTout.parsed.txt.gz \
                                            --pwy_long ${annotation_dir}/${sample}.long.pwy.txt.gz \
                                            --ncbi_tree ${python_resources}/ncbi_taxonomy_tree.txt \
                                            --megan_names ${python_resources}/ncbi.map \
                                            -o ${output_dir}/${sample}.long.pwy.taxa.txt \
                                            --lca \
                                            --rpkm ${annotation_dir}/${sample}.orf_rpkm.txt.gz

done

