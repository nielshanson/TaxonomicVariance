#!/bin/bash

# script to extract long-pathway tables from each LTSP WL pwy
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

for sample in ${samples[@]}
do
  perl extract_pathway_table_from_pgdb.pl -f $sample -out $sample.long.pwy.txt -t long
done
