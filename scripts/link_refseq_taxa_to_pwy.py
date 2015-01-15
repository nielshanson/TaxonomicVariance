#!/usr/bin/python

# import libraries
import sys # for argument vector
import argparse # to parse arguments
import gzip
import re
from python_resources.LCAStar import *

# describe what the script does
what_i_do = "Script to link RefSeq taxonomy to pathway output"

# initialize the parser
parser = argparse.ArgumentParser(description=what_i_do)
parser.add_argument("--parsed_blast", type=str, dest="parsed_blast", default=None,
                   required=True, nargs=1, help='RefSeq parsed B/LAST output file from MetaPathways')
parser.add_argument("--pwy_long", type=str, dest="pwy_long", default=None,
                  required=True, nargs=1, help='Long-table output form extract_pathway_table_from_pgdb.pl (MetaPathways utiltiy script)')
parser.add_argument("--ncbi_tree", type=str, dest="ncbi_tree", default=None,
                  required=True, nargs=1, help='NCBI Tree')
parser.add_argument("--megan_names", type=str, dest="megan_names", default=None,
                    required=False, nargs=1, help='MEGAN preferred names')
parser.add_argument("--lca", dest="lca", action='store_true', default=False,
                    required=False, help='Calcualte LCA')
#parser.add_argument("--lca_star", dest="lca_star", action='store_true', default=False,
#                    required=False, help='Calcualte LCA')
parser.add_argument("--rpkm", type=str, dest="rpkm", default=None,
                    required=False, nargs=1, help='RPKM file.')
parser.add_argument("-o", "--output_file", type=str, dest="output_file", default=None,
                   required=False, nargs=1, help='output file to load into R')

def read_and_clean_parsed_blast(file_name, blast_to_taxonomy, lca=False):
    if lca:
        # calculate LCA
        blast_to_lca = {}
        if ".gz" in file_name:
             with gzip.open(file_name, 'rb') as fh:
                 for line in fh:
                     fields = line.split("\t")
                     hits = taxa_pattern.search(fields[-1])
                     if hits:
                         taxa = hits.group(1)
                         if taxa:
                             if fields[0] not in blast_to_taxonomy:
                                 blast_to_taxonomy[fields[0]] = []
                             blast_to_taxonomy[fields[0]].append(taxa)
             for orf in blast_to_taxonomy:
                temp_list = []
                for t in blast_to_taxonomy[orf]:
                    temp_list.append([t])
                print_id = True
                lca = lcastar.getTaxonomy(temp_list, print_id) # get LCA
                lca_star = lcastar.lca_star(blast_to_taxonomy[orf], True, print_id)
                #print lca_star
                lca_lineage = map(lcastar.translateIdToName, lcastar.get_lineage(lca))
                #lca_star_lineage = map(lcastar.translateIdToName, lcastar.get_lineage(lca_star[0]))
                #print lca_star_lineage
                if orf not in blast_to_lca:
                    blast_to_lca[orf] = ";".join(lca_lineage[::-1])
        blast_to_taxonomy = blast_to_lca
    else:
        if ".gz" in file_name:
            with gzip.open(file_name, 'rb') as fh:
                for line in fh:
                    fields = line.split("\t")
                    hits = taxa_pattern.search(fields[-1])
                    if hits:
                        taxa = hits.group(1)
                        lineage = map(lcastar.translateIdToName, lcastar.get_lineage(lcastar.get_a_Valid_ID([taxa])))
                        if not None in lineage:
                            # only adds the first hit
                            if fields[0] not in blast_to_taxonomy:
                                blast_to_taxonomy[fields[0]] = ";".join(lineage[::-1])
        else:
            with open(file_name, "r") as fh:
                for line in fh:
                    fields = line.split("\t")
                    hits = taxa_pattern.search(fields[-1])
                    if hits:
                        taxa = hits.group(1)
                        lineage = map(lcastar.translateIdToName, lcastar.get_lineage(lcastar.get_a_Valid_ID([taxa])))
                        if not None in lineage:
                            print ";".join(lineage[::-1])
                            # only adds the first hit
                            if fields[0] not in blast_to_taxonomy:
                                blast_to_taxonomy[fields[0]] = lineage
    
    return blast_to_taxonomy

def parse_rpkm_file(fh, orf_id_to_rpkm):
    for line in fh:
        fields = line.split("\t")
        orf_id_to_rpkm[fields[0]] = fields[1].strip("\n")
    return orf_id_to_rpkm

def process_rpkm_file(file_name, orf_id_to_rpkm):
    if ".gz" in file_name:
        with gzip.open(file_name, 'rb') as fh:
            return parse_rpkm_file(fh, orf_id_to_rpkm)
    else:
        with open(file_name, "r") as fh:
            return parse_rpkm_file(fh, orf_id_to_rpkm)

taxa_pattern = re.compile("\[(.*?)\]")

# the main function of the script
def main():
    args = vars(parser.parse_args())
    
    # parse options
    parsed_blast = args["parsed_blast"][0]
    pwy_long = args["pwy_long"][0]
    ncbi_tree = args["ncbi_tree"][0]
    megan_names = args["megan_names"][0]
    output_file = args["output_file"][0]
    lca = args["lca"]
    rpkm_file = args["rpkm"][0]
    
    # parse RPKM value for every ORF
    if rpkm_file:
        orf_id_to_rpkm = {}
        orf_id_to_rpkm = process_rpkm_file(rpkm_file, orf_id_to_rpkm)
    
    # Load NCBI Tree and LCA Star object
    global lcastar
    lcastar = LCAStar(ncbi_tree, megan_names)
    
    # LCA star parameters
    alpha = 0.51
    min_reads = 0
    lcastar.setLCAStarParameters(0, alpha, min_reads)
    print 'Done initializing NCBI Tree'
    
    # capture taxonomic annotation from parsed B/LAST file
    blast_to_taxonomy = {} # map from read to taxonomy
    blast_to_taxonomy = read_and_clean_parsed_blast(parsed_blast, blast_to_taxonomy,lca)
    
    # read through pathway file and append taxonomy field if hit occurs:
    header = ["SAMPLE", "PWY_NAME", "PWY_COMMON_NAME", "RXN_NAME", "RXN_COMMON_NAME", "NUM_REACTIONS",
              "NUM_COVERED_REACTIONS", "ORF_COUNT", "ORF", "TAXONOMY"]
    if orf_id_to_rpkm:
        header.append("RPKM")
    if ".gz" in pwy_long:
        with gzip.open(output_file + ".gz", 'w') as fh_out:
            fh_out.write("\t".join(header) + "\n")
            with gzip.open(pwy_long, 'rb') as fh:
                fh.readline() # read out old header
                for line in fh:
                    fields = line.split("\t")
                    orf_id = fields[-1].strip("\n")
                    fields[-1] = fields[-1].strip("\n")
                    if orf_id in blast_to_taxonomy:
                        fields.append(blast_to_taxonomy[orf_id])
                    else:
                        fields.append("NA")
                    if orf_id_to_rpkm:
                        if orf_id in orf_id_to_rpkm:
                            fields.append(orf_id_to_rpkm[orf_id])
                        else:
                            fields.append("NA")
                    fh_out.write("\t".join(fields) + "\n")
    else:
        with open(output_file + ".gz", 'wb') as fh_out:
            fh_out.write("\t".join(header) + "\n")
            with open(pwy_long, 'rb') as fh:
                fh.readline() # read out old header
                for line in fh:
                    fields = line.split("\t")
                    orf_id = fields[-1].strip("\n")
                    fields[-1] = fields[-1].strip("\n")
                    if orf_id in blast_to_taxonomy:
                        fields.append(blast_to_taxonomy[orf_id])
                    else:
                        fields.append("NA")
                    fh_out.write("\t".join(fields) + "\n")
    
    exit()
    
    

if __name__ == "__main__":
    main()
