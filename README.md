## rsfind_bgen: UK Biobank genotype data access
# Overview
BGEN is the compressed file format used to store UK Biobank Genotype Data, it is indexed (cf Tabix) using stored offsets in sqlite database files. Currently delivered to research projects as one file per chromosome.

# Aims
To understand the *bgen* data format by writing a small 'C' application to access a given bgen file using a modified bgen index database (BGI file) in which the Variant table is indexed (in an SQL sense) on rsid.  

