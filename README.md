# Ion_reporter_Oncology
Several scripts for manipulating and wrangling raw variant files from Thermo Fisher's Ion Reporter Software

# vcf_to_xlsx script
This R script allows the user to extract oncology relevant information (Genomic_Instability, percent_LOH, MSI_Status, etc) from columns in a vcf output and transform it to an excel/csv format. 
In order for the script to work, the input vcf file must be a vcf output after running the "Oncomine Comprehensive Plus - w3.1 - DNA - Single Sample" workflow on the Ion Reporter Software

# Variant_filter_analysis script
This R script automatically reviews the output of variant filters created on Ion Reporter. Sometimes, the created filters on the Ion Reporter miss important variants, so R does the manual checking, review and verification of the filtered results from Ion Reporter to make sure no variants are missed. 
In order for the script to work, Input files must be ".tsv" from Ion Reporter.
