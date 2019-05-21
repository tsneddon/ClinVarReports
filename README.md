# ClinVarReports
Python scripts to generate ClinGen submitter reports from ClinVar FTP files.

## About this project
ClinVar outputs a submission_summary.txt.gz file containing a list of all SCV submissions per variant to its FTP site.
The scripts in this project use this file to generate the following files in the subdirectory ClinVarReports/:

**ClinVarExcelReports.py** - this script takes the argument 'ZeroStar' or 'OneStar' and outputs an Excel file containing each variant in the submission_summary.txt file that a ZeroStar or OneStar submitter needs to update or review. The Excel contains a README with summary stats and 6 structured tabs as detailed below:
  * \#1. Outlier_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is discrepant from the majority (>= 2/3) of 1-star clinical submitters.
  * \#2. Consensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is in the majority (>= 2/3) with 1-star clinical submitters.
  * \#3. NoConsensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] does not have a majority (>= 2/3) over 1-star clinical submitters.
  * \#4. VUSvsLBB: ClinVar variants where the submitter clinical significance is [VUS] vs [LB/B] and where there are no 1-star clinical submitter [P/LP] variants.
  * \#5. IntraLab_discrepancy: ClinVar variants where the submitter has a discrepant clinical significance [P] vs [LP] vs [VUS] vs [LB] vs [B] with themselves.
  * \#6. Lab_vs_EP: ClinVar variants where the submitter clinical significance is discrepant from an Expert Panel (EP) or Practice Guideline.

**ClinVarExcelReports.py** also generates a: 
  * ReportsStats Excel file containing the summary variant counts for each ZeroStar or OneStar submitter
  * Distribution Excel file containing all variants with 2 or more ZeroStar or OneStar submitters

## How to run these scripts
All scripts are run as 'python3 *filename.py* *arg*' where arg = 'ZeroStar' or 'OneStar'
All scripts use FTP to take the most recent ClinVar FTP files as input and to output the files with the date of the FTP submission_summary.txt.gz file appended:

  * ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
  * ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/variation_allele.txt.gz
  * ftp.ncbi.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
  * ftp.ncbi.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz

These ClinVar files are then removed when finished.
