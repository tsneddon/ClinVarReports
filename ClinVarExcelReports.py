#Must run as >python3 ClinVarExcelReports.py arg
#Arg options are 'ZeroStar', 'OneStar'
#The script outputs an Excel file(with the last modfified date of the FTP file) for each submitter:
#1. Outlier_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is discrepant from the majority (>= 2/3) of 1-star clinical submitters.
#2. Consensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is in the majority (>= 2/3) with 1-star clinical submitters.
#3. NoConsensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] does not have a majority (>= 2/3) over 1-star clinical submitters.
#4. VUSvsLBB: ClinVar variants where the submitter clinical significance is [VUS] vs [LB/B] and where there are no 1-star clinical submitter [P/LP] variants.
#5. IntraLab_discrepancy: ClinVar variants where the submitter has a discrepant clinical significance [P] vs [LP] vs [VUS] vs [LB] vs [B] with themselves.
#6. Lab_vs_EP: ClinVar variants where the submitter clinical significance is discrepant from an EP (or Practice Guideline).

from ExcelReportsFunctions import *

def main():

    inputFile1 = 'organization_summary.txt'
    inputFile2 = 'submission_summary.txt.gz'
    inputFile3 = 'variation_allele.txt.gz'
    inputFile4 = 'variant_summary.txt.gz'

    arg = sys.argv[1]
    #Options are 'ZeroStar', 'OneStar'

    dir = 'ClinVar' + arg + 'Reports'

    get_file(inputFile1, 'pub/clinvar/xml/')
    date = get_file(inputFile2, '/pub/clinvar/tab_delimited/')
    get_file(inputFile3, '/pub/clinvar/tab_delimited/')
    get_file(inputFile4, '/pub/clinvar/tab_delimited/')

    ExcelDir = make_directory(dir, date, arg)

    excelFile = arg + 'Report_' + date + '.xlsx'

    statFile = '_' + arg + 'ReportsStats_' + date + '.xlsx'

    create_orgDict(inputFile1)
    create_scvHash(inputFile2, arg)
    create_a2vHash(inputFile3)
    create_HGVSHash(inputFile4)

    create_files(ExcelDir, excelFile, date, statFile, arg)

main()
