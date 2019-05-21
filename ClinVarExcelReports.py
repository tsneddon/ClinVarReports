from ExcelReportsFunctions import *

def main():

    inputFile1 = 'ClinVarVariationRelease_00-latest.xml.gz'
    inputFile2 = 'variation_allele.txt.gz'
    inputFile3 = 'variant_summary.txt.gz'
    inputFile4 = 'submission_summary.txt.gz'

    arg = sys.argv[1]
    #Options are 'ZeroStar', 'OneStar'

    dir = 'ClinVar' + arg + 'Reports'

    get_file(inputFile1, '/pub/clinvar/xml/clinvar_variation/')
    get_file(inputFile2, '/pub/clinvar/tab_delimited/')
    get_file(inputFile3, '/pub/clinvar/tab_delimited/')
    date = get_file(inputFile4, '/pub/clinvar/tab_delimited/')

    create_orgDict(inputFile1)
    create_a2vHash(inputFile2)
    create_HGVSHash(inputFile3)
    create_scvHash(inputFile4, arg)

    ExcelDir = make_directory(dir, date, arg)

    excelFile = arg + 'Report_' + date + '.xlsx'
    distFile = '_' + arg + 'DistributionReport_' + date + '.xlsx'

    statFile = '_' + arg + 'ReportsStats_' + date + '.xlsx'

    create_files(ExcelDir, excelFile, date, statFile, arg)
    create_distFile(ExcelDir, distFile, date, arg)

main()
