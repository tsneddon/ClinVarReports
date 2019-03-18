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
