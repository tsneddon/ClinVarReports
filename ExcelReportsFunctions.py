import xml.etree.ElementTree as ET
from ftplib import FTP
import os
import sys
import datetime
import time
import csv
import gzip
import re
import pprint
import xlsxwriter

orgDict = {}
scvHash = {}
a2vHash = {}
HGVSHash = {}
EPHash = {}
subList = {}
today = datetime.datetime.today().strftime('%Y%m%d') #todays date YYYYMMDD


def get_file(file, path):
    '''This function gets ClinVar files from FTP'''

    domain = 'ftp.ncbi.nih.gov'
    user = 'anonymous'
    password = 'tsneddon@broadinstitute.org'

    ftp = FTP(domain)
    ftp.login(user, password)
    ftp.cwd(path)
    localfile = open(file, 'wb')
    ftp.retrbinary('RETR ' + file, localfile.write)
    raw_date = ftp.sendcmd('MDTM ' + file)
    date = datetime.datetime.strptime(raw_date[4:], "%Y%m%d%H%M%S").strftime("%m-%d-%Y")
    ftp.quit()
    localfile.close()

    return(date)


def make_directory(dir, date, sub):
    '''This function makes a local directory for new files if directory does not already exist'''

    directory = dir + '/' + sub + '_Reports_' + date

    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        sys.exit('Program terminated, ' + directory + ' already exists.')

    return(directory)


def convert_date(date):
    '''This function converts a ClinVar date eg May 02, 2018 -> YYYYMMDD'''

    mon2num = dict(Jan='01', Feb='02', Mar='03', Apr='04', May='05', Jun='06',\
                   Jul='07', Aug='08', Sep='09', Oct='10', Nov='11', Dec='12')

    if '-' not in date:
        newDate = re.split(', | ',date)
        newMonth = mon2num[newDate[0]]
        convertDate = (newDate[2] + newMonth + newDate[1]) #YYYYMMDD, an integer for date comparisons
    else:
        convertDate = date

    return(convertDate)


def print_date(date):
    '''This function converts a date eg YYYYMMDD -> YYYY-MM-DD'''

    printDate = date[0:4] + "-" + date[4:6] + "-" + date[6:8] #MM/DD/YYYY, for printing to file
    return(printDate)


def create_orgDict(gzfile):
    '''This function makes a dictionary from the ClinVar production XML'''

    with gzip.open(gzfile) as input:
        for event, elem in ET.iterparse(input):

            if elem.tag == 'VariationArchive':
                for ClinAss in elem.iter(tag='ClinicalAssertion'):
                    for ClinAcc in ClinAss.iter(tag='ClinVarAccession'):
                        orgID = int(ClinAcc.attrib['OrgID'])
                        accession = ClinAcc.attrib['Accession']

                        if accession not in orgDict:
                            orgDict[accession] = orgID

                elem.clear()

    input.close()
    os.remove(gzfile)
    return(orgDict)


def create_scvHash(gzfile, arg):
    '''This function makes a hash of each SCV in each VarID'''

    global subList
    excludeList = {}

    with gzip.open(gzfile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    varID = int(col[0])
                    clinSig = col[1]
                    rawDate = col[2]
                    dateLastEval = convert_date(rawDate) #convert date eg May 02, 2018 -> YYYYMMDD

                    conditionList = []
                    rawConditionList = col[5].split(';')
                    for item in rawConditionList:
                        if ':'in item:
                            item = item.split(':', 1)[1]
                        conditionList.append(item)

                    condition = '; '.join(sorted(set(conditionList)))

                    revStat = col[6]
                    colMeth = col[7]

                    submitter = col[9]
                    submitter = submitter.rstrip()
                    submitter = re.sub('[^0-9a-zA-Z]+', '_', submitter)

                    submitter = submitter[0:45]

                    SCV = col[10]
                    accession = col[10].split('.', 1)[0]

                    if accession in orgDict:
                        orgID = orgDict[accession]
                    else:
                        orgID = 'None'

                    if (revStat == 'reviewed by expert panel' or revStat == 'practice guideline') and 'PharmGKB' not in submitter: #-- to exclude PharmGKB records
                        EPHash[varID] = {'ClinSig':clinSig, 'Submitter':submitter, 'DateLastEval':dateLastEval, 'OrgID':orgID}

                    if arg == 'OneStar' and submitter not in subList and revStat == 'criteria provided, single submitter' and 'clinical testing' in colMeth:
                        subList[submitter] = orgID

                    if arg == 'ZeroStar':
                        if revStat == 'no assertion criteria provided' and submitter not in subList:
                            subList[submitter] = orgID

                        if revStat == 'criteria provided, single submitter' and 'clinical testing' in colMeth and submitter not in excludeList:
                            excludeList[submitter] = orgID

                    if varID not in scvHash.keys():
                        scvHash[varID] = {}

                    scvHash[varID][SCV] = {'ClinSig':clinSig, 'DateLastEval':dateLastEval, 'Submitter':submitter, 'ReviewStatus':revStat, 'ColMeth':colMeth, 'Condition':condition, 'OrgID':orgID}

    if arg == 'ZeroStar':
        for key in excludeList:
            if key in subList:
                del subList[key]

    input.close()
    os.remove(gzfile)
    return(scvHash, EPHash, subList)


def create_a2vHash(gzfile):
    '''This function makes a dictionary of VarID to AlleleID'''

    with gzip.open(gzfile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()
            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    varID = int(col[0])
                    type = col[1]
                    alleleID = int(col[2])

                    #Ignore rows that are not Variant (simple type)
                    #This excludes Haplotype, CompoundHeterozygote, Complex, Phase unknown, Distinct chromosomes
                    if type == 'Variant':
                        a2vHash[alleleID] = varID

    input.close()
    os.remove(gzfile)
    return(a2vHash)


def create_HGVSHash(gzfile):
    '''This function makes a hash of metadata for each VarID'''

    with gzip.open(gzfile, 'rt') as input:
        line = input.readline()

        while line:
            line = input.readline()

            if not line.startswith('#'): #ignore lines that start with #
                col = re.split(r'\t', line) #split on tabs
                if not col[0] == '': #ignore empty lines
                    alleleID = int(col[0])
                    type = col[1]
                    HGVSname = col[2]
                    geneSym = col[4]
                    phenotype = col[13]
                    guidelines = col[26]

                    if alleleID in a2vHash:
                        HGVSHash[a2vHash[alleleID]] = {'VarType':type, 'HGVSname':HGVSname, 'GeneSym':geneSym,'Phenotype':phenotype,'Guidelines':guidelines}

    input.close()
    os.remove(gzfile)
    return(HGVSHash)


def create_files(ExcelDir, excelFile, date, statFile, arg):
    '''This function creates an Excel file for the stats of each submitter'''

    dir = ExcelDir

    stat_output_file = dir + '/' + statFile

    workbookStat = xlsxwriter.Workbook(stat_output_file)

    worksheetStat = workbookStat.add_worksheet('README')
    worksheetStat.write(0, 0, 'Tab#')
    worksheetStat.write(0, 1, 'Description')
    worksheetStat.write(1, 0, '1. Outlier_PLPvsVLBB')
    worksheetStat.write(1, 1, 'ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is discrepant from the majority (>= 2/3) of 1-star clinical submitters.')
    worksheetStat.write(2, 0, '2. Consensus_PLPvsVLBB')
    worksheetStat.write(2, 1, 'ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is in the majority (>= 2/3) with 1-star clinical submitters.')
    worksheetStat.write(3, 0, '3. NoConsensus_PLPvsVLBB')
    worksheetStat.write(3, 1, 'ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] does not have a majority (>= 2/3) over 1-star clinical submitters.')
    worksheetStat.write(4, 0, '4. VUSvsLBB')
    worksheetStat.write(4, 1, 'ClinVar variants where the submitter clinical significance is [VUS] vs [LB/B] and where there are no 1-star clinical submitter [P/LP] variants.')
    worksheetStat.write(5, 0, '5. IntraLab_discrepancy')
    worksheetStat.write(5, 1, 'ClinVar variants where the submitter has a discrepant clinical significance [P] vs [LP] vs [VUS] vs [LB] vs [B] with themselves.')
    worksheetStat.write(6, 0, '6. Lab_vs_EP')
    worksheetStat.write(6, 1, 'ClinVar variants where the submitter clinical significance is discrepant from an Expert Panel (EP) or Practice Guideline.')

    if arg == 'ZeroStar':
        worksheetStat0 = workbookStat.add_worksheet('0StarStats')
        worksheetStat0.write(0, 0, '0-Star submitter')

    if arg == 'OneStar':
        worksheetStat0 = workbookStat.add_worksheet('1StarStats')
        worksheetStat0.write(0, 0, '1-Star submitter')

    worksheetStat0.write(0, 1, 'OrgID')
    worksheetStat0.write(0, 2, '1. Outlier_PLPvsVLBB')
    worksheetStat0.write(0, 3, '2. Consensus_PLPvsVLBB')
    worksheetStat0.write(0, 4, '3. NoConsensus_PLPvsVLBB')
    worksheetStat0.write(0, 5, '4. VUSvsLBB')
    worksheetStat0.write(0, 6, '5. IntraLab_discrepancy')
    worksheetStat0.write(0, 7, '6. Lab_vs_EP')

    create_tabs(ExcelDir, excelFile, date, workbookStat, worksheetStat0, arg)


def create_tabs(ExcelDir, excelFile, date, workbookStat, worksheetStat0, arg):
    '''This function creates an Excel file for each sub in the subList'''

    dir = ExcelDir
    count = 0

    for sub in subList:

        count += 1
        worksheetStat0.write(count, 0, sub)
        worksheetStat0.write(count, 1, subList[sub])

        sub_output_file = dir + '/' + sub + '[' + str(subList[sub]) + ']_' + excelFile

        workbook = xlsxwriter.Workbook(sub_output_file)
        worksheet0 = workbook.add_worksheet('README')

        worksheet0.write(0, 0, "Date of ClinVar FTP file: " + date)
        worksheet0.write(2, 0, "Clinical submitter: " + sub)
        worksheet0.write(4, 0, "This Excel file is the output of a script that takes the most recent submission_summary.txt file from the ClinVar FTP site and outputs all the variants that are discrepant for " + sub)
        worksheet0.write(5, 0, 'Each tab is the result of a different set of parameters as outlined below:')
        worksheet0.write(6, 0, '#Variants:')
        worksheet0.write(7, 1, '1. Outlier_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is discrepant from the majority (>= 2/3) of 1-star clinical submitters.')
        worksheet0.write(8, 1, '2. Consensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] is in the majority (>= 2/3) with 1-star clinical submitters.')
        worksheet0.write(9, 1, '3. NoConsensus_PLPvsVLBB: ClinVar variants where the submitter clinical significance [P/LP] vs [VUS/LB/B] does not have a majority (>= 2/3) over 1-star clinical submitters.')
        worksheet0.write(10, 1, '4. VUSvsLBB: ClinVar variants where the submitter clinical significance is [VUS] vs [LB/B] compared to 1-star submitters and where there are no 1-star clinical submitter [P/LP] variants.')
        worksheet0.write(11, 1, '5. IntraLab_discrepancy: ClinVar variants where the submitter has a discrepant clinical significance [P] vs [LP] vs [VUS] vs [LB] vs [B] with themselves.')
        worksheet0.write(12, 1, '6. Lab_vs_EP: ClinVar variants where the submitter clinical significance is discrepant from an Expert Panel (EP) or Practice Guideline.')
        worksheet0.write(14, 0, 'Note: Tab classification counts are for unique submissions only i.e. if the same variant is submitted twice as Pathogenic by the same submitter, it will only be counted once')
        worksheet0.write(15, 0, 'Note: A variant can occur in multiple tabs i.e. if the same variant is submitted twice, once as Pathogenic and once as Benign by the same submitter, the variant could be both an outlier and the consensus')

        tabList = [create_tab1, create_tab2, create_tab3, create_tab4, create_tab5, create_tab6]
        for tab in tabList:
            tab(sub, workbook, worksheet0, worksheetStat0, count, arg)

        workbook.close()

    workbookStat.close()


def create_tab1(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#1 (Outlier_PLPvsVLBB) in the Excel file'''

    worksheet1 = workbook.add_worksheet('1.Outlier_PLPvsVLBB')

    tab = 1
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = outlier_medsig(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet1, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet1, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 7, 0, row)
    print_stats2file(worksheetStat0, count, 2, row)


def create_tab2(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#2 (Consensus_PLPvsVLBB) in the Excel file'''

    worksheet2 = workbook.add_worksheet('2.Consensus_PLPvsVLBB')

    tab = 2
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = consensus_medsig(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet2, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet2, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 8, 0, row)
    print_stats2file(worksheetStat0, count, 3, row)


def create_tab3(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#3 (NoConsensus_PLPvsVLBB) in the Excel file'''

    worksheet3 = workbook.add_worksheet('3.NoConsensus_PLPvsVLBB')

    tab = 3
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = noConsensus_medsig(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet3, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet3, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 9, 0, row)
    print_stats2file(worksheetStat0, count, 4, row)


def create_tab4(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#4 (VUSvsLBB) in the Excel file'''

    worksheet4 = workbook.add_worksheet('4.VUSvsLBB')

    tab = 4
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = VUSvsLBB(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet4, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet4, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 10, 0, row)
    print_stats2file(worksheetStat0, count, 5, row)


def create_tab5(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#5 (IntraLab_discrepancy) in the Excel file'''

    worksheet5 = workbook.add_worksheet('5.IntraLab_discrepancy')

    tab = 5
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = IntraLab_discrepancy(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet5, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet5, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 11, 0, row)
    print_stats2file(worksheetStat0, count, 6, row)


def create_tab6(sub, workbook, worksheet0, worksheetStat0, count, arg):
    '''This function creates the Tab#6 (Lab_vs_EP) in the Excel file'''

    worksheet6 = workbook.add_worksheet('6.Lab_vs_EP')

    tab = 6
    row = 0
    p2fileVarIDs = {}
    headerSubs = []

    for varID in scvHash:
        p2fileVarIDs, headerSubs = Outlier_EP(varID, sub, headerSubs, p2fileVarIDs, arg)

    print_header(sub, p2fileVarIDs, headerSubs, worksheet6, tab)

    for varID in p2fileVarIDs:
        varSubs = get_varSubs(sub, varID)
        row = print_variants(sub, worksheet6, row, varID, headerSubs, varSubs, p2fileVarIDs, tab)

    print_stats(worksheet0, 12, 0, row)
    print_stats2file(worksheetStat0, count, 7, row)


def outlier_medsig(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the outlier submitters in a medically significant VarID'''

    subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

    consensus = ''

    if ('PLP' in subSignificance and vlbb/(vlbb+plp) >= 2/3) or ('VLBB' in subSignificance and plp/(vlbb+plp) >= 2/3):
        if varID not in p2fileVarIDs.keys():
            p2fileVarIDs[varID] = {}

        p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other}

        if submitters:
            headerSubs.extend(submitters)

    headerSubs = sorted(set(headerSubs))
    if sub in headerSubs:
        headerSubs.remove(sub)

    if 'PLP' in subSignificance and vlbb/(vlbb+plp) >= 2/3:
        if vus == 0 and lb == 0 and b != 0:
            consensus = 'B'

        if vus == 0 and lb != 0 and b == 0:
            consensus = 'LB'

        if vus == 0 and lb != 0 and b != 0:
            consensus = 'LB/B'

        if vus != 0 and lb == 0 and b == 0:
            consensus = 'VUS'

        if vus != 0 and lb != 0 and b != 0:
            consensus = 'VUS/LB/B'

        if vus != 0 and lb != 0 and b == 0:
            consensus = 'VUS/LB'

        if vus != 0 and lb == 0 and b != 0:
            consensus = 'VUS/B'

    if 'VLBB' in subSignificance and plp/(vlbb+plp) >= 2/3:
        if p != 0 and lp == 0:
            consensus = 'P'

        if p == 0 and lp != 0:
            consensus = 'LP'

        if p != 0 and lp != 0:
            consensus = 'P/LP'

    if consensus != '':
        p2fileVarIDs[varID].update({'Consensus':consensus})

    return(p2fileVarIDs, headerSubs)


def consensus_medsig(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the consensus submitters in a medically significant VarID'''

    subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

    consensus = ''

    if ('PLP' in subSignificance and vlbb != 0 and plp/(vlbb+plp) >= 2/3) or ('VLBB' in subSignificance and plp != 0 and vlbb/(vlbb+plp) >= 2/3):
        if varID not in p2fileVarIDs.keys():
            p2fileVarIDs[varID] = {}

        p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other}

        if submitters:
            headerSubs.extend(submitters)

    headerSubs = sorted(set(headerSubs))
    if sub in headerSubs:
        headerSubs.remove(sub)

    if (('PLP' in subSignificance) and vlbb != 0 and plp/(vlbb+plp) >= 2/3):
        if p != 0 and lp == 0:
            consensus = 'P'

        if p == 0 and lp != 0:
            consensus = 'LP'

        if p != 0 and lp != 0:
            consensus = 'P/LP'

    if ('VLBB' in subSignificance and plp != 0 and vlbb/(vlbb+plp) >= 2/3):
        if vus != 0 and lb == 0 and b == 0:
            consensus = 'VUS'

        if vus == 0 and lb != 0 and b == 0:
            consensus = 'LB'

        if vus == 0 and lb == 0 and b != 0:
            consensus = 'B'

        if vus == 0 and lb != 0 and b != 0:
            consensus = 'LB/B'

        if vus != 0 and lb != 0 and b == 0:
            consensus = 'VUS/LB'

        if vus != 0 and lb == 0 and b != 0:
            consensus = 'VUS/B'

        if vus != 0 and lb != 0 and b != 0:
            consensus = 'VUS/LB/B'

    if consensus != '':
        p2fileVarIDs[varID].update({'Consensus':consensus})

    return(p2fileVarIDs, headerSubs)


def noConsensus_medsig(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the submitters with no consensus in a medically significant VarID'''

    subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

    all = ''

    if ('PLP' in subSignificance or 'VLBB' in subSignificance) and plp != 0 and vlbb != 0 and plp/(vlbb+plp) < 2/3 and vlbb/(vlbb+plp) < 2/3:

        if varID not in p2fileVarIDs.keys():
            p2fileVarIDs[varID] = {}

        p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other}

        if submitters:
            headerSubs.extend(submitters)

        if p != 0 and lp == 0 and vus != 0 and lbb == 0:
            all = 'P vs VUS'

        if p != 0 and lp == 0 and vus == 0 and lb == 0 and b != 0:
            all = 'P vs B'

        if p != 0 and lp == 0 and vus == 0 and lb != 0 and b == 0:
            all = 'P vs LB'

        if p != 0 and lp == 0 and vus != 0 and lb != 0 and b == 0:
            all = 'P vs VUS/LB'

        if p != 0 and lp == 0 and vus != 0 and lb == 0 and b != 0:
            all = 'P vs VUS/B'

        if p != 0 and lp == 0 and vus != 0 and lb != 0 and b != 0:
            all = 'P vs VUS/LB/B'

        if p != 0 and lp == 0 and vus == 0 and lb != 0 and b != 0:
            all = 'P vs LB/B'

        if p == 0 and lp != 0 and vus != 0 and lb == 0 and b == 0:
            all = 'LP vs VUS'

        if p == 0 and lp != 0 and vus == 0 and lb != 0 and b == 0:
            all = 'LP vs LB'

        if p == 0 and lp != 0 and vus == 0 and lb == 0 and b != 0:
            all = 'LP vs B'

        if p == 0 and lp != 0 and vus != 0 and lb != 0 and b == 0:
            all = 'LP vs VUS/LB'

        if p == 0 and lp != 0 and vus != 0 and lb == 0 and b != 0:
            all = 'LP vs VUS/B'

        if p == 0 and lp != 0 and vus != 0 and lb != 0 and b != 0:
            all = 'LP vs VUS/LB/B'

        if p == 0 and lp != 0 and vus == 0 and lb != 0 and b != 0:
            all = 'LP vs LB/B'

        if p != 0 and lp != 0 and vus != 0 and lb == 0 and b == 0:
            all = 'P/LP vs VUS'

        if p != 0 and lp != 0 and vus == 0 and lb != 0 and b == 0:
            all = 'P/LP vs LB'

        if p != 0 and lp != 0 and vus == 0 and lb == 0 and b != 0:
            all = 'P/LP vs B'

        if p != 0 and lp != 0 and vus != 0 and lb != 0 and b == 0:
            all = 'P/LP vs VUS/LB'

        if p != 0 and lp != 0 and vus != 0 and lb == 0 and b != 0:
            all = 'P/LP vs VUS/B'

        if p != 0 and lp != 0 and vus != 0 and lb != 0 and b != 0:
            all = 'P/LP vs VUS/LB/B'

        if p != 0 and lp != 0 and vus == 0 and lb != 0 and b != 0:
            all = 'P/LP vs LB/B'

        if plp == 0 and vus != 0 and lb == 0 and b != 0:
            all = 'VUS vs B'

        if plp == 0 and vus != 0 and lb != 0 and b == 0:
            all = 'VUS vs LB'

        if plp == 0 and vus != 0 and lb != 0 and b != 0:
            all = 'VUS vs LB/B'

        if all != '':
            p2fileVarIDs[varID].update({'All_significances':all})

    headerSubs = sorted(set(headerSubs))
    if sub in headerSubs:
        headerSubs.remove(sub)

    return(p2fileVarIDs, headerSubs)


def VUSvsLBB(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the submitters with VUS vs LBB significant VarID'''

    subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

    consensus = ''

    if ('VLBB' in subSignificance) and plp == 0 and vus != 0 and lbb != 0:

        if varID not in p2fileVarIDs.keys():
            p2fileVarIDs[varID] = {}

        p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other}

        if submitters:
            headerSubs.extend(submitters)

        if vus/vlbb >= 2/3:
            consensus = 'VUS'

        if lb == 0 and b != 0 and b/vlbb >= 2/3:
            consensus = 'B'

        if lb != 0 and b == 0 and lb/vlbb >= 2/3:
            consensus = 'LB'

        if lb != 0 and b != 0 and lbb/vlbb >= 2/3:
            consensus = 'LB/B'

        if consensus != '':
            p2fileVarIDs[varID].update({'Consensus':consensus})
        else:
            p2fileVarIDs[varID].update({'Consensus':'N/A'})

    headerSubs = sorted(set(headerSubs))

    return(p2fileVarIDs, headerSubs)


def IntraLab_discrepancy(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the submitters with a [P/LP] vs [VUS] vs [LB/B] discrepant clinical significances'''

    subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

    duplicate_subs, conflict = get_duplicates(sub, varID)

    if sub in duplicate_subs and len(conflict) > 1:
        if varID not in p2fileVarIDs.keys():
            p2fileVarIDs[varID] = {}

        p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other, 'Conflict':conflict}

        if submitters:
            headerSubs.extend(submitters)

        headerSubs = sorted(set(headerSubs))
        if sub in headerSubs:
            headerSubs.remove(sub)

    return(p2fileVarIDs, headerSubs)


def Outlier_EP(varID, sub, headerSubs, p2fileVarIDs, arg):
    '''This function returns the submitters where the clinical significance is discrepant from an Expert Panel or Practice Guideline'''

    clinSig = ''
    EPconflict = ''

    if varID in EPHash:
        for SCV in scvHash[varID]:
            if scvHash[varID][SCV]['Submitter'] == sub and (scvHash[varID][SCV]['ClinSig'] == 'Pathogenic' or scvHash[varID][SCV]['ClinSig'] == 'Likely pathogenic' or \
               scvHash[varID][SCV]['ClinSig'] == 'Uncertain significance' or scvHash[varID][SCV]['ClinSig'] == 'Likely benign' or scvHash[varID][SCV]['ClinSig'] == 'Benign'):
                clinSig = scvHash[varID][SCV]['ClinSig']

                if ((EPHash[varID]['ClinSig'] == 'Pathogenic' or EPHash[varID]['ClinSig'] == 'Likely pathogenic') and \
                   (scvHash[varID][SCV]['ClinSig'] == 'Uncertain significance' or scvHash[varID][SCV]['ClinSig'] == 'Likely benign' or \
                   scvHash[varID][SCV]['ClinSig'] == 'Benign')) or ((EPHash[varID]['ClinSig'] == 'Uncertain significance' or EPHash[varID]['ClinSig'] == 'Likely benign' or \
                   EPHash[varID]['ClinSig'] == 'Benign') and (scvHash[varID][SCV]['ClinSig'] == 'Pathogenic' or scvHash[varID][SCV]['ClinSig'] == 'Likely pathogenic')):
                    EPconflict = 'P/LP vs VUS/LB/B'
                if ((EPHash[varID]['ClinSig'] == 'Uncertain significance' and (scvHash[varID][SCV]['ClinSig'] == 'Likely benign' or \
                   scvHash[varID][SCV]['ClinSig'] == 'Benign')) or ((EPHash[varID]['ClinSig'] == 'Likely benign' or EPHash[varID]['ClinSig'] == 'Benign') and \
                   scvHash[varID][SCV]['ClinSig'] == 'Uncertain significance')):
                    EPconflict = 'VUS vs LB/B'
                if (EPHash[varID]['ClinSig'] == 'Pathogenic' and scvHash[varID][SCV]['ClinSig'] == 'Likely pathogenic') or (scvHash[varID][SCV]['ClinSig'] == 'Pathogenic' and \
                   EPHash[varID]['ClinSig'] == 'Likely pathogenic'):
                    EPconflict = 'P vs LP'
                if ((EPHash[varID]['ClinSig'] == 'Benign' and scvHash[varID][SCV]['ClinSig'] == 'Likely benign') or (EPHash[varID]['ClinSig'] == 'Likely benign' and \
                   scvHash[varID][SCV]['ClinSig'] == 'Benign')):
                    EPconflict = 'B vs LB'

        if clinSig != '' and clinSig != EPHash[varID]['ClinSig']:
            if EPconflict != '':
                if varID not in p2fileVarIDs.keys():
                    p2fileVarIDs[varID] = {}

                subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other = get_pathCounts(sub, varID, arg)

                p2fileVarIDs[varID] = {'Total':total, 'PLP':plp, 'VUS':vus, 'LBB':lbb, 'Misc':other, 'EP':EPHash[varID]['Submitter'], 'EP_clinSig':EPHash[varID]['ClinSig']}
                p2fileVarIDs[varID].update({'EPConflict':EPconflict})

                if submitters:
                    headerSubs.extend(submitters)

            else:
                p2fileVarIDs[varID].update({'EPConflict':'N/A'})

    headerSubs = sorted(set(headerSubs))
    if sub in headerSubs:
        headerSubs.remove(sub)

    return(p2fileVarIDs, headerSubs)


def get_pathCounts(sub, varID, arg):
    '''This function returns the counts of ACMG pathogenicities for each VarID'''

    submitters = []
    p = 0
    lp = 0
    vus = 0
    lb = 0
    b = 0
    other = 0
    subSignificance = []
    unique_subs = []

    for SCV in scvHash[varID]:

        if scvHash[varID][SCV]['ClinSig'] == 'Pathogenic':
           #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] == sub and current_sub not in unique_subs:
                if arg == 'OneStar' or (arg == 'ZeroStar' and 'no assertion criteria provided' in scvHash[varID][SCV]['ReviewStatus']):
                    unique_subs.append(current_sub)
                    submitters.append(scvHash[varID][SCV]['Submitter'])
                    subSignificance.append('PLP')
                    p += 1

            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                p += 1


        if scvHash[varID][SCV]['ClinSig'] == 'Likely pathogenic':
           #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] == sub and current_sub not in unique_subs:
                if arg == 'OneStar' or (arg == 'ZeroStar' and 'no assertion criteria provided' in scvHash[varID][SCV]['ReviewStatus']):
                    unique_subs.append(current_sub)
                    submitters.append(scvHash[varID][SCV]['Submitter'])
                    subSignificance.append('PLP')
                    lp += 1

            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                lp += 1

        if scvHash[varID][SCV]['ClinSig'] == 'Uncertain significance':
           #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] == sub and current_sub not in unique_subs:
                if arg == 'OneStar' or (arg == 'ZeroStar' and 'no assertion criteria provided' in scvHash[varID][SCV]['ReviewStatus']):
                    unique_subs.append(current_sub)
                    submitters.append(scvHash[varID][SCV]['Submitter'])
                    subSignificance.append('VLBB')
                    vus += 1

            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                vus += 1

        if scvHash[varID][SCV]['ClinSig'] == 'Likely benign':
           #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] == sub and current_sub not in unique_subs:
                if arg == 'OneStar' or (arg == 'ZeroStar' and 'no assertion criteria provided' in scvHash[varID][SCV]['ReviewStatus']):
                    unique_subs.append(current_sub)
                    submitters.append(scvHash[varID][SCV]['Submitter'])
                    subSignificance.append('VLBB')
                    lb += 1

            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                lb += 1

        if scvHash[varID][SCV]['ClinSig'] == 'Benign':
           #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] == sub and current_sub not in unique_subs:
                if arg == 'OneStar' or (arg == 'ZeroStar' and 'no assertion criteria provided' in scvHash[varID][SCV]['ReviewStatus']):
                    unique_subs.append(current_sub)
                    submitters.append(scvHash[varID][SCV]['Submitter'])
                    subSignificance.append('VLBB')
                    b += 1

            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                b += 1

        else:
            #Don't double count (Illumina's) duplicate submissions!!!
            current_sub = scvHash[varID][SCV]['Submitter'] + scvHash[varID][SCV]['ClinSig']
            if scvHash[varID][SCV]['Submitter'] != sub and current_sub not in unique_subs and scvHash[varID][SCV]['ReviewStatus'] == 'criteria provided, single submitter' and 'clinical testing' in scvHash[varID][SCV]['ColMeth']:
                unique_subs.append(current_sub)
                submitters.append(scvHash[varID][SCV]['Submitter'])
                other += 1

    plp = p+lp
    lbb = lb+b
    vlbb = vus+lb+b
    total = plp+vus+lbb+other

    subSignificance = sorted(set(subSignificance))

    return(subSignificance, submitters, p, lp, plp, vus, lb, b, lbb, vlbb, total, other)


def get_duplicates(sub, varID):
    '''This function returns the duplicate submitters for each VarID'''

    count = 0
    clinSigs = []
    duplicate_subs = []
    conflict = []

    for SCV in scvHash[varID]:

        if scvHash[varID][SCV]['Submitter'] == sub:
            if scvHash[varID][SCV]['ClinSig'] not in clinSigs:
                if scvHash[varID][SCV]['ClinSig'] == 'Pathogenic' or scvHash[varID][SCV]['ClinSig'] == 'Likely pathogenic' or \
                   scvHash[varID][SCV]['ClinSig'] == 'Uncertain significance' or \
                   scvHash[varID][SCV]['ClinSig'] == 'Likely benign' or scvHash[varID][SCV]['ClinSig'] == 'Benign':
                    clinSigs.append(scvHash[varID][SCV]['ClinSig'])
                    count += 1

    if count > 1:
        if 'Pathogenic' in clinSigs:
            conflict.append('P')
        if 'Likely pathogenic' in clinSigs:
            conflict.append('LP')
        if 'Uncertain significance' in clinSigs:
            conflict.append('VUS')
        if 'Likely benign' in clinSigs:
            conflict.append('LB')
        if 'Benign' in clinSigs:
            conflict.append('B')

        duplicate_subs.append(sub)

    return(duplicate_subs, conflict)


def get_varSubs(sub, varID):
    '''This function returns the list of variant submitters'''

    varSubs = []
    if varID in scvHash:
        for SCV in scvHash[varID]:
            if scvHash[varID][SCV]['Submitter'] != sub:
                if scvHash[varID][SCV]['DateLastEval'] != '-':
                    #Convert date from YYYYMMDD -> YYYY-MM-DD
                    subPrintDate = print_date(scvHash[varID][SCV]['DateLastEval'])
                    varSubs.append(scvHash[varID][SCV]['Submitter'] + ' [' + scvHash[varID][SCV]['ClinSig'] + ' (' + subPrintDate + ')]')
                else:
                    varSubs.append(scvHash[varID][SCV]['Submitter'] + ' [' + scvHash[varID][SCV]['ClinSig'] + ' (No DLE)]')

    varSubs = sorted(set(varSubs))

    return(varSubs)


def print_header(sub, p2fileVarIDs, headerSubs, worksheet, tab):
    '''This function prints all the header titles to the Excel tabs'''

    k = 0
    if p2fileVarIDs != {}:
        worksheet.write(0, k, 'VarID')
        k+=1
        worksheet.write(0, k, 'Gene_symbol')
        k+=1
        worksheet.write(0, k, 'All_conditions')
        k+=1
        worksheet.write(0, k, 'HGVS_name')
        k+=1

        if tab == 6:
            worksheet.write(0, k, 'EP')
            k+=1
            worksheet.write(0, k, 'EP_significance')
            k+=1
            worksheet.write(0, k, 'EP_conflict')
            k+=1
        else:
            worksheet.write(0, k, 'EP or Practice Guideline')
            k+=1

            if tab == 3:
                worksheet.write(0, k, 'All_significances')
            elif tab == 5:
                worksheet.write(0, k, 'Intralab_conflicts')
            else:
                worksheet.write(0, k, 'Consensus_significance')
            k+=1

        worksheet.write(0, k, sub + '_significance(s)')
        k+=1
        worksheet.write(0, k, sub + '_SCVid(s)')
        k+=1
        worksheet.write(0, k, sub + '_ReviewStatus(s)')
        k+=1
        worksheet.write(0, k, sub + '_most_recent_DateLastEvaluated')
        k+=1
        worksheet.write(0, k, sub + '_condition(s)')
        k+=1

        for head in headerSubs:
            if head != sub:
                worksheet.write(0, k, head)
                k+=1

        worksheet.write(0, k, 'Total_submissions')
        k+=1
        worksheet.write(0, k, 'Total_PLP')
        k+=1
        worksheet.write(0, k, 'Total_VUS')
        k+=1
        worksheet.write(0, k, 'Total_LBB')
        k+=1
        worksheet.write(0, k, 'Total_Misc')
    else:
        worksheet.write(0, 0, 'No variants found')


def print_variants(sub, worksheet, row, varID, headerSubs, varSubs, p2fileVarIDs, tab):
    '''This function prints all the variants to the Excel tabs'''

    row += 1
    k = 0
    worksheet.write(row, k, varID)
    k+=1
    if HGVSHash[varID]['GeneSym']:
        worksheet.write(row, k, HGVSHash[varID]['GeneSym'])
    k+=1

    if HGVSHash[varID]['Phenotype']:
        worksheet.write(row, k, HGVSHash[varID]['Phenotype'])
    k+=1

    if HGVSHash[varID]['HGVSname']:
        worksheet.write(row, k, HGVSHash[varID]['HGVSname'])
    k+=1

    if tab != 6:
       if varID in EPHash.keys():
           worksheet.write(row, k, EPHash[varID]['Submitter'] + ' (' + EPHash[varID]['ClinSig'] + ')')
       else:
           worksheet.write(row, k, 'N/A')
       k+=1

    if tab == 3:
        worksheet.write(row, k, p2fileVarIDs[varID]['All_significances'])
        k+=1
    elif tab == 5:
        sublist = ' vs '. join(p2fileVarIDs[varID]['Conflict'])
        worksheet.write(row, k, sublist)
        k+=1
    elif tab == 6:
        worksheet.write(row, k, EPHash[varID]['Submitter'])
        k+=1
        if EPHash[varID]['DateLastEval'] != '-':
            #Convert date from YYYYMMDD -> YYYY-MM-DD
            subPrintDate = print_date(EPHash[varID]['DateLastEval'])
            worksheet.write(row, k, EPHash[varID]['ClinSig'] + ' (' + subPrintDate + ')')
        else:
            worksheet.write(row, k, EPHash[varID]['ClinSig'] + ' (No DLE)')
        k+=1
        worksheet.write(row, k, p2fileVarIDs[varID]['EPConflict'])
        k+=1
    else:
        worksheet.write(row, k, p2fileVarIDs[varID]['Consensus'])
        k+=1

    clinSig = []
    scvs= []
    revStats = []
    dle = []
    conditions = []

    for scv in scvHash[varID]:
        if (tab == 1 and scvHash[varID][scv]['Submitter'] == sub and \
           (((scvHash[varID][scv]['ClinSig'] == 'Pathogenic' or scvHash[varID][scv]['ClinSig'] == 'Likely pathogenic') and any(x in p2fileVarIDs[varID]['Consensus'] for x in ['VUS','LB','B','LB/B','VUS/LB/B','VUS/LB','VUS/B'])) or \
           ((scvHash[varID][scv]['ClinSig'] == 'Uncertain significance' or scvHash[varID][scv]['ClinSig'] == 'Likely benign' or scvHash[varID][scv]['ClinSig'] == 'Benign') and \
           any(x in p2fileVarIDs[varID]['Consensus'] for x in ['P','LP','P/LP'])))) or \
           (tab == 2 and scvHash[varID][scv]['Submitter'] == sub and \
           (((scvHash[varID][scv]['ClinSig'] == 'Pathogenic' or scvHash[varID][scv]['ClinSig'] == 'Likely pathogenic') and any(x in p2fileVarIDs[varID]['Consensus'] for x in ['P','LP','P/LP'])) or \
           ((scvHash[varID][scv]['ClinSig'] == 'Uncertain significance' or scvHash[varID][scv]['ClinSig'] == 'Likely benign' or scvHash[varID][scv]['ClinSig'] == 'Benign') and \
           any(x in p2fileVarIDs[varID]['Consensus'] for x in ['VUS','LB','B','LB/B','VUS/LB/B','VUS/LB','VUS/B'])))) or \
           ((tab != 1 or tab != 2) and scvHash[varID][scv]['Submitter'] == sub and \
           (scvHash[varID][scv]['ClinSig'] == 'Pathogenic' or scvHash[varID][scv]['ClinSig'] == 'Likely pathogenic' or \
           scvHash[varID][scv]['ClinSig'] == 'Uncertain significance' or scvHash[varID][scv]['ClinSig'] == 'Likely benign' or scvHash[varID][scv]['ClinSig'] == 'Benign')):
            clinSig.append(scvHash[varID][scv]['ClinSig'])
            scvs.append(scv)
            revStats.append(scvHash[varID][scv]['ReviewStatus'])
            dle.append(print_date(scvHash[varID][scv]['DateLastEval']))
            conditions.append(scvHash[varID][scv]['Condition'])

    clinSigList = ' | '. join(sorted(set(clinSig)))
    scvList = ' | '. join(sorted(set(scvs)))
    revList = ' | '. join(sorted(set(revStats)))
    dleList = sorted(set(dle))
    dle = dleList[-1]
    conditionsList = ' | '. join(sorted(set(conditions)))

    worksheet.write(row, k, clinSigList)
    k+=1
    worksheet.write(row, k, scvList)
    k+=1
    worksheet.write(row, k, revList)
    k+=1
    worksheet.write(row, k, dle)
    k+=1
    worksheet.write(row, k, conditionsList)
    k+=1

    for headerSub in headerSubs:
        p2file = 'no'
        for varSub in varSubs:
            varSubList = varSub.split(' [')
            varSubStr = varSubList[0]
            if headerSub == varSubStr:
                p2file = varSub[varSub.find("[")+1:varSub.find("]")]
        if p2file != 'no':
            worksheet.write(row, k, p2file)
            k += 1
        else:
            k += 1

    if varID in p2fileVarIDs:
       worksheet.write(row, k, p2fileVarIDs[varID]['Total'])
       k+=1
       worksheet.write(row, k, p2fileVarIDs[varID]['PLP'])
       k+=1
       worksheet.write(row, k, p2fileVarIDs[varID]['VUS'])
       k+=1
       worksheet.write(row, k, p2fileVarIDs[varID]['LBB'])
       k+=1
       worksheet.write(row, k, p2fileVarIDs[varID]['Misc'])

    return(row)


def print_stats(worksheet0, line, column, row):
    '''This function prints the total variant count to the README Excel tab'''

    worksheet0.write(line, column, row)


def print_stats2file(worksheetStat0, count, column, row):
    '''This function prints the total variant counts to the Stats file'''

    worksheetStat0.write(count, column, row)
