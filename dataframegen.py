import sys
from Bio import SeqIO
import pandas as pd
import re
from datetime import datetime
import time
from tqdm import tqdm
from collections import Counter
import itertools
import argparse
import numpy as np
import check_aa_mut
from check_aa_mut import getAminoacidMut as am
import distutils
from distutils import util
import multiprocessing
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser(description='Generates mutation table sorted by month, along with amino acid mutations and genes involved')
parser.add_argument('--aln', '--alignment', metavar = 'Alignment.fasta', required=True, help="Inputs alignment sequence (.fasta) with reference genome as first sequence")
parser.add_argument('--name', metavar = 'Name_1.feather and Name_2.feather', required=True, help="Outputs table as .feather format")
args = vars(parser.parse_args())

testfile = args['aln']
lineage = args['name']
minall = 0
maxinit = 100
minchange = 0 

ignore_list = set(["-", "R", "Y", "W", "S", "M", "K", "H", "B", "D", "V", "N", "X"])
mutation1 = []
datelist = []
sequencelist = []
     
print("WARNING: Sequence(s) identical to the reference genome OR with no date in format YYYY-MM-DD will be removed\nStage 1 of 3: Finding mutations for each sequence")
sys.stdout.flush()

for ref in SeqIO.parse(testfile, 'fasta'):
    refseq_record = ref.upper()
    refseq = refseq_record.seq
    break

#Loop 1: Get date, name and mutations of each sequence
    
def loop1(seq_record):
    if str(seq_record.upper()) != str(refseq_record) and re.compile(r'20[0-9][0-9]-(0[1-9]|1[012])-(0[1-9]|[12][0-9]|3[01])').search(str(seq_record)):
        date_match = re.search(r'20[0-9][0-9]-(0[1-9]|1[012])-(0[1-9]|[12][0-9]|3[01])', seq_record.description)
        date = datetime.strptime(date_match.group(), "%Y-%m-%d").date()
        mutation1objects = []
        queryseq = seq_record.seq.upper()
        for i in range(0, len(refseq)):
            nt1, nt2 = refseq[i], queryseq[i]
            if nt2 != nt1 and nt2 not in ignore_list:
                mut = (nt1 + str(i+1) + nt2)
                mutation1objects.append(mut)
        return [date, seq_record.description, mutation1objects] 
    
with Pool(multiprocessing.cpu_count()-1) as p:
    datelist, sequencelist, mutation1 = [list(a) for a in zip(*[x for x in list(tqdm(p.imap(loop1, itertools.islice(SeqIO.parse(testfile, 'fasta'), 1, None), chunksize = 3), position=0, leave=True)) if x is not None])]

mutationdata = {"Date": datelist, "Mutations": mutation1}
mutationdf = pd.DataFrame(mutationdata, index = sequencelist)

mutationunique = list(set(itertools.chain.from_iterable(mutation1)))
mutationuniquenums = [re.findall('\d+',temp1) for temp1 in mutationunique] # extracts numbers from strings
mutationnumsint = [int(*n) for n in mutationuniquenums] # returns 0 for the empty list corresponding to the word
mutationlist = [x for y, x in sorted(zip(mutationnumsint, mutationunique))] # sorts mutationunique based on the sorting of nums2

mutationdf['Year'] = mutationdf['Date'].apply(lambda x: str(x.year))
mutationdf['Month'] = mutationdf['Date'].apply(lambda x: str(x.month))
mutationdf['Date'] = pd.to_datetime(mutationdf['Date'], errors='coerce')
mutationdf['Week'] = mutationdf['Date'].dt.to_period('W').apply(lambda r: r.start_time)
mutationdf['Period'] = mutationdf['Year'] + "-" + mutationdf['Month']

tmp = mutationdf.groupby(["Week"])

periodlist = pd.to_datetime(list((mutationdf['Week'].unique()))).sort_values().strftime('%Y-%m-%d')
sequencetotalbyperiod = []
percentagelistlist = []
totallist = []

print("Stage 2 of 3: Generating mutation table")
sys.stdout.flush()

for period in tqdm(range(len(periodlist)), position=0, leave=True):
    periodsequences = tmp.get_group(periodlist[period])
    percentagelist=[]
    groupmutationtotal = list(itertools.chain.from_iterable(list(periodsequences['Mutations'])))
    groupmutationcount = Counter(groupmutationtotal)
    groupmutationcounts = {}
    for asdf in mutationlist:
        if asdf not in groupmutationcount.keys():
            groupmutationcount[asdf] = 0
        groupmutationcounts[asdf] = groupmutationcount.get(asdf)
    percentagelist = [x * 100/ len(periodsequences) for x in groupmutationcounts.values()]
    #percentagelist = [round(num, 3) for num in percentagelist] # for rounding to 3 sig. fig.
    percentagelistlist.append(percentagelist)
    totallist.append(len(periodsequences))
    
mutationtablebytimeraw = pd.DataFrame(percentagelistlist, index = periodlist, columns = mutationlist).transpose()


print("Stage 3 of 3: Filtering out mutations with prevalence rate below variable% and change in prevalence below variable %, and months with sample size below variable")
sys.stdout.flush()

mutationtablebytimefiltered = pd.DataFrame(index = mutationlist)
totallist1 = []
periodlist1 = []

for totals in range(len(totallist)):
    if totallist[totals] >= 10:
    #if totallist[totals] >= float(args['n']):
        mutationtablebytimefiltered = mutationtablebytimefiltered.join(mutationtablebytimeraw.iloc[:,totals])
        totallist1.append(totallist[totals])
        periodlist1.append(periodlist[totals])

def loop2(mutation2):
    aamut = am.nucl2aa(mutation2)
    for b in aamut:
        if "orf1a" not in b[0] and "orf1b" not in b[0] and max(set(mutationtablebytimefiltered.loc[mutation2])) >= minall and mutationtablebytimefiltered.loc[mutation2][-1] - mutationtablebytimefiltered.loc[mutation2][0] >= minchange and mutationtablebytimefiltered.loc[mutation2][0] <= maxinit:
            initlistmid = mutationtablebytimefiltered.loc[mutation2][0]
            percentagelistlistmid = list(mutationtablebytimefiltered.loc[mutation2])
            mutationlistmid = mutationlist[mutationlist.index(mutation2)]
            maxprevalencelistmid = max(set(mutationtablebytimefiltered.loc[mutation2]))
            maxprevalencechangelistmid = mutationtablebytimefiltered.loc[mutation2][-1] - mutationtablebytimefiltered.loc[mutation2][0]
            maxprevalencechangealllistmid = max(set(mutationtablebytimefiltered.loc[mutation2])) - min(set(mutationtablebytimefiltered.loc[mutation2]))
            return [b[0], b[2], initlistmid, percentagelistlistmid, mutationlistmid, maxprevalencechangealllistmid, maxprevalencelistmid, maxprevalencechangelistmid]

with Pool(multiprocessing.cpu_count()-1) as p:
    genelistmid, aamutationlistmid, initlistmid, percentagelistlistmid, mutationlistmid, maxprevalencechangealllistmid, maxprevalencelistmid, maxprevalencechangelistmid = [d for d in zip(*[x for x in list(tqdm(p.imap(loop2, tqdm(mutationlist, position=0, leave=True), chunksize = 50), position=0, leave=True)) if x is not None])]
            
aamutationlistmid = [none.replace('-', 'None') for none in aamutationlistmid]

genelist = []
aamutationlist = []
mutationlist1 = []
percentagelistlist1 =  []
maxprevalencelist = []
maxprevalencechangelist = []
maxprevalencechangealllist = []
initlist = []

#Remove synonymous mutations
print("Removing synonymous mutations...")

for searchna in range(len(aamutationlistmid)):
    if aamutationlistmid[searchna] != 'None':
        genelist.append(genelistmid[searchna])
        aamutationlist.append(aamutationlistmid[searchna])
        mutationlist1.append(mutationlistmid[searchna])
        percentagelistlist1.append(percentagelistlistmid[searchna])
        maxprevalencelist.append(maxprevalencelistmid[searchna])
        maxprevalencechangelist.append(maxprevalencechangelistmid[searchna])
        maxprevalencechangealllist.append(maxprevalencechangealllistmid[searchna])
        initlist.append(initlistmid[searchna])
        
mutlabellist = [s + ' (' for s in mutationlist1]
aamutationlabellist = [': '+ s + ')' for s in aamutationlist]
labellist = [list(a) for a in zip(mutlabellist, genelist, aamutationlabellist)]
labellist = [''.join([str(elem) for elem in labelsublist]) for labelsublist in labellist]

#Create dataframs for Dash app
pxdf1 = pd.DataFrame()
pxdf1['Mutations'] = mutationlist1
pxdf1['Labels'] = labellist
pxdf1['Maximum prevalence (%)'] = [float(a) for a in maxprevalencelist]
pxdf1['Change in prevalence (Fin - Start)(%)'] = [float(a) for a in maxprevalencechangelist]
pxdf1['Change in prevalence (All)(%)'] = [float(a) for a in maxprevalencechangealllist]
pxdf1['Initial prevalence'] = [float(a) for a in initlist]
pxdf1['Percentage by period'] = [[float(b) for b in a] for a in percentagelistlist1]
pxdf1['Gene'] = genelist
pxdf1['AA Mutations'] = aamutationlist
pxdf1['AA Label'] = pxdf1['Gene'] + ": " + pxdf1['AA Mutations']

pxdf2 = pd.DataFrame()
pxdf2['Totals'] = totallist1
pxdf2['Periods'] = periodlist1
pxdf2['Label'] = 'Total Number of Sequences'

pxdf1.to_feather(lineage + '_1.feather')
pxdf2.to_feather(lineage + '_2.feather')

print("Done")