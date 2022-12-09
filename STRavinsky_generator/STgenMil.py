#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from itertools import product
from re import finditer
import urllib.parse
import pandas as pd
import primer3
import re
import os
import sys
import twobitreader



# In[ ]:
reference = sys.argv[1]
chromosome = sys.argv[2]
million = sys.argv[3]

def getSequence(chromosome, position, reference, flank_size):
    start = int(position) - (flank_size * 1000000)
    if start < 1:
        start = 1
    end = int(position) + (flank_size * 1000000) + 1
    limit = len(twobitreader.TwoBitFile('2bits/' + reference + '.2bit')['chr' + chromosome])
    if end > limit:
        end = limit
    genome = twobitreader.TwoBitFile('2bits/' + reference + '.2bit')
    sequence = genome['chr' + str(chromosome)].get_slice(int(position) - (flank_size * 1000000), int(position) + (flank_size * 1000000) + 1).strip().lower()
    for n in range(2,11):
        if n == 2:
            minToReport = 4
        if n < 5:
            minToReport = 3
        else:
            minToReport = 2
        repeats = [''.join(i) for i in product(('a','c','t','g'), repeat = n)]
        repeats = [repeat for repeat in repeats if len(set(repeat)) != 1]
        indeces = []
        for repeat in repeats:
            indeces += [i.start() for i in re.finditer(repeat * minToReport, sequence)]
        for i in indeces:
            try:
                sequence = sequence[:i] + sequence[i : i + (len(repeat) * minToReport)].upper() + sequence[i + (len(repeat) * minToReport):]
            except:
                continue
    return sequence

def findRepetitiveSequence(STR):
    n = 2
    while n < 11:
        if STR.count(STR[:n]) > 2 and STR[:n] == STR[n: 2*n]:
            return STR[:n] + ' x ' + str(STR.count(STR[:n]))
        else:
            n += 1
    return 'Drop'

def getSequenceLength(sequence_repeats):
    try:
        values = sequence_repeats.split(' x ')
        return len(values[0]) * int(values[1])
    except:
        return 'Drop'


# In[ ]:


def PGTail(reference, chromosome, position, flank_size):
    sequence = getSequence(chromosome, position, reference, flank_size)
    indeces = [i.start() for i in re.finditer(r'[ACTG]', sequence)]
    if len(indeces) == 0:
        return pd.DataFrame()
    STRseqs = ['']
    STRindeces = [flank_size * 1000000]
    range_start = indeces[0]
    n = indeces[0]
    for i in indeces[1:]:
        if i == n + 1:
            n = i
        else:
            STRindeces.append(range_start)
            STRseqs.append(sequence[range_start:n])
            range_start = i
            n = i
    df = pd.DataFrame()
    df['Distance'] = STRindeces
    df['Distance'] = df['Distance'] - (flank_size * 1000000)
    df['Sequence_full'] = [i.strip() for i in STRseqs]
    df['Sequence'] = df['Sequence_full'].apply(findRepetitiveSequence)
    df = df[df['Sequence_full'] != 'Drop']
    df['STR length'] = df['Sequence'].apply(getSequenceLength)
    df = df[df['Sequence'] != 'Drop']
    df['Coordinates'] = df.apply(lambda x : str(chromosome) + ':' + str(int(x['Distance'] + position)) + '-' + str(int(x['Distance'] + position + x['STR length'])), axis = 1)
    df['Distance'] = df['Distance'].astype(int)
    df = df.sort_values('Distance')
    del df['Sequence_full']
    del df['Distance']
    return df


# In[ ]:


df = PGTail(reference, chromosome, int(million) * 1000000, 2)
df.to_csv('Results/' + reference + '/' + reference + '_chr' + chromosome + '_' + str(million) + '_STRs.csv', index = False)