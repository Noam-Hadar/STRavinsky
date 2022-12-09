from os import system
import twobitreader
from time import sleep

chromosomes = [str(i) for i in range(1,23)] + ['X','Y', 'M']
for reference in ('hg19', 'hg38', 'CHM13'):
    for chromosome in chromosomes:
        if reference in ('hg38', 'CHM13'):
            continue
        limit = int(len(twobitreader.TwoBitFile('2bits/' + reference + '.2bit')['chr' + chromosome])/1000000)
        for i in range(1, limit):
            system('python3 STgenMil.py ' + reference + ' ' + chromosome + ' ' + str(i))