#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 17:53:54 2021

@author: Elisavet Iliopoulou
"""
import pandas as pd
import os
from multiprocessing import Pool
import numpy as np


qCNEs_dict = {}
qCNEs_dict_sort = {}
zCNEs_sort = pd.read_csv('DanRer11_Astyanax2_CNE_headers.tsv' , " ", names=['zCNEs'], engine='python').sort_values(by=['zCNEs']).reset_index(drop=True)
zCNEs_sort['ID'] = 'zCNE' + zCNEs_sort.index.astype(str)

for filename in os.listdir('CNEs_directory'):
    qCNEs_dict[filename] = pd.read_csv('CNEs_directory/%s' % filename, " ", header=None, engine='python')


                                # Reverse negatives

    qCNEs_dict[filename]['%s_Reverse' % filename] = pd.Series([])
    j=0
    while j < len(qCNEs_dict[filename])-1:

        end = qCNEs_dict[filename].iloc[j][5]
        start = qCNEs_dict[filename].iloc[j][4]

        if int(qCNEs_dict[filename].iloc[j][5]) < int(qCNEs_dict[filename].iloc[j][4]):
            qCNEs_dict[filename]['%s_Reverse' % filename][j] = '-'
            qCNEs_dict[filename][4][j] = qCNEs_dict[filename].iloc[j][5]
            qCNEs_dict[filename][5][j] = start
        else:
            qCNEs_dict[filename]['%s_Reverse' % filename][j] = '+'

        j+=1

                            #  Choose tCNE based on e-value

    qCNEs_dict_sort[filename] = qCNEs_dict[filename].sort_values(by=[1,4,5]).reset_index(drop=True)
    x = qCNEs_dict_sort[filename]
    check_old = len(x)
    check_new = 0

    while check_old != check_new:
        check_old = check_new

        j=0
        while j < len(x) - 1:

            first_target = x.iloc[j]
            second_target = x.iloc[j+1]

            if first_target[1] == second_target[1] and first_target[5] > second_target[4] and first_target[4] <= second_target[4]:
                if first_target[6] >= second_target[6]:
                    x.drop([j], axis=0, inplace=True)
                else:
                    x.drop([j+1], axis=0, inplace=True)

                x = x.reset_index(drop=True)
            j+=1

        check_new = len(x)
    qCNEs_dict_sort[filename] = x

    x_qSort = x.sort_values(by=[0]).reset_index(drop=True)



                        #  Sort by query and Match tCNE - qCNE

    j=0
    i=0
    k=len(x_qSort)
    zCNEs_sort['%s_chrom' % filename] = pd.Series([])
    zCNEs_sort['%s_start' % filename] = pd.Series([], dtype='Int64')
    zCNEs_sort['%s_end' % filename] = pd.Series([], dtype='Int64')
    zCNEs_sort['%s_qStart' % filename] = pd.Series([], dtype='Int64')
    zCNEs_sort['%s_qEnd' % filename] = pd.Series([], dtype='Int64')
    zCNEs_sort['%s_Reverse' % filename] = pd.Series([])

    while j < k:

        if str(x_qSort.iloc[j][0]) == str(zCNEs_sort['zCNEs'][i]):

            zCNEs_sort['%s_chrom' % filename][i] = x_qSort.iloc[j][1]
            zCNEs_sort['%s_start' % filename][i] = x_qSort.iloc[j][4]
            zCNEs_sort['%s_end' % filename][i] = x_qSort.iloc[j][5]
            zCNEs_sort['%s_qStart' % filename][i] = x_qSort.iloc[j][2]
            zCNEs_sort['%s_qEnd' % filename][i] = x_qSort.iloc[j][3]
            zCNEs_sort['%s_Reverse' % filename][i] = x_qSort.iloc[j]['%s_Reverse' % filename]

            j+=1
            i+=1


        else:
            i+=1



                                 #  SuperMatrix presence

    superMatrix_presence = zCNEs_sort.dropna(inplace=False).reset_index(drop=True)



                                #  Remove overlapping qCNEs


    check_old = len(superMatrix_presence)
    check_new = 0
    while check_old != check_new:
        check_old = check_new
        cne = 0
        while cne < len(superMatrix_presence)-1:
            superMatrix_presence = superMatrix_presence.reset_index(drop=True)

            qChromosome = superMatrix_presence.iloc[cne]['zCNEs'].split(':')[0]
            qChromosome_sec = superMatrix_presence.iloc[cne+1]['zCNEs'].split(':')[0]

            first_query_end = int(superMatrix_presence.iloc[cne]['zCNEs'].split(':')[1].split('-')[1])
            first_query_start = int(superMatrix_presence.iloc[cne]['zCNEs'].split(':')[1].split('-')[0])
            sec_query_start = int(superMatrix_presence.iloc[cne+1]['zCNEs'].split(':')[1].split('-')[0])
            sec_query_end = int(superMatrix_presence.iloc[cne+1]['zCNEs'].split(':')[1].split('-')[1])


            dif1 = abs(first_query_start - first_query_end)
            dif2 = abs(sec_query_start - sec_query_end)

            if qChromosome == qChromosome_sec and first_query_end > sec_query_start and first_query_start <= sec_query_start:
                if dif1 >= dif2:
                    superMatrix_presence.drop([cne+1], axis=0, inplace=True)
                else:
                    superMatrix_presence.drop([cne], axis=0, inplace=True)
            cne += 1

        check_new = len(superMatrix_presence)


                                 #  Extend Target


for filename in os.listdir('CNEs_directory'):

    print('About to extend %s...' % filename)

    superMatrix_presence = superMatrix_presence.sort_values(by=['%s_chrom' % filename, '%s_start' % filename, '%s_end' % filename]).reset_index(drop=True)

    i = 0
    while i < len(superMatrix_presence):

        qLength = 1 + int(superMatrix_presence.iloc[i]['zCNEs'].split(':')[1].split('-')[1]) - int(superMatrix_presence.iloc[i]['zCNEs'].split(':')[1].split('-')[0])
        k = qLength - superMatrix_presence['%s_qEnd' % filename][i]



                            # Extend all ... NO MATTER WHAT

        superMatrix_presence['%s_end' % filename][i] += (k + 50)
        superMatrix_presence['%s_start' % filename][i] -= superMatrix_presence['%s_qStart' % filename][i]
        superMatrix_presence['%s_start' % filename][i] -= 50

        i+=1


                                # Now... overlaps matter...

for filename in os.listdir('CNEs_directory'):
    superMatrix_presence = superMatrix_presence.sort_values(by=['%s_chrom' % filename, '%s_start' % filename, '%s_end' % filename]).reset_index(drop=True)

    i = 0
    while i < len(superMatrix_presence)-1:

        if superMatrix_presence['%s_chrom' % filename][i] == superMatrix_presence['%s_chrom' % filename][i+1]:
            if superMatrix_presence['%s_end' % filename][i] >= superMatrix_presence['%s_start' % filename][i+1]:
                print(superMatrix_presence['%s_chrom' % filename][i], superMatrix_presence['%s_start' % filename][i], superMatrix_presence['%s_end' % filename][i],'<---- HERE IT IS')
                superMatrix_presence['%s_end' % filename][i] = superMatrix_presence['%s_start' % filename][i+1] - 1


        i+=1

                            #Adjust negatives after extending

    print('Adjust negatives...')

    i = 0
    while i < len(superMatrix_presence):
        if superMatrix_presence['%s_start' % filename][i] < 0:
            print('s')
            superMatrix_presence['%s_start' % filename][i] = 1

        i+=1

superMatrix_presence[['QChrom']] = superMatrix_presence['zCNEs'].str.split(':',expand=True)[0]
superMatrix_presence[['Q_coord_Start', 'Q_coord_End']]=superMatrix_presence['zCNEs'].str.split(':',expand=True)[1].str.split('-',expand=True)
superMatrix_presence['Q_coord_End'] = superMatrix_presence['Q_coord_End'].astype(int)
superMatrix_presence['Q_coord_Start'] = superMatrix_presence['Q_coord_Start'].astype(int)

superMatrix_presence = superMatrix_presence.sort_values(by=['QChrom', 'Q_coord_Start', 'Q_coord_End']).reset_index(drop=True)


                        # Teleost_CNEs (before extension) for outgroup BLAST

superMatrix_presence[['QChrom_initial']] = superMatrix_presence['QChrom']
superMatrix_presence[['Q_coord_Start_initial']] =superMatrix_presence['Q_coord_Start'].astype(int)
superMatrix_presence[['Q_coord_End_initial']] = superMatrix_presence['Q_coord_End'] .astype(int)





print('About to extend query...')

                                        #  Extend Query

i = 0
while i < len(superMatrix_presence):

    superMatrix_presence['Q_coord_End'][i] += 50
    superMatrix_presence['Q_coord_Start'][i] -= 50

    i+=1

i = 0
while i < len(superMatrix_presence)-1:

    if superMatrix_presence['QChrom'][i] == superMatrix_presence['QChrom'][i+1]:
        if superMatrix_presence['Q_coord_End'][i] > superMatrix_presence['Q_coord_Start'][i+1]:
            superMatrix_presence['Q_coord_End'][i] = superMatrix_presence['Q_coord_Start'][i+1] - 1

    i+=1



print('Finishing... Writing files...')

superMatrix_presence = superMatrix_presence.sort_values(by=['ID']).reset_index(drop=True)

for filename in os.listdir('CNEs_directory'):
    superMatrix_presence[['%s_chrom' % filename, '%s_start' % filename, '%s_end' % filename, 'ID', 'ID', '%s_Reverse' % filename]].to_csv("zCNE_comparisons/%s.bed" % filename, index=False)

superMatrix_presence[['QChrom', 'Q_coord_Start', 'Q_coord_End', 'ID']].to_csv("zCNE_comparisons/danio.bed", index=False)
superMatrix_presence[['QChrom_initial', 'Q_coord_Start_initial', 'Q_coord_End_initial', 'ID']].to_csv("DanRer11.teleost_CNEs.bed", index=False)
