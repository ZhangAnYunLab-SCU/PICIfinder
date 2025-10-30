"""
    PICI_ModelEval.py

    This file is used to test the performance of trained models and validate the PICI prediction pipeline.

    Note:
        1. The detailed model training process is implemented in model_train.py
        2. All functions in this file (except model training) have optimized versions available
           in PICI_Identifier.py

    Main functionalities:
        - Load pre-trained random forest models
        - Execute PICI identification pipeline
        - Output prediction results and performance evaluation
"""

import os
import sys
import pandas as pd
import numpy as np
import subprocess
import model_train as mt

def int_gene(seq_file_name,filepath):
    """
    Coordinates the PICI search process,
    handles integrase and regulator identification.

    Parameter:
    seq_file_name: Sequence file name prefix, used to construct input/output file paths.
    filepath: Working directory path, stores temporary files and result files.

    return: None
    """

    seq_name = 0
    stl_site = 0
    number = 0

    int_list = pd.read_table(seq_file_name + '_total_table.txt',header = 0)
    length_list = pd.read_table(seq_file_name + '_fasta_name.txt',header = 0)
    sequence = pd.read_table(seq_file_name + '.ffa',header = None)

    # Train and load models
    forest = mt.forest_tran(filepath)
    forest1 = mt.forest_tran1(filepath)

    for i in range(len(int_list)):
        cut_dataframe = {'name':[0],'start_x':[0],'end_x':[0],'gene_number':[0],'dir_x':[0],'seq_number_x':[0],'temp_x':[0],'structure_name':[0],'structure_id':[0],'E_value':[0],'gene':[0],'protein':[0],'gene_name':[0],'domain_start':[0],'domain_end':[0],'bac_name':[0],'phage_name':[0]}
        cut_dataframe = pd.DataFrame(cut_dataframe)
        cut_dataframe1 = {'name':[0],'start_x':[0],'end_x':[0],'gene_number':[0],'dir_x':[0],'seq_number_x':[0],'temp_x':[0],'structure_name':[0],'structure_id':[0],'E_value':[0],'gene':[0],'protein':[0],'gene_name':[0],'domain_start':[0],'domain_end':[0],'bac_name':[0],'phage_name':[0]}
        cut_dataframe1 = pd.DataFrame(cut_dataframe1)
        temp_cut_dataframe = {'name':[0],'start_x':[0],'end_x':[0],'gene_number':[0],'dir_x':[0],'seq_number_x':[0],'temp_x':[0],'structure_name':[0],'structure_id':[0],'E_value':[0],'gene':[0],'protein':[0],'gene_name':[0],'domain_start':[0],'domain_end':[0],'bac_name':[0],'phage_name':[0]}
        temp_cut_dataframe = pd.DataFrame(temp_cut_dataframe)

        HTH_flag = 0
        flag = 0
        flag1 = 0
        confirmed = 0
        true_start = 0
        true_end = 0

        if int(seq_name) != int(int_list.iloc[i,5]):
            start = 0
            end = 0
            int_start = 0
            stl_flag = 0
            stl_site = 0
            number = 0
            int_site = 0
            int_name = ''

        # Integrase-Dependent Integration
        if (int_list.iloc[i,7] == 'Phage_integrase%' or int_list.iloc[i,7] == 'Phage_int_SAM_3%') and 'e' in str(int_list.iloc[i,9]):
            if int_list.iloc[i,7] == 'Phage_int_SAM_3%' and seq_name == int_list.iloc[i,5] and int_list.iloc[i,3] == number:
                continue
            number = int_list.iloc[i,3]
            seq_name = int_list.iloc[i,5]
            int_list.iloc[i, 10] = 'int'
            for m in range(len(length_list)):
                if int(length_list.iloc[m,1]) == int(seq_name):
                    sequence_name = length_list.iloc[m,0]
                    length = length_list.iloc[m,2]
                    break

            # Integration Direction: -1 (indicating leftward integration from the attachment site)
            if int(int_list.iloc[i,4]) == -1:
                print('^^^^^^^^^^^^^^^')
                print('int_name_'+str(int_name))
                print('int_name_'+str(int_list.iloc[i, 5]))
                print('int_site_'+str(int_site))
                print('int_site_'+str(int_list.iloc[i, 3]))
                print('int_dire_'+str(int_list.iloc[i, 4]))

                if int_name == int_list.iloc[i, 5] and int(int_site) != 0 and int(int_list.iloc[i,3]) - int(int_site) < 5 and int(int_list.iloc[i,4]) == -1:
                    continue
                int_name = int(int_list.iloc[i,5])
                int_start = int(int_list.iloc[i,1])
                int_site = int(int_list.iloc[i,3])

                start = int(int_list.iloc[i,1])
                end = int(int_list.iloc[i,1]) + 80000
                temp_end = end
                print('int_name')
                print(int_name)
                print('int_site')
                print(int_site)

                for j in range(i,len(int_list)):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,1] >= start and int_list.iloc[j,2] <= end:
                        cut_dataframe = pd.concat([cut_dataframe, int_list[j:j+1]], axis=0)
                    if seq_name != int_list.iloc[j,5]:
                        break

                dire = int(int_list.iloc[i,4])
                cut_dataframe = cut_dataframe.reset_index(drop=True)
                cut_dataframe = cut_dataframe.drop([0])
                cut_dataframe = cut_dataframe_deal(cut_dataframe,sequence_name,dire)
                if dire == -1:
                    end = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                if dire == 1:
                    start = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                table = HTH_pri_judge(cut_dataframe,cut_dataframe1,seq_file_name,start,end,forest,forest1,length,1,dire,filepath)

                att_seq = 'no'
                if table.split('{')[4] != str(0):
                    PICI_feature_gene = table.split('{')[4].split('+')[-1]
                    for p in range(i-20,i):
                        if (int_list.iloc[p+1,7] != 'Cpn60_TCP1' and int_list.iloc[p,7] == 'Cpn10')  or int_list.iloc[p,7] == 'Cpn60_TCP1' or int_list.iloc[p,11] == '30S ribosomal protein S18' or int_list.iloc[p,11] == 'SsrA-binding protein' or 'MetQ' in int_list.iloc[p,11] or 'GMP synthase' in int_list.iloc[p,11] or 'N-acetyltransferase' in int_list.iloc[p,11] or ('HU' in int_list.iloc[p,11] and 'DNA-binding' in int_list.iloc[p,11]) or int_list.iloc[p,11] == 'putative glucose uptake protein GlcU' or int_list.iloc[p,11] == 'putative oxidoreductase':
                            if int_list.iloc[p,7] == 'Cpn60_TCP1' or int_list.iloc[p,7] == 'Cpn10':
                                name = 'GroEL'
                            if int_list.iloc[p,11] == '30S ribosomal protein S18':
                                name = '30S'
                            if int_list.iloc[p,11] == 'SsrA-binding protein':
                                name = 'SmpB'
                            if 'MetQ' in int_list.iloc[p,11]:
                                name = 'metq'
                            if 'GMP synthase' in int_list.iloc[p,11]:
                                name = 'GMP'
                            if 'N-acetyltransferase' in int_list.iloc[p,11]:
                                name = 'GNAT'
                            if 'HU' in int_list.iloc[p,11] and 'DNA-binding' in int_list.iloc[p,11]:
                                name = 'HU'
                            if int_list.iloc[p,11] == 'putative glucose uptake protein GlcU' or int_list.iloc[p,11] == 'putative oxidoreductase':
                                name = 'GlcU'
                            for q in range(len(sequence)):
                                if '>' in sequence.iloc[q,0] and sequence_name in sequence.iloc[q,0]:
                                    seq = sequence.iloc[q+1,0]
                                    int_seq = seq[int_list.iloc[p,1]:start]
                                    back_seq = seq[int(PICI_feature_gene):temp_end-60000]
                                    print('int_seq_'+str(int_list.iloc[p,1]))
                                    print('int_seq_'+str(start))
                                    print('back_seq_'+str(PICI_feature_gene))
                                    print('back_seq_'+str(temp_end-60000))
                                    PICI = att_blast(seq_file_name,sequence_name,int_seq,back_seq,start,temp_end,name,dire,filepath)
                                    if int(PICI.split('+')[0]) != 0 or int(PICI.split('+')[1]) != 0:
                                        confirmed = 1
                                        att_seq = PICI.split('+')[2]
                                    if int(PICI.split('+')[0]) != 0:
                                        true_start = int(PICI.split('+')[0])+int(int_list.iloc[p,1])
                                    if int(PICI.split('+')[1]) != 0:
                                        true_end = int(PICI.split('+')[1])+int(PICI_feature_gene)
                            break
                Information(seq_file_name,table,sequence_name,1,confirmed,true_start,true_end,start,att_seq,cut_dataframe,cut_dataframe1)

            # Integration Direction: 1 (indicating rightward integration from the attachment site)
            # Same as Direction -1
            if int(int_list.iloc[i,4]) == 1:
                for r in range(i,len(int_list)):
                    if seq_name != int_list.iloc[r,5]:
                        break
                    if int_name == int_list.iloc[i, 5] and int(int_list.iloc[i,3]) - int(int_site) < 5 and int(int_list.iloc[i,4]) == 1:
                        flag1 = 1
                if flag1 == 1:
                    continue

                start = int(int_list.iloc[i,2]) - 80000
                temp_start = start
                end = int(int_list.iloc[i,2])
                int_site = int(int_list.iloc[i,3])
                int_name = int(int_list.iloc[i,5])
                print('int_name')
                print(int_name)
                print('int_site')
                print(int_site)

                for j in range(i,-1,-1):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,1] >= start and int_list.iloc[j,2] <= end:
                        cut_dataframe = pd.concat([cut_dataframe, int_list[j:j+1]], axis=0)
                    if seq_name != int_list.iloc[j,5]:
                        break

                dire = int(int_list.iloc[i, 4])
                cut_dataframe = cut_dataframe.reset_index(drop=True)
                cut_dataframe = cut_dataframe.drop([0])
                cut_dataframe = cut_dataframe_deal(cut_dataframe,sequence_name,dire)

                if dire == -1:
                    end = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                if dire == 1:
                    start = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                table = HTH_pri_judge(cut_dataframe,cut_dataframe1,seq_file_name,start,end,forest,forest1,length,1,dire,filepath)

                att_seq = 'no'
                if table.split('{')[4] != str(0):
                    PICI_feature_gene = table.split('{')[4].split('+')[-1]
                    p_end = i + 20
                    if i + 20 >= len(int_list):
                        p_end = len(int_list) - 1
                    for p in range(int(p_end), i, -1):
                        if (int_list.iloc[p - 1, 7] != 'Cpn60_TCP1' and int_list.iloc[p, 7] == 'Cpn10') or int_list.iloc[p,7] == 'Cpn60_TCP1' or int_list.iloc[p,11] == '30S ribosomal protein S18' or int_list.iloc[p,11] == 'SsrA-binding protein' or 'MetQ' in int_list.iloc[p,11] or 'GMP synthase' in int_list.iloc[p,11] or 'N-acetyltransferase' in int_list.iloc[p,11] or ('HU' in int_list.iloc[p,11] and 'DNA-binding' in int_list.iloc[p,11]) or int_list.iloc[p,11] == 'putative glucose uptake protein GlcU':
                            if int_list.iloc[p,7] == 'Cpn60_TCP1' or int_list.iloc[p,7] == 'Cpn10':
                                name = 'GroEL'
                            if int_list.iloc[p,11] == '30S ribosomal protein S18':
                                name = '30S'
                            if int_list.iloc[p,11] == 'SsrA-binding protein':
                                name = 'SmpB'
                            if 'MetQ' in int_list.iloc[p,11]:
                                name = 'metq'
                            if 'GMP synthase' in int_list.iloc[p,11]:
                                name = 'GMP'
                            if 'N-acetyltransferase' in int_list.iloc[p,11]:
                                name = 'GNAT'
                            if 'HU' in int_list.iloc[p,11] and 'DNA-binding' in int_list.iloc[p,11]:
                                name = 'HU'
                            if int_list.iloc[p,11] == 'putative glucose uptake protein GlcU' or int_list.iloc[p,11] == 'putative oxidoreductase':
                                name = 'GlcU'

                            for q in range(len(sequence)):
                                if '>' in sequence.iloc[q, 0] and sequence_name in sequence.iloc[q, 0]:
                                    seq = sequence.iloc[q + 1, 0]
                                    int_seq = seq[end:int_list.iloc[p, 2]]
                                    back_seq = seq[temp_start+60000:int(PICI_feature_gene)]
                                    print('int_seq_' + str(end))
                                    print('int_seq_' + str(int_list.iloc[p, 2]))
                                    print('back_seq_' + str(temp_start+60000))
                                    print('back_seq_' + str(PICI_feature_gene))
                                    PICI = att_blast(seq_file_name, sequence_name, int_seq, back_seq, temp_start, end, name,dire,filepath)
                                    if int(PICI.split('+')[0]) != 0 or int(PICI.split('+')[1]) != 0:
                                        confirmed = 1
                                        att_seq = PICI.split('+')[2]
                                    if int(PICI.split('+')[0]) != 0:
                                        true_start = int(PICI.split('+')[0]) + int(temp_start+60000)
                                    if int(PICI.split('+')[1]) != 0:
                                        true_end = int(PICI.split('+')[1]) + int(end)

                            Information(seq_file_name,table,sequence_name,1,confirmed,true_start,true_end,start,att_seq,cut_dataframe,cut_dataframe1)
                            break

        # Regulator-Mediated Integration
        if (int_list.iloc[i,7] == 'HTH_3' or int_list.iloc[i,7] == 'HTH_19' or int_list.iloc[i,7] == 'HTH_26' or int_list.iloc[i,7] == 'HTH_31' or int_list.iloc[i,7] == 'HTH_17' or (int_list.iloc[i,10] == 'stl_gene' and '^' not in int_list.iloc[i,7])) and 'e' in str(int_list.iloc[i,9]):
            seq_name = int_list.iloc[i,5]
            # Get information of sequances
            for m in range(len(length_list)):
                if int(length_list.iloc[m,1]) == int(seq_name):
                    sequence_name = length_list.iloc[m,0]
                    length = length_list.iloc[m,2]
                    break

            # Direction -1
            if int_start == 0 and int(int_list.iloc[i,4]) == -1 and int(int_list.iloc[i,1]) - 15000 < 0 and stl_flag == 0 and int(int_list.iloc[i,3]) != int(stl_site):
                stl_site = int_list.iloc[i,3]
                for r in range(i+1,len(int_list)):
                    if seq_name != int_list.iloc[r,5]:
                        break
                    if (int_list.iloc[r,7] == 'Phage_integrase%' or int_list.iloc[r,7] == 'Phage_int_SAM_3%') and int(int_list.iloc[r,3]) - int(stl_site) < 20:
                        HTH_flag = 1
                if HTH_flag == 1:
                    continue

                stl_flag = 1
                start = int(int_list.iloc[i,1])
                end = int(int_list.iloc[i,1]) + 80000
                for j in range(i,len(int_list)):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,1] >= start and int_list.iloc[j,2] <= end:
                        cut_dataframe = pd.concat([cut_dataframe, int_list[j:j+1]], axis=0)
                    if seq_name != int_list.iloc[j,5]:
                        break
                if temp_cut_dataframe_deal(cut_dataframe,stl_site) == 0:
                    stl_flag = 0
                    continue
                for j in range(len(int_list)):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,2] <= end:
                        cut_dataframe1 = pd.concat([cut_dataframe1, int_list[j:j+1]], axis=0)

                dire = int(int_list.iloc[i,4])
                cut_dataframe = cut_dataframe.reset_index(drop=True)
                cut_dataframe = cut_dataframe.drop([0])
                cut_dataframe1 = cut_dataframe1.reset_index(drop=True)
                cut_dataframe1 = cut_dataframe1.drop([0])
                cut_dataframe = cut_dataframe_deal(cut_dataframe,sequence_name,dire)
                cut_dataframe1 = cut_dataframe_deal(cut_dataframe1,sequence_name,dire)
                if dire == -1:
                    end = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                if dire == 1:
                    start = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                table = HTH_pri_judge(cut_dataframe,cut_dataframe1,seq_file_name,start,end,forest,forest1,length,2,dire,filepath)
                att_seq = 'no'
                Information(seq_file_name,table,sequence_name,2,confirmed,true_start,true_end,start,att_seq,cut_dataframe,cut_dataframe1)

            # Direction 1 (same as Direction -1)
            if int(int_list.iloc[i,4]) == 1 and int(int_list.iloc[i,1]) + 15000 > length and int(int_list.iloc[i,3]) != int(stl_site):
                stl_site = int_list.iloc[i,3]
                for r in range(i+1,len(int_list)):
                    if seq_name != int_list.iloc[r,5]:
                        break
                    if (int_list.iloc[r,7] == 'Phage_integrase%' or int_list.iloc[r,7] == 'Phage_int_SAM_3%' or ((int_list.iloc[r,7] == 'HTH_3' or int_list.iloc[r,7] == 'HTH_19' or int_list.iloc[r,7] == 'HTH_26' or int_list.iloc[r,7] == 'HTH_31' or int_list.iloc[r,7] == 'HTH_17' or (int_list.iloc[r,10] == 'stl_gene' and '^' not in int_list.iloc[r,7])) and int(int_list.iloc[r,4]) == 1)) and int(int_list.iloc[r,3]) != int(stl_site):
                        start = int(int_list.iloc[r,2]) - 80000
                        end = int(int_list.iloc[r,2])
                        stl_site1 = int_list.iloc[r,3]
                        for j in range(r,-1,-1):
                            if seq_name == int_list.iloc[j,5] and int_list.iloc[j,1] >= start and int_list.iloc[j,2] <= end:
                                temp_cut_dataframe = pd.concat([temp_cut_dataframe, int_list[j:j+1]], axis=0)
                            if seq_name != int_list.iloc[j,5]:
                                break
                        if temp_cut_dataframe_deal(temp_cut_dataframe,stl_site1) == 1:
                            flag = 1
                if flag == 1:
                    continue
                start = int(int_list.iloc[i,2]) - 80000
                end = int(int_list.iloc[i,2])
                for j in range(i,-1,-1):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,1] >= start and int_list.iloc[j,2] <= end:
                        cut_dataframe = pd.concat([cut_dataframe, int_list[j:j+1]], axis=0)
                    if seq_name != int_list.iloc[j,5]:
                        break
                for j in range(len(int_list)-1,-1,-1):
                    if seq_name == int_list.iloc[j,5] and int_list.iloc[j,2] >= start:
                        cut_dataframe1 = pd.concat([cut_dataframe1, int_list[j:j+1]], axis=0)
                dire = int(int_list.iloc[i,4])
                cut_dataframe = cut_dataframe.reset_index(drop=True)
                cut_dataframe = cut_dataframe.drop([0])
                cut_dataframe1 = cut_dataframe1.reset_index(drop=True)
                cut_dataframe1 = cut_dataframe1.drop([0])
                cut_dataframe = cut_dataframe_deal(cut_dataframe,sequence_name,dire)
                cut_dataframe1 = cut_dataframe_deal(cut_dataframe1,sequence_name,dire)
                if dire == -1:
                    end = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                if dire == 1:
                    start = tRNA_end(cut_dataframe,seq_file_name,start,end,dire)
                table = HTH_pri_judge(cut_dataframe,cut_dataframe1,seq_file_name,start,end,forest,forest1,length,2,dire,filepath)
                att_seq = 'no'
                Information(seq_file_name,table,sequence_name,2,confirmed,true_start,true_end,start,att_seq,cut_dataframe,cut_dataframe1)


def temp_cut_dataframe_deal(cut_dataframe, stl_site):
    """
    Temporary dataframe processing, filtering duplicate and low-quality annotations.

    Parameters:
    param cut_dataframe: candidate region gene dataframe.
    param stl_site: position number of the regulatory factor gene.

    return: bool - 1 indicates valid regulatory factor present, 0 indicates invalid
    """

    flag = 0

    # Step 1：Data cleaning and deduplication
    for i in range(1, len(cut_dataframe)):
        for j in range(i + 1, len(cut_dataframe)):
            if cut_dataframe.iloc[i, 3] != cut_dataframe.iloc[j, 3]:
                break
            if cut_dataframe.iloc[j, 7] == 'na':
                continue
            if cut_dataframe.iloc[i, 5] == cut_dataframe.iloc[j, 5]:
                if (int(cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 14])) or (
                        int(cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 14]) and int(
                        cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 13])) or (
                        int(cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 14])) or (
                        int(cut_dataframe.iloc[i, 13]) >= int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 14]) <= int(cut_dataframe.iloc[j, 14])):
                    if float(cut_dataframe.iloc[i, 9].split('_')[0]) > float(cut_dataframe.iloc[j, 9].split('_')[0]):
                        cut_dataframe.iloc[i, 1] = 0
                    if float(cut_dataframe.iloc[i, 9].split('_')[0]) <= float(cut_dataframe.iloc[j, 9].split('_')[0]):
                        cut_dataframe.iloc[j, 1] = 0

        if 'e' not in str(cut_dataframe.iloc[i, 9]).split('_')[0] and cut_dataframe.iloc[i, 7] != 'na':
            cut_dataframe.iloc[i, 1] = 0
        if cut_dataframe.iloc[i, 9] != 'na' and 'e' not in str(cut_dataframe.iloc[i, 9]).split('_')[1]:
            cut_dataframe.iloc[i, 7] = 'na'
            cut_dataframe.iloc[i, 8] = 'na'

    # remove
    cut_dataframe = cut_dataframe[~cut_dataframe['start_x'].isin([0])]

    # Step 2：Validate regulatory factor presence
    for i in range(len(cut_dataframe)):
        if cut_dataframe.iloc[i, 3] == stl_site and (
                cut_dataframe.iloc[i, 7] == 'HTH_3' or cut_dataframe.iloc[i, 7] == 'HTH_19' or cut_dataframe.iloc[
            i, 7] == 'HTH_26' or cut_dataframe.iloc[i, 7] == 'HTH_31' or cut_dataframe.iloc[i, 7] == 'HTH_17' or (
                        cut_dataframe.iloc[i, 12] == 'stl_gene' and '^' not in cut_dataframe.iloc[i, 7]) or
                cut_dataframe.iloc[i, 7] == 'Phage_integrase%' or cut_dataframe.iloc[i, 7] == 'Phage_int_SAM_3%'):
            flag = 1

    if flag == 1:
        return 1
    else:
        return 0


def cut_dataframe_deal(cut_dataframe, sequence_name, dire):
    """
    DataFrame cleaning and standardization, unifying column names and formats

    Parameters:
        cut_dataframe: gene annotation dataframe
        sequence_name: current sequence name
        dire: gene direction (1 or -1)

    Returns: DataFrame - processed and standardized dataframe
    """

    t = 1
    cut_dataframe[['start_x']] = cut_dataframe[['start_x']].astype(int)
    cut_dataframe.to_csv(seq_file_name + '_cut_dataframe.txt', sep='\t', index=False)
    # Sort by start position for forward direction genes
    if dire == 1:
        cut_dataframe = cut_dataframe.sort_values(by=['start_x'], axis=0)

    # DataFrame information optimization
    for i in range(len(cut_dataframe)):
        cut_dataframe.iloc[i, 0] = sequence_name
        # Preferentially use gene name annotations
        if cut_dataframe.iloc[i, 12] != 'na':
            cut_dataframe.iloc[i, 10] = cut_dataframe.iloc[i, 12]

        if 'phage_resistance' in cut_dataframe.iloc[i, 10]:
            cut_dataframe.iloc[i, 7] = 'phage_resistance'
        # Add domain information for hypothetical proteins
        if cut_dataframe.iloc[i, 11] == 'hypothetical protein' and (
                cut_dataframe.iloc[i, 12] == 'tail_gene' or cut_dataframe.iloc[i, 12] == 'head_gene' or
                cut_dataframe.iloc[i, 12] == 'lysis_gene' or cut_dataframe.iloc[i, 12] == 'capsid_gene'):
            cut_dataframe.iloc[i, 11] = cut_dataframe.iloc[i, 11] + '+' + cut_dataframe.iloc[i, 12]
        # Virulence factor and antibiotic resistance gene annotations
        if '#VFDB' in cut_dataframe.iloc[i, 12] or '#SARG' in cut_dataframe.iloc[i, 12]:
            cut_dataframe.iloc[i, 11] = cut_dataframe.iloc[i, 12]

        # Label bacterial-origin and phage-origin genes
        if cut_dataframe.iloc[i, 10] != 'PICI' and cut_dataframe.iloc[i, 15] == 'bac' and cut_dataframe.iloc[
            i, 9] != 'na' and cut_dataframe.iloc[i, 12] != 'stl_gene':
            if 'e' in cut_dataframe.iloc[i, 9].split('_')[-1]:
                cut_dataframe.iloc[i, 10] = 'bac'
        if cut_dataframe.iloc[i, 16] == 'phage':
            cut_dataframe.iloc[i, 12] = 'phage'

        # Assign consecutive gene numbers
        if i == 0:
            cut_dataframe.iloc[i, 6] = t
        if i > 0 and cut_dataframe.iloc[i, 2] == cut_dataframe.iloc[i - 1, 2]:
            cut_dataframe.iloc[i, 6] = t
        if i > 0 and cut_dataframe.iloc[i, 2] != cut_dataframe.iloc[i - 1, 2]:
            t += 1
            cut_dataframe.iloc[i, 6] = t

        if cut_dataframe.iloc[i, 7] == 'na':
            continue

        # Overlapping domain annotations (keep best E-value)
        for j in range(i + 1, len(cut_dataframe)):
            if cut_dataframe.iloc[i, 3] != cut_dataframe.iloc[j, 3]:
                break
            if cut_dataframe.iloc[j, 7] == 'na':
                continue
            if cut_dataframe.iloc[i, 5] == cut_dataframe.iloc[j, 5]:
                if (int(cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 14])) or (
                        int(cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 14]) and int(
                        cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 13])) or (
                        int(cut_dataframe.iloc[i, 13]) < int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 14]) > int(cut_dataframe.iloc[j, 14])) or (
                        int(cut_dataframe.iloc[i, 13]) >= int(cut_dataframe.iloc[j, 13]) and int(
                        cut_dataframe.iloc[i, 14]) <= int(cut_dataframe.iloc[j, 14])):
                    if float(cut_dataframe.iloc[i, 9].split('_')[0]) > float(cut_dataframe.iloc[j, 9].split('_')[0]):
                        cut_dataframe.iloc[i, 1] = 0
                    if float(cut_dataframe.iloc[i, 9].split('_')[0]) <= float(cut_dataframe.iloc[j, 9].split('_')[0]):
                        cut_dataframe.iloc[j, 1] = 0

        if 'e' not in str(cut_dataframe.iloc[i, 9]).split('_')[0] and cut_dataframe.iloc[i, 7] != 'na' and \
                cut_dataframe.iloc[i, 7] != 'phage_resistance':
            cut_dataframe.iloc[i, 1] = 0
        if cut_dataframe.iloc[i, 9] != 'na' and 'e' not in str(cut_dataframe.iloc[i, 9]).split('_')[1]:
            cut_dataframe.iloc[i, 7] = 'na'
            cut_dataframe.iloc[i, 8] = 'na'

    cut_dataframe = cut_dataframe[~cut_dataframe['start_x'].isin([0])]
    cut_dataframe = cut_dataframe.drop(["domain_start", "domain_end", 'bac_name', 'phage_name'], axis=1)
    cut_dataframe.columns = ['name_x', 'start_x', 'end_x', 'number', 'dire_x', 'seq_x', 'new_name', 'structure_name',
                             'structure_id', 'E_value', 'gene', 'protein', 'gene_name']
    # Sort again by start position
    if dire == 1:
        cut_dataframe = cut_dataframe.sort_values(by=['start_x'], ascending=(False), axis=0)
    cut_dataframe.to_csv(seq_file_name + '_cut_dataframe1.txt', sep='\t', index=False)
    return cut_dataframe


def Information(seq_file_name, table, seq_name, n, confirmed, true_start, true_end, start, att_seq, cut_dataframe,
                cut_dataframe1):
    """
    Output result records and generate PICI prediction information table.

    Parameters:
    seq_file_name: str - output filename prefix.
    table: str - feature string returned by HTH_pri_judge.
    seq_name: str - sequence name.
    n: int - PICI type (1: with integrase, 2: with regulatory factor).
    confirmed: int - attachment site confirmation flag.
    true_start: int - confirmed start position.
    true_end: int - confirmed end position.
    start: int - original start position.
    att_seq: str - attachment site sequence.
    cut_dataframe: DataFrame - candidate region gene dataframe.
    cut_dataframe1: DataFrame - extended region gene dataframe.

    return: None
    """

    # Parse feature table string
    table_parts = table.split('{')
    HTH_hmm = table_parts[0]  # HTH domain name
    pri_hmm = table_parts[1]  # Primase domain name
    HTH = int(table_parts[2])  # HTH gene position
    pri_start_end = table_parts[3].split('+')
    pri = int(float(pri_start_end[0]))  # Primase start position
    all_pri = table_parts[3]  # Complete primase position range
    range1 = table_parts[4]  # PICI prediction range
    phage = table_parts[5]  # Phage feature count
    ICE = int(table_parts[6])  # ICE feature count
    ter = int(table_parts[7])  # Terminase flag
    dire = int(table_parts[8])  # Integration direction

    # Adjust prediction range using confirmed attachment sites    if confirmed == 1 and HTH_hmm != str(0) and pri_hmm != str(0):
    if confirmed == 1 and HTH_hmm != str(0) and pri_hmm != str(0):
        if true_start != 0:
            range1 = str(true_start) + '+' + range1.split('+')[1]
        if true_end != 0:
            range1 = range1.split('+')[0] + '+' + str(true_end)

    # Filtering
    condition = ''
    phage_judge = 'PICI'
    all_table = open(seq_file_name + '_Information.txt', 'a')

    if (n == 1 and HTH == 0 and pri == 0 and ter == 0) or (n == 2 and pri == 0):
        return 0
    #  Determine PICI completeness condition using match-case
    if n == 1 and HTH != 0 and pri != 0:
        condition = 'complete'
    if n == 1 and HTH == 0 and pri != 0:
        condition = 'uncomplete-HTH'
    if n == 1 and HTH != 0 and pri == 0:
        condition = 'uncomplete-pri'
    if n == 1 and HTH == 0 and pri == 0 and ter == 1:
        condition = 'short_PICI'
    if n == 2 and HTH != 0 and pri != 0:
        condition = 'uncomplete-int'
    if HTH_hmm != str(0) and pri_hmm != str(0) and (int(phage.split('+')[0]) > 2 or int(phage.split('+')[2]) == 1):
        phage_judge = 'phage'
    if ICE > 2:
        phage_judge = 'ICE'

    # Filter based on terminase features
    if HTH_hmm == pri_hmm and pri_hmm == HTH:
        if ter == 2:
            phage_judge = 'phage'
        if ter == 3 or ter == 4:
            phage_judge = 'bac'

    # Range validity check
    if range1 == str(0) or range1.split('+')[1] == str(0):
        return 0
    # Exclude specific primase types
    if pri_hmm in ['&HTH_36', '&IstB_IS21']:
        return 0
    # Special filtering for short PICIs
    if condition == 'short_PICI' and int(phage.split('+')[1]) >= 2:
        print('test5!!!')
        return 0

    print('HTH_' + str(HTH))
    print('pri_' + str(pri))

    # Confidence score calculation
    true_pri_list = ['PriCT_1', 'DUF5906', 'D5_N', 'Pox_D5', 'VirE', 'Prim-Pol', 'DUF927', 'HTH_56', 'phage_resistance']
    lev = 0
    if n == 1:
        lev = lev + 25
    print('int_lev_' + str(lev))
    if 'HTH_17' in HTH_hmm:
        lev = lev + 20
    elif HTH != 0:
        lev = lev + 25
    print('stl_lev_' + str(lev))
    # Primase feature scorin
    if pri != 0:
        pri_tab1 = pri_hmm.count('&')
        pri_tab2 = pri_hmm.count('|')
        pri_tab = pri_tab1 + pri_tab2
        pri_count = 0
        for line in true_pri_list:
            if pri_hmm.count(line):
                pri_count += 1
        if pri_count / pri_tab >= 0.5:
            lev = lev + 35
        else:
            lev = lev + 30
        print('pri_count_' + str(pri_count / pri_tab))
    print('pri_lev_' + str(lev))

    if confirmed == 1:
        lev = lev + 45
    print('att_lev_' + str(lev))
    if int(phage.split('+')[1]) >= 1:
        lev = lev + 40
    print('phage_lev_' + str(lev))

    # Virulence factor and antibiotic resistance gene annotation extraction
    VFDB = ''
    SARG = ''
    Abi = ''

    print('start!!!!__' + str(int(range1.split('+')[0])))
    print('end!!!!__' + str(int(range1.split('+')[1])))
    cut_dataframe = cut_dataframe.drop_duplicates(subset='new_name', keep='first', inplace=False)

    # Gene annotation extraction
    if n == 1:
        for i in range(len(cut_dataframe)):
            if int(cut_dataframe.iloc[i, 1]) >= int(range1.split('+')[0]) and int(cut_dataframe.iloc[i, 2]) <= int(
                    range1.split('+')[1]) and 'VFDB' in cut_dataframe.iloc[i, 11]:
                VFDB = VFDB + cut_dataframe.iloc[i, 11].split('#')[0] + ','
            if int(cut_dataframe.iloc[i, 1]) >= int(range1.split('+')[0]) and int(cut_dataframe.iloc[i, 2]) <= int(
                    range1.split('+')[1]) and 'SARG' in cut_dataframe.iloc[i, 11]:
                SARG = SARG + cut_dataframe.iloc[i, 11].split('#')[0] + ','
            if int(cut_dataframe.iloc[i, 1]) >= int(range1.split('+')[0]) and int(cut_dataframe.iloc[i, 2]) <= int(
                    range1.split('+')[1]) and '#Abi' in cut_dataframe.iloc[i, 11]:
                Abi = Abi + cut_dataframe.iloc[i, 11].split('#')[0] + ','
    if n == 2:
        for i in range(len(cut_dataframe)):
            if ((dire == -1 and int(cut_dataframe.iloc[i, 2]) <= int(range1.split('+')[1])) or (
                    dire == 1 and int(cut_dataframe.iloc[i, 2]) >= int(range1.split('+')[0]))) and 'VFDB' in \
                    cut_dataframe.iloc[i, 11]:
                VFDB = VFDB + cut_dataframe.iloc[i, 11].split('#')[0] + ','
            if ((dire == -1 and int(cut_dataframe.iloc[i, 2]) <= int(range1.split('+')[1])) or (
                    dire == 1 and int(cut_dataframe.iloc[i, 2]) >= int(range1.split('+')[0]))) and 'SARG' in \
                    cut_dataframe.iloc[i, 11]:
                SARG = SARG + cut_dataframe.iloc[i, 11].split('#')[0] + ','
            if ((dire == -1 and int(cut_dataframe.iloc[i, 2]) <= int(range1.split('+')[1])) or (
                    dire == 1 and int(cut_dataframe.iloc[i, 2]) >= int(range1.split('+')[0]))) and '#Abi' in \
                    cut_dataframe.iloc[i, 11]:
                Abi = Abi + cut_dataframe.iloc[i, 11].split('#')[0] + ','

    VFDB = VFDB.rstrip(',')
    SARG = SARG.rstrip(',')
    Abi = Abi.rstrip(',')

    all_table.write(seq_name + '\t')
    all_table.write(HTH_hmm + '\t')
    all_table.write(pri_hmm + '\t')
    all_table.write(range1 + '\t')
    all_table.write(phage + '\t')
    all_table.write(str(ICE) + '\t')
    all_table.write(str(ter) + '\t')
    all_table.write(phage_judge + '\t')
    all_table.write(condition + '\t')
    all_table.write(str(dire) + '\t')
    all_table.write(str(att_seq) + '\t')
    all_table.write(str(VFDB) + '\t')
    all_table.write(str(SARG) + '\t')
    all_table.write(str(Abi) + '\t')
    all_table.write(str(lev) + '\n')
    all_table.close()


def att_blast(seq_file_name, seq_name, int_seq, back_seq, start, end, name, dire, filepath):
    """
    Attachment site alignment, identifying attL/attR sites.

    Parameters:
    seq_file_name: str - sequence filename prefix.
    seq_name: str - sequence name.
    int_seq: str - integrase-side sequence (for attL detection).
    back_seq: str - background sequence (for attR detection).
    start: int - PICI candidate region start position.
    end: int - PICI candidate region end position.
    name: str - attachment site type name (e.g., GroEL, HU, etc.).
    dire: int - integration direction.
    filepath: str - working directory path.
    PICI_finder_path: str - tool installation path.

    Returns: 'PICI_start+PICI_end+att_seq' format containing site positions and sequence.
    """

    print('att_blast')
    # attL
    att_seq = 'no'
    int_file = open(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_int_cut_seq.txt', 'a')
    int_file.write('>' + seq_name + '\n')
    int_file.write(int_seq + '\n')
    int_file.close()
    # attR
    back_file = open(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_back_cut_seq.txt', 'a')
    back_file.write('>' + seq_name + '\n')
    back_file.write(back_seq + '\n')
    back_file.close()

    # BLAST Database Construction and Alignment
    # attL
    if os.path.exists(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_int_cut_seq.txt') and os.path.getsize(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_int_cut_seq.txt') > 1:
        subprocess.run('makeblastdb -in ' + seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_int_cut_seq.txt -input_type fasta -dbtype nucl -title ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_int_cut_seq -out ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_int_cut_seq', shell=True)
        subprocess.run('blastn -db ' + seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_int_cut_seq' + ' -query ' + filepath + '/databases/att.txt -outfmt 6 -task blastn-short -word_size 7 -evalue 1 -max_hsps 1 -perc_identity 90 -qcov_hsp_perc 60 -out ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_int_blast.txt', shell=True)
    # attR
    if os.path.exists(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_back_cut_seq.txt') and os.path.getsize(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_back_cut_seq.txt') > 1:
        subprocess.run('makeblastdb -in ' + seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_back_cut_seq.txt -input_type fasta -dbtype nucl -title ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_back_cut_seq -out ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_back_cut_seq', shell=True)
        subprocess.run('blastn -db ' + seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(
            end) + '_back_cut_seq' + ' -outfmt 6 -query ' + filepath + '/databases/att.txt -task blastn-short -word_size 7 -evalue 1 -max_hsps 1 -perc_identity 90 -qcov_hsp_perc 60 -out ' + seq_file_name + '_' + seq_name + '_' + str(
            start) + '_' + str(end) + '_back_blast.txt', shell=True)


    # BLAST result parsing
    PICI_start = 0
    PICI_end = 0
    length = 0
    # attL
    if os.path.exists(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_int_blast.txt') and os.path.getsize(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_int_blast.txt') > 1:
        int_att = pd.read_table(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_int_blast.txt',
                                sep='\t', header=None)
        print('name')
        print(name)
        for i in range(len(int_att)):
            if name in int_att.iloc[i, 0]:
                print('int_att')
                print(int_att.iloc[i, 0])
                print(int_att.iloc[i, 3])
                if int(int_att.iloc[i, 3]) >= length:
                    length = int(int_att.iloc[i, 3])
                else:
                    continue
                if int_att.iloc[i, 8] > int_att.iloc[i, 9]:
                    temp = int_att.iloc[i, 9]
                    int_att.iloc[i, 9] = int_att.iloc[i, 8]
                    int_att.iloc[i, 8] = temp
                    att_seq = int_att.iloc[i, 0]

                if dire == -1:
                    PICI_start = int_att.iloc[i, 8]
                if dire == 1:
                    PICI_end = int_att.iloc[i, 9]
    else:
        if dire == -1:
            PICI_start = 0
        elif dire == 1:
            PICI_end = 0

    # attR
    length = 0
    if os.path.exists(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_back_blast.txt') and os.path.getsize(
            seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_back_blast.txt') > 1:
        back_att = pd.read_table(seq_file_name + '_' + seq_name + '_' + str(start) + '_' + str(end) + '_back_blast.txt',
                                 sep='\t', header=None)
        for i in range(len(back_att)):
            if name in back_att.iloc[i, 0]:
                print('back_att')
                print(back_att.iloc[i, 0])
                print(back_att.iloc[i, 3])
                if back_att.iloc[i, 3] >= length:
                    length = int(back_att.iloc[i, 3])
                else:
                    continue

                if back_att.iloc[i, 8] > back_att.iloc[i, 9]:
                    temp = back_att.iloc[i, 9]
                    back_att.iloc[i, 9] = back_att.iloc[i, 8]
                    back_att.iloc[i, 8] = temp
                    att_seq = back_att.iloc[i, 0]

                if dire == -1:
                    PICI_end = back_att.iloc[i, 9]
                if dire == 1:
                    PICI_start = back_att.iloc[i, 8]
    else:
        if dire == -1:
            PICI_end = 0
        if dire == 1:
            PICI_start = 0

    # Attachment site sequence query
    att_table = pd.read_table(filepath + '/databases/att.txt', header=None)
    for i in range(len(att_table)):
        if att_seq in att_table.iloc[i, 0]:
            att_seq = att_table.iloc[i + 1, 0]
            break

    print('PICI_start+PICI_end')
    print(str(PICI_start) + '+' + str(PICI_end))
    return str(PICI_start) + '+' + str(PICI_end) + '+' + att_seq


def tRNA_end(cut_dataframe, seq_file_name, start, end, dire):
    """
    tRNA boundary detection, determining PICI integration boundaries.

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe.
    seq_file_name: str - sequence filename prefix.
    start: int - initial start position.
    end: int - initial end position.
    dire: int - integration direction (1: forward, -1: reverse).

    Returns: int - adjusted boundary position.
    """

    tRNA_end = cut_dataframe
    if dire == -1:
        true_end = end
    if dire == 1:
        true_end = start

    for i in range(len(tRNA_end)):
        if tRNA_end.iloc[i, 1] > start and tRNA_end.iloc[i, 2] < end:
            # Detect tRNA genes or integrases/recombinases in specific directions
            if 'tRNA-' in tRNA_end.iloc[i, 10] or (
                    ('Phage_integrase%' in tRNA_end.iloc[i, 7] or 'Recombinase%' in tRNA_end.iloc[i, 7]) and int(
                    tRNA_end.iloc[i, 4]) == -1 and dire == -1) and abs(
                    int(tRNA_end.iloc[i, 3]) - int(tRNA_end.iloc[0, 3])) > 5:
                true_end = tRNA_end.iloc[i, 2]
                break

            if 'tRNA-' in tRNA_end.iloc[i, 10] or (
                    ('Phage_integrase%' in tRNA_end.iloc[i, 7] or 'Recombinase%' in tRNA_end.iloc[i, 7]) and int(
                    tRNA_end.iloc[i, 4]) == 1 and dire == 1) and abs(
                    int(tRNA_end.iloc[i, 3]) - int(tRNA_end.iloc[0, 3])) > 5:
                true_end = tRNA_end.iloc[i, 1]
                break
    return true_end


def HTH_pri_judge(cut_dataframe, cut_dataframe1, seq_file_name, start, end, forest, forest1, length, n, dire, filepath):
    """
    Regulatory element identification, detecting HTH and primase features.

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe.
    cut_dataframe1: DataFrame - extended region gene dataframe (for background analysis).
    seq_file_name: str - sequence filename prefix.
    start: int - candidate region start position.
    end: int - candidate region end position.
    length: int - total sequence length.
    n: int - PICI type (1: with integrase, 2: with regulatory factor).
    dire: int - integration direction.
    filepath: str - working directory path.
    PICI_finder_path: str - tool installation path.

    Returns: str - feature string.
    """

    # HTH
    # HTH regulatory factor identification
    HTH_list = cut_dataframe

    print(HTH_list)
    print('start_' + str(start))
    print('end_' + str(end))

    HTH_site = 0
    flag = 0
    flag1 = 0
    flag2 = 0
    flag3 = 0
    HTH_site_1 = 0
    HTH_site_2 = 0
    HTH_site_3 = 0
    HTH_name_1 = ''
    HTH_hmm_1 = ''
    HTH_gene_1 = ''

    # Search for valid HTH regulatory factors
    for i in range(len(HTH_list)):
        flag = 0
        flag1 = 0
        flag2 = 0
        flag3 = 0
        if HTH_list.iloc[i, 1] >= start and HTH_list.iloc[i, 2] <= end and (
                HTH_list.iloc[i, 10] == 'stl_gene' or 'HTH' in HTH_list.iloc[i, 7]) and abs(
                int(HTH_list.iloc[i, 2]) - int(HTH_list.iloc[0, 2])) <= 25000 and int(HTH_list.iloc[i, 4]) == dire:
            HTH_site = int(HTH_list.iloc[i, 6])

            # Check neighboring genes to identify HTH regulatory modules
            for j in range(i + 1, len(HTH_list)):
                if int(HTH_list.iloc[j, 6]) == HTH_site:
                    print('break1')
                    continue
                if abs(int(HTH_list.iloc[j, 6]) - HTH_site) == 1:
                    flag3 = 1
                if (abs(int(HTH_list.iloc[j, 6]) - HTH_site) == 1 or abs(
                        int(HTH_list.iloc[j, 6]) - HTH_site) == 2) and int(HTH_list.iloc[j, 4]) == -dire:
                    if flag3 == 1 and abs(int(HTH_list.iloc[j, 6]) - HTH_site) == 2 and HTH_list.iloc[
                        j, 7] != 'HTH_3' and HTH_list.iloc[j, 7] != 'HTH_19' and HTH_list.iloc[j, 7] != 'HTH_26' and \
                            HTH_list.iloc[j, 7] != 'HTH_31' and HTH_list.iloc[j, 7] != 'HTH_17':
                        continue
                    print('break3')

                    if abs(int(HTH_list.iloc[j, 6]) - HTH_site) == 1:
                        if n == 2 and HTH_list.iloc[0, 6] != HTH_site:
                            flag2 = 1
                            HTH_site_1 = 0
                            HTH_site_2 = 0
                            HTH_site_3 = 0
                            break
                        HTH_site_1 = HTH_site
                        flag = 1
                    if flag1 == 1 and flag == 0:
                        flag = 1

                    if abs(int(HTH_list.iloc[j, 6]) - HTH_site) == 2:
                        if n == 2 and HTH_list.iloc[0, 6] != HTH_site:
                            flag2 = 1
                            HTH_site_1 = 0
                            HTH_site_2 = 0
                            HTH_site_3 = 0
                            break
                        HTH_site_1 = HTH_site
                        flag1 = 1
                    print('break4')

                    if 'HTH' in HTH_list.iloc[j, 7] or HTH_list.iloc[j, 10] == 'stl_gene':
                        HTH_site_2 = HTH_list.iloc[j, 6]

                    for k in range(j + 1, len(HTH_list)):
                        if HTH_site_2 != 0 and abs(int(HTH_list.iloc[k, 6]) - int(HTH_site_2)) == 1 and int(
                                HTH_list.iloc[k, 4]) == -dire and (
                                'HTH' in HTH_list.iloc[k, 7] or HTH_list.iloc[k, 10] == 'stl_gene') and int(
                                HTH_site_2) != int(HTH_list.iloc[k, 6]):
                            HTH_site_3 = HTH_list.iloc[k, 6]
                            break
                        if HTH_site_2 == 0 and abs(int(HTH_list.iloc[k, 6]) - int(HTH_site_1)) == 2 and int(
                                HTH_list.iloc[k, 4]) == -dire and (
                                'HTH' in HTH_list.iloc[k, 7] or HTH_list.iloc[k, 10] == 'stl_gene') and int(
                                HTH_site_2) != int(HTH_list.iloc[k, 6]):
                            HTH_site_3 = HTH_list.iloc[k, 6]
                            break
                if flag == 1 or flag2 == 1:
                    break

        if flag == 1 or flag2 == 1:
            break

    # Collect detailed information for HTH sites
    for i in range(len(HTH_list)):
        if HTH_site_1 != 0 and HTH_site_1 == HTH_list.iloc[i, 6]:
            HTH_name_1 = HTH_name_1 + '|' + HTH_list.iloc[i, 7]
            HTH_hmm_1 = HTH_hmm_1 + '|' + HTH_list.iloc[i, 8]
            HTH_gene_1 = HTH_gene_1 + '|' + HTH_list.iloc[i, 10]

    # HTH validation and filtering
    if HTH_site_1 != 0:
        for i in range(len(HTH_list)):
            # Neighboring regulatory elements
            if (dire == -1 and int(HTH_list.iloc[i, 6]) - int(HTH_site_1) <= 5 and int(HTH_list.iloc[i, 6]) - int(
                    HTH_site_1) > 0) or (
                    dire == 1 and int(HTH_list.iloc[i, 6]) - int(HTH_site_1) >= -5 and int(HTH_list.iloc[i, 6]) - int(
                    HTH_site_1) < 0):
                if HTH_list.iloc[i, 7] == 'HTH_3' or HTH_list.iloc[i, 7] == 'HTH_19' or HTH_list.iloc[
                    i, 7] == 'HTH_26' or HTH_list.iloc[i, 7] == 'HTH_31' or HTH_list.iloc[i, 10] == 'stl_gene':
                    for j in range(i + 1, len(HTH_list)):
                        if abs(int(HTH_list.iloc[j, 6]) - int(HTH_site_1)) <= 6 and abs(
                                int(HTH_list.iloc[j, 6]) - int(HTH_list.iloc[i, 6])) == 1 and int(
                                HTH_list.iloc[j, 4]) == -dire and int(HTH_list.iloc[i, 4]) == dire:
                            print(int(HTH_list.iloc[i, 6]))
                            print(int(HTH_list.iloc[j, 6]))
                            print(int(HTH_list.iloc[j, 4]))
                            HTH_site_1 = 0
                            HTH_site_2 = 0
                            HTH_site_3 = 0
                            print('break5')
                            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                                0) + '{' + str(0) + '{' + str(0) + '{' + str(dire)

    if 'PF01381' not in HTH_hmm_1 and 'HTH_19' not in HTH_name_1 and 'HTH_17' not in HTH_name_1 and 'HTH_26' not in HTH_name_1 and 'HTH_31' not in HTH_name_1 and 'stl_gene' not in HTH_gene_1:
        HTH_site_1 = 0
        HTH_site_2 = 0
        HTH_site_3 = 0

    # Label PICI-related genes
    for i in range(len(HTH_list)):
        if HTH_list.iloc[i, 1] >= start and HTH_list.iloc[i, 2] <= end:
            # Label HTH-neighboring genes, integrases, and resistance genes as PICI genes
            if (HTH_site_1 != 0 and ((int(HTH_list.iloc[i, 6]) == HTH_site_1 or int(
                    HTH_list.iloc[i, 6]) == HTH_site_1 + 1 or int(
                    HTH_list.iloc[i, 6]) == HTH_site_1 + 2) and dire == -1) or ((int(
                    HTH_list.iloc[i, 6]) == HTH_site_1 or int(HTH_list.iloc[i, 6]) == HTH_site_1 - 1 or int(
                    HTH_list.iloc[i, 6]) == HTH_site_1 - 2) and dire == 1)) or HTH_list.iloc[i, 10] == 'int' or 'Abi' in \
                    HTH_list.iloc[i, 7] or 'tsst-1' in HTH_list.iloc[i, 7] or 'DUF1474' in HTH_list.iloc[
                i, 7] or 'SAPIS-gp6' in HTH_list.iloc[i, 7] or 'Phage_pRha' in HTH_list.iloc[i, 7] or 'SAPI' in \
                    HTH_list.iloc[i, 7] or 'DUF1492' in HTH_list.iloc[i, 7]:
                HTH_list.iloc[i, 10] = 'PICI_gene'
            # Label initial integration site genes
            if n == 1 and int(HTH_list.iloc[i, 6]) == int(HTH_list.iloc[0, 6]):
                HTH_list.iloc[i, 10] = 'PICI_gene'
            # Exclude integrases near HTH sites
            if HTH_site_1 != 0 and ((dire == -1 and (
                    int(HTH_list.iloc[i, 6]) - HTH_site_1 == 1 or int(HTH_list.iloc[i, 6]) - HTH_site_1 == 2)) or (
                                            dire == 1 and (int(HTH_list.iloc[i, 6]) - HTH_site_1 == -1 or int(
                                            HTH_list.iloc[i, 6]) - HTH_site_1 == -2))):
                print('test1!!!!!')
                if HTH_list.iloc[i, 7] == 'Phage_integrase%' or HTH_list.iloc[i, 7] == 'Phage_int_SAM_3%':
                    print('test2!!!!!')
                    return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                        0) + '{' + str(0) + '{' + str(0) + '{' + str(dire)

    for i in range(len(cut_dataframe1)):
        if n == 2 and cut_dataframe1.iloc[i, 1] >= start and cut_dataframe1.iloc[i, 2] <= end:
            if (HTH_site_1 != 0 and ((int(cut_dataframe1.iloc[i, 6]) == HTH_site_1 or int(
                    cut_dataframe1.iloc[i, 6]) == HTH_site_1 + 1 or int(
                    cut_dataframe1.iloc[i, 6]) == HTH_site_1 + 2) and dire == -1) or ((int(
                    cut_dataframe1.iloc[i, 6]) == HTH_site_1 or int(cut_dataframe1.iloc[i, 6]) == HTH_site_1 - 1 or int(
                    cut_dataframe1.iloc[i, 6]) == HTH_site_1 - 2) and dire == -1)) or cut_dataframe1.iloc[
                i, 10] == 'int' or 'Abi' in cut_dataframe1.iloc[i, 7] or 'tsst-1' in cut_dataframe1.iloc[
                i, 7] or 'DUF1474' in cut_dataframe1.iloc[i, 7] or 'SAPIS-gp6' in cut_dataframe1.iloc[
                i, 7] or 'Phage_pRha' in cut_dataframe1.iloc[i, 7] or 'SAPI' in cut_dataframe1.iloc[
                i, 7] or 'DUF1492' in cut_dataframe1.iloc[i, 7]:
                cut_dataframe1.iloc[i, 10] = 'PICI_gene'

    # Primase
    # Initialize primase-related variables and data structures
    pri_list = pd.DataFrame(columns=(['pri_site']))
    pri_file_deal = HTH_list
    pri_name = pd.read_table(filepath + '/databases/pri_6_7.txt', header=None)

    Prim_site = 0
    Prim_dire = 0
    pri_site = 0
    pri_x = 0
    temp = 0
    primase_site = 0
    IstB_site = 0
    pri_site_temp = 0

    # Identify primase features
    for i in range(len(pri_file_deal)):
        for j in range(len(pri_name)):
            if pri_file_deal.iloc[i, 1] >= start and pri_file_deal.iloc[i, 2] <= end and (
                    pri_file_deal.iloc[i, 8] == pri_name.iloc[j, 0] or 'phage_resistance' in pri_file_deal.iloc[
                i, 7]) and abs(int(pri_file_deal.iloc[i, 3]) - int(pri_file_deal.iloc[0, 3])) <= 25000 and abs(
                    int(pri_file_deal.iloc[i, 6]) - int(pri_file_deal.iloc[0, 6])) > 1 and int(
                    pri_file_deal.iloc[i, 4]) == -dire:
                temp_1 = 0
                temp_2 = 0

                # Filter false positives: exclude specific
                if (pri_file_deal.iloc[i, 7] == 'AAA_25' and pri_file_deal.iloc[
                    i, 11] == 'DNA repair protein RadA') or (pri_file_deal.iloc[i, 7] == 'AAA'):
                    continue

                # IstB_IS21: Check for nearby transposases
                if pri_file_deal.iloc[i, 7] == 'IstB_IS21':
                    IstB_site = pri_file_deal.iloc[i, 6]
                    for k in range(len(pri_file_deal)):
                        if (int(pri_file_deal.iloc[k, 6]) == int(IstB_site) + 1 or int(pri_file_deal.iloc[k, 6]) == int(
                                IstB_site) - 1) and ('rve' in pri_file_deal.iloc[k, 7] or 'Tnp' in pri_file_deal.iloc[
                            k, 7] or 'Transposase' in pri_file_deal.iloc[k, 7]):
                            temp_1 = 1
                            break
                    if temp_1 == 1:
                        continue
                # Prim-Pol: Check for nearby integrases
                if pri_file_deal.iloc[i, 7] == 'Prim-Pol':
                    Prim_site = pri_file_deal.iloc[i, 6]
                    Prim_dire = pri_file_deal.iloc[i, 4]
                    for k in range(len(pri_file_deal)):
                        if (int(pri_file_deal.iloc[k, 6]) <= int(Prim_site) + 5 and int(
                                pri_file_deal.iloc[k, 6]) >= int(Prim_site) - 5) and (
                                pri_file_deal.iloc[k, 7] == 'Phage_int_SAM_3%' or pri_file_deal.iloc[
                            k, 7] == 'Phage_integrase%') and pri_file_deal.iloc[k, 4] == Prim_dire:
                            temp_2 = 1
                            break
                    if temp_2 == 1:
                        continue

                # Elevate phage_resistance annotation to domain name level
                if 'phage_resistance' in pri_file_deal.iloc[i, 10]:
                    pri_file_deal.iloc[i, 7] = pri_file_deal.iloc[i, 10]
                pri_site = pri_file_deal.iloc[i, 6]
                # Ensure reasonable distance between primase candidates
                if pri_site_temp != 0:
                    if abs(int(pri_site) - int(pri_site_temp)) > 5:
                        break
                    else:
                        pri_site_temp = int(pri_site)
                if pri_site_temp == 0:
                    pri_site_temp = int(pri_site)
                pri_list = pd.concat([pri_list, pri_file_deal[i:i + 1]], axis=0)
                break

    pri_list = pri_list.drop(['pri_site'], axis=1)
    # Merge multiple domain annotations for the same primase site and remove duplicates
    for i in range(len(pri_list)):
        if primase_site != pri_list.iloc[i, 6]:
            pri_x = i
            primase_site = pri_list.iloc[i, 6]
            pri_list.iloc[i, 1] = pri_list.iloc[i, 6]
        if temp == pri_list.iloc[i, 6]:
            pri_list.iloc[pri_x, 7] = pri_list.iloc[pri_x, 7] + '|' + pri_list.iloc[i, 7]
            pri_list.iloc[i, 7] = None
        temp = pri_list.iloc[i, 6]
    pri_list = pri_list.dropna()

    if pri_list.empty:
        pri_list = pd.DataFrame({'0': ['na'], '1': [0], '2': [0], '3': [0], '4': [0], '5': [0], '6': [0], '7': [0]})
    else:
        pri_list_new = pd.concat([pri_list.iloc[:, [0]], pri_list.iloc[:, [7]]], axis=1)
        pri_list_new = pd.concat([pri_list_new, pri_list.iloc[:, [6]]], axis=1)

        # Label primase-related genes as PICI genes
        for i in range(len(pri_file_deal)):
            for j in range(len(pri_list_new)):
                if int(pri_file_deal.iloc[i, 6]) == int(pri_list_new.iloc[j, 2]):
                    pri_file_deal.iloc[i, 10] = 'PICI_gene'

                # Exclude specific types of false positive primases
                if 'Resolvase' in pri_list_new.iloc[j, 1] or 'Phage_GPA' in pri_list_new.iloc[j, 1] or 'RepA_N' in \
                        pri_list_new.iloc[j, 1] or 'RecT' in pri_list_new.iloc[j, 1]:
                    return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                        0) + '{' + str(0) + '{' + str(0) + '{' + str(dire)
        # Similarly label primase-related genes in the extended dataframe
        for i in range(len(cut_dataframe1)):
            for j in range(len(pri_list_new)):
                if int(cut_dataframe1.iloc[i, 6]) == int(pri_list_new.iloc[j, 2]):
                    cut_dataframe1.iloc[i, 10] = 'PICI_gene'

                if 'Resolvase' in pri_list_new.iloc[j, 1] or 'Phage_GPA' in pri_list_new.iloc[j, 1] or 'RepA_N' in \
                        pri_list_new.iloc[j, 1] or 'RecT' in pri_list_new.iloc[j, 1]:
                    print('test3!!!!!')
                    return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                        0) + '{' + str(0) + '{' + str(0) + '{' + str(dire)
        print('pri_list_new')
        print(pri_list_new)
    print(pri_list)

    pri_hmm = ''
    for i in range(len(pri_list)):
        pri_hmm = pri_hmm + '&' + str(pri_list.iloc[i, 7])

    # Terminase detection: check for presence of terminase (phage packaging signal)
    ter = ter_judge(HTH_list, start, end)
    if n == 1 and HTH_site_1 == 0 and pri_list.iloc[0, 1] == 0:
        if ter == 0:
            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                0) + '{' + str(0) + '{' + str(dire)
        else:
            if dire == -1:
                end_hp = pri_file_deal.iloc[0, 6] + 11
            if dire == 1:
                end_hp = pri_file_deal.iloc[0, 6] - 11
    if n == 2 and (HTH_site_1 == 0 or pri_list.iloc[0, 1] == 0):
        return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
            0) + '{' + str(0) + '{' + str(dire)

    # PICI boundary heuristic rules
    if pri_list.iloc[0, 1] == 0 and dire == -1:
        end_hp = HTH_site_1 + 10
    if pri_list.iloc[0, 1] == 0 and dire == 1:
        end_hp = HTH_site_1 - 10
    if (HTH_site_1 == 0 or (pri_list.iloc[0, 1] != 0 and HTH_site_1 != 0)) and dire == -1:
        end_hp = pri_list.iloc[len(pri_list) - 1, 1] + 5
    if (HTH_site_1 == 0 or (pri_list.iloc[0, 1] != 0 and HTH_site_1 != 0)) and dire == 1:
        end_hp = pri_list.iloc[len(pri_list) - 1, 1] - 5

    # Exclusion gene filtering: check for incompatible gene types within PICI region
    for i in range(len(pri_file_deal)):
        if (int(pri_file_deal.iloc[i, 6]) <= end_hp and dire == -1) or (
                int(pri_file_deal.iloc[i, 6]) >= end_hp and dire == 1):
            if 'ORF6N' in pri_file_deal.iloc[i, 7] or 'ORF6C' in pri_file_deal.iloc[i, 7] or 'Rep_2' in \
                    pri_file_deal.iloc[i, 7] or 'ResIII' in pri_file_deal.iloc[i, 7] or 'tcpA^' in pri_file_deal.iloc[
                i, 7] or 'Helicase_C' in pri_file_deal.iloc[i, 7] or 'MOBP1^' in pri_file_deal.iloc[i, 7] or 'YqgF' in \
                    pri_file_deal.iloc[i, 7] or 'MOBT^' in pri_file_deal.iloc[i, 7] or 't4cp2^' in pri_file_deal.iloc[
                i, 7] or 'MOBF' in pri_file_deal.iloc[i, 7] or 'FtsK_SpoIIIE' in pri_file_deal.iloc[
                i, 7] or 'DNA_methylase' in pri_file_deal.iloc[i, 7] or 'N6_N4_Mtase' in pri_file_deal.iloc[
                i, 7] or 'Whib' in pri_file_deal.iloc[i, 7] or 'SSB' in pri_file_deal.iloc[i, 7] or 'RadC' in \
                    pri_file_deal.iloc[i, 7] or 'N6_Mtase' in pri_file_deal.iloc[i, 7] or 'SWI2_SNF2' in \
                    pri_file_deal.iloc[i, 7] or 'Replicase' in pri_file_deal.iloc[i, 7] or 'SNF2-rel_dom' in \
                    pri_file_deal.iloc[i, 7] or 'RepB_primase' in pri_file_deal.iloc[i, 7] or 'DNA_pol3_beta' in \
                    pri_file_deal.iloc[i, 7] or 'Myosin_tail_1' in pri_file_deal.iloc[i, 7] or 'ANT' in \
                    pri_file_deal.iloc[i, 7] or 'DNA_pol_A' in pri_file_deal.iloc[i, 7] or 'AntA' in pri_file_deal.iloc[
                i, 7] or 'Rep_trans' in pri_file_deal.iloc[i, 7] or 'FtsK_gamma' in pri_file_deal.iloc[
                i, 7] or 'Myosin_tail_1' in pri_file_deal.iloc[i, 7] or 'RepC' in pri_file_deal.iloc[
                i, 7] or 'Relaxase' in pri_file_deal.iloc[i, 7]:
                return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                    0) + '{' + str(2) + '{' + str(dire)

    no_pri = pd.read_table(filepath + '/databases/pri_7.txt', header=None)
    gene_hp_count = 0
    gene_count = 0
    HTH_list = pri_file_deal.drop_duplicates(subset='start_x', keep='first', inplace=False)
    cut_dataframe1 = cut_dataframe1.drop_duplicates(subset='start_x', keep='first', inplace=False)
    if n == 1:
        for m in range(len(HTH_list)):
            if (int(HTH_list.iloc[m, 6]) <= end_hp and dire == -1) or (
                    int(HTH_list.iloc[m, 6]) >= end_hp and dire == 1):
                gene_count += 1
                if HTH_list.iloc[m, 10] == 'PICI_gene' or HTH_list.iloc[m, 10] != 'bac':
                    gene_hp_count += 1
                for p in range(len(no_pri)):
                    if no_pri.iloc[p, 0] == HTH_list.iloc[m, 8]:
                        return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                            0) + '{' + str(0) + '{' + str(2) + '{' + str(dire)
    # Apply stricter filtering for PICI types lacking integrases
    if n == 2:
        for m in range(len(cut_dataframe1)):
            if (int(cut_dataframe1.iloc[m, 6]) <= end_hp and dire == -1) or (
                    int(cut_dataframe1.iloc[m, 6]) >= end_hp and dire == 1):
                gene_count += 1
                if cut_dataframe1.iloc[m, 10] == 'PICI_gene' or cut_dataframe1.iloc[m, 10] != 'bac':
                    gene_hp_count += 1
                for p in range(len(no_pri)):
                    if no_pri.iloc[p, 0] == cut_dataframe1.iloc[m, 8]:
                        return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                            0) + '{' + str(0) + '{' + str(2) + '{' + str(dire)

    print('gene_count_' + str(gene_count))
    print('gene_hp_count_' + str(gene_hp_count))
    if gene_count != 0 and ter != 1:
        hp_ratio = float(gene_hp_count) / float(gene_count)
        if hp_ratio < 0.6:
            print('hp_ratio' + str(hp_ratio))
            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                0) + '{' + str(3) + '{' + str(dire)

    gene_count = 0
    gene_dir_count = 0
    gene_count1 = 0
    gene_dir_count1 = 0

    # Gene direction analysis downstream of HTH region
    if HTH_site_1 != 0:
        for m in range(len(HTH_list)):
            if (int(HTH_list.iloc[m, 6]) > HTH_site_1 and int(HTH_list.iloc[m, 6]) <= end_hp and dire == -1) or (
                    int(HTH_list.iloc[m, 6]) < HTH_site_1 and int(HTH_list.iloc[m, 6]) >= end_hp and dire == 1):
                gene_count += 1
                if int(HTH_list.iloc[m, 4]) == -dire:
                    gene_dir_count += 1
            if abs(HTH_site_1 - int(HTH_list.iloc[0, 6])) >= 5 and (
                    (int(HTH_list.iloc[m, 6]) <= HTH_site_1 and dire == -1) or (
                    int(HTH_list.iloc[m, 6]) >= HTH_site_1 and dire == 1)):
                gene_count1 += 1
                if int(HTH_list.iloc[m, 4]) == dire:
                    gene_dir_count1 += 1
    print('gene_count_' + str(gene_count))
    print('gene_dir_count_' + str(gene_dir_count))

    if gene_count != 0 and ter != 1:
        dir_ratio = float(gene_dir_count) / float(gene_count)
        if dir_ratio < 0.6:
            print('dir_ratio' + str(dir_ratio))
            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                0) + '{' + str(4) + '{' + str(dire)
    if gene_count1 != 0 and ter != 1:
        dir_ratio1 = float(gene_dir_count1) / float(gene_count1)
        if dir_ratio1 < 0.6:
            print('dir_ratio1' + str(dir_ratio1))
            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                0) + '{' + str(4) + '{' + str(dire)

    print('pri_list')
    print(pri_list)
    print(int(pri_list.iloc[0, 6]))

    # Primase neighboring region analysis
    pri_gene_count = 0
    pri_gene_dire = 0
    if pri_list.iloc[0, 1] != 0:
        for i in range(len(HTH_list)):
            if dire == -1 and int(HTH_list.iloc[i, 6]) > int(pri_list.iloc[0, 6]) and int(HTH_list.iloc[i, 6]) <= int(
                    pri_list.iloc[0, 6]) + 5:
                pri_gene_count += 1
                if int(HTH_list.iloc[i, 4]) == -dire:
                    pri_gene_dire += 1
            if dire == 1 and int(HTH_list.iloc[i, 6]) < int(pri_list.iloc[0, 6]) and int(HTH_list.iloc[i, 6]) >= int(
                    pri_list.iloc[0, 6]) - 5:
                pri_gene_count += 1
                if int(HTH_list.iloc[i, 4]) == -dire:
                    pri_gene_dire += 1
    if pri_gene_count != 0 and ter != 1:
        pri_dir_ratio = float(pri_gene_dire) / float(pri_gene_count)
        if pri_dir_ratio < 0.6:
            print('pri_dir_ratio' + str(pri_dir_ratio))
            return str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(0) + '{' + str(
                0) + '{' + str(4) + '{' + str(dire)

    print('HTH_site_1_' + str(HTH_site_1))
    print('HTH_site_2_' + str(HTH_site_2))
    print('HTH_site_3_' + str(HTH_site_3))
    print(pri_list)
    print('pri_file_deal')
    print(pri_file_deal)
    phage_range = PICI_feature_phage(pri_file_deal, seq_file_name, start, end, forest1, length, dire)
    PICI_true_range = PICI_feature(pri_file_deal, seq_file_name, start, end, forest, length, dire, filepath)
    phage = phage_judge(pri_file_deal, phage_range, dire, PICI_true_range, filepath)
    ICE = ICE_judge(pri_file_deal, PICI_true_range, dire, filepath)
    return HTH_name_1 + '{' + pri_hmm + '{' + str(HTH_site_1) + '{' + str(
        pri_list.iloc[0, 1]) + '{' + PICI_true_range + '{' + str(phage) + '{' + str(ICE) + '{' + str(ter) + '{' + str(
        dire)


def ter_judge(HTH_list, start, end):
    """
    Terminase detection: check for presence of terminase

    Parameters:
    HTH_list: DataFrame - candidate region gene dataframe containing all gene domain annotations

    Returns: int - 1 indicates terminase present, 0 indicates absent
    """

    for i in range(len(HTH_list)):
        if 'Terminase' in HTH_list.iloc[i, 7]:
            if abs(HTH_list.iloc[i, 6] - HTH_list.iloc[0, 6]) < 10:
                return 1
    return 0


def PICI_feature_phage(cut_dataframe, seq_file_name, start, end, forest1, length, dire):
    """
    Phage feature analysis using machine learning model 1

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe
    seq_file_name: str - sequence filename prefix
    start: int - candidate region start position
    end: int - candidate region end position
    length: int - total sequence length
    dire: int - integration direction (1 or -1)
    PICI_finder_path: str - tool installation path
    n: int - PICI type (1: with integrase, 2: with regulatory factor)

    Returns: str - 'start+end' format representing PICI prediction range
    """

    # Data preparation and sliding window feature extraction
    gene_hp_ratio = []
    gene_start = []
    gene_end = []

    PICI_feature = cut_dataframe.drop_duplicates(subset='new_name', keep='first', inplace=False)
    if int(length) < 7000:
        if dire == -1:
            true_start = PICI_feature.iloc[0, 1]
            true_end = PICI_feature.iloc[len(PICI_feature) - 1, 2]
        if dire == 1:
            true_start = PICI_feature.iloc[len(PICI_feature) - 1, 1]
            true_end = PICI_feature.iloc[0, 2]
        return str(true_start) + '+' + str(true_end)
    # Sliding window analysis
    for i in range(0, int(length), 500):
        gene_number = 0
        gene_hp_number = 0
        gene_start_site = 0
        gene_end_site = 0
        if i + 7000 > int(length):
            break
        # Count gene features within current window
        for j in range(len(PICI_feature)):
            if (PICI_feature.iloc[j, 1] >= int(start) + i and PICI_feature.iloc[j, 2] <= int(
                    start) + i + 7000 and dire == -1) or (
                    PICI_feature.iloc[j, 2] <= int(end) - i and PICI_feature.iloc[j, 1] >= int(
                    end) - i - 7000 and dire == 1):
                if dire == -1:
                    if gene_start_site == 0:
                        gene_start_site = PICI_feature.iloc[j, 1]
                    gene_end_site = PICI_feature.iloc[j, 2]
                if dire == 1:
                    if gene_end_site == 0:
                        gene_end_site = PICI_feature.iloc[j, 2]
                    gene_start_site = PICI_feature.iloc[j, 1]
                gene_number += 1
                # Count hypothetical proteins and PICI-related genes
                if '%' in PICI_feature.iloc[j, 7] or PICI_feature.iloc[j, 10] == 'PICI_gene' or '*' in \
                        PICI_feature.iloc[j, 7] or PICI_feature.iloc[j, 10] != 'bac' or (
                        'tail' in PICI_feature.iloc[j, 10] and 'hypothetical protein' in PICI_feature.iloc[j, 11]) or \
                        PICI_feature.iloc[j, 10] == 'capsid_gene' or PICI_feature.iloc[j, 10] == 'head_gene' or \
                        PICI_feature.iloc[j, 10] == 'lysis_gene':
                    gene_hp_number += 1

        # Calculate hypothetical protein ratio for current window and store features
        if gene_number == 0:
            gene_hp_ratio.append(0)
        else:
            gene_hp_ratio.append(gene_hp_number / gene_number)
        gene_start.append(gene_start_site)
        gene_end.append(gene_end_site)

    # Machine learning prediction and feature gene detection
    PICI_feature_table_out = forest1.predict(np.array(gene_hp_ratio).reshape(-1, 1))
    PICI_feature_table = pd.concat([pd.DataFrame(gene_hp_ratio), pd.DataFrame(PICI_feature_table_out)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_start)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_end)], axis=1)

    # Detect boundary positions of PICI feature genes
    feature_gene_end = 0
    for i in range(len(PICI_feature)):
        if ((PICI_feature.iloc[i, 10] == 'PICI_gene' and int(PICI_feature.iloc[i, 2]) <= int(end) and dire == -1) or (
                PICI_feature.iloc[i, 10] == 'PICI_gene' and int(PICI_feature.iloc[i, 1]) >= int(
                start) and dire == 1)) and abs(
                int(PICI_feature.iloc[i, 2]) - int(PICI_feature.iloc[0, 1])) < 25000 and 'Abi' not in PICI_feature.iloc[
            i, 7]:
            print('site_' + str(PICI_feature.iloc[i, 3]))
            print(abs(int(PICI_feature.iloc[i, 2]) - int(PICI_feature.iloc[0, 1])))
            feature_gene_end = PICI_feature.iloc[i, 2]

    print('phage_feature_gene_end')
    print(feature_gene_end)

    # PICI boundary determination logic
    count = 0
    flag = 0
    # Find transition boundaries from PICI to bacterial genome
    for i in range(len(PICI_feature_table)):
        if i == 0 and PICI_feature_table.iloc[i, 1] == 'bac':
            count += 1
        if i == 1 and PICI_feature_table.iloc[i, 1] == 'bac':
            count += 1
        if i == 2 and PICI_feature_table.iloc[i, 1] == 'bac':
            count += 1
        if count == 3:
            true_start = PICI_feature_table.iloc[0, 2]
            true_end = PICI_feature_table.iloc[0, 3]
            break

        if i + 3 < len(PICI_feature_table) and PICI_feature_table.iloc[i, 1] == 'PICI' and PICI_feature_table.iloc[
            i + 3, 1] == 'bac' and PICI_feature_table.iloc[i + 1, 1] != 'PICI':
            if (int(PICI_feature_table.iloc[i, 3]) >= int(feature_gene_end) and dire == -1) or (
                    int(PICI_feature_table.iloc[i, 3]) <= int(feature_gene_end) and dire == 1):
                flag = 1
                if dire == -1:
                    true_start = PICI_feature_table.iloc[0, 2]
                    true_end = PICI_feature_table.iloc[i, 3]
                if dire == 1:
                    true_start = PICI_feature_table.iloc[i, 2]
                    true_end = PICI_feature_table.iloc[0, 3]
                break
    # If no clear transition boundary found, use the entire analysis region
    if flag != 1:
        if dire == -1:
            true_start = PICI_feature_table.iloc[0, 2]
            true_end = PICI_feature_table.iloc[len(PICI_feature_table) - 1, 3]
        if dire == 1:
            true_start = PICI_feature_table.iloc[len(PICI_feature_table) - 1, 2]
            true_end = PICI_feature_table.iloc[0, 3]
    if int(true_end) > int(end) and dire == -1:
        true_end = end
    if int(true_start) < int(start) and dire == 1:
        true_start = start

    print('phage_range2')
    print('true_phage_start')
    print(true_start)
    print(true_end)
    PICI_feature_table.to_csv(
        seq_file_name + '_' + cut_dataframe.iloc[0, 0] + '_' + str(start) + '_' + str(end) + '_feature.txt', sep='\t',
        index=False)

    return str(true_start) + '+' + str(true_end)


def PICI_feature(cut_dataframe, seq_file_name, start, end, forest, length, dire, filepath):
    """
    PICI feature analysis using machine learning model 2 for region classification

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe
    seq_file_name: str - sequence filename prefix
    start: int - candidate region start position
    end: int - candidate region end position
    length: int - total sequence length
    dire: int - integration direction (1 or -1)
    filepath: str - working directory path
    PICI_finder_path: str - tool installation path
    n: int - PICI type (1: with integrase, 2: with regulatory factor)

    Returns: str - 'start+end+feature_gene_position' format
             containing PICI range and key gene information
    """

    # Multi-dimensional feature analysis
    gene_density = []
    gene_hp_ratio = []
    gene_mid_len = []
    gene_dire_ratio = []
    gene_start = []
    gene_end = []
    PICI_feature = cut_dataframe.drop_duplicates(subset='new_name', keep='first', inplace=False)
    ICE = pd.read_table(filepath + '/databases/tail', header=None)

    if int(length) < 7000:
        if dire == -1:
            true_start = PICI_feature.iloc[0, 1]
            true_end = PICI_feature.iloc[len(PICI_feature) - 1, 2]
        if dire == 1:
            true_start = PICI_feature.iloc[len(PICI_feature) - 1, 1]
            true_end = PICI_feature.iloc[0, 2]

    for i in range(0, int(length), 500):
        gene_len = []
        gene_mid = 0
        gene_number = 0
        gene_dire_con = 0
        gene_hp_number = 0
        gene_start_site = 0
        gene_end_site = 0

        if i + 7000 > int(length):
            break

        for j in range(len(PICI_feature)):
            if (PICI_feature.iloc[j, 1] >= int(start) + i and PICI_feature.iloc[j, 2] <= int(
                    start) + i + 7000 and dire == -1) or (
                    PICI_feature.iloc[j, 2] <= int(end) - i and PICI_feature.iloc[j, 1] >= int(
                    end) - i - 7000 and dire == 1):
                if dire == -1:
                    if gene_start_site == 0:
                        gene_start_site = PICI_feature.iloc[j, 1]
                    gene_end_site = PICI_feature.iloc[j, 2]
                if dire == 1:
                    if gene_end_site == 0:
                        gene_end_site = PICI_feature.iloc[j, 2]
                    gene_start_site = PICI_feature.iloc[j, 1]
                gene_number += 1
                gene_len.append(int(PICI_feature.iloc[j, 2]) - int(PICI_feature.iloc[j, 1]))

                # ICE gene filtering: check if gene is conjugation transposon-related
                flag = 0
                for k in range(len(ICE)):
                    if PICI_feature.iloc[j, 7] == ICE.iloc[k, 1]:
                        flag = 1
                        break

                if flag == 1 or '%' in PICI_feature.iloc[j, 7] or PICI_feature.iloc[j, 10] == 'PICI_gene' or '*' in \
                        PICI_feature.iloc[j, 7] or PICI_feature.iloc[j, 10] != 'bac' or (
                        'tail' in PICI_feature.iloc[j, 10] and 'hypothetical protein' in PICI_feature.iloc[j, 11]) or \
                        PICI_feature.iloc[j, 10] == 'capsid_gene' or PICI_feature.iloc[j, 10] == 'head_gene' or \
                        PICI_feature.iloc[j, 10] == 'lysis_gene':
                    gene_hp_number += 1

                # Gene direction continuity statistics
                if j + 1 < len(PICI_feature) and PICI_feature.iloc[j, 4] == PICI_feature.iloc[j + 1, 4]:
                    gene_dire_con += 1
        l = len(gene_len)
        gene_len.sort()

        # Empty window
        if gene_len == []:
            gene_mid = 0
            gene_density.append(0)
            gene_hp_ratio.append(0)
            gene_mid_len.append(0)
            gene_dire_ratio.append(0)
            gene_start.append(gene_start_site)
            gene_end.append(gene_end_site)
        else:
            if l % 2 == 0:
                gene_mid = (gene_len[int(l / 2) - 1] + gene_len[int(l / 2)]) / 2
            else:
                gene_mid = gene_len[int((l - 1) / 2)]
            gene_density.append(gene_number / 7000 * 100)
            gene_hp_ratio.append(gene_hp_number / gene_number)
            gene_mid_len.append(gene_mid)
            gene_dire_ratio.append((gene_number - gene_dire_con) / gene_number)
            gene_start.append(gene_start_site)
            gene_end.append(gene_end_site)

    # Machine learning classification and feature integration
    PICI_feature_table = pd.concat([pd.DataFrame(gene_density), pd.DataFrame(gene_hp_ratio)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_dire_ratio)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_mid_len)], axis=1)
    PICI_feature_table = PICI_feature_table.values
    if PICI_feature_table.size > 0:
        print('TRUE')
    else:
        return str(start) + '+' + str(end)

    PICI_feature_table_out = forest.predict(PICI_feature_table)
    PICI_feature_table = pd.DataFrame(PICI_feature_table)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(PICI_feature_table_out)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_start)], axis=1)
    PICI_feature_table = pd.concat([PICI_feature_table, pd.DataFrame(gene_end)], axis=1)
    print(PICI_feature_table)
    print('start_' + str(start))
    print('end_' + str(end))
    print(dire)

    feature_gene_end = 0
    for i in range(len(PICI_feature)):
        if ((PICI_feature.iloc[i, 10] == 'PICI_gene' and int(PICI_feature.iloc[i, 2]) <= int(end) and dire == -1) or (
                PICI_feature.iloc[i, 10] == 'PICI_gene' and int(PICI_feature.iloc[i, 1]) >= int(
                start) and dire == 1)) and abs(
                int(PICI_feature.iloc[i, 2]) - int(PICI_feature.iloc[0, 1])) < 25000 and 'Abi' not in PICI_feature.iloc[
            i, 7]:
            print('site_' + str(PICI_feature.iloc[i, 3]))
            print(abs(int(PICI_feature.iloc[i, 2]) - int(PICI_feature.iloc[0, 1])))
            feature_gene_end = PICI_feature.iloc[i, 2]

    # Boundary determination logic
    count = 0
    flag = 0
    for i in range(len(PICI_feature_table)):
        if i == 0 and PICI_feature_table.iloc[i, 4] == 'bac':
            count += 1
        if i == 1 and PICI_feature_table.iloc[i, 4] == 'bac':
            count += 1
        if i == 2 and PICI_feature_table.iloc[i, 4] == 'bac':
            count += 1
        if count == 3:
            true_start = PICI_feature_table.iloc[0, 6]
            true_end = PICI_feature_table.iloc[0, 5]
            break
        if i + 1 < len(PICI_feature_table) and PICI_feature_table.iloc[i, 4] == 'PICI' and PICI_feature_table.iloc[
            i + 1, 4] == 'bac':
            if (int(PICI_feature_table.iloc[i, 6]) >= int(feature_gene_end) and dire == -1) or (
                    int(PICI_feature_table.iloc[i, 5]) <= int(feature_gene_end) and dire == 1):
                flag = 1
                if dire == -1:
                    true_start = PICI_feature_table.iloc[0, 5]
                    true_end = PICI_feature_table.iloc[i, 6]
                if dire == 1:
                    true_start = PICI_feature_table.iloc[i, 5]
                    true_end = PICI_feature_table.iloc[0, 6]
                break

    if flag != 1:
        if dire == -1:
            true_start = PICI_feature_table.iloc[0, 5]
            true_end = PICI_feature_table.iloc[len(PICI_feature_table) - 1, 6]
        if dire == 1:
            true_start = PICI_feature_table.iloc[len(PICI_feature_table) - 1, 5]
            true_end = PICI_feature_table.iloc[0, 6]

    # Boundary adjustment and result return
    if int(true_end) > int(end) and dire == -1:
        true_end = end
    if int(true_start) < int(start) and dire == 1:
        true_start = start

    print('feature_gene_end_' + str(feature_gene_end))
    print('true_start')
    print('true_start_' + str(true_start))

    PICI_feature_table.to_csv(
        seq_file_name + '_' + cut_dataframe.iloc[0, 0] + '_' + str(start) + '_' + str(end) + 'true_feature.txt',
        sep='\t', index=False)
    print('TRUR_range')
    print(str(true_start) + '+' + str(true_end))
    return str(true_start) + '+' + str(true_end) + '+' + str(feature_gene_end)


def phage_judge(cut_dataframe, phage_range, dire, PICI_true_range, filepath):
    """
    Phage determination based on tail proteins and other features

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe
    phage_range: str - range obtained from phage feature analysis
    dire: int - integration direction (1 or -1)
    PICI_true_range: str - actual PICI prediction range
    filepath: str - working directory path
    PICI_finder_path: str - tool installation path

    Returns: str - 'tail_count+PICI_tail_count+BLAST_result' format quantifying phage features
    """

    # Phage tail protein detection
    tail = pd.read_table(filepath + '/databases/tail', header=None)
    phage = cut_dataframe
    print('phage')
    print(phage)

    true_start = phage_range.split('+')[0]
    true_end = phage_range.split('+')[1]
    print('true_phage_start')
    print(true_start)
    print(true_end)

    tail_count = 0
    PICI_tail_count = 0
    temp = 0

    for i in range(len(phage)):
        for j in range(len(tail)):
            if ((int(phage.iloc[i, 2]) <= int(true_end) and dire == -1) or (
                    int(phage.iloc[i, 2]) >= int(true_start) and dire == 1)) and '^' not in phage.iloc[
                i, 7] and '%' not in phage.iloc[i, 7] and ((tail.iloc[j, 1] in phage.iloc[i, 7] and 'tail' in tail.iloc[
                j, 2] and 'hypothetical protein' in phage.iloc[i, 11]) or (
                                                                   'tail_gene' in phage.iloc[i, 11] or 'capsid_gene' in
                                                                   phage.iloc[i, 11] or 'head_gene' in phage.iloc[
                                                                       i, 11] or 'lysis_gene' in phage.iloc[i, 11])):
                # Tail proteins within PICI region (within 15kb of integration site)
                if '*' in phage.iloc[i, 7] and abs(int(phage.iloc[i, 2]) - int(phage.iloc[0, 2])) <= 25000:
                    PICI_tail_count += 1
                    phage.iloc[i, 10] = 'PICI_gene'
                    break
                if '*' not in phage.iloc[i, 7] and abs(int(phage.iloc[i, 2]) - int(phage.iloc[0, 2])) <= 15000:
                    break
                if i - 2 >= 0 and '^' in phage.iloc[i - 2, 7]:
                    break
                if i - 1 >= 0 and '^' in phage.iloc[i - 1, 7]:
                    break
                if i + 1 < len(phage) and '^' in phage.iloc[i + 1, 7]:
                    break
                if i + 2 < len(phage) and '^' in phage.iloc[i + 2, 7]:
                    break
                # Other tail proteins (possibly from complete phages)
                if '*' not in phage.iloc[i, 7] or ('*' in phage.iloc[i, 7] and abs(
                        int(phage.iloc[i, 2]) - int(phage.iloc[0, 2])) > 25000) or 'tail_gene' in phage.iloc[
                    i, 11] or 'lysis_gene' in phage.iloc[i, 11]:
                    if temp == int(phage.iloc[i, 6]):
                        break
                    print('phage.iloc[i, 6]')
                    print(phage.iloc[i, 2])
                    print(phage.iloc[i, 6])
                    print(phage.iloc[i, 7])
                    tail_count += 1
                    temp = int(phage.iloc[i, 6])

    # Background region analysis and integration site prediction
    phage = cut_dataframe.drop_duplicates(subset='new_name', keep='first', inplace=False)
    count = 0
    other_count = 0
    for i in range(len(phage)):
        if (phage.iloc[i, 2] > int(true_end) and dire == -1) or (phage.iloc[i, 2] < int(true_start) and dire == 1):
            count += 1
            if count > 7:
                break
            if 'tail_gene' in phage.iloc[i, 11] or 'capsid_gene' in phage.iloc[i, 11] or 'head_gene' in phage.iloc[
                i, 11] or 'lysis_gene' in phage.iloc[i, 11]:
                print('phage.iloc[i, 6]')
                print(phage.iloc[i, 6])
                print(phage.iloc[i, 7])
                other_count += 1

    tail_count = tail_count + other_count
    print('tail_count_' + str(tail_count))
    print('PICI_tail_count_' + str(PICI_tail_count))
    blast_out = blast(phage, PICI_true_range)
    print('blast_out_' + str(blast_out))
    print('other_count_' + str(other_count))
    return str(tail_count) + '+' + str(PICI_tail_count) + '+' + str(blast_out)


def ICE_judge(cut_dataframe, PICI_true_range, dire, filepath):
    """
    ICE element determination, detecting conjugation transfer-related genes

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe
    PICI_true_range: str - actual PICI prediction range string 'start+end'
    dire: int - integration direction (1 or -1)
    filepath: str - working directory path
    PICI_finder_path: str - tool installation path

    Returns: int - ICE feature count, higher values indicate greater likelihood of ICE rather than PICI
    """

    # ICE-related gene detection
    T4SS = pd.read_table(filepath + '/databases/tail', header=None)
    ICE = cut_dataframe
    ICE_count = 0
    true_start = PICI_true_range.split('+')[0]
    true_end = PICI_true_range.split('+')[1]
    temp = 0
    for i in range(len(ICE)):
        if ICE.iloc[i, 7] == 'Terminase_2':
            break
        for j in range(len(T4SS)):
            if ((int(ICE.iloc[i, 2]) <= int(true_end) and dire == -1) or (
                    int(ICE.iloc[i, 2]) >= int(true_start) and dire == 1)) and T4SS.iloc[j, 1] in ICE.iloc[
                i, 7] and 'ICE' in T4SS.iloc[j, 2]:
                if temp == int(ICE.iloc[i, 6]):
                    continue
                ICE_count += 1
                temp = int(ICE.iloc[i, 6])
    print('ICE_count_' + str(ICE_count))
    return ICE_count


def blast(cut_dataframe, PICI_true_range):
    """
    Simple ratio calculation, count phage gene proportion

    Parameters:
    cut_dataframe: DataFrame - candidate region gene dataframe
    PICI_true_range: str - actual PICI prediction range string 'start+end'

    Returns: int - 0 indicates phage gene proportion ≤50%, 1 indicates >50%
    """

    # Data preprocessing: remove duplicate gene annotations
    cut_dataframe = cut_dataframe.drop_duplicates(subset='new_name', keep='first', inplace=False)
    true_start = PICI_true_range.split('+')[0]
    true_end = PICI_true_range.split('+')[1]
    print(cut_dataframe)
    number = 0
    phage_number = 0

    # Count gene composition within PICI region
    for i in range(len(cut_dataframe)):
        if int(cut_dataframe.iloc[i, 1]) >= int(true_start) and int(cut_dataframe.iloc[i, 2]) <= int(true_end):
            number += 1
            if 'phage' in cut_dataframe.iloc[i, 12] and cut_dataframe.iloc[i, 10] != 'PICI_gene':
                phage_number += 1
    print(str(number) + '+' + str(phage_number))
    if number == 0 or phage_number / number <= 0.5:
        return 0
    else:
        return 1


if __name__ == "__main__":
    seq_file_name = sys.argv[1]
    filepath = sys.argv[2]
    print('seq_file_name')
    print(seq_file_name)
    int_gene(seq_file_name, filepath)
    print('final_part3')
    subprocess.run('python3 scripts/PICI_Result_Integration.py ' + seq_file_name, shell=True)

