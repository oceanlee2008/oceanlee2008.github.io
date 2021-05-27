#!/usr/bin/env python

import sys
import argparse
import numpy as np
import time
import zipfile
from collections import defaultdict


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            Data_all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                Data_all_reads.append(ends)
        return Data_all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need
"""

def Input_data_reads(input_reads):
    trinity_reads = []
    Data_all_reads = [item for list_data in input_reads for item in list_data]
    for i in Data_all_reads:
        a = i[:16]
        b = i[16:32]
        c = i[32:48]
        trinity_reads.append([a, b, c])

    all_index_Refer = defaultdict(list)
    for i in range(0, len(reference) - 16):
        kmer = reference[i:i + 16]
        all_index_Refer[kmer].append(i)
    return all_index_Refer,trinity_reads

def Make_genome_hash(reference, key_length):
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash

def Appeared_base_snp(refer_data_Dict):
    snp_appeared = []
    mutation_index = {}
    # local the mutation position
    for i, j in refer_data_Dict.items():
        for k in j:
            if k in snp_appeared:
                for i in all_index_Refer[i]:
                    genomeIndex = k[2] + i
                    mut = [k[1], k[0], genomeIndex]

                    if mutation_index.get(genomeIndex) != None:
                        mutation_index[genomeIndex].append((k[1], k[0]))
                    else:
                        mutation_index[genomeIndex] = [(k[1], k[0])]
            else:
                snp_appeared.append(k)
    return mutation_index

def SNP_position_checking(mutation_index):
    for i, j in mutation_index.items():
        # filter out less common variation
        if len(j) <= 2:  # note: this seems to be the max
            continue
        else:
            # there's several variation, make sure they're consistent!
            consistent = True
            inconsistent = 0
            for k in j:
                if k != j[0]:
                    consistent = False
                    inconsistent += 1
                    break
            # consistent = snp
            if consistent:
                snp = [j[0][0], j[0][1], i]
                variation.append(snp)

            if consistent == False:
                number_SNP = {}
                snp_Count_max = -1
                Base_SNPs = ()
                for k in range(0, len(j)):
                    snp = (j[k][0], j[k][1], i)
                    if number_SNP.get(snp) != None:
                        number_SNP[snp] += 1
                        if number_SNP[snp] > snp_Count_max:
                            snp_Count_max = number_SNP[snp]
                            Base_SNPs = snp
                    else:
                        number_SNP[snp] = 1

                snp_probability = float(snp_Count_max) / float(len(j))
                print('snp',Base_SNPs,snp_probability)
                # majority by over 2/3rds
                if snp_probability > 0.67:
                    snp = list(Base_SNPs)
                    variation.append(snp)
                    continue
                else:
                    pass

    return number_SNP,variation

def Insertion_deletion(indels,all_index_Refer,reference,read):
    # read errors
    snp_appearedInsertions = []
    snp_appearedDeletions = []
    insertions = []
    deletions = []
    for ind in indels:
        first, second, third = ind
        read = first + second + third
        startIndexs = all_index_Refer[first]

        for index in startIndexs:
            i = 0
            # read diverges from reference
            try:
                while read[i] == reference[index + i]:
                    i += 1
            except:
                pass

            maxIndelLen = 5
            for len in range(1, maxIndelLen + 1):
                # to checking the insertion
                if read[i + len:i + 5 + len] == reference[index + i:index + i + 5]:
                    insertTup = (read[i:i + len], index + i)
                    if insertTup in snp_appearedInsertions:

                        insertions.append(insertTup)
                    else:
                        snp_appearedInsertions.append(insertTup)

                # check for deletion
                if read[i:i + 5] == reference[index + i + len:index + i + 5 + len]:
                    deleteTup = (reference[index + i:index + i + len], index + i)
                    if deleteTup in snp_appearedDeletions:
                        deletions.append(deleteTup)
                    else:
                        snp_appearedDeletions.append(deleteTup)
    insertions = list(set(insertions))
    deletions = list(set(deletions))
    return insertions,deletions

def HammingDistance(a, b):
    if len(a) != len(b): return float('inf')
    return sum(1 if a[i] != b[i] else 0 for i in range(len(a)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                                 'of a genome and a set of reads and aligns the reads to the reference genome, '
                                                 'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)
    """
        TODO: Use this space to implement any additional functions you might need
    """
        
    refer_data_Dict = {}
    indels = []
    # input and modify reads, as independent
    all_index_Refer,trinity_reads = Input_data_reads(input_reads)
    Make_genome_hash(reference, 50)
    # Place the mapped index where each part of the triple matches to find all variation and indels!
    for t in trinity_reads:
        read = ""
        random_index = []
        Mutation_Refer_Index = 100
        Mutation_offset = 0
        
        trinity_index = {1: [], 2: [], 3: []}
        Number_mismathes = 0
        Match_Table = False
        for i in [1, 2, 3]:
            read = t[i - 1]
            if all_index_Refer.get(read) != None:
                Match_Table = True
                trinity_index[i].extend(all_index_Refer[read])
                continue
            else:
                Number_mismathes += 1

        if Number_mismathes != 1:
            continue

        if len(trinity_index[1]) == 0:
            read = t[0]
            Mutation_offset = -16
            Mutation_Refer_Index = 2
        elif len(trinity_index[2]) == 0:
            read = t[1]
            Mutation_offset = 16
            Mutation_Refer_Index = 1
        elif len(trinity_index[3]) == 0:
            read = t[2]
            Mutation_offset = 16
            Mutation_Refer_Index = 2

        for i in trinity_index[Mutation_Refer_Index]:
            random_index.append(i + Mutation_offset)

        # hamming distance to find snps
        min_Ham_Dis = float('inf')
        Refer_minHam = ""
        for i in random_index:
            ham = HammingDistance(read, reference[i:i + 16])
            if ham < min_Ham_Dis:
                min_Ham_Dis = ham
                Refer_minHam = reference[i:i + 16]

        if Refer_minHam != "" and min_Ham_Dis == 1:
            ref = Refer_minHam
            #  the position of snp
            for i in range(0, len(ref) - 1):
                if (read[:i] == ref[:i]) and (read[i + 1:] == ref[i + 1:]):
                    if refer_data_Dict.get(ref) != None:
                        refer_data_Dict[ref].append((read[i], ref[i], i))
                    else:
                        refer_data_Dict[ref] = [(read[i], ref[i], i)]
                    break

        else:
            indels.append(t)
    # SNP_position_checking and find read errors
    variation = []
    mutation_index = Appeared_base_snp(refer_data_Dict)
    number_SNP,variation = SNP_position_checking(mutation_index)
    snps = variation
    
    insertions,deletions = Insertion_deletion(indels,all_index_Refer,reference,read)


    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)

