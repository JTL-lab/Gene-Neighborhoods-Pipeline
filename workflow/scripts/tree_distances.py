#!/usr/bin/env python

"""
Script to generate distance matrixes (Robinson - Foulds and Boot Split \
Distance) for all bootstrapped trees contained within a directory of .treefiles
"""

import os
import sys
import argparse
import logging
import csv
import dendropy as dd
from dendropy.calculate import treecompare

#Function to obtain bipartitions as objects and as Newick strings
def get_bipartitions(tree):

    bipartitions = []
    bip_strings = []

    for bip in tree:

        tree_bip = [split.split_as_newick_string(tree.taxon_namespace) \
                    for split in tree.encode_bipartitions()]

        bip_strings.append(tree_bip)
        bipartitions.append(bip)

    return bipartitions, bip_strings

#Get the common set of leaves between the two trees,
#i.e the shared leaf-set of species -> must have at least 4 species in common!

'''Input: List of leaves from two respective trees
    Output: Returns shared list of nodes for trees with >= 4 nodes'''
def getShared(firstTreeLeaves, secondTreeLeaves):

    in_common = list(set(firstTreeNodes)&set(secondTreeNodes))

    assert (len(in_common) < 4), "Trees have less than 4 species in common: \
    cannot be compared!"

    return in_common

#Returns list of leaves not present in both trees
def getDiff(firstTreeLeaves, secondTreeLeaves):

    return list(set(firstTreeLeaves)^set(secondTreeLeaves))

#Step 3: Prune all splits to the common leaf-set of species
'''Input: List of bipartitions as Nodes and as strings for two trees
   Output: Returns list of nodes not present in both trees, \
   as Nodes and as Newick strings'''
def getPrunedBips(bips1, bipStrings1, bips2, bipStrings2):
    bipsDiff = []
    bipsEqual = []

    for i in range(len(bips1)):
        found = False

        if any(x.split_bitmask == bips1[i] for x in bips2):
            found = True

        if found:
            bipsEqual.append(bips1[i])
        else:
            bipsDiff.append(bips1[i])

    for j in range(len(bips2)):
        found = False

        if any(x.split_bitmask == bips2[j] for x in bips1):
            found = True

        if found:
            if not any(x.split_bitmask == bips2[j] for x in bipsEqual):
                bipsEqual.append(bips2[j])
        else:
            if not any(x.split_bitmask == bips2[j] for x in bipsDiff):
                bipsDiff.append(bips2[j])

    return bipsDiff, bipsEqual

def getSum(bipartitions):

    nums = []

    for i in range(len(bipartitions)):
        if (bipartitions[i].label is not None):
            num = bipartitions[i].label
            nums.append(num)

    #cast to int
    BS_vals = list(map(int, nums))
    for i in range(len(BS_vals)):
        BS_vals[i] = float(BS_vals[i]/100)

    sum_BS = sum(BS_vals)


    return sum_BS, nums


#def getSumDiffSplits(bootstrapTree):
'''Input: Outputs from getSum()
   Output: Mean of bootstrap support values'''
def getMeanBSD(sumBS, listBS):

    if len(listBS) == 0:
        return 0

    else:
        return float(sumBS/len(listBS))


def eBSD(sum_total_BS, sum_mutual_split_BS, mean_BS_mutual):

    return (1 - ((sum_mutual_split_BS/sum_total_BS)*mean_BS_mutual))

def dBSD(sum_total_BS, sum_diff_split_BS, mean_BS_diff):

    return ((sum_diff_split_BS/sum_total_BS)*mean_BS_diff)

def BSD(sum_total_BS, sum_mutual_split_BS, sum_diff_split_BS, \
mean_BS_mutual, mean_BS_diff):
    equalBSD = eBSD(sum_total_BS, sum_mutual_split_BS, mean_BS_mutual)
    diffBSD = dBSD(sum_total_BS, sum_diff_split_BS, mean_BS_diff)
    bootSplitDist = ((equalBSD+diffBSD)/2)

    return bootSplitDist, equalBSD, diffBSD

def get_bsd(tree_1, tree_2):

    #Get bipartition and string reps of bipartitions for both trees
    bips_1, names_1 = get_bipartitions(tree_1)
    bips_2, names_2 = get_bipartitions(tree_2)

    #Obtain list of bipartitions shared by both trees and list of bipartitions
    #that differ between trees
    bipsDiff, bipsEqual = getPrunedBips(bips_1,names_1,bips_2,names_2)

    #Obtain total sum(s) of BS supports for bipartitions for both trees
    sum_1, toDiscard1 = getSum(bips_1)
    sum_2, toDiscard2 = getSum(bips_2)

    #Tally the total BS support between both trees
    sum_total_BS = sum_1 + sum_2

    sum_mutual_BS, list_mutual_BS = getSum(bipsEqual)
    sum_diff_BS, list_diff_BS = getSum(bipsDiff)

    mean_BS_mutual = getMeanBSD(sum_mutual_BS, list_mutual_BS)
    mean_BS_diff = getMeanBSD(sum_diff_BS, list_diff_BS)

    #Calculate BSD
    BSD_val, eBSD, dBSD = BSD(sum_total_BS, sum_mutual_BS, sum_diff_BS, \
    mean_BS_mutual, mean_BS_diff)

    return BSD_val, eBSD, dBSD, sum_total_BS, sum_mutual_BS, sum_diff_BS, \
    mean_BS_mutual, mean_BS_diff


def make_distance_matrix(tree_dict, tree_list, dist_type):

    if (dist_type == 'rf'):
        file_name = "rf_matrix.csv"

    else:
        file_name = 'bsd_matrix.csv'

    with open(file_name, 'w') as csvfile:

        csv_writer = csv.writer(csvfile)

        #Write header line
        header = []
        header.append(',')
        for tree in tree_dict:
            header.append(tree_dict[tree])

        csv_writer.writerow(header)


        for tree_ref in tree_list:

            #Append row tree header
            distances = []
            distances.append(tree_dict[tree_ref])

            for i in range(len(tree_list)):

                if trees[i] is not tree_ref:

                    #Robinson-Foulds distances
                    if (dist_type == 'rf'):
                        dist_val = treecompare.weighted_robinson_foulds_distance\
                        (tree_ref,tree_list[i], is_bipartitions_updated=True)

                    #Boot-Split distances
                    else:
                        dist_val, eBSD, dBSD, sum_total, sum_mutual, sum_diff,\
                    mean_mutual, mean_diff = get_bsd(tree_ref, tree_list[i])

                    distances.append(dist_val)

            csv_writer.writerow(distances)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate two distance \
                                     matrixes (RF, BSD) from directory of \
                                     bootstrapped phylogenetic trees stored\
                                     as .treefile files. Please ensure that\
                                     all files are properly labeled (e.g tree_\
                                     1, tree_2, etc.)')

    parser.add_argument('-tree_path','-tp', type=str,
                        help='Path to directory containing all .treefile trees\
                        that will be considered in the distance matrices.')

    args = parser.parse_args()

    #Catch any invalid path!
    path_invalid = False

    try:
        os.path.exists(str(sys.argv[2]))
        logging.debug("Path is valid")
        os.chdir(args.tree_path)

    except:
        logging.error("Path is invalid!")
        path_invalid = True

    if path_invalid:
        logging.error("Path to .treefiles was invalid. Please use valid path.")
        sys.exit(1)

    #Will become dictionary used to inform distance matrix calculations
    files = []
    trees = []

    #Open user inputted directory
    for filename in os.listdir(os.getcwd()):

        #Create tree objects for each and store in an array
        if filename.endswith('.treefile'):

            #Store name of file
            name = filename.split('.')
            tree_name = "tree_" + name[0]
            files.append(tree_name)

            with open(filename, 'r') as file:

                #Obtain Newick tree
                tree_temp = file.read().replace('\n', '')
                trees.append(tree_temp)

    #Create dendropy trees for each tree stored in trees
    taxa = dd.TaxonNamespace(label='Salmonella')
    dd_trees = []

    for i in range(len(trees)):
        tree = dd.Tree.get(data=trees[i], schema="newick",
                           taxon_namespace=taxa)
        #tree.reconstruct_taxon_namespace
        tree.encode_bipartitions()

        dd_trees.append(tree)

    #Link the files with their respective trees
    tree_dict = dict(zip(dd_trees, files))

    #Call matrix functions to create .csv files
    make_distance_matrix(tree_dict, dd_trees, 'rf')
    make_distance_matrix(tree_dict, dd_trees, 'bsd')
