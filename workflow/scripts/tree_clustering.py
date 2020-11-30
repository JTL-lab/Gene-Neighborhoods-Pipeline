#!/usr/bin/env python
"""
Script to generate two distance matrices from a directory of .treefiles rendered
with IQ-TREE respectively using:
a) Robinson-Foulds distances
b) Boot Split Distances (distances weighted according to boot split support)
"""

import os
import sys
import argparse
import logging
import ete3 as et
import dendropy as dd
from dendropy.utility.textprocessing import StringIO
from make_ML_tree_vis import make_tree_vis

def clean_string(dd_cluster):
    """
    Tidy up taxon strings
    """
    clust = dd_cluster.split(" ")
    clust_s = clust[1].replace("'", "")
    clust_str = clust_s.replace("\n","")

    return clust[1]

def get_cluster_vis(RF_matrix, BSD_matrix):
    """
    Create UPGMA and NJ trees based on RF and BSD distances respectively
    """
    pdm_rf = dd.PhylogeneticDistanceMatrix.from_csv(src=StringIO(RF_matrix),\
                                                   delimiter=",")

    pdm_bsd = dd.PhylogeneticDistanceMatrix.from_csv(src=StringIO(BSD_matrix),\
                                                   delimiter=",")

    #UPGMA trees
    upgma_tree_rf = pdm_rf.upgma_tree()
    upgma_rf_str = clean_string(upgma_tree_rf.as_string("newick"))

    upgma_tree_bsd = pdm_bsd.upgma_tree()
    upgma_bsd_str = clean_string(upgma_tree_bsd.as_string("newick"))

    #Neighbor Joining trees
    nj_tree_rf = pdm_rf.nj_tree()
    nj_rf_str = clean_string(nj_tree_rf.as_string("newick"))

    nj_tree_bsd = pdm_bsd.nj_tree()
    nj_bsd_str = clean_string(nj_tree_bsd.as_string("newick"))

    #Create visualizations
    upgma_rf_tree = make_tree_vis(upgma_rf_str, "UPGMA RF", "UPGMA_RF")
    upgma_bsd_tree = make_tree_vis(upgma_bsd_str, "UPGMA BSD", "UPGMA_BSD")
    nj_rf_tree = make_tree_vis(nj_rf_str, "NJ RF","NJ_RF")
    nj_bsd_tree = make_tree_vis(nj_bsd_str, "NJ BSD", "NJ_BSD")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Script to create UPGMA and\
                                     NJ tree cluster visualizations from .csv\
                                     distance matrixes.")

    parser.add_argument('-matrix_path', '-mp', type=str,
                    help='Path to .csv files containing distance matrices.')

    args = parser.parse_args()

    path_invalid = False
    try:
        os.path.exists(str(sys.argv[2]))
        logging.debug("Path to matrices is valid")
        os.chdir(args.matrix_path)

    except:
        logging.error("Path is invalid!")
        path_invalid = True

    if path_invalid:
        logging.error("Path to distance matrix was invalid.\
                    #  Please use valid path.")
        sys.exit(1)

    with open('rf_matrix.csv', 'r') as rf_matrix:
        RF_data = rf_matrix.read()
        print(RF_data)

    with open('bsd_matrix.csv', 'r') as bsd_matrix:
        BSD_data = bsd_matrix.read()
        print(BSD_data)

    get_cluster_vis(RF_data, BSD_data)
