#!/usr/bin/env python
import os
import sys
import subprocess
import argparse
import ete3 as et
import dendropy as dd
from dendropy.utility.textprocessing import StringIO
from make_ML_tree_vis import make_tree_vis

def clean_string(dd_cluster):
    
    clust = dd_cluster.split(" ")
    clust[1].replace("'", "")
    clust.pop(0)
    
    return clust
    
def get_cluster_vis(RF_matrix, BSD_matrix): 
    
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
    nj_rf_str = clean_string((nj_tree_rf.as_string("newick"))

    nj_tree_bsd = pdm_bsd.nj_tree()
    nj_bsd_str = clean_string(nj_tree_bsd.as_string("newick"))
    
    #Create visualizations
    upgma_rf_tree = make_tree_vis(upgma_rf_str, "UPGMA RF")
    upgma_bsd_tree = make_tree_vis(upgma_bsd_str, "UPGMA BSD")
    nj_rf_tree = make_tree_vis(upgma_rf_str, "NJ RF")
    nj_rf_tree = make_tree_vis(upgma_rf_str, "NJ BSD")
    


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Script to create UPGMA and \
                                     NJ clusters from .csv distance matrixes.")
    
    parser.add_argument('-rf_matrix', '-rf', type=str,
                        help='Name of .csv containing RF distance matrix.')
    
    parser.add_argument('-bsd_matrix', '-bsd', type=str,
                        help='Name of .csv containing BSD distance matrix.')
    
    args = parser.parse_args()
    
    with open(str(argsv[2]), 'r') as file_1:
        RF_data = file_1.read().replace('\n', '')
    
    with open(str(argsv[4]), 'r') as file_2:
        BSD_data = file_2.read().replace('\n', '')
        
    get_cluster_vis(RF_data, BSD_data)
        
        

