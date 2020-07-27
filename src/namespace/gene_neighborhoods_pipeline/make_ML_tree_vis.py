#!/usr/bin/env python

#Script to ouput ete3 powered phylogenetic tree visualizations given 
#a homologous sequence file

import sys
import subprocess
import argparse 
import ete3 as et
      
def make_tree_vis(newick_tree, gene_name, file_name):
    
    tree = et.Tree(newick_tree)
    
    #Customizations
    nifty = et.TreeStyle()
    nifty.show_branch_support = True
    nifty.branch_vertical_margin = 15
    nifty.scale = 10
    
    #Leaf node customizations
    style_L = et.NodeStyle()
    style_L['shape'] = 'sphere'
    style_L['size'] = 10
    style_L['fgcolor'] = '#6BC245'
    style_L['hz_line_type'] = 1
    style_L['hz_line_color'] = '#cccccc'

    #Internal node customizations
    style_I = et.NodeStyle()
    style_I['shape'] = 'sphere'
    style_I['size'] = 10
    style_I['fgcolor'] = '#5271FF'
    style_I['hz_line_type'] = 1
    style_I['hz_line_color'] = '#cccccc'

    #Apply stylizations
    for node in tree.traverse():
        if node.is_leaf():
            node.set_style(style_L)
        else:
            node.set_style(style_I)
    
    #Annotate root of tree with gene name
    gene_label = et.TextFace(gene_name, ftype="Helvetica", fsize=5,\
                             fgcolor='green')
    tree.add_face(gene_label, column=0, position = "branch-top")
    
    tree.render(file_name + "_tree.png", tree_style=nifty)
    
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Generate 1000 BS replicate \
                                     ML tree with visualizations from \
                                     homologous protein seq FASTA file') 
                                     
    #User input: name of homologous sequence file, gene name                
    parser.add_argument('-h_seq','-seq', type=str,
                        help='Name file containing homologous sequences:\
                        must be in same directory')
         
    parser.add_argument('-gene', '-g', type=str, required=False,
                        help='Name of gene from homologous sequence file')
    
    args = vars(parser.parse_args())
    
    #Call on MUSCLE alignment tool: MUSCLE binary file must be in directory
    subprocess.run(['muscle','-in',str(sys.argv[2]),'-out',str(sys.argv[2])+\
                    '_align.fasta'])
    
    #Call IQ-TREE: automatically detect best model and perform 1000 BS
    subprocess.run(['iqtree','-s',str(sys.argv[2])+'_align.fasta','-m','MFP',\
                    '-B','1000'])
                     
    #Create visualization from .treefile 
    with open(str(sys.argv[2])+'_align.fasta.treefile','r') as file:
        
    	newick_str = file.read().replace('\n', '')
    
    make_tree_vis(newick_str, str(sys.argv[4]), str(sys.argv[2]))
    
    
    
