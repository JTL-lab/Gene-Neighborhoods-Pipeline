#!/bin/bash

echo "Beginning phylo pipeline: converting homologous sequences to trees..."

#Store names of .treefiles created
declare -a files

#Create .treefiles from each homologous sequence file in directory
for file in ./*.fasta; do
    #Make tree visualizations
    echo "Generating .treefile..."
    output=$(python3 make_ML_tree_vis.py -h_seq $file)
    files+=($output)
done

#Show list of .treefiles
echo ".treefiles created:"
for filename in ${files[@]}; do echo $filename; done

#Make directory to move .treefiles into
if [ ! -d treefiles ]; then
  mkdir treefiles;
  chmod -R o+rw treefiles;
fi;

#Move .treefiles
echo "Storing .treefiles in treefiles directory"
for file in ./*.treefile; do
  mv $file treefiles;
done

#Make distance matrices
echo "Creating RF distance and BSD matrices from .treefiles..."
python3 tree_distances.py -tree_path treefiles

#Make UPGMA and Neighbor-Joining tree clusters
echo "Clustering trees from distance matrices..."
python3 tree_clustering.py -matrix_path treefiles
