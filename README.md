
<h1> Project Description </h1>
A bioinformatics software pipeline to simplify the analysis of gene neighborhoods utilizing tools such as Prokka, RGI, Bioconda, BLAST, MUSCLE, IQ-TREE, ete3, and DendroPy. Priority is placed on creating cluster and phylogenetic visualizations from the identified gene neighborhoods. The pipeline is currently tailored to examine 5 genes upstream and 5 genes downstream from an annotated AMR gene.

<h2> Dependencies </h2>
Currently, the pipeline has been designed and tested on Linux.


All external software/package dependencies are checked at the beginning of the workflow, and if given user permission, will download automatically prior to running. Please note that Python must be version 3.5 or higher to avoid the need for backwards compatibility.

Links to repositories used:

prokka: https://github.com/tseemann/prokka

RGI: https://github.com/arpcard/rgi

Bioconda: https://github.com/bioconda

IQ-TREE: https://github.com/Cibiv/IQ-TREE

ete3: https://github.com/etetoolkit/ete

DendroPy: https://github.com/jeetsukumaran/DendroPy

<h2> Structure </h2> 

```bash
.
├── config
├── docs
│   └── pipeline_visualization.png
├── README.md
├── requirements.txt
├── scripts
│   ├── pipeline_prototype2.sh
│   ├── pipeline_protoype1.sh
│   └── pipeline.sh
├── setup.py
├── test
└── workflow
    ├── envs
    │   └── gene_neighborhoods_pipeline.yml
    └── scripts
        ├── make_ML_tree_vis.py
        ├── __pycache__
        │   └── make_ML_tree_vis.cpython-38.pyc
        ├── tree_clustering.py
        └── tree_distances.py
```

<h2> Authors </h2>
Chandana @ (chandana277)

Julia Lewandowski @ (JTL-lab)
