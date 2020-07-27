from setuptools import setup, find_packages

"""Exception block to modify later...
try:
    from setuptools import setup, find_packages
except ImportError:
    from disutils.core import setup"""

long_description="A bioinformatics pipeline designed to simplify the \
      analysis of gene neighborhoods. Examines 5 genes \
      upstream and 5 genes downstream from an annotated AMR gene, and provides\
      cluster and phylogenetic visualization outputs from the identified \
      gene neighborhoods."
      
setup(
      name="gene_neighborhoods_pipeline",
      version="0.1.0",
      package_dir={"": "src"},
      packages=find_namespace_packages(where="src")),
      description="Bioinformatics pipeline to simplify analysis and \
      visualization of gene neighborhoods."
      long_description=long_description,
      install_requires=["prokka>=1.14.5","rgi>=5.1.1","numpy>=1.19.0",
                        "pandas>=1.0.5","biopython>=1.77","ete3>=3.1.1",
                        "dendropy>=4.4.0"],
      python_requires=">=3.5",
      
      #Things to modify later
      author= "Chandana & Julia",
      scripts="",
      
      
)
