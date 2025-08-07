# PanGene-O-Meter
Quantifying bacterial genomes' gene-content similarity and gene-content-based de-duplication

Introduction
============
Bacterial genome evolution is shaped to a great extent by horizontal gene transfer, detectable as genes with a presence-absence pattern of variation that does not follow phylogenetic relationships across the majority of the genome. While average nucleotide identity (ANI) is the most common measure to measure bacterial genome similarity, its ability to capture gene-content differences is limited, and is use therefore often misses a major factor for functional variability. **PanGene-O-Meter provide an orthogonal gene-content-based method to quantify similarity between bacterial genomes, leveraging knowledge of pan-genome structure and relationships between orthology groups.** PanGene-O-Meter, provides a higher level of granularity than ANI when classifying bacterial isolates of the same species. Furthermore, PanGene-O-Meter is useful for efficient bacterial genomes clustering, allowing for facile selection of one or several 'representative' genomes, effectively de-duplicating a set of genomes by their gene content. 

## Table of contents
  * [Pipeline overview](#pipeline-overview)
  * [Quick start](#quick-start)
  * [How to run](#how-to-run)
  * [Directory structure and analysis output](#directory-structure-and-analysis-output)
  * [Command line arguments](#command-line-arguments)

## Pipeline overview
1.	Orthology group assignment: Assign predicted genes to orthology groups (pangenes) to constructing a pan-genome (different methods can be used for this step).
2.	Phyletic pattern representation: Represent each genome as a binary presence/absence (P/A) vector, where each element corresponds to a pangene and the value is either 1 (present) or 0 (absent).
3.	Gene-content similarity (GCS) calculation: Calculate GCS based on Jaccard similarity \($GCSj$\) between the P/A vectors of two bacterial genomes \($pS1$ and $pS2$\):

$$GCSj(pS1,pS2)=\frac{M11}{(M01+M10+M11)}$$

Where $M11$ represents the total number of elements (orthology-groups) for which both $pS1$ and $pS2$ have a value of 1; $M01$ represents the total number of elements for which the value of $pS1$ is $0$ and $pS1$ is $1$; and $M10$ represents the total number of elements for which the value of $pS1$ is $1$ and the value of $pS2$ is $0$.</span>

Alternatively, the GCS can be defined by the shared pangenes, based on their overlap coefficient \($GCSo$\).  This approach reduces the impact of partial genomes on the GCS metric.

$$GCSo(pS1,pS2)=\frac{M11}{(min(M01+M11,M10+M11))}$$

4. **GeneContRep**, a greedy incremental clustering algorithm to deduplicate a list of genomes. GeneContRep, allows the selection of representative genomes according to specific gene-content similarity cutoffs.

## Quick start
### Installing
0. If not already installed, install miniconda on your system
1. Create a new environment for PanGeneOmeter
```
conda create -n PanGeneOmeter python==2.7.*
```
**OR** [replace /some/where/ with your desired conda environment path]
```
conda create -p /some/where/PanGeneOmeter python==2.7.*
```
2. Activate the conda environmet
```
conda activate PanGeneOmeter
## OR ##
conda activate -p /some/where/PanGeneOmeter
```
3. Downlaod and install PanGeneOmeter and dependencies
```
git clone https://github.com/HaimAshk/PanGene-O-Meter.git $CONDA_PREFIX/PanGeneOmeter_github
git clone https://github.com/neherlab/pan-genome-analysis.git $CONDA_PREFIX/panx_github
conda env update --file $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter-environment.yml  --prune 
mv $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter.pl $CONDA_PREFIX/bin
cp $CONDA_PREFIX/PanGeneOmeter_github/PanX_Patch/sf_extract_sequences.py $CONDA_PREFIX/panx_github/scripts/
```


Support
=======
PanGene-O-Meter is actively supported and developed software. Please use the [issue tracker](https://github.com/HaimAshk/PanGene-O-Meter/issues) for malfunctions and the [GitHub discussions](https://github.com/HaimAshk/PanGene-O-Meter/discussions) for questions, comments, feature requests, etc.

<img src="https://www.seekpng.com/png/detail/137-1379498_work-in-progress.png" alt="Work In Progress@seekpng.com" height=60 width=60> <br>
**THE CODE IS STILL UNDER DEVELOPMENT -- PLEASE VISIT AGAIN SOON**
