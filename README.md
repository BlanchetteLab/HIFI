# HIFI
A tool for high-resolution estimation of DNA-DNA interaction frequency from Hi-C data

## Overview
HIFI is a set of tools to infer true intra-chromosomal interaction frequencies at restriction fragment (RF) resolution from Hi-C data. 
It uses adaptive kernel density estimation and Markov Random Field approaches to provide accurate estimates of interaction frequency matrices from sparse read-count matrices.

Input: A bam file corresponding to Hi-C data from a given chromosome (or part thereof). 

Step 1: Produce a RF-resolution read count matrix, using the BAMtoSparseMatrix.py program.

Step 2: Produce the RF-resolution true IF estimation, using the HIFI program.

Step 3: (optional) Visualize the inferred true IF matrix, using the parseHIFIoutput.py and plotHIFIoutput.py program. 

Output: Estimated interaction frequency matrix at RF resolution

## Software requirements
1) Linux with g++ compiler
2) Samtools: http://www.htslib.org/download/
3) Matplotlib (should be installed with most versions of Python; or can be installed using Anaconda (https://conda.io/docs/user-guide/install/download.html)).
4) NumPy (should be installed with most versions of Python; or can be installed using Anaconda (https://conda.io/docs/user-guide/install/download.html)).

## Installation guide
Installation is expected to take a few minutes:

1) Either download the package by clicking the "Clone or download" button, unziping file in desired location, and renaming the directory "HIFI"   OR   Use the command line ``` git clone https://github.com/BlanchetteLab/HIFI ```.
2) ``` cd HIFI ```
3) ``` make HIFI ```
4) If Matplotlib or NumPy are not installed, install it using Anaconda (https://conda.io/docs/user-guide/install/download.html).
5) If samtools is not installed, install it (http://www.htslib.org/download/) and make sure the samtools executable is in your path.
6) Download into the examples directory the following test Hi-C data set from Rao et al. 2014: https://www.cs.mcgill.ca/~/HIFI/Rao_GM12878.hg19.chr9_example.bam

## Example data set
The bam file linked just above comes from Rao et al. (2014), limited to intrachromosomal contacts in region chr9:122000000-132000000. This is the data we are going to use as example. The bam file was produced from fastq files using HiCUP's standard pipeline to map read pairs to hg19 and perform read-pair quality filtering. Additional filtering (MAPQ value >= 30) ensures unique mappability. 

## Quick start
1) Process bam file to produce input to HIFI (expected run time: 1 minute):

```src/BAMtoSparseMatrix.py examples/Rao_GM12878.hg19.chr9_example.bam examples/hg19.HindIII_fragments.bed ./examples_output```

2) Run HIFI with default parameters (expected run time: 15 minutes):

``` src/HIFI examples_output/Rao_GM12878.hg19.chr9_example.chr9_chr9.RF.tsv examples_output/Rao_GM12878.hg19.chr9_example.chr9_chr9.RF.HIFI_MRF.tsv -method=mrf```

3) Extract a subset of the IF matrix for visualization (positions 125000000-129000000) (expected run time: 1 minute)

```python src/parseHIFIoutput.py examples_output/Rao_GM12878.hg19.chr9_example.chr9_chr9.RF.HIFI_MRF.tsv examples/hg19.HindIII_fragments.bed 125000000 129000000 examples_output```

4) Visualize the output IF matrix (expected run time: 15 minutes):

```python src/plotHIFIoutput.py examples_output/Rao_GM12878.hg19.chr9_example.chr9_chr9.RF.HIFI_MRF.125000000_129000000.tsv examples/hg19.HindIII_fragments.bed 0.0 1.5 examples_output```

## Command line details:
1) Processing of bam file

```
BAMtoSparseMatrix.py [-h] bam_filepath digest_filepath output_dir

positional arguments:
  bam_filepath     input BAM file path
  digest_filepath  expected restriction fragment digest file (BED) path
  output_dir       file path to output directory

optional arguments:
  -h, --help       show this help message and exit
```

2) HIFI
```
HIFI <sparseInputMatrix.tsv> <sparseOutputMatrix.tsv> <options>

Required arguments:
Choice of inference method
	-method=[fixed|kde|mrf<default>]
		fixed: fixed resolution (binning) approach; use with -fragmentsPerBin
		kde: kernel density estimation; use with -kdeBandWidth
		akde: adaptive kernel density estimation; use with -kdeMinCount and -kdeMaxBandwidth
		mrf: Markov Random Field estimation
	
Optional arguments:
Estimation method:
	-method=[fixed|kde|akde|mrf]
		fixed: fixed resolution (binning) approach; use with -fragmentsPerBin
		kde: kernel density estimation; use with -kdeBandWidth
		akde: adaptive kernel density estimation; use with -kdeMinCount and -kdeMaxBandwidth
		mrf: Markov Random Field estimation
		Default: mrf
	
Fixed binning options
	-fragmentsPerBin=<INTEGER> : Number of fragments per bin. Default: 1

KDE options
	-kdeBandwidth=<FLOAT> : Bandwidth of the KDE. Default: 20.000000

Adaptive KDE options
	-kdeMinCount=<INTEGER> : Minimum read pair count for KDE. Default: 50
	-kdeMaxBandwidth=<FLOAT> : Maximum bandwidth to be used in AKDE. Default: 25.000000

MRF options
	-boundaryKS=<FLOAT> : Minimum value of the KS test statistic to call a boundary. Default: 2.000000
	-mrfMaxIter=<INTEGER> : Number of iterations of MRF. Default: 5
	-sigmaMultiplier=<FLOAT> : Sigma multiplier term in MRF (see Methods). Default: 0.040000
	-minSigma=<FLOAT> : Minimum value of sigma to be considered (see Methods). Default: 0.001000

Output options
	-minOutput=<FLOAT> : Smallest IF value to be outputed. IF values below this number are not outputed. Defaul: 0.000010
	-outputNotNormalized : Does not output bias-corrected IFs, but instead produces non-bias-corrected IFs.

Subsetting options
	-firstRow=<INTEGER> : First row to be analyzed
	-lastRow=<INTEGER> : Last row to be analyzed
	-firstCol=<INTEGER> : First column to be analyzed
	-lastCol=<INTEGER> : Last column to be analyzed

Optimizations
	-bandSize=<INTEGER> : Limits the analysis to a band of the given width along the main diagonal. Useful when analyzing very large matrices, to limit the analysis to short range contacts. Default: inactive
```

3) Extraction of a subset of the IF matrix for visualization
```
parseHIFIoutput.py [-h]
                          HIFI_output digest_filepath bp_start bp_end
                          output_dir

positional arguments:
  HIFI_output      file path to HIFI sparse matrix
  digest_filepath  expected restriction fragment digest (BED) filepath
  bp_start         starting base pair for region of interest
  bp_end           final base pair for region of interest
  output_dir       file path to output directory

optional arguments:
  -h, --help       show this help message and exit
```

4) Visualize the output IF matrix 
```
plotHIFIoutput.py [-h] HIFI_output digest_filepath vmin vmax output_dir

positional arguments:
  HIFI_output      file path to HIFI sparse matrix (recommend filtering before
                   plotting - see 'parseHIFIoutput.py')
  digest_filepath  expected restriction fragment digest (BED) filepath
  vmin             minimum value to be used for colour palette range
  vmax             maximum value to be used for colour palette range
  output_dir       file path to output directory

optional arguments:
  -h, --help       show this help message and exit
```

## Testing
All software was tested on Linux Ubuntu 12.04.5 LTS (GNU/Linux 3.2.0-86-generic x86_64).

## License:
HIFI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.

HIFI is distributed in the hopes that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
