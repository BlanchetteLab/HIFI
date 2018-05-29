# HIFI
A tool for high-resolution estimation of DNA-DNA interaction frequency from Hi-C data

## Overview
HIFI is a set of tools to infer true intra-chromosomal interaction frequencies at restriction fragment resolution from Hi-C data. 
It uses adaptive kernel density estimation and Markov Random Field approaches to provide accurate estimates of interaction frequency matrices from sparse read-count matrices.

Input: A bam file corresponding to Hi-C data from a given chromosome. 

Step 1: Produce a restriction fragment resolution read count matrix, using the ?? program.

Step 2: Produce the restriction fragment resolution true IF estimation, using the HIFI program itself.

Step 3: (optional) Visualize the inferred true IF matrix, using the ?? program. 

Output: Estimated interaction frequency matrix at restriction-fragment resolution

## Software requirements
1) Linux with gcc compiler
2) ???CHRIS???

## Installation guide
Installation is expected to take a few minutes:
1) Download the package by clicking the "Clone or download" button
2) Click "download ZIP", download file, and unzip file in desired location, and rename directory "HIFI"
OR
Copy the git link provided, and use "git clone https://github.com/BlanchetteLab/HIFI" from the command line.
3) cd HIFI
4) make HIFI
?? CHRIS ??

## Quick start
The package includes the bam file corresponding to a subset of ??AUTHORS?? Hi-C data (REFERENCE), limited to region chr?:1234-4567. This is the data we are going to use as example.
1) Process bam file (expected run time: ??):
COMMAND LINE?

2) Run HIFI with default parameters (expected run time: ??):
./HIFI example/example.tsv example/HIFI_output.tsv

3) Visualize the output IF matrix (expected run time: ??):
COMMAND LINE?

## Details:
1) Processing of bam file

2) HIFI

HIFI <sparseInputMatrix.tsv> <sparseOutputMatrix> <options>
Options:
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

Subsetting options
	-firstRow=<INTEGER> : First row to be analyzed
	-lastRow=<INTEGER> : Last row to be analyzed
	-firstCol=<INTEGER> : First column to be analyzed
	-lastCol=<INTEGER> : Last column to be analyzed

Optimizations
	-bandSize=<INTEGER> : Limits the analysis to a band of the given width along the main diagonal. Useful when analyzing very large matrices, to limit the analysis to short range contacts. Default: inactive


3) True-size RF-resolution IF matrix visualization

## Testing
All software was tested on Linux Ubuntu 12.04.5 LTS (GNU/Linux 3.2.0-86-generic x86_64).

## License:
RobusTAD is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.

RobusTAD is distributed in the hopes that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
