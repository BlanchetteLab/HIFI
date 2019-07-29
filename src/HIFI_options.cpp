#include "HIFI_options.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void HIFI_options::read(int argc, char *argv[]) {
  for (int i=3;i<argc;i++) {

    if (strstr(argv[i],"-method=")) {
      int x=sscanf(argv[i],"-method=%s",method);
      if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);}
      if (!strcmp(method,"mrf")) {method_mrf=true; continue;}
      if (!strcmp(method,"kde")) {method_kde=true; method_mrf=false; continue;}
      if (!strcmp(method,"akde")) {method_akde=true; method_mrf=false; continue;}
      if (!strcmp(method,"fixed")) {method_fixed=true; method_mrf=false; continue;}
      fprintf(stderr,"Unknown method: %s\nChoices are: mrf, kde, fixed\n",method);
    }

    if (strstr(argv[i],"-fragmentsPerBin=")) {int x=sscanf(argv[i],"-fragmentsPerBin=%d",&fragmentsPerBin);if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} if (!method_fixed) {fprintf(stderr,"-fragmentsPerBin option can only be used with -method=kde\n"); exit(1);} continue;}    
    if (strstr(argv[i],"-kdeBandwidth=")) {int x=sscanf(argv[i],"-kdeBandwidth=%lf",&kdeBandwidth);if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} if (!method_kde) {fprintf(stderr,"-kdeBandwidth option can only be used with -method=kde\n"); exit(1);} continue;}
    if (strstr(argv[i],"-kdeMaxBandwidth=")) {int x=sscanf(argv[i],"-kdeMaxBandwidth=%lf",&kdeMaxBandwidth);if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} continue;}
    if (strstr(argv[i],"-kdeMinCount=")) {int x=sscanf(argv[i],"-kdeMinCount=%d",&kdeMinCount);if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} continue;}
    if (strstr(argv[i],"-boundaryKS=")) {int x=sscanf(argv[i],"-boundaryKS=%lf",&boundaryKS);if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} continue;}

    if (!strcmp(argv[i],"-outputNotNormalized")) {outputNormalized=false; continue;}
    if (!strcmp(argv[i],"-noBias")) {noBias=true; continue;}
    if (strstr(argv[i],"-minOutput=")) {int x=sscanf(argv[i],"-minOutput=%lf",&minOutput); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-mrfMaxIter=")) {int x=sscanf(argv[i],"-mrfMaxIter=%d",&mrfMaxIter); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-firstRow=")) {int x=sscanf(argv[i],"-firstRow=%d",&firstRow); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-lastRow=")) {int x=sscanf(argv[i],"-lastRow=%d",&lastRow); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}

    if (strstr(argv[i],"-firstColumn=")) {int x=sscanf(argv[i],"-firstColumn=%d",&firstCol); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-lastColumn=")) {int x=sscanf(argv[i],"-lastColumn=%d",&lastCol); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-boundaryOutput=")) {int x=sscanf(argv[i],"-boundaryOutput=%s",boundaryOutput); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-minSigma=")) {int x=sscanf(argv[i],"-minSigma=%lf",&minSigma); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-sigmaMultiplier=")) {int x=sscanf(argv[i],"-sigmaMultiplier=%lf",&sigmaMultiplier); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-bandSize=")) {int x=sscanf(argv[i],"-bandSize=%d",&bandSize); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}
    if (strstr(argv[i],"-normalizationMatrix=")) {int x=sscanf(argv[i],"-normalizationMatrix=%s",normalizationMatrix); if (x!=1) {fprintf(stderr,"Error: Can't read %s\n",argv[i]);} ; continue;}

    fprintf(stderr,"Unknown option: %s\n",argv[i]);
    exit(1);
  }
}


void HIFI_options::print() {
  printf("HIFI version 1.0\n");
  printf("option -method= %s\n",method);
  printf("option -fragmentsPerBin= %d\n",fragmentsPerBin);
  printf("option -kdeMinCount = %d\n",kdeMinCount);
  printf("option -kdeBandwidth = %lf\n",kdeBandwidth);
  printf("option -kdeMaxBandwidth = %lf\n",kdeMaxBandwidth);
  printf("option -mrfMaxIter = %d\n",mrfMaxIter);
  printf("option -minSigma = %f\n",minSigma);
  printf("option -sigmaMultiplier = %f\n",sigmaMultiplier);
  printf("option -minSigma = %f\n",minSigma);
  printf("option -firstRow = %d\n",firstRow);
  printf("option -lastRow = %d\n",lastRow);
  printf("option -firstCol = %d\n",firstCol);
  printf("option -lastCol = %d\n",lastCol);
  printf("option -boundaryKS = %f\n",boundaryKS);
  printf("option -minOutput = %lf\n",minOutput);
  printf("option -outputNotNormalized = %d\n",(int)(1-outputNormalized));

  printf("option -bandSize = %d\n",bandSize);
}


void HIFI_options::help() {
  fprintf(stderr,"HIFI <sparseInputMatrix.tsv> <sparseOutputMatrix> <options>\n");
  printf("Required arguments:\n");
  printf("\tsparseInputMatrix.tsv: A tab-separated or space-separated file, where\n\t\t- the first line specifies the index of the restriction fragment \n\t\tcorresponding to the first and last row, and first and last column\n\t\t of the input read count matrix. \n\t\tExample: \"#3 7577 3 7577\"\n");
  printf("\t\t- each subsequent row corresponds to one entry of the matrix, using the format:\n\t\tReadcount_ij Chr_i RF_i Chr_j RF_j.\n\t\tExample: \"3 22 91 22 5005\" means that there are 3 read pairs connecting\n\t\trestriction fragments 91 and 5005 on chromosome 22.\n");
  printf("\tsparseOutputFile.tsv: Name of file where output will be written, in the same format as the input file.\n");


  printf("Options:\n");
  printf("Choice of inference method\n");
  printf("\t-method=[fixed|kde|mrf<default>]\n");
  printf("\t\tfixed: fixed resolution (binning) approach; use with -fragmentsPerBin\n\t\tkde: kernel density estimation; use with -kdeBandWidth\n\t\takde: adaptive kernel density estimation; use with -kdeMinCount and -kdeMaxBandwidth\n\t\tmrf: Markov Random Field estimation\n"); 
  printf("\nFixed binning options\n");
  printf("\t-fragmentsPerBin=<INTEGER> : Number of fragments per bin. Default: %d\n",fragmentsPerBin);

  printf("\nKDE options\n");
  printf("\t-kdeBandwidth=<FLOAT> : Bandwidth of the KDE. Default: %lf\n",kdeBandwidth);

  printf("\nAdaptive KDE options\n");
  printf("\t-kdeMinCount=<INTEGER> : Minimum read pair count for KDE. Default: %d\n",kdeMinCount);
  printf("\t-kdeMaxBandwidth=<FLOAT> : Maximum bandwidth to be used in AKDE. Default: %lf\n",kdeMaxBandwidth);

  printf("\nMRF options\n");
  printf("\t-boundaryKS=<FLOAT> : Minimum value of the KS test statistic to call a boundary. Default: %f\n",boundaryKS);
  printf("\t-mrfMaxIter=<INTEGER> : Number of iterations of MRF. Default: %d\n",mrfMaxIter);
  printf("\t-sigmaMultiplier=<FLOAT> : Sigma multiplier term in MRF (see Methods). Default: %f\n",sigmaMultiplier);
  printf("\t-minSigma=<FLOAT> : Minimum value of sigma to be considered (see Methods). Default: %f\n",minSigma);

  printf("Output options\n");
  printf("\t-minOutput=<FLOAT> : Smallest IF value to be outputed. IF values below this number are not outputed. Defaul: %lf\n",minOutput);

  printf("\nSubsetting options\n");
  printf("\t-firstRow=<INTEGER> : First row to be analyzed\n");
  printf("\t-lastRow=<INTEGER> : Last row to be analyzed\n");
  printf("\t-firstCol=<INTEGER> : First column to be analyzed\n");
  printf("\t-lastCol=<INTEGER> : Last column to be analyzed\n");

  printf("\nOptimizations\n");
  printf("\t-bandSize=<INTEGER> : Limits the analysis to a band of the given width along the main diagonal. Useful when analyzing very large matrices, to limit the analysis to short range contacts. Default: inactive\n");
}

