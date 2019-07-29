class HIFI_options {
 public:
  char method[100];
  bool method_mrf;
  bool method_fixed;
  bool method_kde;
  bool method_akde;
  int bandSize;

  // fixed binning
  int fragmentsPerBin;

  // KDE
  double kdeBandwidth;

  // adaptive KDE
  double kdeMaxBandwidth;
  int kdeMinCount;

  // MRF
  double minSigma;
  double sigmaMultiplier;
  int mrfMaxIter;

  // normalization
  bool noBias;
  bool outputNormalized;

  // boundaries
  double boundaryKS;
  char boundaryOutput[1000];

  // output
  double minOutput;

  // internal
  bool diagRenormalization;
  int firstRow;
  int firstCol;
  int lastRow;
  int lastCol;

  // user-define normalization matrix
  char normalizationMatrix[1000];

  
  HIFI_options() {
    method_mrf=true;
    method_fixed=false;
    method_kde=false;
    fragmentsPerBin=1;
    kdeBandwidth=20;
    kdeMaxBandwidth=25;
    minSigma=0.001;
    sigmaMultiplier=0.04;
    kdeMinCount=50;
    noBias=false;
    boundaryKS=2.0;
    diagRenormalization=true;
    outputNormalized=true;
    minOutput=0.00001;
    mrfMaxIter=5;
    firstRow=-1;
    firstCol=-1;
    lastRow=-1;
    lastCol=-1;
    bandSize=99999999;
    normalizationMatrix[0]=0;
  }

  void read(int argc, char *argv[]);
  void print();
  void help();
};
