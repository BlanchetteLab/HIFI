// Mathieu Blanchette Copyright (2019)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Just large enough to handle HindIII data in human
#define MATRIX_SIZE 72000
#define MAX_DISTANCE 1500
#define DIAGONAL_BIN_SIZE 10000
double min(double a, double b){
  return a<b?a:b;
}

double min4(double a, double b, double c, double d){
  return min(min(a,b),min(c,d));
}

int compare_floats (const float *a, const float *b)
{
  return (int) (*a - *b);
}

double **cum;
double **cumE;
int ** nbGood;

class BandedMatrixSym {
 public:
  float *mat[MATRIX_SIZE];
  BandedMatrixSym(int size, int b) {
    for (int i=0;i<size;i++) {
      mat[i]=(float*)malloc(b*sizeof(int));
    }
  }

  float get(int i, int j) {
    if (j<i) {
      fprintf(stderr,"ERROR get %d %d\n",i,j);
      exit(1);
    }
      
    if (j-i<MAX_DISTANCE)   return mat[i][j-i];    
    else return 0;
  }

  float set(int i, int j, float x) {
    if (j<i) {
      fprintf(stderr,"ERROR set\n");
      exit(1);
    }
      
    if (j-i<MAX_DISTANCE) mat[i][j-i]=x;
    else {
      fprintf(stderr,"ERROR2 set\n");
      exit(1);
      
    }
  }
  
};


void rect(int i1, int j1, int i2, int j2, double &sum, int &n, double **c) {
  sum=c[i2][j2]-c[i1-1][j2]-c[i2][j1-1]+c[i1-1][j1-1];
  n=nbGood[i2][j2]-nbGood[i1-1][j2]-nbGood[i2][j1-1]+nbGood[i1-1][j1-1];
}


int main(int argc, char *argv[]) {

  if (argc<7) {
    fprintf(stderr,"callPeaks <interactionFrequencyFile.tsv> <peakSize> <windowSize> <chromosome> <restrictionEnzymeCutSites.bed> <mergeDistance>\n");
  }
  
  FILE *f=fopen(argv[1],"r");
  if (!f) {fprintf(stderr,"Can't open interaction frequency file %s",argv[1]);exit(1);}
  int p=atoi(argv[2]);
  int w=atoi(argv[3]);
  char  chr[100];  
  sprintf(chr,"%s",argv[4]);
  FILE *rf=fopen(argv[5],"r");
  if (!rf) {fprintf(stderr,"Can't open restriction enzyme cut sites file %s",argv[5]);exit(1);}
  double  mergeDist=atof(argv[6]);
  float multFactor=1;
  if (argc>7) multFactor=atof(argv[7]);
  float minScore=1;
  int localDistanceForBadRows=100;
  float minRowZscore=-1.5;
  



  fprintf(stderr,"Reading RF file");
  int RF_end[MATRIX_SIZE];
  int rf_count=0;
  char line[100];
  int RF_size_sum=0;
  float mean_RF_size=0;

  while (fscanf(rf,"%[^\n]\n",line)!=EOF) {
    char c[100];
    int s,e;
    sscanf(line,"%s %d %d",c,&s,&e);
    if (!strcmp(c,chr)){
      RF_end[rf_count]=e;
      if (e-s<100000) {
	RF_size_sum+=e-s; // exclude telomerese, centromeres...
      }
      rf_count++;}
  }
  mean_RF_size=(float)RF_size_sum/rf_count;
  fclose(rf);

    fprintf(stderr,"Allocating memory\n");
  float **m=(float**)malloc(MATRIX_SIZE*sizeof(float*));
  for (int i=0;i<MATRIX_SIZE;i++) {m[i]=(float*)malloc(MATRIX_SIZE*sizeof(float));memset(m[i],0,MATRIX_SIZE*sizeof(float));}

  int size = rf_count+100;

  BandedMatrixSym mE=BandedMatrixSym(size,MAX_DISTANCE);
  
  char **valid=(char**)malloc((size)*sizeof(char*));
  for (int i=0;i<(size);i++) {valid[i]=(char*)malloc((size)*sizeof(char));memset(valid[i],0,(size)*sizeof(char));}

  char ***output=(char***)malloc((size)*sizeof(char**));
  for (int i=0;i<(size);i++) {output[i]=(char**)malloc((size)*sizeof(char*));memset(output[i],0,(size)*sizeof(char*));}

  cum=(double**)malloc((size)*sizeof(double*));
  for (int i=0;i<(size);i++) cum[i]=(double*)malloc((size)*sizeof(double));

  cumE=(double**)malloc((size)*sizeof(double*));
  for (int i=0;i<(size);i++) {cumE[i]=(double*)malloc((size)*sizeof(double));memset(cumE[i],0,(size)*sizeof(double));}

  nbGood=(int**)malloc((size)*sizeof(int*));
  for (int i=0;i<(size);i++) nbGood[i]=(int*)malloc((size)*sizeof(int));


  BandedMatrixSym scores=BandedMatrixSym(size,MAX_DISTANCE);
  fprintf(stderr,"Done allocating memory\n");


  fprintf(stderr,"Reading matrix...");
  float sumRow[MATRIX_SIZE];
  int bad[MATRIX_SIZE];


  float sumDiagTrue[100000];
  int nDiagTrue[100000];
  memset(sumDiagTrue,0,100000*sizeof(float));
  memset(nDiagTrue,0,100000*sizeof(int));
  
  
  while (fscanf(f,"%[^\n]\n",line)!=EOF) {
    if (line[0]!='#') {
    int i,j,foo;
    double x;
    sscanf(line,"%lf %d %d %d %d",&x,&foo,&i,&foo,&j);
    if (i>=MATRIX_SIZE-localDistanceForBadRows || j>=MATRIX_SIZE-localDistanceForBadRows) continue;
    float xx=x/1000000;
    xx*=multFactor;
    m[i][j]=xx;
    
    int diag = (RF_end[j]-RF_end[i])/DIAGONAL_BIN_SIZE;
    sumDiagTrue[diag]+=xx;
    nDiagTrue[diag]++;
    
    if (j-i<localDistanceForBadRows) {
      sumRow[i]+=xx;
      sumRow[j]+=xx;
    }
    
    if (i>rf_count) rf_count=i;
    if (j>rf_count) rf_count=j;
    }
  }
  fclose(f);

  fprintf(stderr,"Looking for badrows...\n");

  // Get the mean and variance of sumRows
  double sum=0;
  double sum2=0;
  for (int i=0;i<rf_count;i++) {sum+=sumRow[i];sum2+=sumRow[i]*sumRow[i];}
  float mean = sum/rf_count;
  float stddev=sqrt(sum2/rf_count - sum*sum/rf_count/rf_count);
  //  fprintf(stderr,"mean = %lf, stddev = %lf\n",mean,stddev);
  int nbBad=0;  
  for (int i=0;i<rf_count;i++) {
    float z = (sumRow[i]-mean)/stddev;
    //fprintf(stderr,"z=%f\n",z);
    if (z < minRowZscore) {
      //      fprintf(stderr,"bad %f %lf\n",sumRow[i],z);
      bad[i]=1;
      nbBad++;
    }
  }
  fprintf(stderr,"Nb bad = %d\n",nbBad);

  // calculate expected per diagonal
  fprintf(stderr,"Calculating E matrix\n");
  for (int i=0;i<rf_count;i++) {
    for (int j=i+1;j<i+MAX_DISTANCE;j++) {
      int diag=(RF_end[j]-RF_end[i])/DIAGONAL_BIN_SIZE;
      float x=sumDiagTrue[diag]/(nDiagTrue[diag]+0.001);
      mE.set(i,j,x);
    }
  }
  
  fprintf(stderr,"Calculating cumulative matrix\n");

  cum[0][0]=m[0][0];
  cumE[0][0]=m[0][0];
  for (int i=1;i<rf_count+100;i++) {
    if (!bad[i] && !bad[0]) {
      cum[i][0]=0+cum[i-1][0];
      cumE[i][0]=0+cumE[i-1][0];
      nbGood[i][0]=nbGood[i-1][0];
    }
    else {
      cum[i][0]=cum[i-1][0];
      cumE[i][0]=cumE[i-1][0];
      nbGood[i][0]=nbGood[i-1][0];
    }
  }
  
  for (int i=1;i<rf_count+100;i++) {
    if (!bad[0] && !bad[i]) {
      cum[0][i]=m[0][i]+cum[0][i-1];
      cumE[0][i]=mE.get(0,i)+cumE[0][i-1];
      nbGood[0][i]=1+nbGood[0][i-1];
    }
    else {
      cum[0][i]=cum[0][i-1];
      cumE[0][i]=cumE[0][i-1];
      nbGood[0][i]=nbGood[0][i-1];
    }
  }

  for (int i=1;i<rf_count+100;i++) {
    for (int j=1;j<rf_count+100;j++) {
      if(j>i && !bad[i] && !bad[j]) {
	cum[i][j]=m[i][j]+cum[i-1][j]+cum[i][j-1]-cum[i-1][j-1];
	cumE[i][j]=mE.get(i,j)+cumE[i-1][j]+cumE[i][j-1]-cumE[i-1][j-1];
	nbGood[i][j]=1+nbGood[i-1][j]+nbGood[i][j-1]-nbGood[i-1][j-1];
	
      }
      else {
	cum[i][j]=cum[i-1][j]+cum[i][j-1]-cum[i-1][j-1];
	cumE[i][j]=cumE[i-1][j]+cumE[i][j-1]-cumE[i-1][j-1];
	nbGood[i][j]=nbGood[i-1][j]+nbGood[i][j-1]-nbGood[i-1][j-1];
      }
    }
  }

  int peakSize=0.25*(2*p/mean_RF_size+1)*(2*p/mean_RF_size+1);
  if (peakSize<10) peakSize=10;
  int donutSize=peakSize;
  
  printf("#chr1\tx1\tx2\tchr2\ty1\ty2\tpeak\tscore\tstrand1\tstrand2\tcolor\tobserved\texpected\tBL\tDonut\tH\tV\tratioBL\tratioDonut\tratioH\tratioV\n");
  printf("# juicer_tools version 1.11.09\n");
  fprintf(stderr,"Computing scores...\n");
  for (int i=1;i<rf_count;i++) {

    if (i>1000) {
      free(cum[i-1000]);
      free(cumE[i-1000]);
      free(nbGood[i-1000]);
    }
    

    for (int j=i+1;j<rf_count;j++) {
      if (m[i][j]>0) {
	int l_peak,l_window,r_peak,r_window;
	int t_peak,t_window,b_peak,b_window;
	int cur_i=RF_end[i];
	int  cur_j=RF_end[j];
	l_peak=j;
	while (l_peak>1 && cur_j-RF_end[l_peak]<p) l_peak--;
	r_peak=j;
	while (r_peak<rf_count &&  RF_end[r_peak]-cur_j<p) r_peak++;
	l_window=j;
	while (l_window>1 && cur_j-RF_end[l_window]<w) l_window--;
	r_window=j;
	while (r_window<rf_count &&  RF_end[r_window]-cur_j<w) r_window++;

	t_peak=i;
	while (t_peak>1 && cur_i-RF_end[t_peak]<p) t_peak--;
	b_peak=i;
	while (b_peak<rf_count &&  RF_end[b_peak]-cur_i<p) b_peak++;
	t_window=i;
	while (t_window>1 && cur_i-RF_end[t_window]<w) t_window--;
	b_window=i;
	while (b_window<rf_count &&  RF_end[b_window]-cur_i<w) b_window++;

	double full, h, v, bl, peak, donut;
	int full_n, h_n, v_n, bl_n, peak_n,donut_n;
	double aboveslice,belowslice,leftslice,rightslice;
	int aboveslice_n,belowslice_n,leftslice_n,rightslice_n;

	double efull, eh, ev,ebl, epeak, edonut;
	int efull_n, eh_n, ev_n, ebl_n, epeak_n,edonut_n;
	double eaboveslice,ebelowslice,eleftslice,erightslice;
	int eaboveslice_n,ebelowslice_n,eleftslice_n,erightslice_n;

	rect(t_window,l_window,b_window,r_window,full,full_n,cum);rect(t_window,l_window,b_window,r_window,efull,efull_n,cumE);
	rect(t_peak,l_peak,b_peak,r_peak,peak,peak_n,cum);	rect(t_peak,l_peak,b_peak,r_peak,epeak,epeak_n,cumE);

	if (peak_n==0) continue;
	
	rect(i,l_window,i,l_peak-1,leftslice,leftslice_n,cum);	rect(i,l_window,i,l_peak-1,eleftslice,eleftslice_n,cumE);
	rect(i,r_peak+1,i,r_window,rightslice,rightslice_n,cum);	rect(i,r_peak+1,i,r_window,erightslice,erightslice_n,cumE);

	rect(t_window,j,t_peak-1,j,aboveslice,aboveslice_n,cum);	rect(t_window,j,t_peak-1,j,eaboveslice,eaboveslice_n,cumE);
	rect(b_peak+1,j,b_window,j,belowslice,belowslice_n,cum);	rect(b_peak+1,j,b_window,j,ebelowslice,ebelowslice_n,cumE);

	donut=(full-peak-leftslice-rightslice-aboveslice-belowslice);	edonut=(efull-epeak-eleftslice-erightslice-eaboveslice-ebelowslice);
	donut_n=(full_n-peak_n-leftslice_n-rightslice_n-aboveslice_n-belowslice_n);	edonut_n=(efull_n-epeak_n-eleftslice_n-erightslice_n-eaboveslice_n-ebelowslice_n);

	// horizontal
	rect(t_window,l_peak,b_window,r_peak,h,h_n,cum);	rect(t_window,l_peak,b_window,r_peak,eh,eh_n,cumE);
	double hor=(h-peak);	double ehor=(eh-epeak);
	int hor_n=h_n-peak_n; 	int ehor_n=eh_n-epeak_n;
	
	// vertical
	rect(t_peak,l_window,b_peak,r_window,v,v_n,cum);	rect(t_peak,l_window,b_peak,r_window,ev,ev_n,cumE);
	double ver=(v-peak);	double ever=(ev-epeak);
	int ver_n=v_n-peak_n;	int ever_n=ev_n-epeak_n;
	
	// bottom left
	double bl1,bl2;	double ebl1,ebl2;
	int bl1_n,bl2_n;	int ebl1_n,ebl2_n;
	rect(i+1,l_window,b_window,j-1,bl1,bl1_n,cum);	rect(i+1,l_window,b_window,j-1,ebl1,ebl1_n,cumE);
	rect(i+1,l_peak,b_peak,j-1,bl2,bl2_n,cum);	rect(i+1,l_peak,b_peak,j-1,ebl2,ebl2_n,cumE);
	bl=bl1-bl2;	ebl=ebl1-ebl2;
	bl_n=bl1_n-bl2_n;	ebl_n=ebl1_n-ebl2_n;
	
	double ratio_donut=donut/edonut;
	double ratio_hor=hor/ehor;
	double ratio_ver=ver/ever;
	double ratio_bl=bl/ebl;
	double ratio_peak=peak/epeak;

	double score=min4(ratio_peak/ratio_donut,ratio_peak/ratio_hor,ratio_peak/ratio_ver,ratio_peak/ratio_bl);
	scores.set(i,j,score);

	if (score>=minScore && ratio_peak>=1 &&  donut_n>=donutSize && hor_n>=1 && ver_n>=1 && bl_n>=1 && peak_n>=peakSize ) {
	  
	  valid[i][j]=1;
	  char *out=(char*) malloc(400);
	  
	  float maxi=-100;
	  int maxi_i,maxi_j;
	  for (int a=t_peak;a<=b_peak;a++) {
	    for (int b=l_peak;b<=r_peak;b++) {
	      if (b>a && !bad[a] && !bad[b] && m[a][b]/mE.get(a,b)>maxi) {
		maxi=m[a][b]/mE.get(a,b);
		maxi_i=a;
		maxi_j=b;
	      }	  
	    }
	  }
//	  sprintf(out,"%s\t%d\t%d\t%s\t%d\t%d\t%f\t%f\t.\t.\t0,0,0\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",chr,RF_end[maxi_i],RF_end[maxi_i]+1,chr,RF_end[maxi_j], RF_end[maxi_j]+1,ratio_peak,scores.get(i,j),m[maxi_i][maxi_j],mE.get(maxi_i,maxi_j),ratio_bl,ratio_donut,ratio_hor,ratio_ver,ratio_peak/ratio_bl,ratio_peak/ratio_donut,ratio_peak/ratio_hor,ratio_peak/ratio_ver, i, j, l_window, r_window,t_window,b_window,l_peak,r_peak,t_peak,b_peak,peak_n,donut_n);
	  sprintf(out,"%s\t%d\t%d\t%s\t%d\t%d\t%f\n",chr,RF_end[maxi_i],RF_end[maxi_i]+1,chr,RF_end[maxi_j], RF_end[maxi_j]+1,scores.get(i,j));
	  output[i][j]=out;
	}
	else {
	  valid[i][j]=0;
	}
      }
    }
  }
  
  // Merge and output peaks
  for (int i=0;i<rf_count;i++) {
    for (int j=i+1;j<rf_count;j++) {
      if (valid[i][j]) {
	
	// check if there isn't a better peak with mergeDist
	int left=i-1;
	int right=i+1;
	int top=j-1;
	int bottom=j+1;
	while (left>0 && fabs(RF_end[i]-RF_end[left])<mergeDist) left--;
	while (right<rf_count && fabs(RF_end[i]-RF_end[right])<mergeDist) right++;
	while (top>0 && fabs(RF_end[j]-RF_end[top])<mergeDist) top--;
	while (bottom<rf_count && fabs(RF_end[j]-RF_end[bottom])<mergeDist) bottom++;

	float maxscore=-9999;
	for (int ii=left;ii<=right;ii++)  {
	  for (int jj=top;jj<=bottom;jj++)  {
	    if (jj>ii && (ii!=i || jj!=j)) {
	      float x=scores.get(ii,jj);
	      if (x>maxscore) maxscore=x;
	      
	    }
	  }
	}

	if (scores.get(i,j)>maxscore) {
	  printf("%s",output[i][j]);
	}
      }
    }
  }
}
