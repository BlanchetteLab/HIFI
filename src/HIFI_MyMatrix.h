
template <class T> class MyMatrix {
 public:
  int firstRow,firstCol,lastRow,lastCol; // inclusive                                                                                  
  int bandSize; // number of fragments for banded analysis
   
  T **data;

  MyMatrix() {};
  MyMatrix(  int fr,int fc, int lr,int lc) {};
  int startRow() {return firstRow;}
  int endRow() {return lastRow+1;}
  int startCol(int row) {return -999;}
  int endCol(int row) {return -999;}
  T get(int i, int j) {};
  void set(int i, int j, T x) {};
  void clear() {};
};


template <class T>
class UpperDiag: public MyMatrix<T> {
 public:

  UpperDiag() {
  }

  int startCol(int row) {
    return max(row,this->firstCol);
  }
  int endCol(int row) {
    return min(this->lastCol+1,row+this->bandSize);
  }

  UpperDiag(int fr,int fc, int lr,int lc, int bs) {
    this->firstRow=fr;                                                                                                                     
    this->firstCol=fc;                                                                                                                     
    this->lastRow=lr;                                                                                                                      
    this->lastCol=lc;                           
    this->bandSize=bs;

    this->data = (T**) malloc((this->lastRow-this->firstRow+1)*sizeof(T*));                                                                 
    for (int i=this->firstRow;i<=this->lastRow;i++) {
      int s=(endCol(i)-startCol(i));
      this->data[i-this->firstRow]=(T*)malloc(s*sizeof(T));
      memset(this->data[i-this->firstRow],0,s*sizeof(T));
    }
  }

  T get(int i, int j) {
    if (i>j || i<this->startRow() || i>=this->endRow() || j<this->startCol(i) || j>=this->endCol(i) || j-startCol(i)>this->bandSize) {
      fprintf(stderr,"ERROR: Getting element outside matrix i=%d j=%d %d %d %d %d\n",i,j,this->startRow() , this->endRow(), this->startCol(i), this->endCol(i));exit(1);}
    return this->data[i-this->firstRow][j-startCol(i)];

  }

  void set(int i, int j, T x) {
    if (i>j || i<this->startRow() || i>=this->endRow() || j<this->startCol(i) || j>=this->endCol(i) || j-startCol(i)>this->bandSize) {
      fprintf(stderr,"ERROR: Setting element outside matrix i=%d j=%d %d %d %d %d\n",i,j,this->startRow() , this->endRow(), this->startCol(i), this->endCol(i));exit(1);}
    this->data[i-this->firstRow][j-startCol(i)]=x;
  }                                                                                                                                         
  void clear() {                                                                                                                           for (int i=this->firstRow;i<=this->lastRow;i++) {                                                                                        free(this->data[i-this->firstRow]);
    }
    free(this->data);
  }
};

class HIFI_matrixProperties {
 public:
  int fullMatrix_firstRow,fullMatrix_firstCol,fullMatrix_lastRow, fullMatrix_lastCol, fullMatrix_size;
  int chr1,chr2;

};

