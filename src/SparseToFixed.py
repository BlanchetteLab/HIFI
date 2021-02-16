#! /usr/bin/env python2.7
#
#   SparseToFixed.py
#   Produces fixed-binning sparse matrix file from sparse matrix containing either read counts or interaction frequencies (TSV format).
#
#   Required arguments:
#       1) sparse matrix file path
#       2) expected digest file path
#       3) fixed bin size (e.g., 1, 5, 10, or 25 kb) --- must be an integer
#       4) output file path (short with score format: https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format)
#   Optional arguments:
#       --start     bin interactions by fragment start instend of end (i.e., fend)
#
#   Requires Python NumPy module.
#
#   Expected TSV format for digest file:
#       chromosome  start   end
#       *start and end expected to be one-based
#       **ensure chromosome names match those found in the sparse matrix file
#       ***file may include more than three columns, but only the first three will be considered
#
#   Steps:
#       1) determine sparse matrix dimensions
#       2) collect required fragment start or end loci
#       3) bin interactions by desired fixed intervals
#       4) write fixed binning interaction frequencies to storage in short with score format
#
#   Author: Christopher JF Cameron
#

import argparse,gzip,os,re,sys
import numpy as np
np.set_printoptions(precision=15)

def get_matrix_dimensions(filepath):
    """returns sparse matrix dimensions"""
    row_min,row_max,col_min,col_max = None,None,None,None
    with open(filepath,'rt') as f:
        for line in f:
            try: freq,row_chrom,row,col_chrom,col = line.rstrip().split()
            except ValueError: print "Error - sparse matrix file is not in the expected format"
            row,col = int(row),int(col)
            row_min = row if row_min == None or row < row_min else row_min
            row_max = row if row_max == None or row > row_max else row_max
            col_min = col if col_min == None or col < col_min else col_min
            col_max = col if col_max == None or col > col_max else col_max
            
    return row_min,row_max,col_min,col_max

def get_chrom_fends(filepath,chrom,frag_start):
    """returns expected fragment ends (or fends) for provided chromosome"""
    fends = []
    pattern = ''.join(['^',chrom,'\s'])
    with open(filepath,'rt') as f:
        for line in f:
            if re.match(pattern,line):   #   only retain chromosome of interest
                chrom,start,end = line.rstrip().split()[:3]
                fends.append(int(start) if frag_start else int(end))

    return sorted(fends,key=int)

parser=argparse.ArgumentParser()
parser.add_argument("sparse_matrix_filepath",help="input sparse matrix file path",type=str)
parser.add_argument("digest_filepath",help="expected restriction fragment digest file (BED-like) path",type=str)
parser.add_argument("resolution",help="fixed binning resolution (must be an integer)",type=str)
parser.add_argument("output_filepath",help="output (fixed bins) sparse matrix file path",type=str)
parser.add_argument("--start",help="use fragment start instead of end when binning",default=False, action='store_true')
parser.add_argument("--no_score",help="output Juicebox 'short format'",default=False, action='store_true')
args=parser.parse_args()

assert(int(args.resolution)),''.join(["Error - provided resolution '",str(args.resolution),"' is not an integer"])
args.resolution = int(args.resolution)
#   ensure output directory exists
out_dir = os.path.dirname(args.output_filepath)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
del out_dir
                                                                                                                                                
##
#   parse sparse matrix file
##

fend_dict = {}
FILE = gzip.open(args.sparse_matrix_filepath,'rt') if args.sparse_matrix_filepath.endswith(".gz") else open(args.sparse_matrix_filepath,'rt')
#   parse header line
header = FILE.readline().rstrip()
if header.startswith('#'):
    row_min,row_max,col_min,col_max = [int(val) for val in header.replace('#','').split()]
else:
    #   no header line present in file
    print "Warning - header not found in sparse matrix file. Parsing file to determine dimensions."
    row_min,row_max,col_min,col_max = get_matrix_dimensions(args.sparse_matrix_filepath)
    FILE.seek(0)

cis = False
float_16 = False
for i,line in enumerate(FILE):
    try: freq,row_chrom,row,col_chrom,col = line.rstrip().split()
    except ValueError: print "Error - sparse matrix file is not in the expected format"
    row_chrom = row_chrom if row_chrom.startswith("chr") else ''.join(["chr",row_chrom])
    col_chrom = col_chrom if col_chrom.startswith("chr") else ''.join(["chr",col_chrom])

    #   get required fragment ends (or fends) for observed chromosomes
    try: fend_dict[row_chrom]
    except KeyError:
        unknown_chrom = False
        for chrom in set([row_chrom,col_chrom]):
            chrom = chrom if chrom.startswith("chr") else ''.join(["chr",chrom])
            print >> sys.stderr,(''.join(["Importing '",chrom,"' expected fends ... "])),
            fend_dict[chrom] = get_chrom_fends(args.digest_filepath,chrom,args.start)
            if len(fend_dict[chrom]) == 0:
                print >> sys.stderr,(''.join(["Warning - unknown chromosome '",chrom,"' encountered and interaction ignored"])),
                del fend_dict[chrom]
                unknown_chrom = True
            print >> sys.stderr,("done")
        if unknown_chrom:
            continue
        if row_chrom == col_chrom:
            cis = True
            
        print >> sys.stderr,(''.join(["Initializing frequncy matrix ... "])),
        #   initialize fixed binned matrix
        max_val = max(fend_dict[row_chrom][row_max]//args.resolution,
            fend_dict[col_chrom][col_max]//args.resolution)
        n = m = max_val+1
        #row_min,row_max = (fend_dict[row_chrom][val]//args.resolution for val in [row_min,row_max])
        #col_min,col_max = (fend_dict[col_chrom][val]//args.resolution for val in [col_min,col_max])
        #n = abs(row_max-row_min)+1
        #m = abs(col_max-col_min)+1
        try: 
            matrix = np.zeros((n,m),dtype=np.float64)
            total = np.float64(0.0)
        except MemoryError:
            print >> sys.stderr,("Warning - not enough available memory to allocate for a matrix of datatype 'np.float64'. Trying 'np.float16' instead. Expect a loss of precision in matrix values."),
            matrix = np.zeros((n,m),dtype=np.float16)
            float_16 = True
            total = np.float16(0.0)
        print >> sys.stderr,("done")

        print >> sys.stderr,(''.join(["Parsing sparse matrix file ... "])),
    
    freq = np.float16(freq) if float_16 else np.float64(freq)
    total += freq
    row,col = int(row),int(col)
    if col < row and cis:
        print "Warning - lower triangle entries of the matrix have been encountered and will be ignored"
        total -= freq
        continue 
    elif col == row and cis:
        print "Warning - main diagonal entries of the matrix have been encountered and will be ignored"
        total -= freq
        continue
    #   bin interaction values
    row = (fend_dict[row_chrom][row]//args.resolution)
    col = (fend_dict[col_chrom][col]//args.resolution)
    matrix[row,col] += freq
FILE.close()
del fend_dict,float_16

assert(np.isclose(total,np.sum(matrix,dtype=np.float64))),"Error - matrix sum does not equal input frequency sum"
print >> sys.stderr,("done")

### print number of reads lost by ignoring the main diagonal
#print np.sum(matrix.diagonal(),dtype=np.float64)

##
#   write binned frequency values to storage
##

total = 0
print >> sys.stderr,(''.join(["Writing binned values to storage ... "])),
with open(args.output_filepath,'wt') as o:
    #   write interaction lines
    indices = zip(*matrix.nonzero())
    for (row,col) in indices:
        if cis and row > col:
            continue
        freq = matrix[row,col]
        if args.no_score:
            #   short format: https://github.com/aidenlab/juicer/wiki/Pre#short-format
            count = freq
            while count > 0.0:
                o.write('\t'.join(['0',row_chrom,str((row*args.resolution)+1),str(row),'0',col_chrom,str((col*args.resolution)+1),str(col)])+'\n')
                count -= 1.
        else:
            #   short with score format: https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
            o.write('\t'.join(['0',row_chrom,str((row*args.resolution)+1),str(row),'0',col_chrom,str((col*args.resolution)+1),str(col),str(freq)])+'\n')
        total += freq
print >> sys.stderr,("done")
assert(np.isclose(total,np.sum(matrix,dtype=np.float64))),"Error - entries in the binned matrix not written to storage"
