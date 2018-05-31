#! /usr/bin/env python2.7
#
#   parseHIFIoutput.py
#   Parses sparse matrix output file of HIFI and filters for region of interest
#   Requires standard Python modules.
#   Steps:
#       1) determine first and last restriction fragment (RF) to be observed in filtered sparse matrix file - determineRangeFragments()
#       2) parse HIFI output sparse matrix file and filter for desired locus - parseHIFIoutput()
#   Author: Christopher JF Cameron
#

import sys,os,argparse,re

def determineRangeFragments(fend_filepath,target_chrom,bp_start,bp_end):
    """Parse restriction fragment BED file for fragments within desired range"""
    min_RF,max_RF=None,None
    print >> sys.stderr,("Parsing restriction fragment file ... "),
    with open(fend_filepath,"rU") as f:
        i=-1  #   fragment id
        for line in f:
            chrom,start,end=line.split()
            if(chrom==target_chrom):
                i+=1
                #determine restriction fragments (RF) within desired range
                if(bp_start<=int(end)-1<=bp_end):
                    min_RF=i if min_RF==None or i<min_RF else min_RF
                    max_RF=i if max_RF==None or i>max_RF else max_RF
    print >> sys.stderr,("Done")
    
    return min_RF,max_RF

def parseHIFIoutput(input_filepath,output_filepath,RF_start,RF_end,target_chrom):
    """Parses HIFI sparse matrix file for desired restriction fragment (RF) pairwise interactions within range""" 
    print >> sys.stderr,("Parsing HIFI sparse matrix file ... "),
    with open(output_filepath,"w") as o:
        o.write("#"+" ".join([str(val) for val in [RF_start,RF_end,RF_start,RF_end]])+"\n")
        f=gzip.open(input_filepath,"rb") if input_filepath.endswith(".gz") else open(input_filepath,"rU")
        for line in f:
            if(not line.startswith("#")):
                try: freq,row_chrom,row_RF,col_chrom,col_RF=line.rstrip().split()
                except:
                    print "Error - sparse matrix format not recognized"
                    sys.exit(-2)
                row_RF,col_RF=int(row_RF),int(col_RF)
                row_chrom="chr"+row_chrom if not row_chrom.startswith("chr") else row_chrom
                col_chrom="chr"+col_chrom if not col_chrom.startswith("chr") else col_chrom
                if(re.match("chr([0-9]+|[XY])",row_chrom) and re.match("chr([0-9]+|[XY])",col_chrom)):
                    assert(row_chrom==col_chrom),"Error - script is only meant to be applied to cis-contact maps"
                    assert(row_RF<=col_RF),"Error - lower triangle pair-wise interactions have been found in HIFI output file and should not be present"
                    #   check if RF is located within region of interest
                    if((RF_start<=row_RF<=RF_end)and(RF_start<=col_RF<=RF_end)and(not row_RF==col_RF)):
                        o.write("\t".join([freq,row_chrom,str(row_RF),col_chrom,str(col_RF)])+"\n")
                    if(row_RF>RF_end): 
                        #   assuming HIFI sparse matrix is ordered by RF integer id to save on processing time
                        #   stop iterating over lines in file when current row is downstream of region of interest
                        break
        f.close()
    print >> sys.stderr,("Done")
    
#   parse comamand line arguments
parser=argparse.ArgumentParser()
parser.add_argument("HIFI_output",help="file path to HIFI sparse matrix")
parser.add_argument("digest_filepath",help="expected restriction fragment digest (BED) filepath")
parser.add_argument("bp_start",help="starting base pair for region of interest")
parser.add_argument("bp_end",help="final base pair for region of interest")
parser.add_argument("output_dir",help="file path to output directory")
args=parser.parse_args()

#   ensure output directory exists
output_directory=os.path.join(args.output_dir,"")
if(not os.path.exists(output_directory)):
    os.makedirs(output_directory)

#   determine chromosome described within HIFI output
with open(args.HIFI_output,"r") as f:
    header=f.readline()
    chrom=f.readline().rstrip().split()[1]
#   determine first and last restriction fragment (RF) to be observed in filtered sparse matrix file
RF_start,RF_end=determineRangeFragments(args.digest_filepath,chrom,int(args.bp_start),int(args.bp_end))

#   parse HIFI output sparse matrix file and filter for desired locus
output_filepath=os.path.join(output_directory,".".join(os.path.basename(args.HIFI_output).split(".")[:-1]+[args.bp_start+"_"+args.bp_end,"tsv"]))
parseHIFIoutput(args.HIFI_output,output_filepath,RF_start,RF_end,chrom)