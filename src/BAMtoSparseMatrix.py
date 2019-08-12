#! /usr/bin/env python2.7
#
#   BAMtoSparseMatrix.py
#   Produces sparse Hi-C matrix format from filtered BAM and enzyme digest file (TSV format).
#   Requires samtools (http://samtools.sourceforge.net/) and Python NumPy module.
#   Expected TSV format for digest file:
#       chromosome  start   end
#       *start and end expected to be one-based
#       **ensure chromosome names match sequence names in BAM and include "chr" prefix
#       ***file may include more than three columns, but only the first three will be considered
#   Steps:
#       1) convert BAM to SAM and split by chr_A and chr_B interactions - BASH code, see awk command in subprocess calls
#       2) parse restriction enzyme digest file and build chromosome arrays - parseDigestFile() and buildChromosomeArray()
#       3) iterate over chromosome SAM files and bin interactions by restriction fragments (RF) - iterateOverSAMFiles()
#   Author: Christopher JF Cameron
#

import os,sys,argparse,subprocess,glob,re,datetime,distutils.spawn
import numpy as np

#   if samtools is not present in your PATH environment variable, please adjust the following line to include the full filepath
samtools_exe=distutils.spawn.find_executable("samtools")
if(samtools_exe==None):
    print("Error - 'samtools' cannot be found in PATH. Please install/adjust PATH variable before conitnuing or set within Python script")
    sys.exit(-2)

def buildChromosomeArray(cut_sites):
    """Builds numpy int32 array for expected chromosome digest"""
    chrom_array=np.empty(cut_sites[-1],dtype=np.int32)
    #   label array elements by restriction fragment ID
    for i,start in enumerate(cut_sites[:-1]):
        chrom_array[start:cut_sites[i+1]]=i
    
    return chrom_array

def parseDigestFile(digest_filepath,obs_chroms):
    """Parse formatted enzyme digest file"""
    array_dict={}
    cut_sites,prev_chrom=None,None
    print >> sys.stderr,("Building chromosome digest arrays ... "),
    with open(digest_filepath,"r") as FILE:
        for line in FILE:
            cur_chrom,start,end=line.rstrip().split()[:3]
            start,end=int(start)-1,int(end)
            if(not cur_chrom==prev_chrom):
                if(not cut_sites==None and prev_chrom in obs_chroms):
                    #   build chromosome array
                    array_dict[prev_chrom]=buildChromosomeArray(cut_sites)
                if(cur_chrom in obs_chroms):
                    print >> sys.stderr, (chrom),
                cut_sites=[start,end]
            else:
                #   append to existing cut site list
                cut_sites.append(end)
            prev_chrom=cur_chrom
        if(prev_chrom in obs_chroms):
            #   add last chromosome encountered
            array_dict[prev_chrom]=buildChromosomeArray(cut_sites)
    print >> sys.stderr, ("Done")

    return array_dict

def iterateOverSAMFiles(array_dict,LOG_filepath,output_prefix):
    """Iterate over chromosome temporary files and bin by expected restriction fragments in sparse matrix format"""
    print >> sys.stderr, ("Binning temporary files by expected enzyme digest ..."),
    with open(LOG_filepath,"a") as LOG:
        for filepath in glob.glob(output_prefix+"*.tmp"):
            basename=os.path.basename(filepath)
            LOG.write("Binning '"+basename+"' by expected fends ... ")
            freq_dict={}
            count=0
            file_problem=False
            min_row,max_row,min_col,max_col=None,None,None,None
            with open(filepath,"rU") as f:
                for line in f:
                    ID,chrom_A,start_A,chrom_B,start_B=line.rstrip().split()
                    start_A,start_B=int(start_A)-1,int(start_B)-1
                    #   index chromosomes based on digest arrays
                    try: array_dict[chrom_A]
                    except KeyError:
                        LOG.write("Error - '"+chrom_A+"' not found in the expected enzyme digest (skipping file)\n")
                        file_problem=True
                        break
                    try: array_dict[chrom_B]
                    except KeyError:
                        LOG.write("Error - '"+chrom_B+"' not found in the expected enzyme digest (skipping file)\n")
                        file_problem=True
                        break
                    #   ensure consistent ordering of chromosomes, based on digest arrays
                    (chrom_A,start_A),(chrom_B,start_B)=sorted([[chrom_A,start_A],[chrom_B,start_B]],key=lambda x:str(x[0]))
                    #   bin reads by restriction fragment ends (or fends)
                    RF_A=array_dict[chrom_A][start_A]
                    RF_B=array_dict[chrom_B][start_B]
                    #ensure fragment ordering for cis-contact maps
                    if(chrom_A==chrom_B and RF_A>RF_B):
                        RF_A,RF_B=RF_B,RF_A
                    #   store frequency in Python dictionary
                    try:
                        freq_dict[RF_A][RF_B]+=1
                    except KeyError:
                        try:
                            freq_dict[RF_A][RF_B]=1
                        except KeyError:
                            freq_dict[RF_A]={RF_B:1}
                    #track max and min row and cols
                    min_col=RF_B if RF_B<min_col or min_col==None else min_col
                    max_col=RF_B if RF_B>max_col or max_col==None else max_col
                    min_row=RF_A if RF_A<min_row or min_row==None else min_row
                    max_row=RF_A if RF_A>max_row or max_row==None else max_row
            if(file_problem):
                os.system("rm "+filepath) #   remove SAM file and conserve disk space
                continue
            #   write binned data to disk
            with open(output_prefix+".".join(["",chrom_A+"_"+chrom_B,"RF","tsv"]),"w") as o:
                o.write(" ".join(["#",str(min_row),str(max_row),str(min_col),str(max_col)])+"\n")   #header, defines boundaries of matrix
                for row in sorted(freq_dict.keys(),key=int):
                    for col in sorted(freq_dict[row].keys(),key=int):
                        o.write("\t".join([str(freq_dict[row][col]),chrom_A,str(row),chrom_B,str(col)])+"\n")
            LOG.write("Done\n")
            os.system("rm "+filepath)
        LOG.write("HI-C BINNING COMPLETE\n")
    print >> sys.stderr, ("Done")
    
#   parse comamand line arguments
parser=argparse.ArgumentParser()
parser.add_argument("bam_filepath",help="input BAM file path")
parser.add_argument("digest_filepath",help="expected restriction fragment digest file (BED) path")
parser.add_argument("output_dir",help="file path to output directory")
args=parser.parse_args()
    
#   define output directory
basename=".".join(os.path.basename(args.bam_filepath).split(".")[:-1])
output_directory=os.path.join(args.output_dir,"")
if(not os.path.exists(output_directory)):
    os.makedirs(output_directory)
    
output = ["b;ahb;ahs"]
output_prefix=os.path.join(args.output_dir,basename)
date_stamp="_".join(re.split("\W+|:|\.",str(datetime.datetime.now()))[:-1])
LOG_filepath=output_directory+basename+".BAMtoSparseMatrix."+date_stamp+".log"
with open(LOG_filepath,"w") as LOG:
    LOG.write("#  Log file for BAMtoSparseMatrix v1.03 (author: Christopher JF Cameron)\n\n")
    print >> sys.stderr,("Converting '"+basename+"' to temporary chromosome interaction files ..."),
    #   remove old tmp files from output directory
    for filepath in glob.glob(''.join([output_prefix,"*.tmp"])):
        os.remove(filepath)
    #   convert BAM to SAM and split by chr_A and chr_B interactions
    convert_process=subprocess.Popen([samtools_exe,"view",args.bam_filepath],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    split_process=subprocess.Popen(["awk","-F","\t","{ID_A=$1; chr_A=$3; start_A=$4; getline; ID_B=$1; chr_B=$3; start_B=$4; if(ID_A!=ID_B){print \"Error - BAM read pairs are not sequential (Stopping)\"; exit} if(chr_A>chr_B){tmp_chr=chr_A; tmp_start=start_A; chr_A=chr_B; start_A=start_B; chr_B=tmp_chr; start_B=tmp_start} filename=\""+output_prefix+".\"chr_A\"_\"chr_B\".tmp\"; print ID_A,chr_A,start_A,chr_B,start_B > filename}"],stdin=convert_process.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output=split_process.communicate()[0].rstrip().split("\n") #   wait for samtools and awk to finish before proceeding with Python script
    if(len(glob.glob(''.join([output_prefix,"*.tmp"]))) > 0):  
        print >> sys.stderr, ("Done")
    else:   #   problem encountered while parsing BAM file
        print " Error - BAM file was not processed correctly. Please check that SAMtools and awk are installed correctly (and work). See examples of called SAMtools and awk commands below:"
        print '\n'.join([' '.join(["\nSAMtools command (should print the first ten lines of the BAM to stdout):\n\t",samtools_exe,"view",args.bam_filepath,'|',"head"]),' '.join(["\nSAMtools+awk command (should create sparse matrix files):\n\t",samtools_exe,"view",args.bam_filepath,'|',"awk","-F","\'\\t\'","\'{ID_A=$1; chr_A=$3; start_A=$4; getline; ID_B=$1; chr_B=$3; start_B=$4; if(ID_A!=ID_B){print \"Error - BAM read pairs are not sequential (Stopping)\"; exit} if(chr_A>chr_B){tmp_chr=chr_A; tmp_start=start_A; chr_A=chr_B; start_A=start_B; chr_B=tmp_chr; start_B=tmp_start} filename=\""+output_prefix+".\"chr_A\"_\"chr_B\".tmp\"; print ID_A,chr_A,start_A,chr_B,start_B > filename}\'"])])
        sys.exit()
    LOG.write("BAM TO SAM CONVERSION COMPLETE\n\n")

    #   get list of temporary chromosome files from
    chroms=set([])
    for filepath in glob.glob(os.path.join(args.output_dir,"*.tmp")):
        #   parse filename
        for chrom in set(filepath.split(".")[-2].split("_")):
            chroms.add(chrom)
    #   parse digest file
    array_dict=parseDigestFile(args.digest_filepath,chroms)  #   returned 'chroms' variable has ordered chromosomes with respect to 'chrom_arrays'
    LOG.write("CHROMOSOME ARRAYS CREATED CORRECTLY\n\n")

#   iterate over chromosome SAM files and bin interactions by restriction fragments (RF)
iterateOverSAMFiles(array_dict,LOG_filepath,output_prefix)
