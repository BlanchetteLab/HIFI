#! /usr/bin/env python2.7
#
#   plotHIFIoutput.py
#   Plots sparse matrix output file of HIFI
#   Requires Python's NumPy and Matplotlib modules
#   Steps:
#       1) build full representation of sparse matrix in memory - buildIFmatrix()
#       2) determine restriction fragment sizes to create true-size representation - getRFdata()
#       3) plot sparse matrix in true-size representation using Matplotlib - plotSparseMatrix()
#   Author: Christopher JF Cameron
#

import argparse,os,sys,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

cmap="jet"  #   colour map to be used by Matplotlib, see others here:https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html

def buildIFmatrix(input_filepath,vmax):
    """Build full Hi-C IF matrix"""
    f=gzip.open(input_filepath,"rb") if(input_filepath.endswith(".gz")) else open(input_filepath,"rU")
    #process header
    header=f.readline().rstrip()
    if(not header.startswith("#")): #first line in HIFI sparse matrix file is not a header
        print "Error - unspecified range for plotting and no header line present in HIFI output file"
        sys.exit(2)    
    first_row,last_row,first_col,last_col=[int(val) for val in header[1:].split()]

    #create empty matrix
    m=(last_col-first_col)+1
    n=(last_row-first_row)+1
    matrix=np.zeros(shape=(n,m), dtype=float)

    #   fill in matrix entries
    total=0.0
    for line in f:
        if(not line.startswith(("#","0"))):   #incase header is present and user defines ranges
            freq,chrom_A,row_RF,chrom_B,col_RF=line.strip().split()
            freq=float(freq)
            #   convert HIFI restriction fragment (RF) ids to indices
            row_RF=int(row_RF)-first_row
            col_RF=int(col_RF)-first_col
            if(0<=row_RF<=(last_row-first_row) and 0<=col_RF<=(last_col-first_col)):  #contact is within desired matrix
                matrix[row_RF][col_RF]=freq
                total+=freq
                try:
                    matrix[col_RF][row_RF]=freq   #   create symmetrical Hi-C matrix
                except ValueError: continue #   non-square matrix
                total+=freq
            elif(row_RF>(last_row-first_row)): #outside matrix boundaries, skip remaining file
                break
    f.close()
    
    #   read count normalize IF matrix
    for i,row in enumerate(matrix):
        for j,val in enumerate(row):
            freq=np.log10(((val/total)*1000000.0)+1.0)
            matrix[i][j]=freq
    
    #   set diagonal to largest value in matrix
    for i,j in zip(xrange(0,n,1),xrange(0,m,1)):
        matrix[i][j]=vmax

    return matrix,chrom_A,[first_row,last_row,first_col,last_col]

def getRFdata(digest_filepath,target_chrom,(min_y,max_y,min_x,max_x)):
    "Parse expected restriction fragment (RF) BED file and determine RF sizes"
    #   determine fragment sizes
    i=-1
    RF_dict={"x":{"sizes":[],"start":None,"end":None},"y":{"sizes":[],"start":None,"end":None}}
    with open(digest_filepath,"rU") as f:
        for line in f:
            cur_chrom,start,end=line.rstrip().split()
            if(cur_chrom==target_chrom):
                i+=1
                start,end=int(start)-1,int(end)-1
                if(min_y<=i<=max_y):
                    RF_dict["y"]["start"]=start if RF_dict["y"]["start"]==None else RF_dict["y"]["start"]
                    RF_dict["y"]["sizes"].append(end-start)
                if(min_x<=i<=max_x):
                    RF_dict["x"]["start"]=start if RF_dict["x"]["start"]==None else RF_dict["x"]["start"]
                    RF_dict["x"]["sizes"].append(end-start)
                if(i>max_y):
                    RF_dict["y"]["end"]=start-1 if RF_dict["y"]["end"]==None else RF_dict["y"]["end"]
                if(i>max_x):
                    RF_dict["x"]["end"]=start-1 if RF_dict["x"]["end"]==None else RF_dict["x"]["end"]
                if(i>max_y and i>max_x):    #   remaing fragments for chromosome are not within observed matrix
                    break
            elif(not RF_dict["y"]["start"]==None):  #   chromosome of interest already observed
                break
    
    for axis in ["x","y"]:
        #   normalize RF lengths
        total=float(sum(RF_dict[axis]["sizes"]))
        RF_dict[axis]["sizes"]=[val/total for val in RF_dict[axis]["sizes"]]
        #   define axis labels
        RF_dict[axis]["label"]="".join([target_chrom.replace("chr","Chromosome "),": ",str(RF_dict[axis]["start"])," - ",str(RF_dict[axis]["end"])])
        
    return RF_dict

def plotSparseMatrix(matrix,chrom,vmin,vmax,RF_dict,output_filepath):
    "Plots HIFI sparse matrix using Python's Matplotlib module"
    plt.ioff()  #prevent view window from openning over ssh connection
    fig,ax=plt.subplots(figsize=(8,8))
    #create heat map using 'PatchCollection' rectangles and associated coordinates
    x=0.0
    patches,colours=[],[]
    y=sum(RF_dict["y"]["sizes"])
    for i,row_height in enumerate(RF_dict["y"]["sizes"]):
        y-=row_height
        for j,col_width in enumerate(RF_dict["x"]["sizes"]):
            colours.append(matrix[i][j])
            patch=Rectangle((x,y),col_width,row_height)
            patches.append(patch)
            x+=col_width
        x=0.0
    #plot 'PatchCollection' of rectangles
    pc=PatchCollection(patches,cmap=cmap)
    pc.set_clim([vmin,vmax])
    pc.set_array(np.array(colours))
    pc.set_edgecolor("face")
    pc.set_lw(0.1)
    ax.add_collection(pc)
    #remove outline from color bar
    cbar=plt.colorbar(pc,fraction=0.04575,pad=0.04)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=16,length=4,width=0.5,color="grey")
    ax.tick_params(axis="x",which="both",bottom=False,top=False,direction="out",length=4,width=0.5,pad=10,color="grey",labelsize=16)
    ax.xaxis.set_label_position("top") 
    ax.set_xlabel(RF_dict["x"]["label"],fontsize=16,labelpad=8.25)
    ax.xaxis.tick_top()
    ax.tick_params(axis="y",which="both",left=False,right=False,direction="out",length=4,width=0.5,color="grey",labelsize=16)
    ax.set_ylabel(RF_dict["y"]["label"],fontsize=16,labelpad=None)
    ax.grid(False)
    #   remove x and y ticks  
    plt.xticks([],[])
    plt.yticks([],[])
    #   remove axis edges
    for axis in ["bottom","top","left","right"]:    
        plt.gca().spines[axis].set_visible(False)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(output_filepath+".png",format="png",dpi=600,transparent=False,bbox_inches="tight")
    plt.savefig(output_filepath+".pdf",format="pdf",transparent=False,bbox_inches="tight")
    plt.close()
    
    #save matrix
    np.savetxt(output_filepath+".csv",matrix,delimiter=",")

#   parse comamand line arguments
parser=argparse.ArgumentParser()
parser.add_argument("HIFI_output",help="file path to HIFI sparse matrix (recommend filtering before plotting - see 'parseHIFIoutput.py')")
parser.add_argument("digest_filepath",help="expected restriction fragment digest (BED) filepath")
parser.add_argument("vmin",help="minimum value to be used for colour palette range")
parser.add_argument("vmax",help="maximum value to be used for colour palette range")
parser.add_argument("output_dir",help="file path to output directory")
args=parser.parse_args()

#   ensure output directory exists
output_directory=os.path.join(args.output_dir,"")
if(not os.path.exists(output_directory)):
    os.makedirs(output_directory)

print >> sys.stderr,("Building HIFI sparse matrix ... "),
matrix,chrom,matrix_boundaries=buildIFmatrix(args.HIFI_output,args.vmax)
print >> sys.stderr, ("Done")

print >> sys.stderr,("Importing restriction fragment information ... "),
RF_dict=getRFdata(args.digest_filepath,chrom,matrix_boundaries)
print >> sys.stderr, ("Done")

print >> sys.stderr,("Plotting full matrix ... "),
output_filename=".".join(os.path.basename(args.HIFI_output).split(".")[:-1])
plotSparseMatrix(matrix,chrom,args.vmin,args.vmax,RF_dict,output_directory+output_filename)
print >> sys.stderr, ("Done")