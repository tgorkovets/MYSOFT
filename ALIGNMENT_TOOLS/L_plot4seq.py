# -*- coding: utf-8 -*-
"""
This library plots profiles on top of sequences.
Using combination of R and python.
It can generate visual representation of sequences by itself.

"""

from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import cPickle as pickle
from Bio import SeqIO

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

import numpy as np
import pandas as pd
from StringIO import StringIO

from hist_ss import get_hist_ss
from hist_ss import get_hist_ss_in_aln, get_hist_ss_in_aln_for_shade
from Bio.Align.AlignInfo import SummaryInfo
import L_aln_tools
from L_shade_aln import shade_aln2png


TEMP_DIR='/Users/alexeyshaytan/junk/temptemp'

def plot_mat4seq(filename='default',data=[],seq1=[],seq1lab='Chain 1',seq1offset=0,features1=[],seq2=[],seq2lab='Chain 2',seq2offset=0,features2=[],title=''):
	"""
	will plot a 2D matrix for intersection of two sequences
	data is a dataframe with three columns Resid1, Resid2, Value.
	Resids - 0 based numbering, or offset specified
	you have to have zeros for non interacting residue(!)
	"""
	tempdf=TEMP_DIR+'/temp.csv'
	temppng1=TEMP_DIR+'/tempprofseq1.png'
	temppng2=TEMP_DIR+'/tempprofseq2.png'

	lenseq1=len(seq1[0])
	lenseq2=len(seq2[0])
	# print profile
	#convert matrix to a dataframe
	# df=pd.DataFrame(np.array([profile,range(len(profile))]).T,columns=[axis,'Resid'])
	data.to_csv(tempdf,index=False)
	shade_aln2png(seq1,filename=temppng1,shading_modes=['charge_functional'], legend=False, features=features1,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
	shade_aln2png(seq2,filename=temppng2,shading_modes=['charge_functional'], legend=False, features=features2,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False,rotate=True)


	#let's write an R-script
	a=open(TEMP_DIR+'/mat4seq.r','w')

	a.write(r"""
	library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img1 <- readPNG("%s")
img2 <- readPNG("%s")


"""%(tempdf,temppng1,temppng2))

	a.write("""
	seqimg1 <- rasterGrob(img1, interpolate=TRUE,width=1)
	seqimg2 <- rasterGrob(img2, interpolate=TRUE,height=1)


theme_set(theme_bw(base_size=24)+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=df,aes(x=Resid1+1-%d,y=Resid2+1-%d,color=Value))+

xlab("%s")+
ylab("%s")+
xlim(%f,%f)+
ylim(%f,%f)+
ggtitle("%s")+
# geom_hline(yintercept = c_c$PROT2_resid, colour="green", linetype = "longdash",size=0.5)+
# geom_vline(xintercept = h2azimp, colour="green", linetype = "longdash",size=0.5)+
geom_point(size=5)+scale_colour_gradient(low="blue", high="red",name="Value")+

annotation_custom(seqimg1, ymin=%f, ymax=0, xmin=0.5,xmax=%f)+
annotation_custom(seqimg2, xmin=%f,xmax=0,ymin=0.5, ymax=%f)
"""%(seq1offset,seq2offset,seq1lab,seq2lab,\
	-lenseq1*0.05,lenseq1*1.01,-lenseq2*0.05,lenseq2*1.01,\
	title[-lenseq1:len(title)],\
	-lenseq2*0.05,lenseq1+0.5,\
	-lenseq1*0.05,lenseq2+0.5))

	a.write("""
	ggsave("%s",plot=a,height=%f,width=%f)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),12./60.*lenseq2+3,12./60.*lenseq1+3))

	a.close()
	os.system('/usr/bin/R --vanilla --slave < '+TEMP_DIR+'/mat4seq.r')



def plot_prof4seq(filename='default',profile=[],seqmsa=[],features=[],axis='X',title=''):
	"""
	will plot a profile for every position in a sequence or small msa
	profile is a list of values.
	"""
	tempdf=TEMP_DIR+'/temp.csv'
	temppng=TEMP_DIR+'/tempprofseq.png'
	# print profile
	#convert profile to a dataframe
	df=pd.DataFrame(np.array([profile,range(len(profile))]).T,columns=[axis,'Resid'])
	df.to_csv(tempdf,index=False)
	shade_aln2png(seqmsa,filename=temppng,shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)

	#let's write an R-script
	a=open(TEMP_DIR+'/prof4seq.r','w')

	a.write(r"""
	library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img <- readPNG("%s")

"""%(tempdf,temppng))

	a.write("""
	seqimg <- rasterGrob(img, interpolate=TRUE,width=1)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=df,aes(x=Resid+1,y=%s))+
geom_bar(stat='identity',position='identity')+#scale_y_continuous(limits=c(-5,6),breaks=seq(0,6),labels=seq(0,6,by=1))+
# scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_color_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
xlab('Resid')+
ylim(%f,%f)+
# geom_point(data=h3data_amore,aes(color=color),fill='red',size=2)+
#ylab("RMSD, Å")+
ggtitle("%s")+
annotation_custom(seqimg, ymin=%f, ymax=0, xmin=0.5,xmax=%f)
"""%(axis,-max(profile)*1.01,max(profile)*1.01,title[-len(profile):len(title)],-max(profile),len(profile)+0.5))

	a.write("""
	ggsave("%s",plot=a,height=4,width=%f)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),12./60.*len(profile)))

	a.close()
	os.system('/usr/bin/R --vanilla --slave < '+TEMP_DIR+'/prof4seq.r')


if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A')])
    features=get_hist_ss_in_aln_for_shade(msa,below=True)
    plot_prof4seq('default',map(np.abs,map(np.sin,range(len(msa[0])))),msa,features,axis='conservation')
    # print features
    # shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    


