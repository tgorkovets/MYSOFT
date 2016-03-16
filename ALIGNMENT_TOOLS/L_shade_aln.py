# -*- coding: utf-8 -*-
"""
This is a library that makes good images of shaded alignments 
through TeXShade.

Input:
1) alignment (might be a single sequence).
2) Shading options
3) Features list
Output:
pdf or image (needed for automaitc plotting)


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

import networkx as nx
from StringIO import StringIO

from hist_ss import get_hist_ss
from hist_ss import get_hist_ss_in_aln, get_hist_ss_in_aln_for_shade
from Bio.Align.AlignInfo import SummaryInfo
import L_aln_tools


TEMP_DIR='/Users/alexeyshaytan/junk/temptemp'

def shade_aln2png(msa,filename='default',shading_modes=['similar'],features={},title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,rotate=False):
    intf=TEMP_DIR+'/tempshade.pdf'
    shade_aln2pdf(msa,intf,shading_modes,features,title,legend, logo,hideseqs,splitN,setends,ruler,show_seq_names,show_seq_length)
    #let's use imagemagic
    if rotate:
        os.system('convert -density 150 '+intf+' -trim -rotate -90 %s'%(filename if filename[-3:]=='png' else filename+'.png'))
    else:
        os.system('convert -density 150 '+intf+' -trim %s'%(filename if filename[-3:]=='png' else filename+'.png'))


def shade_aln2pdf(msa,filename='default',shading_modes=['similar'],features={},title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True):
    """
will convert msa to a shaded pdf.
shading_modes: similar, ... see write_texshade code
features - a list of dictionaries:
{'style':'block','helix','loop','-->','--','<--',',-,' - all texshade features types of features + frameblock
'position':'top','bottom','ttop','bbottom', etc. if no - automatic
'seqref':number of sequence for selection - default consensus
'sel':[begin,end] - region selection, this is in 0 based numbering (as in msa - we override here the TexShade 1-based numbering)
'text':'text'
'color':'Red' this is for frame block
'thickness':1.5 -1.5pt
}

    """
    ns='consensus'
    hideseqs_by_name=[]
    if(len(msa)==1):
        msa=msa[:]
        msa.extend([SeqRecord(msa[0].seq,id='dum',name='dum')])
        hideseqs_by_name.append('dum')
    #####if we are splitting the alignment into blocks - get number of blocks

    a_len=len(msa)
    num=int(a_len/splitN)+1
    while ((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    print "Chosen splitting parameters"    
    print a_len, splitN

    ####iterate over blocks and create alignment fasta files
    for i in range(num):
        t_aln=msa[(i*splitN):((i+1)*splitN)]
        AlignIO.write(t_aln,open(TEMP_DIR+'/alignment%d.fasta'%i,'w'), 'fasta')

    res_per_line=len(msa[0])

    #prepare feature section

    #alias dict
    aliasf={'alpha':'helix','beta':'-->','domain':'loop'}


    features_code=''
    for i in features:
        if(i['style']=='block'):
            features_code+="\\frameblock{%s}{%d..%d}{%s[%.1fpt]}"%(i.get('position','top'),i['sel'][0]+1,i['sel'][1]+1,i.get('color','Red'),i.get('thickness',1.5))
        else:
            features_code+="\\feature{%s}{%s}{%d..%d}{%s}{%s}"%(i.get('position','top'),str(i.get('seqref','consensus')),i['sel'][0]+1,i['sel'][1]+1,aliasf.get(i.get('style','loop'),i.get('style','loop')),i.get('text',''))
            
    a=open(TEMP_DIR+'/align.tex','w')

    a.write(r"""\documentclass[11pt,landscape]{article}
%\documentclass{standalone}
%\usepackage[a0paper]{geometry}
\usepackage{hyperref}
""")

    h=((a_len/30.*18 + (2.5 if legend else 0.0)) if (a_len/30.*18 + (2.5 if legend else 0.0) <18.0) else 18)
    w=(22/200.*res_per_line+2.5)

    if title:
        h+=1
        if w<len(title)*0.4:
            w=len(title)*0.4
    a.write("""
\\usepackage[paperwidth=%fin, paperheight=%fin]{geometry}
        """%(w,h))

    a.write(r"""
\usepackage{texshade}

\begin{document}""")
    a.write("""
\\pagenumbering{gobble}
\\centering
\\Huge{%s}
\\vspace{-0.5in}"""%title)

    for i in range(num-1):
        write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%i , features_code, res_per_line,False,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name)
    write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%(num-1) , features_code, res_per_line,legend,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name)

    a.write(r"""
\end{document} """)
    a.close()

    command='/usr/bin/env pdflatex --file-line-error --synctex=1 -output-directory=%s --save-size=10000  %s/align.tex > /dev/null'%(TEMP_DIR,TEMP_DIR)

    print('Launcning Latex:')
    print(command)
    os.system(command)
    os.system('mv '+TEMP_DIR+'/align.pdf %s'%(filename if filename[-3:]=='pdf' else (filename+'.pdf')))




def write_texshade(file_handle,aln_fname,features,res_per_line=120,showlegend=True,shading_modes=['similar'],logo=False,hideseqs=False,setends=[],ruler=False,numbering_seq='consensus',hide_ns=False,show_seq_names=True,show_seq_length=True,hideseqs_by_name=[]):

    for shading in shading_modes:
        shading=str(shading)
        file_handle.write("""
    \\begin{texshade}{%s}
    \\residuesperline*{%d}
    """%(aln_fname,res_per_line))
       

        file_handle.write(r"""
    \seqtype{P}
    \defconsensus{{}}{*}{upper}
    """)
        #a very dirty hack
        if(setends):
            if(numbering_seq=='consensus'):
                file_handle.write("""
    \\setends[%d]{%s}{%d..%d}
    """%(setends[0]+1,numbering_seq,setends[0]+setends[0],setends[1]+setends[0]+setends[0]-1))
            else:
                    file_handle.write("""
    \\setends{%s}{%d..%d}
    """%(numbering_seq,setends[0],setends[1]))
           
        if(ruler):
            file_handle.write("""
    \\showruler{top}{%s}
    """%numbering_seq)
        if(hide_ns):
            file_handle.write("""
    \\hideseq{%s}
    """%numbering_seq)

        file_handle.write(r"""
    \seqtype{P}
    """)
        if(not show_seq_names):
            file_handle.write(r"""
    \hidenames
    """) 
        if(not show_seq_length):
            file_handle.write(r"""
    \hidenumbering
    """) 
        if((shading=='similar')|(shading=='0')):
            file_handle.write(r"""
    \shadingmode{similar}
    \threshold[80]{50}
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo{top} \showlogoscale{leftright}
    \namesequencelogo{logo}
    """)

        if((shading=='hydropathy_functional')|(shading=='1')):
            file_handle.write(r"""
    \shadingmode[hydropathy]{functional}
    \shadeallresidues
    \threshold[80]{50}
    
    """)     
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[hydropathy]{top} \showlogoscale{leftright}

    """)
        if((shading=='chemical_functional')|(shading=='2')):
            file_handle.write(r"""

    \shadingmode[chemical]{functional}
    \shadeallresidues
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[chemical]{top} \showlogoscale{leftright}

    """)
        if((shading=='structure_functional')|(shading=='3')):
            file_handle.write(r"""
    \shadingmode[structure]{functional}
    \shadeallresidues
    """)

        if((shading=='charge_functional')|(shading=='4')):
            file_handle.write(r"""
    \shadingmode[charge]{functional}
    \shadeallresidues
    """)

        if((shading=='diverse')|(shading=='5')):
            file_handle.write(r"""
    \shadingmode{diverse}

    """)
       
        if(hideseqs):
            file_handle.write(r"""
    \hideseqs

    """)
        for s in hideseqs_by_name:
            file_handle.write("""
    \\hideseq{%s}
    """%s)

        file_handle.write(r"""

    %\setends{consensus}{1..160}
    %\setends{consensus}{1..160}


    %\feature{ttop}{1}{1..160}{bar:conservation}{}
    %\showfeaturestylename{ttop}{conserv.}
    \ttopspace{-\baselineskip}

    %\feature{top}{1}{1..160}{color:charge}{}
    %\showfeaturestylename{top}{charge}

    %\feature{bottom}{1}{1..160}{color:molweight[ColdHot]}{}
    %\showfeaturestylename{bottom}{weight}

    %\bbottomspace{-\baselineskip}
    %\feature{bbottom}{2}{1..160}{bar:hydrophobicity[Red,Gray10]}{}
    %\showfeaturestylename{bbottom}{hydrophob.}

    %\bargraphstretch{3}
    %\featurestylenamescolor{Red}
    %\featurestylenamesrm  \featurestylenamesit

    %\showsequencelogo{top}


    %\showconsensus[ColdHot]{bottom}
    \showconsensus[black]{top}

    %\defconsensus{.}{lower}{upper}
    %\defconsensus{{}}{lower}{upper}
    %\defconsensus{{}}{*}{upper}

    """)
        file_handle.write(features)
        file_handle.write(r"""
    \featurerule{1mm}""")
        if(showlegend):
            file_handle.write(r"""
    \showlegend""")

        align = AlignIO.read(open(aln_fname,'r'), "fasta")
        for a,i in zip(align,range(len(align))):
            # print a.id.replace('|',' | ')
            file_handle.write("""
    \\nameseq{%d}{%s}"""%(i+1,a.id.replace('|',' | ')))

        file_handle.write(r"""

    \end{texshade}
    """)

def feature_str2dict(featurestring,position='top'):
    """converts string of secondary structure annotation (like in VMD) to our type of dict"""
    #numbering should be 0 based
    #HHHHHHEEEEEBBBBBBCCCCCbTGI
    features=[]
    style=''
    begin=-1
    for i,s in enumerate(featurestring):
        if(s in ['H','G','I']): #helices
            if style!='helix':
                style='helix'
                begin=i
        else:
            if style=='helix':
                style=''
                end=i-1
                features.append({'style':'helix','sel':[begin,end],'position':position})

    for i,s in enumerate(featurestring):
        if(s in ['E','B','b']): #helices
            if style!='-->':
                style='-->'
                begin=i
        else:
            if style=='-->':
                style=''
                end=i-1
                features.append({'style':'-->','sel':[begin,end],'position':position})
    return features




    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A')])
    features=get_hist_ss_in_aln_for_shade(msa,below=True)
    # features=[{'style':'fill:$\uparrow$','sel':[5,10],'text':'test'}]
    print features
    shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    



    
            