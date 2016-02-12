#Example of how to change DNA sequence in nucleosome
# via simple library


import pandas as pd
import numpy as np
from datetime import datetime
import sys
import os

import pandas as pd

#sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/X3DNA_analysis_modeling')
#sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')


from dna_tools_simple import change_dna_seq_in_pdb,X3DNA_find_pair,X3DNA_analyze_bp_step,build_dna

data_frame=X3DNA_analyze_bp_step('1kx5.pdb',X3DNA_find_pair('1kx5.pdb'))

#BPname       Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist  BPnum

s=pd.DataFrame([['A-T',0,0,0,0,0,0,0,0,3.4,0,0,36,1]],columns=['BPname','Shear','Stretch','Stagger','Buckle','Prop-Tw','Opening','Shift','Slide','Rise','Tilt','Roll','Twist','BPnum'])

edf=data_frame
for i in range(30):
	edf=edf.append(s,ignore_index=True)


print edf

seq='TAGGGACGAGATGGTACTTTGTGTCTCCTGCTCTGTCAGCAGGGCACTGTACTTGCTGATACCAGGGAATGTTTGTTCTTAAATACCATCATTCCGGACGTGTTTGCCTTGGCCAGTTTTCCATGTACATGCAGAAAGAAGTTTGGACTGATCAATACAGTCCTCTGCCTTTAAAGC'
seq=map(lambda x: x,seq)
print seq

build_dna(edf,'1kx5_new.pdb',seq)

# change_dna_seq_in_pdb('1kx5.pdb','1kx5_new.pdb',['A']*147)
# change_dna_seq_in_pdb('1kx5.pdb','1kx5_new.pdb',None)


# exit
exit()


