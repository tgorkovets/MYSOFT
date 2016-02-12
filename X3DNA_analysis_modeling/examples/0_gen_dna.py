#This script homology models yeast cenp-a nucleosome....

#Try1
#based on 1kx5 nucleosome structure
#we take 
#yeast HTA1 gene (gi 398366187) - as H2A,
# HTB1 (gi 398366183) as H2B,
#HHT1 (gi 6319482) - this should be H3
# but we will take CSE4 gi 27808712
# HHF1 (gi 6319481) - as H4
# we put the 601 DNA sequence
#we truncate the tails according to our NCP model simulations (see Shaytan et al. NAR(?) 2015)
#

#DNA is following
#aaGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAG T TTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAaaga

#Alexey Shaytan 2015

import os
import sys
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS')
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/X3DNA_analysis_modeling')


sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/structure_analysis')
sys.path.append('/Library/modeller-9.14/modlib')
sys.path.append('/Library/modeller-9.14/lib/mac10v4')
# sys.path.append(';'.join(['', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/setuptools-0.9.6-py2.7.egg', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages/PROPKA-3.1-py2.7.egg', '/Library/modeller-9.14/examples/atom_files', '/Applications/Bioinf/biana', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python27.zip', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-darwin', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/plat-mac/lib-scriptpackages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-tk', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-old', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/readline', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/lib-dynload', '/Users/alexeyshaytan/Library/Python/2.7/lib/python/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PIL', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyObjC', '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/gtk-2.0']))
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import L_fasta2pir
from L_aln_tools import muscle_aln

from dna_tools_simple import change_dna_seq_in_pdb


#let's generate a new DNA sequence
#aaGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAG
#T
#TTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAaaga
#len 162 - need to add to 147, i.e. 12 to each side
cen3_seq='AAGTCACATGATGATATTTGATTTTATTATATTTTTAAAAAAAGTAAAAAATAAAAAGTAGTTTATTTTTAAAAAATAAAATTTAAAATATTAGTGTATTTGATTTCCGAAAGTTAAAAAAGA'
# cen3_seq_flanked=['A']*12+cen3_seq+['A']*12
#we will add the 1kx5 sequence to the sides.
cen3_seq_flanked='ATCAATATCCAC'+cen3_seq+'GTGGATATTGAT'
cen3_seq_flanked=map(lambda x: x, cen3_seq_flanked)


change_dna_seq_in_pdb('1kx5.pdb','1kx5_cen3.pdb',cen3_seq_flanked)



