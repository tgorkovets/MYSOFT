#Example of how to change DNA sequence in nucleosome
# via simple library


import pandas as pd
import numpy as np
from datetime import datetime
import sys
import os

import pandas as pd

sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/X3DNA_analysis_modeling')
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')


from dna_tools_simple import change_dna_seq_in_pdb



change_dna_seq_in_pdb('1kx5.pdb','1kx5_new.pdb',['A']*147)
# change_dna_seq_in_pdb('1kx5.pdb','1kx5_new.pdb',None)


# exit
exit()


