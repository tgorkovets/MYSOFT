#!/usr/bin/env python2.7
"""
This is a module provides analysis of simulated SAXS
within VMD through external calls to CRYSOL
Both programs should be installed
and configured beforehand

"""
import os
import subprocess
from VMD import *
from Molecule import *
from atomsel import *

import pandas as pd
import re

import uuid

# from scipy.spatial import KDTree
# from scipy.spatial import cKDTree
import numpy as np
from collections import OrderedDict

__author__="Alexey Shaytan"

TEMP='/Users/alexeyshaytan/junk/tmp2/'
P_ATSAS_DIR='/Applications/ATSAS'
P_ATSAS_CRYSOL=P_ATSAS_DIR+'/bin/crysol'

def crysol(atomsel):
	""" Runs crysol program from ATSAS and returns a dataframe

	See http://www.embl-hamburg.de/biosaxs/manuals/crysol.html
	for columns description

	"""

	#At first we need to makup a couple of unique file names
	unique=str(uuid.uuid4())
	pdb = unique+'.pdb'
	outf = unique+'00.int'

	print("Writing coords to "+pdb)
	atomsel.write('pdb',TEMP+'/'+pdb)
	cmd=P_ATSAS_CRYSOL+' '+pdb
	p = subprocess.Popen(cmd,shell=True,cwd=TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# wait for the process to terminate
	out, err = p.communicate()
	errcode = p.returncode
	print('OUT:'+out)
	print('ERR:'+err)
	#let's read int file.
	df=pd.read_csv(TEMP+'/'+outf, sep=u'  ',skiprows=1,header=None,names=['Vector','Int_in_sol','Int_in_vac','Solv_int','Border_int'])
	return(df)






if __name__ == '__main__':
	print "Kuku"
	X3DNA_find_pair('x')
	help(atomsel)

