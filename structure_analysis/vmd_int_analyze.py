#!/usr/bin/env python2.7
"""
This is a module for VMD provides basic analysis routines
for interactions in molecular systems, it is dependent 
on int_analyze library.
"""


from VMD import *
from Molecule import *
from atomsel import *
import pandas as pd
import numpy as np
from datetime import datetime
import itertools
import sys
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/structure_analysis')

from int_analyze import find_contacts
from int_analyze import find_hbonds
from int_analyze import find_watermed_bonds
from int_analyze import find_watermed_bonds2
from int_analyze import find_ionmed_bonds

def get_contacts(SEL1,SEL2,d_threshold=3.9,exclude_bonded=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1'],code='C',chainfield='segname',half_matrix=True):
	"""Getting contacts"""

	SEL1_ind=SEL1.get('index')
	if(len(SEL1_ind)==0):
		return None
	SEL1_name=SEL1.get('name')
	SEL1_resid=SEL1.get('resid')
	SEL1_chain=SEL1.get(chainfield)
	SEL1_xyz=zip(SEL1.get('x'),SEL1.get('y'),SEL1.get('z'))
	SEL1_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL1_ind,SEL1_chain,SEL1_resid,SEL1_name))
	#DNA_xyz=zip(DNA.get('x'),DNA.get('y'),DNA.get('z'))
	
	# print SEL1_ind.index(3559)

	SEL2_ind=SEL2.get('index')
	if(len(SEL2_ind)==0):
		return None
	SEL2_name=SEL2.get('name')
	SEL2_resid=SEL2.get('resid')
	SEL2_chain=SEL2.get(chainfield)
	SEL2_xyz=zip(SEL2.get('x'),SEL2.get('y'),SEL2.get('z'))
	SEL2_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL2_ind,SEL2_chain,SEL2_resid,SEL2_name))
	
	
	## Now we can call a contact kernel
	# Input: coordinates and indices of group A, coordinates ans indices of group B, threshold
	# if(exclude_bonded):
		# contacts=find_contacts(SEL1_xyz,SEL1_ind,SEL1_xyz,SEL1_ind,d_threshold,exclude_bonded)
	# else:
	contacts=find_contacts(SEL1_xyz,SEL1_ind,SEL2_xyz,SEL2_ind,d_threshold,exclude_bonded,half_matrix=half_matrix)

	# print contacts
	num_cont=len(contacts)
	print "SEL1-SEL2 contacts found: ", num_cont
	# print(contacts)
	#let's form numpy arrays and convert them to data frame
	df_cont=pd.DataFrame()
	if(num_cont!=0):
		cont_np=np.array(contacts)
		
		#get SEL1 chain list
		SEL1_chain_list=np.array([SEL1_params[int(x)][0] for x in cont_np[:,0] ],ndmin=2)
		SEL1_resid_list=np.array([SEL1_params[int(x)][1] for x in cont_np[:,0] ],ndmin=2)
		SEL1_atom_list=np.array([SEL1_params[int(x)][2] for x in cont_np[:,0] ],ndmin=2)
		# print SEL1_chain_list
		# print chains
		
		SEL2_chain_list=np.array([SEL2_params[int(x)][0] for x in cont_np[:,1] ],ndmin=2)
		SEL2_resid_list=np.array([SEL2_params[int(x)][1] for x in cont_np[:,1] ],ndmin=2)
		SEL2_atom_list=np.array([SEL2_params[int(x)][2] for x in cont_np[:,1] ],ndmin=2)
		
		#let's make a numpy array for data frame
		np.set_printoptions(precision=5)
		
		df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array([code]*num_cont,ndmin=2).T,np.expand_dims(cont_np[:,2],axis=1).astype('|S7')))

		# print df_array
		df_cont=pd.DataFrame(df_array,columns=columns)

	return df_cont

def get_hbonds(SEL1_heavy,SEL1_H,SEL2_heavy,SEL2_H,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1','param2','param3'],code='HB',chainfield='segname'):
	"""Get H-bonds"""

	SEL1_heavy_ind=SEL1_heavy.get('index')
	if(len(SEL1_heavy_ind)==0):
		return None
	SEL1_heavy_name=SEL1_heavy.get('name')
	SEL1_heavy_resid=SEL1_heavy.get('resid')
	SEL1_heavy_chain=SEL1_heavy.get(chainfield)
	SEL1_heavy_xyz=zip(SEL1_heavy.get('x'),SEL1_heavy.get('y'),SEL1_heavy.get('z'))
	SEL1_H_xyz=zip(SEL1_H.get('x'),SEL1_H.get('y'),SEL1_H.get('z'))
	SEL1_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL1_heavy_ind,SEL1_heavy_chain,SEL1_heavy_resid,SEL1_heavy_name))
	
	# print SEL1_heavy_ind[1711]
	
	SEL2_heavy_ind=SEL2_heavy.get('index')
	if(len(SEL2_heavy_ind)==0):
		return None
	SEL2_heavy_name=SEL2_heavy.get('name')
	SEL2_heavy_resid=SEL2_heavy.get('resid')
	SEL2_heavy_chain=SEL2_heavy.get(chainfield)
	SEL2_heavy_xyz=zip(SEL2_heavy.get('x'),SEL2_heavy.get('y'),SEL2_heavy.get('z'))
	SEL2_H_xyz=zip(SEL2_H.get('x'),SEL2_H.get('y'),SEL2_H.get('z'))
	SEL2_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL2_heavy_ind,SEL2_heavy_chain,SEL2_heavy_resid,SEL2_heavy_name))
	
	## Now we can call a hbonds kernel
	# Input: list of donors coordinates and indices, list of coordinates of all H, list of aceptor coordinates and indices
	
	#NOTE: for crystal we set D-H-A - any, since Hydrogen positions are arbitrary!
	hbonds=find_hbonds(SEL1_heavy_xyz,SEL1_heavy_ind,SEL1_H_xyz,SEL2_heavy_xyz,SEL2_heavy_ind,SEL2_H_xyz,D_A_thresh,D_H_A_ang_thresh,same)
	
	num_hbonds=len(hbonds)
	print "SEL1-SEL2 H-bonds found: ", num_hbonds
	
	# print hbonds
	df_hbonds=pd.DataFrame()
	if(num_hbonds!=0):
		hbonds_np=np.array(hbonds)
		#get SEL1 chain list
		# print type(SEL1_heavy_params.keys()[1726])
		# print type(hbonds_np[0,0])
		SEL1_chain_list=np.array([SEL1_heavy_params[int(x)][0] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_resid_list=np.array([SEL1_heavy_params[int(x)][1] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_atom_list=np.array([SEL1_heavy_params[int(x)][2] for x in hbonds_np[:,0] ],ndmin=2)
		# print SEL1_chain_list
		# print chains
		
		SEL2_chain_list=np.array([SEL2_heavy_params[int(x)][0] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_resid_list=np.array([SEL2_heavy_params[int(x)][1] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_atom_list=np.array([SEL2_heavy_params[int(x)][2] for x in hbonds_np[:,1] ],ndmin=2)
		
		
		
		#let's make a numpy array for data frame
		df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array([code]*num_hbonds,ndmin=2).T,np.expand_dims(hbonds_np[:,3],axis=1).astype('|S7'),np.expand_dims(hbonds_np[:,4],axis=1).astype('|S7'),np.expand_dims(hbonds_np[:,2],axis=1).astype('|S7')))
		# print df_array
		df_hbonds=pd.DataFrame(df_array,columns=columns)
	return df_hbonds
	

def get_hbonds_imp(SEL1,SEL2,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1','param2','param3'],code='HB'):
	"""Get H-bonds implicitly, i.e. determine donors, acceptors  and Hydrogens automatically"""
	#nucleic and noh and (name \"N.*\" or name \"O.*\")
	#protein and noh and (name \"N.*\" or name \"O.*\" or name \"S.*\")
	SEL1_heavy=atomsel("(noh and (name \"N.*\" or name \"O.*\" or name \"S.*\")) and index %s"%" ".join(map(str,SEL1.get('index'))))
	SEL2_heavy=atomsel("(noh and (name \"N.*\" or name \"O.*\" or name \"S.*\")) and index %s"%" ".join(map(str,SEL2.get('index'))))

	
	SEL1_H=atomsel("hydrogen and within 1.4 of index %s" % " ".join(map(str,SEL1_heavy.get('index'))))
	SEL2_H=atomsel("hydrogen and within 1.4 of index %s" % " ".join(map(str,SEL2_heavy.get('index'))))
	
	df_hbonds=get_hbonds(SEL1_heavy,SEL1_H,SEL2_heavy,SEL2_H,D_A_thresh,D_H_A_ang_thresh,same,columns=columns,code=code)
	return df_hbonds
	
def get_salt_bridges_imp(SEL1,SEL2,d_threshold=3.9,exclude_bonded=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1'],code='SB'):
	"""Get salt-bridges implicitly, i.e. determine all SB forming atoms automatically"""
	#################################################
	## Now let's calculate salt bridges
	## less than 3.9A distance between chared atoms
	## 
	## Charged in  protein, atom names: ARG - NH1, NH2, LYS - NZ, ASP - OD1, OD2, GLU - OE1, OE2
	#but also terminal N and OT1 OT2, terminal N is unfortunately just N, it has type NH3 in CHARMM - together with LYS NH3.
	#O1P O2P is charged in DNA, check amber needed!!! as I remember O names are different there.
	#############
# acidic - basic
	df_ionpairs1=pd.DataFrame()
	df_ionpairs2=pd.DataFrame()
	#Still need to fix it for amber!!!
	#C-term is OXT and O - how to add O???
	#N-term also problematic - just N???
	SELt1=atomsel("((resname ASP GLU and (name OD1 OD2 OE1 OE2)) or (name OT1 OT2 OXT) or (name OP1 OP2 O1P O2P)) and index %s"%" ".join(map(str,SEL1.get('index'))),frame=SEL1.frame)
	SELt2=atomsel("(resname ARG LYS and (name NH1 NH2 NZ) or (backbone and type NH3)) and index %s"%" ".join(map(str,SEL2.get('index'))),frame=SEL2.frame)	
	if(len(SELt1)&len(SELt2)):
		df_ionpairs1=get_contacts(SELt1,SELt2,d_threshold,exclude_bonded,columns=columns,code=code)

	# Basic - acidic
	SELt2=atomsel("((resname ASP GLU and (name OD1 OD2 OE1 OE2)) or (name OT1 OT2 OXT) or (name OP1 OP2 O1P O2P)) and index %s"%" ".join(map(str,SEL2.get('index'))),frame=SEL2.frame)
	SELt1=atomsel("(resname ARG LYS and (name NH1 NH2 NZ) or (backbone and type NH3)) and index %s"%" ".join(map(str,SEL1.get('index'))),frame=SEL1.frame)
	if(len(SELt1)&len(SELt2)):
		df_ionpairs2=get_contacts(SELt1,SELt2,d_threshold,exclude_bonded,columns=columns,code=code)

	return pd.concat([df_ionpairs1,df_ionpairs2])
	
def get_watermed_bonds(SEL1_heavy,SEL2_heavy,WAT_heavy,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param3'],code='WM'):
	"""Find water mediated bonds between SEL1 and SEL2, through WAT, only heavy atoms should be supplied"""

	#The bonds to Hydrogen in charmm are all less than 1.2 A
	#Let's get the corresponding hydrogens
	SEL1_H=atomsel("hydrogen and within 1.4 of index %s" % " ".join(map(str,SEL1_heavy.get('index'))))
	SEL2_H=atomsel("hydrogen and within 1.4 of index %s" % " ".join(map(str,SEL2_heavy.get('index'))))
	WAT_H=atomsel("hydrogen and within 1.4 of index %s" % " ".join(map(str,WAT_heavy.get('index'))))

	SEL1_heavy_ind=SEL1_heavy.get('index')
	SEL1_heavy_name=SEL1_heavy.get('name')
	SEL1_heavy_resid=SEL1_heavy.get('resid')
	SEL1_heavy_chain=SEL1_heavy.get('segname')
	SEL1_heavy_xyz=zip(SEL1_heavy.get('x'),SEL1_heavy.get('y'),SEL1_heavy.get('z'))
	SEL1_H_xyz=zip(SEL1_H.get('x'),SEL1_H.get('y'),SEL1_H.get('z'))
	SEL1_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL1_heavy_ind,SEL1_heavy_chain,SEL1_heavy_resid,SEL1_heavy_name))
	
	# print SEL1_heavy_ind[1711]
	
	SEL2_heavy_ind=SEL2_heavy.get('index')
	SEL2_heavy_name=SEL2_heavy.get('name')
	SEL2_heavy_resid=SEL2_heavy.get('resid')
	SEL2_heavy_chain=SEL2_heavy.get('segname')
	SEL2_heavy_xyz=zip(SEL2_heavy.get('x'),SEL2_heavy.get('y'),SEL2_heavy.get('z'))
	SEL2_H_xyz=zip(SEL2_H.get('x'),SEL2_H.get('y'),SEL2_H.get('z'))
	SEL2_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL2_heavy_ind,SEL2_heavy_chain,SEL2_heavy_resid,SEL2_heavy_name))
	
	WAT_heavy_ind=WAT_heavy.get('index')
	WAT_heavy_name=WAT_heavy.get('name')
	WAT_heavy_resid=WAT_heavy.get('resid')
	WAT_heavy_chain=WAT_heavy.get('segname')
	WAT_heavy_xyz=zip(WAT_heavy.get('x'),WAT_heavy.get('y'),WAT_heavy.get('z'))
	WAT_H_xyz=zip(WAT_H.get('x'),WAT_H.get('y'),WAT_H.get('z'))
	WAT_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(WAT_heavy_ind,WAT_heavy_chain,WAT_heavy_resid,WAT_heavy_name))
	
	
	#NOTE: for crystal we set D-H-A - any, since Hydrogen positions are arbitrary!
	hbonds=find_watermed_bonds(SEL1_heavy_xyz,SEL1_heavy_ind,SEL1_H_xyz,SEL2_heavy_xyz,SEL2_heavy_ind,SEL2_H_xyz,WAT_heavy_xyz,WAT_heavy_ind,WAT_H_xyz,D_A_thresh,D_H_A_ang_thresh,same)
	
	num_hbonds=len(hbonds)
	print "SEL1-SEL2 water mediated bonds found: ", num_hbonds
	
	# print hbonds
	df_watermed_bonds=pd.DataFrame()
	
	if(num_hbonds!=0):
		hbonds_np=np.array(hbonds)
		#get SEL1 chain list
		# print type(SEL1_heavy_params.keys()[1726])
		# print type(hbonds_np[0,0])
		SEL1_chain_list=np.array([SEL1_heavy_params[int(x)][0] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_resid_list=np.array([SEL1_heavy_params[int(x)][1] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_atom_list=np.array([SEL1_heavy_params[int(x)][2] for x in hbonds_np[:,0] ],ndmin=2)
		# print SEL1_chain_list
		# print chains
		
		SEL2_chain_list=np.array([SEL2_heavy_params[int(x)][0] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_resid_list=np.array([SEL2_heavy_params[int(x)][1] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_atom_list=np.array([SEL2_heavy_params[int(x)][2] for x in hbonds_np[:,1] ],ndmin=2)
		
		
		
		#let's make a numpy array for data frame
		df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array([code]*num_hbonds,ndmin=2).T,np.expand_dims(hbonds_np[:,2],axis=1).astype(np.str)))
		# df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array(['WM']*num_hbonds,ndmin=2).T,np.expand_dims(hbonds_np[:,2],axis=1)))
		# print df_array
		df_watermed_bonds=pd.DataFrame(df_array,columns=columns)
	return df_watermed_bonds
	

def get_watermed_bonds_imp(SEL1,SEL2,D_A_thresh=3.5,D_H_A_ang_thresh=30,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param3'],code='WM'):
	"""Get watermed bonds implicitly determining all atom groups, water should be present around"""
	SEL1_heavy=atomsel("(noh and (name \"N.*\" or name \"O.*\" or name \"S.*\")) and index %s"%" ".join(map(str,SEL1.get('index'))))
	SEL2_heavy=atomsel("(noh and (name \"N.*\" or name \"O.*\" or name \"S.*\")) and index %s"%" ".join(map(str,SEL2.get('index'))))
	
	#This should be checked for AMBER, amber name is O resname WAT
	WAT_heavy=atomsel("((name OH2) or ((resname WAT) and (name O))) and within 3.9 of index %s" % " ".join(map(str,SEL1_heavy.get('index'))))
	return get_watermed_bonds(SEL1_heavy,SEL2_heavy,WAT_heavy,D_A_thresh,D_H_A_ang_thresh,same,columns=columns,code=code)
	

def get_ionmed_bonds(SEL1_heavy,SEL2_heavy,ION_heavy,d_threshold=3.9,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1','param3'],code='IM'):
	"""Get ion mediated bonds between SEL1 and SEL2 through ION, input heavy atoms only"""

	SEL1_heavy_ind=SEL1_heavy.get('index')
	SEL1_heavy_name=SEL1_heavy.get('name')
	SEL1_heavy_resid=SEL1_heavy.get('resid')
	SEL1_heavy_chain=SEL1_heavy.get('segname')
	SEL1_heavy_xyz=zip(SEL1_heavy.get('x'),SEL1_heavy.get('y'),SEL1_heavy.get('z'))
	SEL1_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL1_heavy_ind,SEL1_heavy_chain,SEL1_heavy_resid,SEL1_heavy_name))
	
	
	SEL2_heavy_ind=SEL2_heavy.get('index')
	SEL2_heavy_name=SEL2_heavy.get('name')
	SEL2_heavy_resid=SEL2_heavy.get('resid')
	SEL2_heavy_chain=SEL2_heavy.get('segname')
	SEL2_heavy_xyz=zip(SEL2_heavy.get('x'),SEL2_heavy.get('y'),SEL2_heavy.get('z'))
	SEL2_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(SEL2_heavy_ind,SEL2_heavy_chain,SEL2_heavy_resid,SEL2_heavy_name))
	
	ION_heavy_ind=ION_heavy.get('index')
	ION_heavy_name=ION_heavy.get('name')
	ION_heavy_resid=ION_heavy.get('resid')
	ION_heavy_chain=ION_heavy.get('segname')
	ION_heavy_xyz=zip(ION_heavy.get('x'),ION_heavy.get('y'),ION_heavy.get('z'))
	ION_heavy_params=dict((key, (chain,resid,name)) for (key,chain,resid,name) in zip(ION_heavy_ind,ION_heavy_chain,ION_heavy_resid,ION_heavy_name))
	
	# print ION_heavy_ind
	## Now we can call a hbonds kernel
	
	hbonds=find_ionmed_bonds(SEL1_heavy_xyz,SEL1_heavy_ind,SEL2_heavy_xyz,SEL2_heavy_ind,ION_heavy_xyz,ION_heavy_ind,d_threshold,same)
	
	num_hbonds=len(hbonds)
	print "SEL1-SEL2 ion mediated bonds found: ", num_hbonds
	
	df_ionmed_bonds=pd.DataFrame()
	
	if(num_hbonds!=0):
	# print hbonds
	
		hbonds_np=np.array(hbonds)
		# print hbonds_np
		#get SEL1 chain list
		# print type(SEL1_heavy_params.keys()[1726])
		# print type(hbonds_np[0,0])
		SEL1_chain_list=np.array([SEL1_heavy_params[int(x)][0] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_resid_list=np.array([SEL1_heavy_params[int(x)][1] for x in hbonds_np[:,0] ],ndmin=2)
		SEL1_atom_list=np.array([SEL1_heavy_params[int(x)][2] for x in hbonds_np[:,0] ],ndmin=2)
		# print SEL1_chain_list
		# print chains
		
		SEL2_chain_list=np.array([SEL2_heavy_params[int(x)][0] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_resid_list=np.array([SEL2_heavy_params[int(x)][1] for x in hbonds_np[:,1] ],ndmin=2)
		SEL2_atom_list=np.array([SEL2_heavy_params[int(x)][2] for x in hbonds_np[:,1] ],ndmin=2)
		
		ION_atom_list=np.array([ION_heavy_params[int(x)][2] for x in hbonds_np[:,2] ],ndmin=2)
		
		
		#let's make a numpy array for data frame
		df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array([code]*num_hbonds,ndmin=2).T,ION_atom_list.T,np.expand_dims(hbonds_np[:,2],axis=1).astype(np.str)))
	
		# df_array=np.hstack((SEL1_chain_list.T,SEL1_resid_list.T,SEL1_atom_list.T,SEL2_chain_list.T,SEL2_resid_list.T,SEL2_atom_list.T,np.array(['IM']*num_hbonds,ndmin=2).T,ION_atom_list.T,np.expand_dims(hbonds_np[:,2],axis=1)))
		# print df_array
		df_ionmed_bonds=pd.DataFrame(df_array,columns=columns)
	return df_ionmed_bonds

def get_ionmed_bonds_imp(SEL1,SEL2,d_threshold=3.9,same=False,columns=['SEL1_chain','SEL1_resid','SEL1_atom','SEL2_chain','SEL2_resid','SEL2_atom','type','param1','param3'],code='IM'):
	"""Get ion mediated bonds between SEL1 and SEL2, implicitly defining polar bondable groups through ions from environment"""

	SEL1_heavy=atomsel("(noh and (name \"N.*\" or name \"O.*\" or name \"S.*\" or name P)) and index %s"%" ".join(map(str,SEL1.get('index'))))
	SEL2_heavy=atomsel("( noh and (name \"N.*\" or name \"O.*\" or name \"S.*\" or name P)) and index %s"%" ".join(map(str,SEL2.get('index'))))
	
	#This needs to be rechecked for AMBER and other crystal ions, like RB, etc.
	ION_heavy=atomsel("name SOD CLA MG POT 'Na+' 'Cl-'")
	return get_ionmed_bonds(SEL1_heavy,SEL2_heavy,ION_heavy,d_threshold,same,columns=columns,code=code)



if __name__ == '__main__':
	
	#Understanding KDTree
	import numpy as np                                                                                                                                                                                             
	import scipy.spatial as ss                                                                                                                                                                                     
	from itertools import combinations
	from itertools import permutations
	import time
	d_threshold = 1.5
	x=[1,2,3,4,5,6,7,8,9,10]
	y=[1,2,3,4,5,6,7,8,9,10]
	points = zip(x, y)
	x2=[1,3,4,5,6,7,8,8,10,16]
	y2=[1,2,3,4,5,6,7,8,9,16.1]
	points2 = zip(x, y)
	start = time.clock()
	A_tree=KDTree(points)
	B_tree=KDTree(points2)
	neighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	# print neighbors #it is a list of lists
	distances=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# print distances
	A_tree=cKDTree(points)
	B_tree=cKDTree(points2)
	# print points[0]
	cneighbors=A_tree.query_ball_tree(B_tree,d_threshold)
	# print neighbors #it is a list of lists
	cdistances=A_tree.sparse_distance_matrix(B_tree,d_threshold)
	# print distances
	# print neighbors
	# print cneighbors
	#print distances
	# if(neighbors==cneighbors):
		# print "Yes!!!"

	# di=distances.items()
	# cdi=cdistances.items()
	# if(cdi==di):
		# print "Yes!!!"
	# for i in distances.iterkeys():
		# print i
	# print distances
	# print cdistances
	D=(-2,0,0)
	H=(0,0,0)
	A=(-1,0.1,0)
	print get_D_H_A(D,H,A)

	elapsed = (time.clock() - start)
	print "Time elapsed: %f s" % elapsed



