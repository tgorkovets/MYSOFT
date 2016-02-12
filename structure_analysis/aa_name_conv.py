#!/usr/bin/env python2.7
"""
Converts three letter names to one letter
"""

conv_table={
'ALA':'A',
'CYS':'C',
'ASP':'D',
'GLU':'E',
'PHE':'F',
'GLY':'G',
'HIS':'H',
'ILE':'I',
'LYS':'K',
'LEU':'L',
'MET':'M',
'ASN':'N',
'PRO':'P',
'GLN':'Q',
'ARG':'R',
'SER':'S',
'THR':'T',
'VAL':'V',
'TRP':'W',
'TYR':'Y',
'MSE':'X',
'UNK':'X',
'ADE':'A',
'CYT':'C',
'GUA':'G',
'THY':'T',
'A':'A',
'T':'T',
'G':'G',
'C':'C',
'DA':'A',
'DT':'T',
'DG':'G',
'DC':'C',
'U':'U'

}


def three2one(list):
	"""This function converts three letter codes to one
	"""
	return [conv_table[i] if i in conv_table else 'X' for i in list]


