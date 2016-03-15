#!/usr/bin/env python

"""
Program for plotting dat files with legends and names

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv


def myplot(file):
	"Plot dat file"

	
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['ps.useafm'] = True
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['pdf.fonttype'] = 42

	ox=[]
	y=np.array([[1,1],[1,1]])
	with open(file,'r') as f:
		f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title=r.next().pop()
		description=r.next().pop()
		xlabel=r.next().pop()
		ylabel=r.next().pop()
		legends=r.next()
		#print legends
		first_row=r.next()
		ox.append(first_row[0])
		y=np.array(first_row[1:])
		for row in r:
				ox.append(row[0])
				y=np.vstack((y,row[1:]))
 		

 	plt.figure()
 	for i in range(len(legends)-1):
		plt.plot(ox,y[:,i], linewidth=2, label=legends[i+1])

	plt.ylabel(ylabel)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.legend(loc="upper left")
	

	plt.savefig(file.replace('dat','png'),dpi=(600))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	myplot("rmsd.dat")
	#myplot(args.file_name)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 