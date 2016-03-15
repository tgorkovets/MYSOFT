#!/usr/bin/env python

"""
Program for plotting dat files with legends and names
Modified for foldX

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
import itertools


def myplot(file1, file2):
  "Plot dat file"


  rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
  plt.rcParams['ps.useafm'] = True
  rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
  plt.rcParams['pdf.fonttype'] = 42

  ox=[]
  y=np.array([[1,1],[1,1]])
  with open(file1,'r') as f:
    f=open(file1,'r')
    r=csv.reader(f,delimiter='\t')
    for i in range(6):
      print r.next()

    title1=r.next().pop()
    r.next()
		
    #description=r.next().pop()
		#xlabel=r.next().pop()
		#ylabel=r.next().pop()
    legends1=r.next()
		#print legends
    first_row=r.next()
    print first_row[1:]
    print legends1[1:]
    ox.append(first_row[0])
		# y=np.array(first_row[1:])
		# for row in r:
		# 		ox.append(row[0])
		# 		y=np.vstack((y,row[1:]))
  with open(file2,'r') as f:
    f=open(file2,'r')
    r=csv.reader(f,delimiter='\t')
    for i in range(6):
      print r.next()

    title2=r.next().pop()
    r.next()
    
    #description=r.next().pop()
    #xlabel=r.next().pop()
    #ylabel=r.next().pop()
    legends2=r.next()
    #print legends
    first_row2=r.next()
    print first_row2[1:]
    print legends2[1:]
    ox.append(first_row[0])

  data=list(itertools.chain.from_iterable(zip(first_row[1:],first_row2[1:])))
  legend=list(itertools.chain.from_iterable(zip(legends1[1:],['']*len(legends1[1:]))))

  plt.figure()
  p1=plt.bar(range(len(legends1[1:])),map(float,first_row[1:]),width=0.4,align='center',color="yellow",label=legends1[0])
  p2=plt.bar(np.arange(0.5,len(legends1[1:])+0.5,1),map(float,first_row2[1:]),width=0.4,align='center',color="green",label=legends2[0])
  plt.title(title1+' '+title2)
  plt.ylabel("Energy, kcal/mol")
    #p1.margins(y=.1, x=.1)
  plt.autoscale(enable=True, axis='both', tight="True")
  plt.xticks(range(len(legends1[1:])), legends1[1:], size='large',rotation='vertical')
  plt.xlabel("Terms")
  plt.tight_layout()

  plt.savefig((file1.replace('.txt','')+'_'+file2.replace('txt','png')),dpi=(600))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	myplot("stability_3lz0.txt","stability_3afa.txt")
	#myplot(args.file_name)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 