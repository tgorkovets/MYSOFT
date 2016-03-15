#!/usr/bin/env python
# a stacked bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

EL_inner_body=93876
PhiKZ_inner_body=47520
capsid_vol=744676

A=280./(capsid_vol-PhiKZ_inner_body)*2.8*2.8


def genome_size(d,inner_body):
	s=(A/(d*d))*(capsid_vol-inner_body)
	return s

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)


print A

dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

plt.figure()

p2 = plt.plot(dd,genome_size(dd,EL_inner_body),'g-', linewidth=3)
p1 = plt.plot(dd,genome_size(dd,PhiKZ_inner_body),'b-',linewidth=3)

plt.axvline(x=2.8,color='k',ls='dashed',linewidth=2)
#horiz line
plt.axhline(y=280.,color='k',ls='dashed',linewidth=2)

plt.axvline(x=3.1,color='k',ls='dashed',linewidth=2)
#horiz line
plt.axhline(y=211.,color='k',ls='dashed',linewidth=2)

plt.ylabel('Genome size, kb')
#plt.title('Scores by group and gender')
#plt.xticks(ind+width/2., ('G1', 'G2', 'G3', 'G4', 'G5') )
#plt.yticks(np.arange(0,81,10))
plt.xlabel('Distance between DNA rods, nm')
plt.legend( (p1[0], p2[0]), ('PhiKZ', 'EL') )

#plt.show()
plt.savefig('gen_packing.png',dpi=(600))