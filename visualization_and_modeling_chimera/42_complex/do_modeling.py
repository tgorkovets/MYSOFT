import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

os.environ['AMBERHOME'] = '/Users/alexeyshaytan/soft/amber12'

# change to folder with data files
seq_i='ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAGAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT'

#Position of the last base on I strand that is curled in DNA
for str_post in range(-15,70):
# for str_post in range(5,6):

#str_post=-1



	str_seq=seq_i[73-18:str_post+73+1]
	os.system('./gen_dna.x %s dna.pdb'%str_seq)

	rc('open 1kx5.pdb')
	rc('delete #0:HOH')
	rc('delete #0:CL')
	rc('delete #0:-73-%d.I,%d-73.J'%(str_post-1,(str_post-1)*(-1)))
	rc('delete #0:MN')


	rc('open dna.pdb')
	rc('resrenumber -18 #1:.I')
	rc('resrenumber %d #1:.J'%(str_post*(-1)))
	rc('write format pdb #1 dna.pdb')

	##Now let's do superpositions

	rc('match #1:%d.I,%d.J@N1,N2,N3,N4,N7,N9 #0:%d.I,%d.J@N1,N2,N3,N4,N7,N9'%(str_post,(-1)*str_post,str_post,str_post*(-1)))

	rc('open 2o5i_renum.pdb')
	rc('delete #2:HOH')
	rc('delete #2:MN')
	rc('delete #2:ZN.J')

	rc('changechains A,B,C,D,E,H K,L,M,N,O,P #2')

	rc('match #2:-18.I,18.J@N1,N2,N3,N4,N7,N9 #1:-18.I,18.J@N1,N2,N3,N4,N7,N9')

	rc('delete #1:-18.I,18.J')
	rc('delete #1:%d.I,%d.J'%(str_post,(-1)*str_post))


	rc('combine #0,1,2 name temp newchainids false close false')

	#Let's add bonds
	rc('bond #3:-18.I@O3\':-17.I@P')
	rc('bond #3:%d.I@O3\':%d.I@P'%(str_post-1,str_post))
	rc('bond #3:18.J@P:17.J@O3\'')
	rc('bond #3:%d.J@P:%d.J@O3\''%((-1)*str_post+1,(-1)*str_post))

	rc('write format pdb #3 complex_tt.pdb')

	os.system('cp complex_tt.pdb models_tt/complex_%d.pdb'%str_post)

	rc('close session')

#at this point we have model of T. thermophilus complex with nucleosome as model #3
#Let's load E.coli superimpose it and replace it
	rc('open complex_tt.pdb')
	rc('open 4jkr.pdb')
	rc('delete #1:HOH')
	rc('delete #1:SR')
	rc('delete #1:.F')
	rc('changechains A,B,C,D,E Q,R,S,T,U #1')
	#align using beta chain
	rc('mmaker #0:.N #1:.T pair ss alg sw')
	rc('combine #0,1 newchainids false close false name temp2')
	rc('delete #2:.K:.L:.M:.N:.O')
	

	rc('write format pdb #2 complex_ecoli.pdb')

	os.system('cp complex_ecoli.pdb models_ecoli/complex_%d.pdb'%str_post)

	rc('close session')
# # gather the names of .pdb files in the folder
# file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

# # loop through the files, opening, processing, and closing each in turn
# for fn in file_names:
# 	replyobj.status("Processing " + fn) # show what file we're working on
# 	rc("open " + fn)
# 	rc("align ligand ~ligand") # put ligand in front of remainder of molecule
# 	rc("focus ligand") # center/zoom ligand
# 	rc("surf") # surface receptor
# 	rc("preset apply publication 1") # make everything look nice
# 	rc("surftransp 15") # make the surface a little bit see-through
# 	# save image to a file that ends in .png rather than .pdb
# 	png_name = fn[:-3] + "png"
# 	rc("copy file " + png_name + " supersample 3")
# 	rc("close all")
# uncommenting the line below will cause Chimera to exit when the script is done
#rc("stop now")
# note that indentation is significant in Python; the fact that
# the above command is exdented means that it is executed after
# the loop completes, whereas the indented commands that 
# preceded it are executed as part of the loop.
