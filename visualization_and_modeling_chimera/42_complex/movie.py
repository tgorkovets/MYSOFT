import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages


# rc('open complex.pdb')




for str_post,i in zip(range(69,-16,-1),range(1,86)):
	rc('open models/complex_%d.pdb'%str_post)
	rc('color blue :.A,.E')
	rc('color green :.B,.F')
	rc('color yellow :.C,.G')
	rc('color red :.D,.H')
	rc('color grey :.I,.J')
	rc('color orange :.P')
	rc('color cyan :.K:.L:.M:.N:.O')
	rc('matrixset ~/junk/m.dat')

	rc('2dlabel create title text \'+42 complex: T. thermophilus RNAP and nucleosome\' color white size 50 xpos 0.01 ypos 0.96')
	rc('2dlabel create h3 text \'Histone H3\' color blue size 50 xpos 0.01 ypos 0.90')
	rc('2dlabel create h4 text \'Histone H4\' color green size 50 xpos 0.01 ypos 0.86')
	rc('2dlabel create h2a text \'Histone H2A\' color yellow size 50 xpos 0.01 ypos 0.82')
	rc('2dlabel create h2b text \'Histone H2B\' color red size 50 xpos 0.01 ypos 0.78')
	rc('2dlabel create dna text \'DNA\' color grey size 50 xpos 0.15 ypos 0.90')
	rc('2dlabel create rna text \'RNA Polymerase\' color cyan size 50 xpos 0.15 ypos 0.86')
	rc('2dlabel create rnap text \'RNA strand\' color orange size 50 xpos 0.15 ypos 0.82')

	rc('2dlabel create num text \'Protected by nucleosome: %d bp\' color white size 50 xpos 0.32 ypos 0.86'%(73-str_post+1))

	rc('copy file img/complex_%d.png'%i)
	rc('close session')



#rc('copy file image.png')

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
