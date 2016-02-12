#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import subprocess
try:
    from VMD import *
    from Molecule import *
    from atomsel import *
    NOVMD=False
except:
    print "No VMD support found"
    NOVMD=True
    
from force_constants import get_consts
from scipy.optimize import minimize
import uuid

import numpy as np

__author__="Armeev Grigoriy,Alexey Shaytan"

class py3dna():
    def __init__(self,pdbfilename=None,VMD_atomsel=None,tempdir=None,path='/home/armeev/Software/Source/x3dna-v2.1'):
        if (tempdir==None):
            os.mkdir('temp')
            self.TEMP='temp'
        else:
            self.TEMP=tempdir
        
        self.set_3dna_path(path)
                
        self.AVERAGE,self.FORCE_CONST,self.DISP=get_consts()
        
        #copying some files for backbone reconstruction
        cmd=self.P_X3DNA_utils + ' cp_std BDNA'   
        p = subprocess.Popen(cmd,shell=True,cwd=self.TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.pairs_id=''
        self.distCoef=0.1
        
        
        if VMD_atomsel!=None:
            if NOVMD:
                print "VMD Libs aren't avaiable"
                return
            DNA=VMD_atomsel
            pairparam=self.X3DNA_find_pair(DNA_atomsel=DNA)
            
        elif pdbfilename!=None:            
            pairparam=self.X3DNA_find_pair(pdbfile=pdbfilename)
        else:
            print "You must provide VMD atomsel or pdb filename"
            return
        
        self.par_header,self.par_pairs,self.par_frame=self.X3DNA_analyze(pairparam)
        self.num_of_res=self.par_frame.shape[0]
        self.set_movable_bp([[1,self.num_of_res]])

    def set_3dna_path(self,path='/home/armeev/Software/Source/x3dna-v2.1'):
        self.P_X3DNA_DIR=path
        self.P_X3DNA_analyze=self.P_X3DNA_DIR+'/bin/analyze'
        self.P_X3DNA_find_pair=self.P_X3DNA_DIR+'/bin/find_pair'
        self.P_X3DNA_utils=self.P_X3DNA_DIR+'/bin/x3dna_utils'
        self.P_X3DNA_rebuild=self.P_X3DNA_DIR+'/bin/rebuild'
        os.environ['X3DNA']=self.P_X3DNA_DIR
        
    def set_pairs_list_and_dist(self,pairs_list,dist_list):
        '''
        This function sets constraints for distances between b.p.
        provide pairs and length list like
        [[1,147],[3,18],[6,22]],[20,30,40]
        for 3 distances and their lengths 
        WARNING - one bp can't be in two distance pairs
        '''
        pairs_id=[]
        for pair in pairs_list:                    
            bp1=[pair[0],self.num_of_res*2-pair[0]+1]
            bp2=[pair[1],self.num_of_res*2-pair[1]+1]
            pairs_id.append([bp1,bp2])
        self.pairs_id=np.array(pairs_id)
        self.pairs_dist=dist_list
    
    def set_movable_bp(self,listofpairs):
        '''
        This function sets b.p. as movable e.g. only them will
        be affected during minimization
        If you want bp number 1 to 8, 15 to 18, 118 to 147 (inclusive)
        to be movable provide list like
        [[1,8],[15,18],[118,147]]
        other pairs's variables will not change during minimisation
        NOTE! - numeration starts from 1 ALL the time.
        if you aren't shure about numbers create PDB with frame_to_pdb()
        and check them        
        '''
        movable_bp=[]
        for bp in listofpairs:
            movable_bp.append(range(bp[0],bp[1]))
        self.movable_bp=np.hstack(movable_bp)
    
    def X3DNA_find_pair(self,DNA_atomsel=None,pdbfile=None):
        """Performs the analysis using X3DNA

        Parameters
        ----------
        DNA_atomsel - DNA segments selected by atomsel command in VMD.
        (NOT AtomSel!)
        
        Return
        --------
        string from sdoutput of find_pair
        """
        if (DNA_atomsel != None):
            #At first we need to makup a couple of unique file names
            unique=str(uuid.uuid4())
            pdb = unique+'.pdb'
            outf = unique

            print("Writing coords to "+pdb)
            DNA_atomsel.write('pdb',self.TEMP+'/'+pdb)
        if pdbfile != None:
            pdb = os.path.abspath(pdbfile)
            
        cmd=self.P_X3DNA_find_pair+' '+pdb+' stdout'
        #print(cmd)
        p = subprocess.Popen(cmd,shell=True,cwd=self.TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for the process to terminate
        stdout, err = p.communicate()
        errcode = p.returncode
        #print('OUT:'+out)
        print('ERR:'+err)
        return(stdout)

    def X3DNA_analyze(self,ref_fp_id):
        """Performs the analysis using X3DNA

        Parameters
        ----------
        DNA_atomsel - DNA segments selected by atomsel command in VMD.
        (NOT AtomSel!)
        ref_fp_id - this is output id from X3DNA_find_pair function,
        
        Return
        --------
        header - list with header of .par file from analyse
        pairtypes - list with b.p
        par_frame - numpy array with pb params bp in rows params in columns
        """

        #Now we can run X3DNA_analyze
        cmd=self.P_X3DNA_analyze + ' stdin'
        p = subprocess.Popen(cmd,shell=True,cwd=self.TEMP,stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate(input=ref_fp_id)
        print('OUT:'+out+err)

        header,pairtypes,par_frame=self.par_to_frame(self.TEMP+'/'+'bp_step.par')

        return header,pairtypes,par_frame

    def par_to_frame(self,file):
        """
        Parse bp_step.par file as output by X3DNA
        Get a data frame with base and base-pair step
        parameters

        Note that in each line of the data frame there are parameters
        for some base pair and base pair step that preceeds (!) this base-pair.
        I.e. in the first line no base pair step parameters are specified!

        offest - offset for DNA numbering.
        """
        print "Processing ", file
        params=[]
        print file
        with open(file,'r') as f:
	        for line in f:
		        params.append(line.split())			    
        header=params[0:3]
        pairtypes=zip(*params[3:])[0]
        par_frame=np.array(zip(*params[3:])[1:]).transpose()
        return header,pairtypes,par_frame
        
    def frame_to_par(self,frame,filename):
        """
        this fuction writes analog of bp_step.par file from frame to filename.
    
        """
        string=''
        for line in self.par_header:
            for word in line[:-1]:
                string+=(str(word) + '\t')
            string+='\n'

        data1=np.append(self.par_pairs,frame.transpose().astype(str)).reshape(frame.shape[1]+1,-1).transpose()
       
        np.savetxt(filename,data1,fmt='%6s',header=string[:-1],comments='  ',delimiter='\t')

    def frame_to_pdb(self,frame,filename):
        """
        this fuction Builds pdb from frame to filename
        """
        #creating par file
        self.frame_to_par(frame,self.TEMP+'/tempfile.par')
        #rebuilding
        cmd=self.P_X3DNA_rebuild+' -atomic tempfile.par ' + filename
        p = subprocess.Popen(cmd,shell=True,cwd=self.TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        print out,err

    def frame_to_dist(self,frame):
        '''
        returns np.array of distances between pairs set by self.set_pairs_list()
        '''
        #creating par file
        self.frame_to_par(frame,self.TEMP+'/tempfile.par')
        cmd=self.P_X3DNA_rebuild+' -atomic tempfile.par stdout'
        p = subprocess.Popen(cmd,shell=True,cwd=self.TEMP,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pdb_str, err = p.communicate()
            
        coords=np.zeros((len(self.pairs_id),2,2,3))
        res=np.array([line.split() for line in pdb_str.split('\n')])
        for line in np.array(res):
            try:
                if (line[2]=='P') and (line[0]=='ATOM'):
                    if int(line[5]) in self.pairs_id:   
                        coords[np.where(self.pairs_id==int(line[5]))]=np.array([line[6],line[7],line[8]]).astype(float)
            except:
               continue
        bp_centers=coords.sum(axis=2)/2.0       
                
        return np.linalg.norm(bp_centers[:,0]-bp_centers[:,1],axis=1)
        
    def frame_to_energy(self,frame,usepairs=False):
        '''
        returns energy of dna conformation, calculated with force params
        usepairs flag enables calculation of distances between pairs set
        by self.set_pairs_list()
        '''
        self.copy_frame=self.par_frame.copy()
        energy=self.calc_energy(frame[self.movable_bp,6:12],usepairs=usepairs,verbose=false)
        return energy
        
    def run_minimize(self,frame=None,usepairs=False,verbose=True,maxiter=100,maxfev=100):
        '''
        Preformes minimization of energy
        NOTE! it only works with b.p. steps selected with set_movable_bp()
        usepairs flag enables calculation of distances between pairs set
        by self.set_pairs_list()
        maxiter- amount of steps for CG descent
        maxfev- amount of calculations of gradient?
        '''
        if frame==None:
            frame=self.par_frame

        self.copy_frame=self.par_frame.copy()
        x0=frame[self.movable_bp,6:12].flatten().astype(float)
        res = minimize(self.calc_energy, x0, (usepairs,verbose), method='Powell',options={'maxiter':maxiter,'maxfev': maxfev,'xtol' :2.0, 'disp': True})
        local_frame=self.par_frame.copy()
            
        local_frame[self.movable_bp,6:12]=np.array(res.x).reshape(-1,6)
        return local_frame,res
        
    def calc_energy(self,frame,usepairs=False,verbose=True):
        '''
        !!!DO NOT USE THIS FUNCTION DIRECTLY!!!
        instead use frame_to_energy()
        returns energy of dna conformation, calculated with force params
        usepairs flag enables calculation of distances between pairs set
        by self.set_pairs_list()
        '''
        local_frame=self.copy_frame
        local_frame[self.movable_bp,6:12]=np.array(frame).reshape(-1,6)
        subframe=local_frame[1:,6:12]
        dif_matrix=[]
        force_matrix=[]
        for i in range(len(self.par_pairs)-1):
            step=str(self.par_pairs[i][0]+self.par_pairs[i+1][0])
            a=subframe[i].astype(float)-self.AVERAGE[step].astype(float)
            dif_matrix.append(np.tile(a,(6,1)).transpose()*np.tile(a,(6,1)))
            force_matrix.append(self.FORCE_CONST[step])
        
        if not(usepairs):
            result = float(np.multiply(dif_matrix, force_matrix).sum()/2.0)
        else:
            result = (float(np.multiply(dif_matrix, force_matrix).sum()/2.0)            
             + self.distCoef*np.power(self.pairs_dist-self.frame_to_dist(local_frame),2).sum())
        print result
        return result
