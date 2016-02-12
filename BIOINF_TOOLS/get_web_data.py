# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan
Library with different tools, that helps to download NCBI data via E-utils

"""
# from __future__ import division
# from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import csv
import collections
from Bio import Entrez
from Bio import SeqIO
import subprocess
import uuid

# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML
# import itertools

#sysg.path.append('/Users/alexeyshaytan/work_HD/histones_work/13_hist_seq_analysis/sec_str')

Entrez.email = "alexey.shaytan@nih.gov" 

def get_fasta_records(gi_list):
    """
    Downloads fasta records returns as dict, with gi=>(id,sequence)
    """
    print("Getting FASTA ")
    num=len(gi_list)
    records_fasta=dict()
    for i in range(int(num/1000)+1):
        print("Fetching %d th thousands from %d"%(i,num))
        strn = ",".join(map(str,gi_list)[i*1000:(i+1)*1000])
        # print strn
        request=Entrez.epost(db="protein",id=strn)
        result=Entrez.read(request)
        webEnv=result["WebEnv"]
        queryKey=result["QueryKey"]
        handle=Entrez.efetch(db="protein",rettype='fasta',retmode='text',webenv=webEnv, query_key=queryKey)
        for r in SeqIO.parse(handle,'fasta'):
            records_fasta[r.id.split('|')[1]]=(r.id,r.seq)
    print("FASTA Records downloaded for proteins:")
    print(len(records_fasta))
    return records_fasta

def get_pdb_assemblies(pdb_code_list,full_path):
    """
    Downloads a list of macromolecular assemblies for PDB codes and saves it to path
    rsync -rlpt -v -z --delete --port=33444 \
    rsync.wwpdb.org::ftp_data/structures/divided/pdb/ ./pdb
    """
    fn=str(uuid.uuid4())
    with(open(fn,'wb')) as f:
        for i in pdb_code_list:
            f.write('{0}.pdb1.gz\n'.format(i.lower()))
    subprocess.call(['rsync', '-rlpt','-z','-v','-L','--port=33444','rsync.wwpdb.org::ftp_data/biounit/PDB/all/','--files-from={0}'.format(fn),'{0}/'.format(full_path)])
    os.system('rm %s'%fn)

def get_pdbs(pdb_code_list,full_path):
    """
    Downloads a list of plain pdb files for PDB codes and saves it to path
    rsync -rlpt -v -z --delete --port=33444 \
    rsync.wwpdb.org::ftp_data/structures/divided/pdb/ ./pdb
    """
    fn=str(uuid.uuid4())
    with(open(fn,'wb')) as f:
        for i in pdb_code_list:
            f.write('pdb{0}.ent.gz\n'.format(i.lower()))
    subprocess.call(['rsync', '-rlpt','-z','-v','-L','--port=33444','rsync.wwpdb.org::ftp_data/structures/all/pdb/','--files-from={0}'.format(fn),'{0}/'.format(full_path)])
    os.system('rm %s'%fn)

# Here is some other example
# def download():
#     gis = list_gis()
#     base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
#     db = 'protein'

#     for gi in gis:
#         try:
#             with open("download_gi/{}.faa".format(gi), 'r'):
#                 pass
#         except:
#             print "Loading", gi
#             try:
#                 with open("download_gi/{}.faa".format(gi), 'w') as o:
#                     url = "{}efetch.fcgi?db={}&id={}&rettype=fasta&retmode=text".format(base, db, gi)
#                     print url
#                     response = urllib2.urlopen(url)
#                     fasta = response.read()
#                     o.write(fasta)
#             except Exception, e:
#                 print e
#                 print "Failed", gi
#                 # raise



if __name__ == '__main__':
    get_pdb_assemblies(['1aoi','1kx5'],'/Users/alexeyshaytan/junk/kuku/')
    
            