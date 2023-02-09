#!/home/kydoh/.conda/envs/emto/bin/python
import os
import sys
import time
import random
import subprocess
from itertools import permutations # for L12, D019
from itertools import combinations # for bcc, fcc, hcp, B2
import numpy as np
import pyemto
'''
Update
20230126 Add check_error, replace_efgs and modify run_emto
20230127 Add convergence test in check_error
20230130 kgrn process termination problem, add kgrn.kill()
20230130 load previous calculations
20230131 shell script combine into this file
'''

def gen_input(material,lat,jobname,latpath,emtopath,atoms,concs,splts,iqs,its,itas,**kwargs):
   material.bulk_new(lat=lat,jobname=jobname,latpath=latpath,atoms=atoms,concs=concs,splts=splts,iqs=iqs,its=its,itas=itas,**kwargs)
   #for key, value in kwargs.items():
   #   material.emto.set_values(key,value)
   material.emto.kgrn.write_input_file(folder=emtopath)
   material.emto.kfcd.write_input_file(folder=emtopath)

def run_kgrn(jobid,EMTOdir='/usr/local/bin'):
   with open(jobid+'.kgrn','r') as infile:
      kgrn = subprocess.Popen(EMTOdir+'/kgrn_cpa',stdin=infile)
   ### running check ###
   while kgrn.poll() == None:
      time.sleep(10)
      ### Too many iterations ###
      error = 'Too many iterations'
      cmd = f"grep {error} kgrn/{jobid}.prn|wc -l"
      p = subprocess.run(cmd,shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
      if int(p.stdout) > 0:
         kgrn.terminate()
         time.sleep(1)
         if kgrn.poll() == None:
            kgrn.kill()
            time.sleep(1)
         print(f'ERROR: {error}')
         return 1
         
   ### ending check ###
   ### Fermi level not found ###
#   error = 'Fermi level not found'
   # 여러가지 에러 추가하기

   ### No error ###
   ### Convergence test ###
   cmd = f"grep Converged kgrn/{jobid}.prn|wc -l"
   p = subprocess.run(cmd,shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      universal_newlines=True)
   if int(p.stdout) > 0:
      print(f'{jobid} is converged.')
      return 0
   else:
      print(f'{jobid} is NOT converged.')
      return 1

def run_kfcd(jobid,EMTOdir='/usr/local/bin'):
   with open(jobid+'.kfcd','r') as infile:
      subprocess.run(EMTOdir+'/kfcd_cpa',stdin=infile)

def replace_efgs(jobid,efgs):
   cmd = f"sed -i '26 c\EFGS...=  {efgs:.3f} HX....=  0.100 NX...=  5 NZ0..=  6 STMP..= N' {jobid}.kgrn"
   subprocess.run(cmd,shell=True)

def run_emto(jobid,folder,emtopath,EMTOdir='/usr/local/bin'):
   '''
   if 'Too many iterations' issue raise than kill process and change efgs and rerun
   if fermi level issue than change efgs and rerun
   '''
   os.chdir(emtopath)
   check_error = run_kgrn(jobid)
   for count in range(1,401):
      if check_error == 1:
         replace_efgs(jobid,count/100)
         check_error = run_kgrn(jobid)
      else:
         run_kfcd(jobid)
         break

   os.chdir(folder)

def count_data(folder, jobname):
   cmd = "awk 'BEGIN{n=0}{if(NF == 2){n+=1}}END{print n}' emto/fit/"+jobname+".dat"
   p = subprocess.run(cmd,shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      universal_newlines=True)
   return int(p.stdout)
#   child = subprocess.Popen([folder+'/count_data.sh',jobname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#   (stdout,stderr) = child.communicate()
#   return int(stdout)

def get_dE(folder, jobname):
   cmd = "touch emto/fit/"+jobname+".dat;awk 'BEGIN{E_old=-1;E=0}{if(NF == 2){E_old=E;E=$2}}END{print E-E_old}' emto/fit/"+jobname+".dat"
   p = subprocess.run(cmd,shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      universal_newlines=True)
   return float(p.stdout)
#   child = subprocess.Popen([folder+'/get_dE.sh',jobname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#   (stdout,stderr) = child.communicate()
#   return float(stdout)

def get_fp(folder, jobname):
   cmd = "awk 'BEGIN{s2=0;s1=0;s=0;E2=0;E1=0;E=0}{if(NF == 2){s2=s1;s1=s;s=$1;E2=E1;E1=E;E=$2}}END{print E*(2*s-s1-s2)/(s-s1)/(s-s2)+E1*(s-s2)/(s1-s)/(s1-s2)+E2*(s-s1)/(s2-s)/(s2-s1)}' emto/fit/"+jobname+".dat"
   p = subprocess.run(cmd,shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      universal_newlines=True)
   return float(p.stdout)
#   child = subprocess.Popen([folder+'/get_fp.sh',jobname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#   (stdout,stderr) = child.communicate()
#   return float(stdout)

def get_fpp(folder, jobname):
   cmd = "awk 'BEGIN{s2=0;s1=0;s=0;E2=0;E1=0;E=0}{if(NF == 2){s2=s1;s1=s;s=$1;E2=E1;E1=E;E=$2}}END{print 2*E/(s-s1)/(s-s2)+2*E1/(s1-s)/(s1-s2)+2*E2/(s2-s)/(s2-s1)}' emto/fit/"+jobname+".dat"
   p = subprocess.run(cmd,shell=True,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      universal_newlines=True)
   return float(p.stdout)
#   child = subprocess.Popen([folder+'/get_fpp.sh',jobname],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#   (stdout,stderr) = child.communicate()
#   return float(stdout)

def next_sws(folder,jobname,sws,damp):
   fp = get_fp(folder,jobname)
   fpp = get_fpp(folder,jobname)
   sws = sws - fp/(abs(fpp)+damp)
   if sws > 2 and sws < 4:
      return sws
   else:
      return random.random()*2 + 2

def update_damp(dE,damp,damp_lowering,damp_raising):
   if dE < 0:
      return damp/damp_lowering
   else :
      return 1
      #return damp*damp_raising

def check_done(folder,jobname):
   if abs(get_dE(folder, jobname)) < 0.00001:
      print(f'{jobname} done')
      return True
   else:
      print(f'{jobname} NOT yet')
      return False

def generate_job_list(lat,element_list):
   job_list = list()
   # For bcc or fcc
   if lat == 'bcc' or lat == 'fcc':
      nky   = 21
      nkz   = 0
      # pure metal
      for x in element_list:
         jobname = x
         atoms = [x]
         concs = [1]
         splts = [1]
         iqs   = [1]
         its   = [1]
         itas  = [1]
         ncpa  = 1
         job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
      # binary RSS
      element_pairs = list(combinations(element_list,2))
      for x in element_pairs:
         for n in range(1,10):
            jobname = f'{x[0]}{n}{x[1]}{10-n}'
            atoms   = [x[0],x[1]]
            concs   = [n/10,1-n/10]
            splts   = [1,1]
            iqs     = [1,1]
            its     = [1,1]
            itas    = [1,2]
            ncpa    = 11
            job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
   # For hcp
   elif lat == 'hcp':
      nky   = 21
      nkz   = 15
      # pure metal
      for x in element_list:
         jobname = x
         atoms = [x,x]
         concs = [1,1]
         splts = [1,1]
         iqs   = [1,2]
         its   = [1,2]
         itas  = [1,1]
         ncpa  = 1
         job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
      # binary RSS
      element_pairs = list(combinations(element_list,2))
      for x in element_pairs:
         for n in range(1,10):
            jobname = f'{x[0]}{n}{x[1]}{10-n}'
            atoms   = [ x[0],x[1], x[0],x[1] ]
            concs   = [ n/10,1-n/10, n/10,1-n/10 ]
            splts   = [ 1,1, 1,1 ]
            iqs     = [ 1,1, 2,2 ]
            its     = [ 1,1, 2,2 ]
            itas    = [ 1,2, 1,2 ]
            ncpa    = 11
            job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
   # For B2
   elif lat == 'B2':
      nky = 11
      nkz = 0
      ncpa = 1
      # binary ordered phase
      element_pairs = list(combinations(element_list,2))
      for x in element_pairs:
         jobname = f'{x[0]}{x[1]}'
         atoms   = [x[0],x[1]]
         concs   = [1,1]
         splts   = [1,1]
         iqs     = [1,2]
         its     = [1,2]
         itas    = [1,1]
         job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
   # For L12
   elif lat == 'L12':
      nky = 11
      nkz = 0
      ncpa = 1
      # binary ordered phase
      element_pairs = list(permutations(element_list,2))
      for x in element_pairs:
         jobname = f'{x[0]}{x[1]}3'
         atoms   = [x[0],x[1],x[1],x[1]]
         concs   = [1,1,1,1]
         splts   = [1,1,1,1]
         iqs     = [1,2,3,4]
         its     = [1,2,3,4]
         itas    = [1,1,1,1]
         job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
   elif lat == 'D019':
      nky = 11
      nkz = 7
      ncpa = 1
      # binary ordered phase
      element_pairs = list(permutations(element_list,2))
      for x in element_pairs:
         jobname = f'{x[0]}{x[1]}3'
         atoms   = [x[0],x[1],x[1],x[1],x[1],x[1],x[1],x[0]]
         concs   = [1,1,1,1,1,1,1,1]
         splts   = [1,1,1,1,1,1,1,1]
         iqs     = [1,2,3,4,5,6,7,8]
         its     = [1,2,3,4,5,6,7,8]
         itas    = [1,1,1,1,1,1,1,1]
         job_list.append([jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz])
   else:
      sys.exit(f'{lat} is NOT supported yet!!!')
   return job_list

################## MAIN ######################
# set parameters of materials
## path ##
folder = os.getcwd()
latpath = folder+"/lat"
emtopath = folder+"/emto"
EMTOdir = "/usr/local/bin"
material = pyemto.System(folder=emtopath)
## USER INPUTS ##
lat = sys.argv[1]   
element_list = ['Cr','Fe','Ni','Co','Ta','Al','V','Nb']
## common ##
niter = 300
afm = 'F'    
sofc = 'Y'
nz2 = 16
depth = 0.94
amix = 0.01
efmix = 0.9
tole = 1.e-8
iex = 7
dirac_np = 1001
nes = 50
dirac_niter = 500  # relaxation condition

# set hyperparameters of optimization
Nmax = 50
Ediff = 0.000001
damp_lowering = 3
damp_raising = 2

# generate job list
job_list = generate_job_list(lat,element_list)

for jobname,atoms,concs,splts,iqs,its,itas,ncpa,nky,nkz in job_list:
   # check done job
   if check_done(folder,jobname):
      continue

   # initialize data
   subprocess.run(['rm',folder+'/emto/fit/'+jobname+'.dat'])
   subprocess.run(['touch',folder+'/emto/fit/'+jobname+'.dat'])

   # load previous calculations
   cmd = f"grep TOT-PBE {folder}/emto/kfcd/{jobname}*|awk '{{print \$8,\$5}}'|sort -nk2 -r > {folder}/emto/fit/{jobname}.dat"
   p = subprocess.run(cmd,shell=True)
   
   # initialize variables of materials
   sws = random.random()*0.3 + 2.8     # initialize for each phase 

   # initialize variables of optimization
   N=0
   dE = 1
   damp = 1
   loop = 0

   while( loop < Nmax and abs(dE) > Ediff ):
      loop+=1
      # generate input
      gen_input(material,lat,jobname,
             latpath,emtopath,
             atoms,concs,splts,
             iqs=iqs,its=its,itas=itas,
             sws=sws,afm=afm,nky=nky,nkz=nkz,depth=depth,
             niter=niter,ncpa=ncpa,sofc=sofc,nz2=nz2,
             amix=amix,efmix=efmix,tole=tole,iex=iex,
             dirac_np=dirac_np,nes=nes,dirac_niter=dirac_niter)
      # run kgrn & kfcd
      jobid = f'{jobname}_{sws:8f}'
      run_emto(jobid,folder,emtopath)
      # collect data
      cmd = "grep TOT-PBE emto/kfcd/"+jobid+".prn |awk '{print $7,$4}' >> emto/fit/"+jobname+".dat"
      subprocess.run(cmd,shell=True)
#      subprocess.run([folder+'/kfcd_prn2sws_e.sh',jobid,jobname])
      # set next sws
      N = count_data(folder,jobname)
      if N < 3 :
         sws += 0.015
         print(f'Now_sws= {sws}')
      else :
         sws = next_sws(folder,jobname,sws,damp)
         #print(f'Now_sws= {sws}, fp= {fp}, fpp= {fpp}')
         dE = get_dE(folder,jobname)
         print(f'dE= {dE}')
         damp = update_damp(dE,damp,damp_lowering,damp_raising)


