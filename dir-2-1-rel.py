from decimal import *
getcontext().prec = 7

import sys
inp = sys.argv[1]

import numpy as np
from numpy import linalg

import datetime  

import os.path

import subprocess

#create matrix for atomic positions from .xyz file

with open(inp+'.xyz') as f:
  list = f.readlines()
  
matrix = [[pos[0:2],Decimal(pos[10:20]),Decimal(pos[26:36]),Decimal(pos[42:52])] for pos in list]

nat = len(matrix)

atoms = [(row[0]) for row in matrix]

atom_count = dict()

for i in atoms:
	atom_count[i] = atom_count.get(i, 0) + 1


atomic_numbers = dict([('H ', 1), ('He', 2), ('Li', 3), ('Be', 4), ('B ', 5), ('C ', 6), 
('N ', 7), ('O ', 8), ('F ', 9), ('Ne', 10), ('Na', 11), ('Mg', 12), ('Al', 13), ('Si', 14),
('P ', 15), ('S ', 16), ('Cl', 17), ('Ar', 18), ('K ', 19), ('Ca', 20), ('Sc', 21), ('Ti', 22),
('V ', 23), ('Cr', 24), ('Mn', 25), ('Fe', 26), ('Co', 27), ('Ni', 28), ('Cu', 29), ('Zn', 30),
('Ga', 31), ('Ge', 32), ('As', 33), ('Se', 34), ('Br', 35), ('Kr', 36), ('Rb', 37), ('Sr', 38),
('Y ', 39), ('Zr', 40), ('Nb', 41), ('Mo', 42), ('Tc', 43), ('Ru', 44), ('Rh', 45), ('Pd', 46),
('Ag', 47), ('Cd', 48), ('In', 49), ('Sn', 50), ('Sb', 51), ('Te', 52), ('I ', 53), ('Xe', 54),
('Cs', 55), ('Ba', 56), ('La', 57), ('Ce', 58), ('Pr', 59), ('Nd', 60), ('Pm', 61), ('Sm', 62),
('Eu', 63), ('Gd', 64), ('Tb', 65), ('Dy', 66), ('Ho', 67), ('Er', 68), ('Tm', 69), ('Yb', 70),
('Lu', 71), ('Hf', 72), ('Ta', 73), ('W ', 74), ('Re', 75), ('Os', 76), ('Ir', 77), ('Pt', 78),
('Au', 79), ('Hg', 80), ('Tl', 81), ('Pb', 82), ('Bi', 83), ('Po', 84), ('At', 85), ('Rn', 86),
('Fr', 87), ('Ra', 88), ('Ac', 89), ('Th', 90), ('Pa', 91), ('U ', 92), ('Np', 93), ('Pu', 94),
('Am', 95), ('Cm', 96), ('Bk', 97), ('Cf', 98), ('Es', 99)])


#read in basis sets, functional and displacements from .inp file

with open(inp+'.inp') as f:
	basis = {}
	isotopes = {}
	
	for line in f:
		if line == 'METHOD\n':
			continue
		elif line == 'BASIS\n':
			break
		else:
			method = line.replace("\n","")

	for line in f:
		if line == 'BASIS\n':
			continue
		elif line == 'FUNCTIONAL\n':
			break
		else:
			basis[line[:2]] = line.replace("\n","")[6:]
	
	for line in f:
		if line == 'FUNCTIONAL\n':
			continue
		elif line == 'DISPLACEMENT_1\n':
			break
		else:
			functional = line.replace("\n","")
		
	for line in f:
		if line == 'DISPLACEMENT_1\n':
			continue
		elif line == 'DISPLACEMENT_2\n':
			break
		else:
			dis = Decimal(line[0:].replace("\n",""))
	
	for line in f:
		if line == 'DISPLACEMENT_2\n':
			continue
		elif line == 'ISOTOPES\n':
			break
		else:
			dis2 = Decimal(line[0:].replace("\n",""))
			
	for line in f:
		if line == 'ISOTOPES\n':
			continue
		elif line[0:3] == 'END':
			break
		else:
			isotopes[line[:2].replace(" ","")] = float(int(line.replace("\n","")[6:]))


		
#define functions to create dirac mol and run files		
		
def create_mol(name, coord_matrix):
		
	with open(name+'.mol', 'w+') as f:
		f.write('INTGRL\n')
		f.write('\n')
		f.write('\n')
		f.write('C   '+str(len(atom_count))+'    0         A\n')

		j = 0
		for i in range(len(atom_count)):
			el = atoms[j]
			f.write(str(atomic_numbers[el]).rjust(8)+'.'+str(atom_count[el]).rjust(5)+'\n')
			for a in range(nat):
				row = matrix[a]
				row2 = coord_matrix[a]
				if row[0] == el:
					f.write(str(row[0]).ljust(2)+str(round(row2[0],6)).rjust(17)+str(round(row2[1],6)).rjust(15)+str(round(row2[2],6) ).rjust(15)+'\n')
			f.write('LARGE BASIS '+ basis[el]+'\n')
			j = j + atom_count[el]
		f.write('FINISH')


def create_run1(name):

	with open('run-'+name+'_4c-dir', 'w+') as f:
		f.write('#!/bin/bash\n')
		f.write('\n')
		f.write('#PBS -M aaa@aaa.pl\n')
		f.write('#PBS -m abe\n')
		f.write('#PBS -l nodes=1:ppn=1\n')
		f.write('#PBS -l mem=10gb\n')
		f.write('#PBS -l walltime=360:00:00\n')
		f.write('#PBS -o /STORAGE/DATA/kjakubow/stdout_${PBS_JOBID}\n')
		f.write('#PBS -e /STORAGE/DATA/kjakubow/stderr_${PBS_JOBID}\n')
		f.write('\n')
		f.write('cd $PBS_O_WORKDIR\n')
		f.write('current=`pwd`\n')
		f.write('\n')
		f.write('export GFORTRAN_UNBUFFERED_ALL=Y\n')
		f.write('\n')
		f.write('scratch=/STORAGE/DATA/$USER/$PBS_JOBID\n')
		f.write('mkdir -p $scratch\n')
		f.write('\n')
		f.write('export pam=/home/apps/dirac/dirac-17.0/share/dirac/pam\n')
		f.write('export inp='+inp+'-J-4c\n')
		f.write('export mol='+name+'\n')
		f.write('export out=$inp\_$mol.out\n')
		f.write('\n')
		f.write('$pam --noarch --scratch=$scratch --mw=2100 --mol=$mol.mol --inp=$inp.inp --put="DFCOEF" --get="AOPROPER" --outcmo\n')
		f.write('\n')
		f.write('exit 0\n')
		
def create_run2(name):

	with open('run-'+name+'_4c-dir', 'w+') as f:
		f.write('#!/bin/bash\n')
		f.write('\n')
		f.write('#PBS -M aaa@aaa.pl\n')
		f.write('#PBS -m abe\n')
		f.write('#PBS -l nodes=1:ppn=1\n')
		f.write('#PBS -l mem=10gb\n')
		f.write('#PBS -l walltime=360:00:00\n')
		f.write('#PBS -o /STORAGE/DATA/kjakubow/stdout_${PBS_JOBID}\n')
		f.write('#PBS -e /STORAGE/DATA/kjakubow/stderr_${PBS_JOBID}\n')
		f.write('\n')
		f.write('cd $PBS_O_WORKDIR\n')
		f.write('current=`pwd`\n')
		f.write('\n')
		f.write('export GFORTRAN_UNBUFFERED_ALL=Y\n')
		f.write('\n')
		f.write('scratch=/STORAGE/DATA/$USER/$PBS_JOBID\n')
		f.write('mkdir -p $scratch\n')
		f.write('\n')
		f.write('export pam=/home/apps/dirac/dirac-17.0/share/dirac/pam\n')
		f.write('export inp='+inp+'-E-4c\n')
		f.write('export mol='+name+'\n')
		f.write('export out=$inp\_$mol.out\n')
		f.write('\n')
		f.write('$pam --noarch --scratch=$scratch --mw=2100 --mol=$mol.mol --inp=$inp.inp --put="DFCOEF" --get="AOPROPER" --outcmo\n')
		f.write('\n')
		f.write('exit 0\n')
		
#create dirac inp files

if method == 'rel':

	with open(inp+'-J-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('.PROPERTIES F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*TWOINT\n')
		f.write('.SCREEN\n')
		f.write(' 1.0D-18\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.LVCORR\n')
		f.write('.URKBAL\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('**PROPERTIES\n')
		f.write('.SPIN-SPIN COUPLING\n')
		f.write('*END OF\n')
	
	with open(inp+'-E-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*TWOINT\n')
		f.write('.SCREEN\n')
		f.write(' 1.0D-18\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.LVCORR\n')
		f.write('.URKBAL\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('*END OF\n')
		
elif method == 'nrel':

	with open(inp+'-J-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('.PROPERTIES F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*TWOINT\n')
		f.write('.SCREEN\n')
		f.write(' 1.0D-18\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.LVCORR\n')
		f.write('.URKBAL\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**GENERAL\n')
		f.write('.CVALUE\n')
		f.write('2000.0\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('**PROPERTIES\n')
		f.write('.SPIN-SPIN COUPLING\n')
		f.write('*END OF\n')
	
	with open(inp+'-E-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*TWOINT\n')
		f.write('.SCREEN\n')
		f.write(' 1.0D-17\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.LVCORR\n')
		f.write('.URKBAL\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**GENERAL\n')
		f.write('.CVALUE\n')
		f.write('2000.0\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('*END OF\n')
		
elif method == 'noso':

	with open(inp+'-J-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('.PROPERTIES F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.NOSPIN\n')
		f.write('.LVCORR\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('**PROPERTIES\n')
		f.write('.SPIN-SPIN COUPLING\n')
		f.write('*END OF\n')
	
	with open(inp+'-E-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
		f.write('*READIN\n')
		f.write('.UNCONTRACT\n')
		f.write('**HAMILTONIAN\n')
		f.write('.NOSPIN\n')
		f.write('.LVCORR\n')
		f.write('.DFT\n')
		f.write(' '+functional+'\n')
		f.write('**WAVE FUNCTIONS\n')
		f.write('.SCF\n')
		f.write('*SCF\n')
		f.write('.EVCCNV\n')
		f.write('1.0E-8  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('50\n')
		f.write('*END OF\n')
		
#sum up of the input

with open(inp+'-'+method+'.out', 'w+') as f:
	now = datetime.datetime.today()
	f.write('Date and time: '+now.strftime("%Y-%m-%d %H:%M:%S")+'\n\n')
	f.write('***   Input summary   ***\n\n')
	if method == 'rel':
		f.write('Method: relativistic\n\n')
	elif method == 'nrel':
		f.write('Method: nonrelativistic\n\n')
	elif method == 'noso':
		f.write('Method: relativistic with spin-orbit effects switched off\n\n')
	f.write('Displacement 1: '+str(dis)+' Angstrom\n')
	f.write('Displacement 2: '+str(dis2)+' bohr * amu ** 1/2\n\n')
	f.write('Geometry:\n')
	for k in range(nat):
		lline = matrix[k]
		f.write(lline[0].rjust(5)+str(lline[1]).rjust(12)+str(lline[2]).rjust(12)+str(lline[3]).rjust(12)+'\n')
	f.write('\n')
	f.write('Basis set:\n')
	for record in basis:
		space = "        "
		f.write(space+record+": "+basis[record]+'\n')
	f.write('\n')
	f.write('Isotopes:\n')
	for record in isotopes:
		space = "        "
		f.write(space+record.ljust(2)+": "+str(int(isotopes[record]))+'\n')
	f.write('\n')
	f.write('DFT functional: '+functional+'\n\n\n')

#check if all energy files exist

a = 0

#for the molecule

if os.path.isfile(inp+'-4c_'+inp+'.out') == False:
	subprocess.call(['qsub run-'+inp+'-4c-dir'], shell = True)
	a = a + 1
			
#all singular displacements

for i in range(nat):
	
	for j in range(1,4):
		if os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+.out') == False:
			subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+_4c-dir'], shell = True)
			a = a + 1

			
		elif os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-.out') == False:
			subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-_4c-dir'], shell = True)
			a = a + 1
		else:
			continue		
			
#for all double displacements 

#for ++ and --
for i in range(nat):	
	for j in range(1,4):
		for k in range(nat):
			for l in range (1,4):
				if k > i:
					if os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+.out') == False:
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
						a = a + 1

						
					elif os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-.out') == False:
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
						a = a + 1
					else:
						continue

	
				elif k == i:
					if l > j:
						if os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+.out') == False:
							subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
							a = a + 1
					
						elif os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-.out') == False:
							subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
							a = a + 1
						else:
							continue
					
		
#for +-

for i in range(nat):	
	for j in range(1,4):
		for k in range(nat):	
			for l in range (1,4):	
				if k != i:
					if os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-.out') == False:
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
						a = a + 1
					else:
						continue				
								
				elif k == i and l != j:
					if os.path.isfile(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-.out') == False:
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
						a = a + 1
					else:
						continue									
						
if a != 0:
	print('Not all energies calculated!')
	exit()
	


#create dictionary to store energies

energy = {}
en_null = 0

			
#read in energy for the molecule

with open(inp+'-4c_'+inp+'.out') as f:
        for line in f:
                if line[:45] == '   Total energy                             :':
                        energy['0'] = float(line[45:])
                        break
if not '0' in energy:
        subprocess.call(['qsub run-'+inp+'-4c-dir'], shell = True)
        en_null = en_null + 1

#read in energies for all singular displacements

for i in range(nat):

        for j in range(1,4):
                with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+.out') as f:
                        for line in f:
                                if line[:45] == '   Total energy                             :':
                                        energy[str(i)+str(j)+'+'] = float(line[45:])
                                        break
                if not str(i)+str(j)+'+' in energy:
                        subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+_4c-dir'], shell = True)
                        en_null = en_null + 1

                with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-.out') as f:
                        for line in f:
                                if line[:45] == '   Total energy                             :':
                                        energy[str(i)+str(j)+'-'] = float(line[45:])
                                        break
                if not str(i)+str(j)+'-' in energy:
                        subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-_4c-dir'], shell = True)
                        en_null = en_null + 1

#for ++ and --
for i in range(nat):
        for j in range(1,4):
                for k in range(nat):
                        for l in range (1,4):
                                if k > i:
                                        with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+.out') as f:
                                                for line in f:
                                                        if line[:45] == '   Total energy                             :':
                                                                energy[str(i)+str(j)+'+'+str(k)+str(l)+'+'] = float(line[45:])
                                                                break
                                        if not str(i)+str(j)+'+'+str(k)+str(l)+'+' in energy:
                                                subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
                                                en_null = en_null + 1

                                        with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-.out') as f:
                                                for line in f:
                                                        if line[:45] == '   Total energy                             :':
                                                                energy[str(i)+str(j)+'-'+str(k)+str(l)+'-'] = float(line[45:])
                                                                break
                                        if not str(i)+str(j)+'-'+str(k)+str(l)+'-' in energy:
                                                subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
                                                en_null = en_null + 1

                                elif k == i:
                                        if l > j:
                                                with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+.out') as f:
                                                        for line in f:
                                                                if line[:45] == '   Total energy                             :':
                                                                        energy[str(i)+str(j)+'+'+str(k)+str(l)+'+'] = float(line[45:])
                                                                        break
                                                if not str(i)+str(j)+'+'+str(k)+str(l)+'+' in energy:
                                                        subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
                                                        en_null = en_null + 1

                                                with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-.out') as f:
                                                        for line in f:
                                                                if line[:45] == '   Total energy                             :':
                                                                        energy[str(i)+str(j)+'-'+str(k)+str(l)+'-'] = float(line[45:])
                                                                        break
                                                if not str(i)+str(j)+'-'+str(k)+str(l)+'-' in energy:
                                                        subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
                                                        en_null = en_null + 1

#for +-

for i in range(nat):
        for j in range(1,4):
                for k in range(nat):
                        for l in range (1,4):
                                if k != i:
                                        with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-.out') as f:
                                                for line in f:
                                                        if line[:45] == '   Total energy                             :':
                                                                energy[str(i)+str(j)+'+'+str(k)+str(l)+'-'] = float(line[45:])
                                                                break
                                        if not str(i)+str(j)+'+'+str(k)+str(l)+'-' in energy:
                                                subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
                                                en_null = en_null + 1

                                elif k == i and l != j:
                                        with open(inp+'-4c_'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-.out') as f:
                                                for line in f:
                                                        if line[:45] == '   Total energy                             :':
                                                                energy[str(i)+str(j)+'+'+str(k)+str(l)+'-'] = float(line[45:])
                                                                break
                                        if not str(i)+str(j)+'+'+str(k)+str(l)+'-' in energy:
                                                subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
                                                en_null = en_null + 1


if en_null != 0:
        print('Not all energies calculated!')
        exit()
	
	
#create empty hessian matrix

hessian = np.zeros((3*nat,3*nat))

	
#calculate diagonal of the hessian

for i in range(nat):
	for j in range(1,4):
		var = (energy[str(i)+str(j)+'+'] + energy[str(i)+str(j)+'-'] - (2 * energy['0'])) / (float(dis) ** 2)
		hessian[i*3+(j-1), i*3+(j-1)] = var

			
#calculate non-diagonal part of the hessian

for i in range(nat):
	for j in range(1,4):
		for k in range(nat):
			for l in range(1,4):
				if k > i or (k == i and l > j):
			
					var = (energy[str(i)+str(j)+'+'+str(k)+str(l)+'+'] + energy[str(i)+str(j)+'-'+str(k)+str(l)+'-'] - energy[str(i)+str(j)+'+'+str(k)+str(l)+'-'] - energy[str(k)+str(l)+'+'+str(i)+str(j)+'-']) / (4 * (float(dis) ** 2))
					
#over diagonal
					hessian[i*3+j-1, k*3+l-1] = var
					
#under diagonal
					hessian[k*3+l-1, i*3+j-1] = var
					
					
#hessian in a.u.
bohr_to_angstrom =  0.529177211
hessian_bohr = hessian * (bohr_to_angstrom ** 2)

cart = ['x', 'y', 'z']
 
#write hessian into a file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Hessian (a.u.)   ***\n\n')
	f.write('        ')
	for j in range (3*nat):
		f.write((atoms[int(j/3)]+cart[j%3]+'   ').rjust(12))
	f.write('\n')
	ll = 0
	for line in hessian_bohr:
		f.write((atoms[int(ll/3)]+cart[ll%3]+' ').rjust(6))
		for row in line:
			f.write(str(round(Decimal(row),4)).rjust(12))
		ll = ll+1
		f.write('\n')
	f.write('\n\n')


#calculate mass-weighted hessian

atomic_masses = dict([('H ', 1.008), ('He', 4.0026), ('Li', 6.94), ('Be', 9.0122), ('B ', 10.81), ('C ', 12.011), 
('N ', 14.007), ('O ', 15.999), ('F', 18.998), ('Ne', 20.180), ('Na', 22.990), ('Mg', 24.305), ('Al', 26.982), ('Si', 28.085),
('P ', 30.974), ('S ', 32.06), ('Cl', 35.45), ('Ar', 39.948), ('K ', 39.098), ('Ca', 40.078), ('Sc', 44.956), ('Ti', 47.867),
('V ', 50.942), ('Cr', 51.996), ('Mn', 54.938), ('Fe', 55.845), ('Co', 58.993), ('Ni', 58.693), ('Cu', 63.546), ('Zn', 65.38),
('Ga', 69.723), ('Ge', 72.630), ('As', 74.922), ('Se', 78.971), ('Br', 79.904), ('Kr', 83.798), ('Rb', 85.468), ('Sr', 87.62),
('Y ', 88.906), ('Zr', 91.224), ('Nb', 92.906), ('Mo', 95.95), ('Tc', 98), ('Ru', 101.07), ('Rh', 102.91), ('Pd', 106.42),
('Ag', 107.87), ('Cd', 112.41), ('In', 114.82), ('Sn', 118.71), ('Sb', 121.76), ('Te', 127.60), ('I ', 126.90), ('Xe', 131.29),
('Cs', 132.91), ('Ba', 137.33), ('La', 138.91), ('Ce', 140.12), ('Pr', 140.91), ('Nd', 144.24), ('Pm', 145), ('Sm', 150.36),
('Eu', 151.96), ('Gd', 157.25), ('Tb', 158.93), ('Dy', 162.50), ('Ho', 164.93), ('Er', 167.26), ('Tm', 168.93), ('Yb', 173.05),
('Lu', 174.97), ('Hf', 178.49), ('Ta', 180.95), ('W ', 183.84), ('Re', 186.21), ('Os', 190.23), ('Ir', 192.22), ('Pt', 195.08),
('Au', 196.97), ('Hg', 200.59), ('Tl', 204.38), ('Pb', 207.2), ('Bi', 208.98), ('Po', 209.0), ('At', 210), ('Rn', 222),
('Fr', 223), ('Ra', 226), ('Ac', 227), ('Th', 232.04), ('Pa', 231.04), ('U ', 238.03), ('Np', 237), ('Pu', 244),
('Am', 243), ('Cm', 247), ('Bk', 247), ('Cf', 251), ('Es', 252)])



at_mass_sqrt = [0 for i in range(3*nat)]
for i in range(nat):
	for j in range(3):
		at_mass_sqrt[i*3 + j] = isotopes[atoms[i].rstrip()] ** -0.5
	
mass_matrix = np.diag(at_mass_sqrt)

hessian_mass = np.dot(mass_matrix, np.dot(hessian_bohr, mass_matrix))


#diagonilise mass-weighted hessian

(eig_val, eig_vec) = linalg.eig(hessian_mass)
hessian_diag = np.diag(eig_val)


#calculate vibrational frequencies

vib_freq = [0 for i in range(3*nat) ]

amu_to_au = 1822.88848
hartree_to_cm_rev = 219474.62934

for i in range(3*nat):
	if eig_val[i] > 0.00050:
		vib_freq[i] = (eig_val[i] / amu_to_au ) ** 0.5
		

vib_freq_hartree = []
num_freq = 0
for i in range(3*nat):
	if vib_freq[i] != 0:
		vib_freq_hartree.append(vib_freq[i])
		num_freq = num_freq + 1
		
vib_freq_cm = [hartree_to_cm_rev * vib_freq_hartree[i] for i in range(num_freq)]


#write number of modes into a file
with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Number of vibrational modes   ***\n')
	f.write(str(num_freq)+'    \n')
	f.write('\n\n')

#write vibrational frequancies into a file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Vibrational frequencies   ***\n\n')
	f.write(' mode          cm-1          hartree\n')
	f.write('--------------------------------------\n')
	for i in range(num_freq):
		f.write(str(i+1).rjust(4)+str(round(Decimal(vib_freq_cm[i]),2)).rjust(16)+str(round(Decimal(vib_freq_hartree[i]),7)).rjust(17)+'\n')
	f.write('\n\n')


#get normal coordinates

norm_coord = np.zeros((num_freq,3*nat))
num_freq_2 = 0
for i in range(3*nat):
	if vib_freq[i] != 0:
		norm_coord[num_freq_2] = np.transpose(np.dot(mass_matrix, eig_vec))[i]
		num_freq_2 = num_freq_2 + 1
		

#write normal coordinates into a file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Normal Coordinates (bohrs*amu**(1/2))   ***\n\n')
	f.write('mode   '.rjust(17))
	for j in range (3*nat):
		f.write((atoms[int(j/3)]+cart[j%3]).rjust(10))
	f.write('\n')
	for i in range(num_freq):
		f.write(str(i+1).rjust(3)+str(round(Decimal(vib_freq_cm[i]),2)).rjust(10)+' cm-1')
		for j in range(3*nat):
			f.write(str(round(Decimal(norm_coord[i,j]),4)).rjust(10))
		f.write('\n')
	f.write('\n')
	
	
#write mass matrix into a file

#with open(inp+'-rel.out', 'a') as f:
#	f.write('***   Mass matrix   ***\n\n')
#	for i in range(nat):
#		for k in range(0,3):
#			for j in range(nat):
#				for l in range(0,3):
#					f.write(str(mass_matrix[i*3+k, j*3+l])+'   ')
#			f.write('\n')
		

#write energies into a file

#with open(inp+'-rel.out', 'a') as f:
#		f.write('***  Energies   ***\n\n')
#		
#		#for the molecule 
#		f.write('0       '+str(energy['0'])+'\n')
#		
#		#for singular displacements
#		for i in range(nat):
#			for j in range(1,4):
#				f.write(str(i)+str(j)+'+     '+str(energy[str(i)+str(j)+'+'])+'\n')
#				f.write(str(i)+str(j)+'-     '+str(energy[str(i)+str(j)+'-'])+'\n')
		
		#for double displacements
		#for ++ and --
#		for i in range(nat):	
#			for j in range(1,4):
#				for k in range(nat):
#					for l in range (1,4):
#						if k > i:
#							f.write(str(i)+str(j)+'+'+str(k)+str(l)+'+  '+str(energy[str(i)+str(j)+'+'+str(k)+str(l)+'+'])+'\n')
#							f.write(str(i)+str(j)+'-'+str(k)+str(l)+'-  '+str(energy[str(i)+str(j)+'-'+str(k)+str(l)+'-'])+'\n')
#						
#						elif k == i:
#							if l > j:
#								f.write(str(i)+str(j)+'+'+str(k)+str(l)+'+  '+str(energy[str(i)+str(j)+'+'+str(k)+str(l)+'+'])+'\n')
#								f.write(str(i)+str(j)+'-'+str(k)+str(l)+'-  '+str(energy[str(i)+str(j)+'-'+str(k)+str(l)+'-'])+'\n')
#		#for +-
#		for i in range(nat):	
#			for j in range(1,4):
#				for k in range(nat):	
#					for l in range (1,4):	
#						if k != i:
#							f.write(str(i)+str(j)+'+'+str(k)+str(l)+'-  '+str(energy[str(i)+str(j)+'+'+str(k)+str(l)+'-'])+'\n')
#				
#						elif k == i and l != j:
#							f.write(str(i)+str(j)+'+'+str(k)+str(l)+'-  '+str(energy[str(i)+str(j)+'+'+str(k)+str(l)+'-'])+'\n')
							

#		f.write('\n')


# calculate transformation matrix

mass_matrix_amu = mass_matrix * (amu_to_au ** (-0.5))

#print(mass_matrix_amu)

norm_coord_A = norm_coord * bohr_to_angstrom

#norm_coord_inv = np.linalg.pinv(norm_coord_A)

# create matrix with atoms xyz coordinates

atoms_xyz = np.delete(np.array(matrix), 0,1)

for row in atoms_xyz:
	for j in range(3):
		row[j] = float(row[j])
		
#print(norm_coord)
# create displacement vectors for derivative calculations

for i in range (num_freq):
	
# transform the vector into not mass-weighted xyz coordinates

	vect = float(dis2) * norm_coord_A[i]
	vect_xyz = np.dot(mass_matrix_amu, vect)
	dis_xyz = vect_xyz.reshape((nat,3))
#	print(dis_xyz)
	
	
# create singular displaced geometry and run files

	atoms_dis_p = atoms_xyz + dis_xyz
	atoms_dis_m = atoms_xyz - dis_xyz
	

	create_mol(inp+'-J_'+str(i)+'+', atoms_dis_p)
	create_run1(inp+'-J_'+str(i)+'+')
	
	
	create_mol(inp+'-J_'+str(i)+'-', atoms_dis_m)
	create_run1(inp+'-J_'+str(i)+'-')
	
# run the scripts
	
	subprocess.call(['chmod u+x run-'+inp+'-J_'+str(i)+'+'+'_4c-dir'], shell = True)
	subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'+'+'_4c-dir'], shell = True)
	
	subprocess.call(['chmod u+x run-'+inp+'-J_'+str(i)+'-'+'_4c-dir'], shell = True)
	subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'-'+'_4c-dir'], shell = True)
	


# create singular displaced geometry and run files for third derivatives of the energy

	atoms_dis_pp = atoms_xyz + 2 * dis_xyz
	atoms_dis_mm = atoms_xyz - 2 * dis_xyz
	

	create_mol(inp+'-E_'+str(i)+'++', atoms_dis_pp)
	create_run2(inp+'-E_'+str(i)+'++')
	
	

	create_mol(inp+'-E_'+str(i)+'--', atoms_dis_mm)
	create_run2(inp+'-E_'+str(i)+'--')
	
# run the scripts
	
	subprocess.call(['chmod u+x run-'+inp+'-E_'+str(i)+'++'+'_4c-dir'], shell = True)
	subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'++'+'_4c-dir'], shell = True)
	
	subprocess.call(['chmod u+x run-'+inp+'-E_'+str(i)+'--'+'_4c-dir'], shell = True)
	subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'--'+'_4c-dir'], shell = True)
	

create_run1(inp+'-J')

subprocess.call(['cp '+inp+'.mol '+inp+'-J.mol'], shell = True)
subprocess.call(['chmod u+x run-'+inp+'-J_4c-dir'], shell = True)
subprocess.call(['qsub run-'+inp+'-J_4c-dir'], shell = True)



#create dirac .mol and run files for all the double displacements

#for ++ and --

for i in range(num_freq):
	
	for j in range(num_freq):
	
		if j > i:
	
			vect = float(dis2) * norm_coord_A[i] + float(dis2) * norm_coord_A[j]
			vect_xyz = np.dot(mass_matrix_amu, vect)
			dis_xyz = vect_xyz.reshape((nat,3))
	

			atoms_dis_p = atoms_xyz + dis_xyz
			atoms_dis_m = atoms_xyz - dis_xyz
	

			create_mol(inp+'-E_'+str(i)+'+'+str(j)+'+', atoms_dis_p)
			create_run2(inp+'-E_'+str(i)+'+'+str(j)+'+')
	
	
			create_mol(inp+'-E_'+str(i)+'-'+str(j)+'-', atoms_dis_m)
			create_run2(inp+'-E_'+str(i)+'-'+str(j)+'-')
	
# run the scripts
	
			subprocess.call(['chmod u+x run-'+inp+'-E_'+str(i)+'+'+str(j)+'+_4c-dir'], shell = True)
			subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'+_4c-dir'], shell = True)
	
			subprocess.call(['chmod u+x run-'+inp+'-E_'+str(i)+'-'+str(j)+'-'+'_4c-dir'], shell = True)
			subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'-'+str(j)+'-'+'_4c-dir'], shell = True)

#for +-

for i in range(num_freq):
	
	for j in range(num_freq):
		
		if j != i:
		
			vect = float(dis2) * norm_coord_A[i] - float(dis2) * norm_coord_A[j]
			vect_xyz = np.dot(mass_matrix_amu, vect)
			dis_xyz = vect_xyz.reshape((nat,3))
	

			atoms_dis = atoms_xyz + dis_xyz

	
			create_mol(inp+'-E_'+str(i)+'+'+str(j)+'-', atoms_dis)
			create_run2(inp+'-E_'+str(i)+'+'+str(j)+'-')
	
	
# run the scripts
	
			subprocess.call(['chmod u+x run-'+inp+'-E_'+str(i)+'+'+str(j)+'-_4c-dir'], shell = True)
			subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'-_4c-dir'], shell = True)
			
print(mass_matrix)
			

	
	
		
	


	
	
	

	
	


