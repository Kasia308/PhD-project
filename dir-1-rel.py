from decimal import getcontext, Decimal
getcontext().prec = 7

import sys
inp = sys.argv[1]
#dis = Decimal(sys.argv[2])

#create matrix for atomic positions from .xyz file

with open(inp+'.xyz') as f:
  list = f.readlines()
  
matrix = [[pos[0:2].replace(" ",""),Decimal(pos[10:20]),Decimal(pos[26:36]),Decimal(pos[42:52])] for pos in list]

nat = len(matrix)

atoms = [(row[0]) for row in matrix]

atom_count = dict()

for i in atoms:
	atom_count[i] = atom_count.get(i, 0) + 1


atomic_numbers = dict([('H', 1), ('He', 2), ('Li', 3), ('Be', 4), ('B', 5), ('C', 6), 
('N', 7), ('O', 8), ('F', 9), ('Ne', 10), ('Na', 11), ('Mg', 12), ('Al', 13), ('Si', 14),
('P', 15), ('S', 16), ('Cl', 17), ('Ar', 18), ('K', 19), ('Ca', 20), ('Sc', 21), ('Ti', 22),
('V', 23), ('Cr', 24), ('Mn', 25), ('Fe', 26), ('Co', 27), ('Ni', 28), ('Cu', 29), ('Zn', 30),
('Ga', 31), ('Ge', 32), ('As', 33), ('Se', 34), ('Br', 35), ('Kr', 36), ('Rb', 37), ('Sr', 38),
('Y', 39), ('Zr', 40), ('Nb', 41), ('Mo', 42), ('Tc', 43), ('Ru', 44), ('Rh', 45), ('Pd', 46),
('Ag', 47), ('Cd', 48), ('In', 49), ('Sn', 50), ('Sb', 51), ('Te', 52), ('I', 53), ('Xe', 54),
('Cs', 55), ('Ba', 56), ('La', 57), ('Ce', 58), ('Pr', 59), ('Nd', 60), ('Pm', 61), ('Sm', 62),
('Eu', 63), ('Gd', 64), ('Tb', 65), ('Dy', 66), ('Ho', 67), ('Er', 68), ('Tm', 69), ('Yb', 70),
('Lu', 71), ('Hf', 72), ('Ta', 73), ('W ', 74), ('Re', 75), ('Os', 76), ('Ir', 77), ('Pt', 78),
('Au', 79), ('Hg', 80), ('Tl', 81), ('Pb', 82), ('Bi', 83), ('Po', 84), ('At', 85), ('Rn', 86),
('Fr', 87), ('Ra', 88), ('Ac', 89), ('Th', 90), ('Pa', 91), ('U ', 92), ('Np', 93), ('Pu', 94),
('Am', 95), ('Cm', 96), ('Bk', 97), ('Cf', 98), ('Es', 99)])


#read in method, basis sets, functional and displacements from .inp file

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
			basis[line[:2].replace(" ","")] = line.replace("\n","")[6:]
	
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

		

#define function to create dirac .mol files

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
					f.write(str(row[0]).ljust(2)+str(round(row2[1],6)).rjust(17)+str(round(row2[2],6)).rjust(15)+str(round(row2[3],6) ).rjust(15)+'\n')
			f.write('LARGE BASIS '+ basis[el]+'\n')
			j = j + atom_count[el]
		f.write('FINISH')

def create_run(name):

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
		f.write('export inp='+inp+'-4c\n')
		f.write('export mol='+name+'\n')
		f.write('export out=$inp\_$mol.out\n')
		f.write('\n')
		f.write('$pam --noarch --scratch=$scratch --mw=2100 --mol=$mol.mol --inp=$inp.inp --put="DFCOEF" --get="AOPROPER" --outcmo\n')
		f.write('\n')
		f.write('exit 0\n')
							

#create dirac .inp file

if method == 'rel':
	
	with open(inp+'-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
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
		f.write('1.0E-6  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('200\n')
		f.write('*END OF\n')

elif method == 'nrel':

	with open(inp+'-4c.inp', 'w+') as f:
		f.write('**DIRAC\n')
		f.write('.WAVE F\n')
		f.write('**INTEGRALS\n')
		f.write('.NUCMOD\n')
		f.write('2\n')
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
		f.write('1.0E-6  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('200\n')
		f.write('*END OF\n')

elif method == 'noso':

	with open(inp+'-4c.inp', 'w+') as f:
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
		f.write('1.0E-6  1.0D-6\n')
		f.write('.MAXITR\n')
		f.write('200\n')
		f.write('*END OF\n')


	
#create dirac .mol file for the molecule			
		
create_mol(inp, matrix)

#create dirac .mol files for all the singular displacements	

for i in range(nat):
	
	for j in range(1,4):
		newmatrix = matrix
		line = matrix[i]
		line[j] = line[j] + dis
		newmatrix[i] = line
		
		create_mol(inp+'_'+str(i)+str(j)+'+', newmatrix)

		line[j] = line[j] - dis
		
for i in range(nat):
	
	for j in range(1,4):
		newmatrix = matrix
		line = matrix[i]
		line[j] = line[j] - dis
		newmatrix[i] = line
		
		create_mol(inp+'_'+str(i)+str(j)+'-', newmatrix)

		line[j] = line[j] + dis
		
				
#create dirac .mol files for all the double displacements

#for + +

for i in range(nat):
	
	for j in range(1,4):
		newmatrix = matrix
		line = matrix[i]
		line[j] = line[j] + dis
		newmatrix[i] = line
		
		for k in range(nat):
			
			for l in range (1,4):
			
				if k > i:
						
					line2 = matrix[k]
					line2[l] = line2[l] + dis
					newmatrix[k] = line2
	
					create_mol(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+', newmatrix)
					
					line2[l] = line2[l] - dis
					
				elif k == i:
					if l > j:
						line[l] = line[l] + dis
						newmatrix[k] = line
						
						create_mol(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+', newmatrix)
					
						line[l] = line[l] - dis						

		line[j] = line[j] - dis
		
#for + -

for i in range(nat):
	
	for j in range(1,4):
		newmatrix = matrix
		line = matrix[i]
		line[j] = line[j] + dis
		newmatrix[i] = line
		
		for k in range(nat):
			
			for l in range (1,4):
			
				if k != i:
						
					line2 = matrix[k]
					line2[l] = line2[l] - dis
					newmatrix[k] = line2
					
					create_mol(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-', newmatrix)
					
					line2[l] = line2[l] + dis
					
				elif k == i and l != j:
					line[l] = line[l] - dis
					newmatrix[k] = line
						
					create_mol(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-', newmatrix)
					
					line[l] = line[l] + dis						

		line[j] = line[j] - dis

		
#for - - 
		
for i in range(nat):
	
	for j in range(1,4):
		newmatrix = matrix
		line = matrix[i]
		line[j] = line[j] - dis
		newmatrix[i] = line
		
		for k in range(nat):
			
			for l in range (1,4):
			
				if k > i:
						
					line2 = matrix[k]
					line2[l] = line2[l] - dis
					newmatrix[k] = line2
					
					create_mol(inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-', newmatrix)
	
					line2[l] = line2[l] + dis
					
				elif k == i:
					if l > j:
						line[l] = line[l] - dis
						newmatrix[k] = line
						
						create_mol(inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-', newmatrix)

						line[l] = line[l] + dis						

		line[j] = line[j] + dis
		
#create 1c and 4c run scripts for the molecule

create_run(inp)

#create run scripts for all the singular displacements	

for i in range(nat):
	
	for j in range(1,4):
		create_run(inp+'_'+str(i)+str(j)+'+')
		
for i in range(nat):

	for j in range(1,4):
		create_run(inp+'_'+str(i)+str(j)+'-')

#create run scripts for all the double displacements	

#for + +

for i in range(nat):
	for j in range(1,4):
		for k in range(nat):
			for l in range (1,4):
				if k > i:
					create_run(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+')	
				elif k == i:
					if l > j:
						create_run(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+')

		
#for + -

for i in range(nat):
	for j in range(1,4):
		for k in range(nat):		
			for l in range (1,4):		
				if k != i:
					create_run(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-')	
				elif k == i and l != j:	
					create_run(inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-')
					

		
#for - - 
		
for i in range(nat):
	for j in range(1,4):
		for k in range(nat):
			for l in range (1,4):
				if k > i:
					create_run(inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-')	
				elif k == i:
					if l > j:
						create_run(inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-')



#start 4c respect jobs

import subprocess

#for the molecule

subprocess.call(['chmod u+x run-'+inp+'_4c-dir'], shell = True)
subprocess.call(['qsub run-'+inp+'_4c-dir'], shell = True)

#for all singular displacements

for i in range(nat):
	
	for j in range(1,4):
		subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'+_4c-dir'], shell = True)
		subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+_4c-dir'], shell = True)
		
for i in range(nat):
	for j in range(1,4):
		subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'-_4c-dir'], shell = True)
		subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-_4c-dir'], shell = True)

		

#for all the double displacements	

#for + +

for i in range(nat):
	for j in range(1,4):
		for k in range(nat):
			for l in range (1,4):
				if k > i:
					subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
					subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
				elif k == i:
					if l > j:
						subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'+_4c-dir'], shell = True)

		
#for + -

for i in range(nat):
	for j in range(1,4):
		for k in range(nat):		
			for l in range (1,4):		
				if k != i:
					subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
					subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
				elif k == i and l != j:	
						subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'+'+str(k)+str(l)+'-_4c-dir'], shell = True)
					

		
#for - - 
		
for i in range(nat):
	for j in range(1,4):
		for k in range(nat):
			for l in range (1,4):
				if k > i:
					subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
					subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
				elif k == i:
					if l > j:
						subprocess.call(['chmod u+x run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
						subprocess.call(['qsub run-'+inp+'_'+str(i)+str(j)+'-'+str(k)+str(l)+'-_4c-dir'], shell = True)
