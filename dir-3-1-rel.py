from decimal import *
getcontext().prec = 7

import sys
inp = sys.argv[1]
#dis = Decimal(sys.argv[2])
#dis2 = Decimal(sys.argv[3])

import numpy as np
from numpy import linalg

import datetime  

import os.path

import subprocess

import math

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

		
#read in vibrational modes from -rel.out file

with open(inp+'-'+method+'.out') as f:

	for lnum, line in enumerate(f, 0):
			if '***   Number of vibrational modes   ***' in line:
				break
	for line in f:
		num_freq = int(line[0:2])
		break
		
	for lnum, line in enumerate(f, 0):
			if '--------------------------------------' in line:
				break
	
	vib_freq = [ 0 for i in range(num_freq)]
	vib_freq_hartree = [ 0 for i in range(num_freq)]
	
	a = int(0)
	for line in f:
		if a < num_freq:
			vib_freq[a] = float(line[11:22])
			vib_freq_hartree[a] = float(line[25:])
			a = a+1
		else:
			break
			
			
#check if all energy and sscc files exist

a = 0

#for the molecule

if os.path.isfile(inp+'-J-4c_'+inp+'-J.out') == False:
	subprocess.call(['qsub run-'+inp+'-J_4c-dir'], shell = True)
	a = a + 1
			
#all singular displacements

for i in range(num_freq):
	
	if os.path.isfile(inp+'-J-4c_'+inp+'-J_'+str(i)+'+.out') == False:
		subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'+'+'_4c-dir'], shell = True)
		a = a + 1
		
	elif os.path.isfile(inp+'-J-4c_'+inp+'-J_'+str(i)+'-.out') == False:
		subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'-'+'_4c-dir'], shell = True)
		a = a + 1
			
	elif os.path.isfile(inp+'-E-4c_'+inp+'-E_'+str(i)+'++.out') == False:
		subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'++'+'_4c-dir'], shell = True)
		a = a + 1
		
	elif os.path.isfile(inp+'-E-4c_'+inp+'-E_'+str(i)+'--.out') == False:
		subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'--'+'_4c-dir'], shell = True)
		a = a + 1
	else:
		continue

#all double displacements

#for ++ and --

for i in range(num_freq):
	
	for j in range(num_freq):
	
		if j > i:
		
			if os.path.isfile(inp+'-E-4c_'+inp+'-E_'+str(i)+'+'+str(j)+'+.out') == False:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'+_4c-dir'], shell = True)
				a = a + 1
		
			elif os.path.isfile(inp+'-E-4c_'+inp+'-E_'+str(i)+'-'+str(j)+'-.out') == False:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'-'+str(j)+'-_4c-dir'], shell = True)
				a = a + 1
			
			else:
				continue

#for +-

for i in range(num_freq):
	
	for j in range(num_freq):
		
		if j != i:
		
			if os.path.isfile(inp+'-E-4c_'+inp+'-E_'+str(i)+'+'+str(j)+'-.out') == False:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'-_4c-dir'], shell = True)
				a = a + 1
								
if a != 0:
	print('Not all energies and ssccs calculated!')
	exit()			


		
#read in isotopes from optimized geometry sscc calculations

#isotopes = {}

#with open(inp+'-J-4c_'+inp+'-J.out') as f:
#	for line in f:
#		if line[22:23] == '1' and line[0:1]== '@':
#			isotopes[line[2:4]] = float(line[25:36])
			
#get the number of sscc's to read in

sscc_num = math.factorial(nat)/math.factorial(2)/math.factorial(nat-2)

#read in the sscc's

sscc = [[[] for i in range(int(sscc_num+1))] for i in range (2*num_freq+1)]

#for the molecule

sscc1 = sscc[0]

sscc1[0] = '0'

iso_num = {}

with open(inp+'-J-4c_'+inp+'-J.out') as f:
	for lnum, line in enumerate(f, 0):
			if '@ Nuclear species:' in line:
				break
	
	for line in f:
		for element in isotopes:
			if line[2:4].rstrip() == element and float(line[26:32].replace(' ', '')) == isotopes[element]:
				iso_num[element] = int(line[22:23])
		
		if '                        ! Final spin-spin-couplings J (Hz) !' in line:
			break
			

with open(inp+'-J-4c_'+inp+'-J.out') as f:
	for lnum, line in enumerate(f, 0):
		if '                        ! Final spin-spin-couplings J (Hz) !' in line:
			break
				
	for line in f:	
		if '+----------------------------------+' in line:
			break
	for line in f:	
		if '' in line:
			break
	for line in f:	
		if 'At1  #  iso  : At2  #  iso  :   Isotropic     Anisotropic   Asymmetry     S parameter   A parameter' in line:
			break
	for line in f:	
		if '----------------------------------------------------------------------------------------------------' in line:
			break
	
	nh = 1
	for line in f:
		
		if not '----------------------------------------------------------------------------------------------------' in line:
			if line[10:11] == str(iso_num[line[0:2].rstrip()]) and line[25:26] == str(iso_num[line[15:17].rstrip()]):
				sscc1[nh] = [str(line[0:2]), str(line[15:17]), float(line[33:44])]
				nh = nh + 1
		
		if '----------------------------------------------------------------------------------------------------' in line:
			break

# for the single displacements

ap = 1

for i in range(num_freq):
	sscc1 = sscc[ap]
	sscc1[0] = str(i)+'+'
	with open(inp+'-J-4c_'+inp+'-J_'+str(i)+'+.out') as f:
		for lnum, line in enumerate(f, 0):
			if '                        ! Final spin-spin-couplings J (Hz) !' in line:
				break
				
		for line in f:	
			if '+----------------------------------+' in line:
				break
		for line in f:	
			if '' in line:
				break
		for line in f:	
			if 'At1  #  iso  : At2  #  iso  :   Isotropic     Anisotropic   Asymmetry     S parameter   A parameter' in line:
				break
		for line in f:	
			if '----------------------------------------------------------------------------------------------------' in line:
				break
	
		nh = 1
		for line in f:
		
			if not '----------------------------------------------------------------------------------------------------' in line:
				if line[10:11] == str(iso_num[line[0:2].rstrip()]) and line[25:26] == str(iso_num[line[15:17].rstrip()]):
					sscc1[nh] = [str(line[0:2]), str(line[15:17]), float(line[33:44])]
					nh = nh + 1
		
			if '----------------------------------------------------------------------------------------------------' in line:
				break
	ap = ap + 1
		
	sscc1 = sscc[ap]
	sscc1[0] = str(i)+'-'
	with open(inp+'-J-4c_'+inp+'-J_'+str(i)+'-.out') as f:
		for lnum, line in enumerate(f, 0):
			if '                        ! Final spin-spin-couplings J (Hz) !' in line:
				break
				
		for line in f:	
			if '+----------------------------------+' in line:
				break
		for line in f:	
			if '' in line:
				break
		for line in f:	
			if 'At1  #  iso  : At2  #  iso  :   Isotropic     Anisotropic   Asymmetry     S parameter   A parameter' in line:
				break
		for line in f:	
			if '----------------------------------------------------------------------------------------------------' in line:
				break
	
		nh = 1
		for line in f:
			
			if not '----------------------------------------------------------------------------------------------------' in line:
				if line[10:11] == str(iso_num[line[0:2].rstrip()]) and line[25:26] == str(iso_num[line[15:17].rstrip()]):
					sscc1[nh] = [str(line[0:2]), str(line[15:17]), float(line[33:44])]
					nh = nh + 1
		
			if '----------------------------------------------------------------------------------------------------' in line:
				break
	ap = ap + 1


#write SSCC's for the molecule into the file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Spin-spin-couplings J (Hz) for '+inp+'   ***\n')
	
	sscc_pr = sscc[0]
	
	for n in range(1,int(sscc_num)+1):
		sscc_pr2 = sscc_pr[n]
		f.write('J('+sscc_pr2[0].rstrip()+'-'+sscc_pr2[1].rstrip()+'):      '+str(sscc_pr2[2]).rjust(12)+'\n')
	f.write('\n')
	
	
#calculate first-order derivatives of SSCC's

first_derivative = [[[] for i in range(int(num_freq+2))] for i in range(int(sscc_num))]
second_derivative = [[[] for i in range(int(num_freq+2))] for i in range(int(sscc_num))]



for n in range(int(sscc_num)):
	der1 = first_derivative[n]
	der2 = second_derivative[n]
	sscc_pr1 = sscc_pr[n+1]
	der1[0] = sscc_pr1[0]
	der1[1] = sscc_pr1[1]
	der2[0] = sscc_pr1[0]
	der2[1] = sscc_pr1[1]
	
	i = 1
	for j in range(num_freq):
		sscc_p = sscc[i]
		sscc_pp = sscc_p[n+1]
		sscc_m = sscc[i+1]
		sscc_mm = sscc_m[n+1]
		
		der1[j+2] = (sscc_pp[2] - sscc_mm[2])/(2 * float(dis2))
		der2[j+2] = (sscc_pp[2] + sscc_mm[2] - 2 * sscc_pr1[2])/(float(dis2) ** 2)
		
		i = i + 2

#write SSCCs' derivatives into the file		

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Spin-spin-couplings J first derivatives with repect to normal coordinates [Hz / (bohr * amu ** 1/2)] for '+inp+'   ***\n')
	
	
	for n in range(0,int(sscc_num)):
		f_der_pr = first_derivative[n]
		f.write('dJ('+f_der_pr[0].rstrip()+'-'+f_der_pr[1].rstrip()+')/dQ:\n')
		f.write('[')
		
		for i in range(num_freq):
			f.write(str(round(f_der_pr[i+2], 5)).rjust(15))
		
		f.write(']\n')
	f.write('\n')
	
	
	f.write('***   Spin-spin-couplings J second derivatives with repect to normal coordinates [Hz / (bohr^2 * amu)] for '+inp+'   ***\n')
	
	
	for n in range(0,int(sscc_num)):
		s_der_pr = second_derivative[n]
		f.write('d^2J('+s_der_pr[0].rstrip()+'-'+s_der_pr[1].rstrip()+')/dQ^2:\n')
		f.write('[')
		
		for i in range(num_freq):
			f.write(str(round(s_der_pr[i+2],5)).rjust(20))
		
		f.write(']\n')
	f.write('\n')


#create dictionary to store energies

en3 = {}

en_null = 0

#read in energy for the molecule

with open(inp+'-J-4c_'+inp+'-J.out') as f:
    for line in f:
        if line[:45] == '   Total energy                             :':
            en3['0'] = float(line[45:])
            break
if not '0' in en3:
    subprocess.call(['qsub run-'+inp+'-J_4c-dir'], shell = True)
    en_null = en_null + 1
		
#read in energy for all singular displacements

for i in range(num_freq):

#+
	with open(inp+'-J-4c_'+inp+'-J_'+str(i)+'+.out') as f:
		for line in f:
			if line[:45] == '   Total energy                             :':
				en3[str(i)+'+'] = float(line[45:])
				break
	if not str(i)+'+' in en3:
		subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'+'+'_4c-dir'], shell = True)
		en_null = en_null + 1
	
#++
	with open(inp+'-E-4c_'+inp+'-E_'+str(i)+'++.out') as f:
		for line in f:
			if line[:45] == '   Total energy                             :':
				en3[str(i)+'++'] = float(line[45:])
				break
	if not str(i)+'++' in en3:
		subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'++'+'_4c-dir'], shell = True)
		en_null = en_null + 1
			
#-
	with open(inp+'-J-4c_'+inp+'-J_'+str(i)+'-.out') as f:
		for line in f:
			if line[:45] == '   Total energy                             :':
				en3[str(i)+'-'] = float(line[45:])
				break
	if not str(i)+'-' in en3:
		subprocess.call(['qsub run-'+inp+'-J_'+str(i)+'-'+'_4c-dir'], shell = True)
		en_null = en_null + 1

#--
	with open(inp+'-E-4c_'+inp+'-E_'+str(i)+'--.out') as f:
		for line in f:
			if line[:45] == '   Total energy                             :':
				en3[str(i)+'--'] = float(line[45:])
				break
	if not str(i)+'--' in en3:
		subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'--'+'_4c-dir'], shell = True)
		en_null = en_null + 1
			
#read in energy for all double displacements

#for ++ and --

for i in range(num_freq):
	
	for j in range(num_freq):
	
		if j > i:
		
			with open(inp+'-E-4c_'+inp+'-E_'+str(i)+'+'+str(j)+'+.out') as f:
				for line in f:
					if line[:45] == '   Total energy                             :':
						en3[str(i)+'+'+str(j)+'+'] = float(line[45:])
						break
			if not str(i)+'+'+str(j)+'+' in en3:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'+_4c-dir'], shell = True)
				en_null = en_null + 1
				
			with open(inp+'-E-4c_'+inp+'-E_'+str(i)+'-'+str(j)+'-.out') as f:
				for line in f:
					if line[:45] == '   Total energy                             :':
						en3[str(i)+'-'+str(j)+'-'] = float(line[45:])
						break
			if not str(i)+'-'+str(j)+'-' in en3:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'-'+str(j)+'-_4c-dir'], shell = True)
				en_null = en_null + 1
		
#for +-

for i in range(num_freq):
	
	for j in range(num_freq):
		
		if j != i:
		
			with open(inp+'-E-4c_'+inp+'-E_'+str(i)+'+'+str(j)+'-.out') as f:
				for line in f:
					if line[:45] == '   Total energy                             :':
						en3[str(i)+'+'+str(j)+'-'] = float(line[45:])
						break
			if not str(i)+'+'+str(j)+'-' in en3:
				subprocess.call(['qsub run-'+inp+'-E_'+str(i)+'+'+str(j)+'-_4c-dir'], shell = True)
				en_null = en_null + 1
				
if en_null != 0:
        print('Not all energies calculated!')
        exit()
				

#create empty third-order derivative of energy matrix (only xyy elements)

third_derivative_en = np.zeros((num_freq,num_freq))

#calculate third-order derivative of energy

#diagonal

for i in range(num_freq):
	third_derivative_en[i,i] = (en3[str(i)+'++'] - (2 * en3[str(i)+'+']) + (2 * en3[str(i)+'-']) - en3[str(i)+'--']) / (2 * (float(dis2)**3))

	for j in range(num_freq):
		
#over diagonal
		if i < j:
		
			third_derivative_en[i,j] = (en3[str(i)+'+'+str(j)+'+'] - (2 * en3[str(j)+'+']) + en3[str(i)+'+'+str(j)+'-'] - en3[str(j)+'+'+str(i)+'-'] + (2 * en3[str(j)+'-']) - en3[str(i)+'-'+str(j)+'-']) / (2 * (float(dis2)**3))
				

#under diagonal
		elif i>j:
			third_derivative_en[i,j] = (en3[str(j)+'+'+str(i)+'+'] - (2 * en3[str(j)+'+']) + en3[str(i)+'+'+str(j)+'-'] - en3[str(j)+'+'+str(i)+'-'] + (2 * en3[str(j)+'-']) - en3[str(j)+'-'+str(i)+'-']) / (2 * (float(dis2)**3))
		

#print(third_derivative_en)

#write third-derivative of energy into the out file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('***   Molecular energy third derivative with respect to normal coordinates - only xyy elements [hartree / (bohr**3 * amu ** 3/2)]    ***\n')
	for line in third_derivative_en:
		for row in line:
			f.write(str(round(row,5)).rjust(20))
		f.write('\n')
		

#define constants and converters

#Planck constant [kg * m^2 * s^(-1)]

h_red = 1.054571817 * 10 ** (-34)

hartree_to_Hz = 6.579683920502 * (10 ** (15))

amu_to_kg = 9.1093837015 * (10 ** (-31))

m_to_bohr = 18897259885.789


#convert frequencies from hartree to Hz [s ^ (-1)]

vib_freq_Hz = [hartree_to_Hz * vib_freq_hartree[i] for i in range(int(num_freq))]


#convert everything to bohr and amu

h_red_conv = h_red * (m_to_bohr ** 2) * (1 / amu_to_kg)



#calculate the vibrational corrections

vib_corr = [[] for i in range(int(sscc_num))]


#for each SSCC

for sscc_n in range(int(sscc_num)):

	
	first_derivative_1 = first_derivative[sscc_n]
	first_derivative_u = first_derivative_1[2:num_freq+2]
	
	atoms_sscc = first_derivative_1[0:2]
	
	vib_corr_u = []
	vib_corr_u.append(atoms_sscc[0])
	vib_corr_u.append(atoms_sscc[1])
	
	second_derivative_1 = second_derivative[sscc_n]
	second_derivative_u = second_derivative_1[2:num_freq+2]
	
	element_1 = 0
	element_2 = 0
	
	for i in range(num_freq):
		
#		element_1 = element_1 + ((h_red_conv/ vib_freq_Hz[i]) * second_derivative_u[i])
		
		element_1 = element_1 + ((1/ vib_freq_hartree[i])  * second_derivative_u[i])
	
		third_derivative_en_u = third_derivative_en[i]
		
		element_3 = 0
		
		for j in range(num_freq):

			element_3 = element_3 + (third_derivative_en_u[j] / (vib_freq_hartree[j]))
#			element_3 = element_3 + (third_derivative_en_u[j])

		
#		element_2 = element_2 + (((h_red_conv ** 2)/( vib_freq_Hz[i] ** 2)) * first_derivative_u[i] * element_3)
		element_2 = element_2 + ((1/( vib_freq_hartree[i] ** 2)) * first_derivative_u[i] * element_3)

		
	
	corr = 0.25 * (element_1 -  element_2)
	
	corr_harm = 0.25 * element_1
	corr_anharm = 0.25 * element_2
	
	vib_corr_u.append(corr)
	vib_corr_u.append(corr_harm)
	vib_corr_u.append(corr_anharm)



	
	vib_corr[sscc_n] = vib_corr_u
	
#write vibrational corrections for SSCCs for the molecule into the file

with open(inp+'-'+method+'.out', 'a') as f:
	f.write('\n')
	f.write('***   Vibrational corrrections to spin-spin-couplings J (Hz) for '+inp+'   ***\n')
	
	for n in range(int(sscc_num)):
		vib_corr_2 = vib_corr[n]
		f.write('J('+vib_corr_2[0].rstrip()+'-'+vib_corr_2[1].rstrip()+'):          '+str(round(vib_corr_2[2],2)).rjust(12)+'\n'+ 'harmonic term:   '+str(round(vib_corr_2[3],2)).rjust(12)+'\n'+ 'anharmonic term:  '+str(round(vib_corr_2[4],2)).rjust(12)+'\n'+'\n')
	f.write('\n')
		
	

	












				
		 
