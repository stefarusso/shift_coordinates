#!/Users/stefano/anaconda3/envs/p4env/bin/python
import sys
import subprocess
import xyz2gro

#./shift.py fle.xyz
#arc require 2 line before coordinates
#require topol.top and all the .itp files outside the folder 
#
#	shift.py
#	amoeba.key
#	amoeba.prm
#	file.arc
#	folder
#		../shift.py ../file.arc



N1=int(input('Number of atoms for molecule 1 \n'))
charge_multiplicity1=input('Charge and Multiplicity for molecule 1 [0 1] \n')
if len(charge_multiplicity1) == 0:
    charge_multiplicity1 = "0 1"

N2=int(input('Number of atoms for molecule 2 \n'))
charge_multiplicity2=input('Charge and Multiplicity for molecule 1 [-1 1] \n')
if len(charge_multiplicity2) == 0:
    charge_multiplicity2 = "-1 1"


A1=int(input('Atom index for molecule 1 (start from 1)\n'))
A2=int(input('Atom index for molecule 2 (start from 1)\n'))
s= float(input('Max displacment\n'))
step= int(input('Number of steps\n'))
s=s/step

box_length=float(input('Length of the box [A] \n'))


if(not s==0):
	#origin
	step = step + 1
	
is_arc=False



A1=A1-1
A2=A2-1

if len(sys.argv) > 1:
	filename= sys.argv[1]
else:
	sys.exit("Formato richiesto: python shift.py file.xyz")




#save coordinates
#NEED ARC FILE
if(is_arc):
	atoms = []
	coordinates = []
	topo = []
	xyz = open(filename)
	n_atoms = int(xyz.readline().split()[0])
	title = xyz.readline()
	for line in xyz:
		line = line.strip()
		idx,atom , x , y , z , top = line.split(maxsplit=5)
		if(len(atom)>1):
			atom=atom[0]
		atoms.append(atom)
		coordinates.append([float(x), float(y), float(z)])
		topo.append(top)
	xyz.close()
else:
	atoms = []
	coordinates = []
	xyz = open(filename)
	n_atoms = int(xyz.readline().split()[0])
	title = xyz.readline()
	for line in xyz:
	    atom,x,y,z = line.split()
	    atoms.append(atom)
	    coordinates.append([float(x), float(y), float(z)])
	xyz.close()


#M2 positive shift
#M1 negative shift
delta = [  coordinates[A2][0]-coordinates[A1][0] , coordinates[A2][1]-coordinates[A1][1] , coordinates[A2][2]-coordinates[A1][2]  ]
norm = (delta[0]**2+delta[1]**2+delta[2]**2)**0.5

#Normalized versor
delta = [i/norm  for i in delta]

D=norm
d=0
#LOOP OVER N STEPS + INITIAL FRAME
for j in range(step):
	n_coordinates = []
	print("d:"+str(d))
	print("D:"+str(D))
		
	#DISPLACEMENT
	#FIRST MOLECULE, NEGATIVE SHIFT
	for xyz in coordinates[0:N1]:
		n_coordinates.append([xyz[0]-delta[0]*d,xyz[1]-delta[1]*d,xyz[2]-delta[2]*d])
	
	# SECOND MOLECULE, POSITIVE SHIFT
	for xyz in coordinates[N1:]:
		n_coordinates.append([xyz[0]+delta[0]*d,xyz[1]+delta[1]*d,xyz[2]+delta[2]*d])
	
	
	#XYZ_AMOEBA
	filename=f"{D:.1f}.xyz"
	with open(filename, 'w') as f:
	        f.write("%s\n" %n_atoms)
	        f.write("\n")
	
	with open(filename, 'a') as f:
		for i in range(len(n_coordinates)):
			if(is_arc):
				f.write("%d\t%s\t%.6f\t%.6f\t%.6f\t%s\n" %(i+1,atoms[i],n_coordinates[i][0],n_coordinates[i][1],n_coordinates[i][2],topo[i]))
			else:
				f.write("%s\t%.6f\t%.6f\t%.6f\n" %(atoms[i],n_coordinates[i][0],n_coordinates[i][1],n_coordinates[i][2]))
        

	#INPUT SAPT
	filename_sapt=f'{D:.1f}.in'
	with open(filename_sapt, 'w') as f:
                f.write("memory 500gb\nset_num_threads(2)\nmolecule dimer {\n%s\n" %(charge_multiplicity1))

	with open(filename_sapt, 'a') as f:
		for i in range(N1):
			f.write("%s\t%.6f\t%.6f\t%.6f\n" %(atoms[i],n_coordinates[i][0],n_coordinates[i][1],n_coordinates[i][2]))

	with open(filename_sapt, 'a') as f:
		f.write("--\n%s\n" %(charge_multiplicity2))
	
	with open(filename_sapt, 'a') as f:
		for i in range(N1,N1+N2):
			f.write("%s\t%.6f\t%.6f\t%.6f\n" %(atoms[i],n_coordinates[i][0],n_coordinates[i][1],n_coordinates[i][2]))
	with open(filename_sapt, 'a') as f:
		f.write("units angstrom\nno_reorient\nsymmetry c1\n}\nset basis 6-311G\nenergy('sapt0')\n")
	
#	# XYZ -> GRO
#	import glob
#	for filename in glob.glob('./*.xyz'):
#                f=filename.split(".")
#                f.pop()
#                f='.'.join(str(i) for i in f)
#                xyz2gro.convert(box_length,"topol.top",f+".xyz",f+".gro")
	

	#START SAPT
	#cmd=["psi4",filename_sapt,f'{filename_sapt}.psi']
	#psi_out = subprocess.run(cmd, stdout=subprocess.PIPE)
	
	
	
	
	d=d+s/2
	D=D+s

import glob
for filename in glob.glob('./*.xyz'):
        f=filename.split(".")
        f.pop()
        f='.'.join(str(i) for i in f)
        xyz2gro.convert(box_length,"topol.top",f+".xyz",f+".gro")

cmd=["../gromacs.sh"]
out = subprocess.run(cmd, stdout=subprocess.PIPE)
