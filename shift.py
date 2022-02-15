import sys

#modify with your atoms
#-----------------------#
shifted_atoms=[0,21]	#
total_n=34		#
#-----------------------#


filename = sys.argv[1]
if len(sys.argv) > 3:
	s= float(sys.argv[2])
	step= int(sys.argv[3])
else:
	sys.exit("Formato reqeuired: python shift.py coordinates.xyz total_shift n_steps")


s=s/step
atoms = []
coordinates = []
xyz = open(filename)
n_atoms = int(xyz.readline())
title = xyz.readline()
for line in xyz:
    atom,x,y,z = line.split()
    atoms.append(atom)
    coordinates.append([float(x), float(y), float(z)])
xyz.close()


deltax=coordinates[shifted_atoms[0]][0]-coordinates[shifted_atoms[1]][0]
deltay=coordinates[shifted_atoms[0]][1]-coordinates[shifted_atoms[1]][1]
deltaz=coordinates[shifted_atoms[0]][2]-coordinates[shifted_atoms[1]][2]
norm=(deltax**2+deltay**2+deltaz**2)**0.5
deltax=deltax/norm
deltay=deltay/norm
deltaz=deltaz/norm

d=norm
for j in range(step):
	n_coordinates = []
	for i in range(shifted_atoms[1]):
	     n_coordinates.append([coordinates[i][0]+deltax*s/2,coordinates[i][1]+deltay*s/2,coordinates[i][2]+deltaz*s/2])
	
	for i in range(shifted_atoms[1],total_n):
	     n_coordinates.append([coordinates[i][0]-deltax*s/2,coordinates[i][1]-deltay*s/2,coordinates[i][2]-deltaz*s/2])
	
	d=d+s
	
	with open('%.1f.xyz' %d, 'a') as f:
	        f.write("%s\n" %n_atoms)
	        f.write("%s" %title)
	
	with open('%.1f.xyz' %d, 'a') as f:
		for i in range(len(n_coordinates)):
			f.write("%s\t%.6f\t%.6f\t%.6f\n" %(atoms[i],n_coordinates[i][0],n_coordinates[i][1],n_coordinates[i][2]))
	coordinates=n_coordinates
