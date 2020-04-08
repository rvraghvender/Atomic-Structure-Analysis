#!/usr/bin/python3.7
#program to calculate the partial distribution funtion and plot it                    ##
from typing import Any

import platform                                                                       ##
import sys                                                                            ##
																					  ##
print("System Information : ",sys.version, sys.platform)                              ##
print("Build Info : ",platform.machine(),platform.system())                           ##
print()                                                                               ##
print()                                                                               ##
																					  ##
import os                                                                             ##
cwd = os.getcwd()  # Get the current working directory (cwd)                          ##
files = os.listdir(cwd)  # Get all the files in that directory                        ##
print("Available Files in directory %r : " % (cwd))                                   ##
print()                                                                               ##
for counter, file_items in enumerate(files,1):                                        ##
	print(counter,"--" ,file_items)                                                   ##
print()							                              ##
print("--------------------------------------------------------------------------")   ##
																					  ##
print()                                                                               ##
print()                                                                               ##
																					  ##
########################################################################################


import matplotlib.pyplot as plt
import numpy as np


print("This program calculate and plot the Pair distribution function")
print()

print("Note : This program is only valid for orthorhombic cell as of now ")
print()







A = float(input("Enter lattice parameter a : ")) #Entering cell parameters A,B,C
B = float(input("Enter lattice parameter b : "))
C = float(input("Enter lattice parameter c : "))

print()


def distance(a, b):
	dx = a[0] - b[0]
	abs_x = abs(dx)
	if abs_x > A:
		dx = abs_x % A
	x = min(abs_x, (A - abs_x))

	dy = a[1] - b[1]
	abs_y = abs(dy)
	if abs_y > B:
		dy = abs_y % B
	y = min(abs_y, (B - abs_y))

	dz = a[2] - b[2]
	abs_z = abs(dz)
	if abs_z > C:
		dz = abs_z % C
	z = min(abs_z, (C - abs_z))

	return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def volume_sphere(r):
	volume = 4.0 / 3.0 * np.pi * r**3    #  volume of sphere

	if r > min(A,B,C) /2.0 :             # this will remove the spherical cap if we go more than defined.
		volume -= 6 * np.pi * ( 2 * r +  (min(A,B,C) / 2.0) ) * ( r - ( min(A,B,C) / 2.0 ))**2 / 3 # formula from wiki
	return volume


class Trajectory:
	def __init__(self,filename, skip, interface_offset = 4.0, resolution = 500):
		self.filename = filename
		self.skip = skip
		self.resolution = resolution
		self.parse_input()
		self.compute_volume_per_TeO2()

	def parse_input(self):  # for getting coordinates from the trajectory file
		with open(self.filename, 'r') as fo:  #opening file
			data = fo.readlines()              # reading each line

		self.n_atoms = int(data[0].split()[0])  # reading first line as number of atoms which is expected to be integer
		self.n_steps_total = int(len(data) / (self.n_atoms +2))  # we then find the total timestep by dividing the total
		# no. of lines in trajectory files by n_atoms (no. of atoms) +2 header (2 is for additional two lines (1. no.
		# of atoms 2. comment or blank line))
		self.n_steps = self.n_steps_total // self.skip # this allows to check on every interval timestep as atoms does
		#not move appreciably in each step. Of course setting skip = 1, will give 100% accuracy.
		self.atom_list = [line.split()[0] for line in data[2:self.n_atoms + 2]] # here i am trying to capture each atoms
		# i will start reading from 3 rows and will go upto total no. of atoms + 2 that's why data[2:self.n_atoms + 2]]
		# also line.np.it()[0] will read the first element of the list which is indeed atom name
		self.coordinates = np.zeros((self.n_steps, self.n_atoms, 3))
		# this will create n_steps 2D matrix with each (n_atoms rows) and (3 columns) all filled with zeros.
		# n_atoms rows can be regarded as total no. of atoms and 3 columns will capture x,y,z coordinates details.
		# similarty n_step matrices will act as timestep .
		for step in range(self.n_steps): # here i am creating a loop with nsteps
			coords = np.zeros((self.n_atoms, 3)) # for each step i am creating zero matrix with n_atoms rows and 3 colu.
			i = step * self.skip * (self.n_atoms + 2)  # now this i will iterate each atom index in increasing timestep.
			# and thus this i represents move forward in timestep starting from 0.
			for j, line in enumerate(data[i + 2 : i + self.n_atoms + 2]): # here i am enumerating in each current
				# timestep "i". Thus this loop scans whole atoms in a single timestep.
				coords[j, :] = [float(value) for value in line.split()[1:]] # here data is a list where number of
				# element is equals to no. of lines in trajectory file. Also, line is temporary variable name
				# assigned to each element of data. And this line is alias for atom name, x ,y and z of each atom.
				#Therefore the list has started from 1 i.e. [1:] so that it ignores the atom name.
			self.coordinates[step] = coords # this will assign the current value of configuration of timestep
			# to coordinates list. In this way we have n_step matrices with atomic coordinates filled up.

	def compute_volume_per_TeO2(self):
		global compute_volume_per_TeO2
		n_desired_atom = self.n_steps * self.n_atoms
		volume = A * B * C

		average_n_desired_atom = (n_desired_atom - self.n_steps ) / self.n_steps
		compute_volume_per_TeO2 = volume / average_n_desired_atom


	def compute_radial_distribution(self):
		r_cutoff = 9.995
		dr = r_cutoff / self.resolution # giving the resolution means how small (thick) i need my bin size.
		self.radii = np.linspace(0.010, r_cutoff, self.resolution)
		volumes = np.zeros(self.resolution) # this will create an 1-D zero array with resolution no. of elements.
		self.g_of_r = np.zeros(self.resolution) # this will create an 1-D zero array with resolution no. of elements.


		for step in range(self.n_steps):
			print('{:4d} : {:4d}'.format(step, self.n_steps)) #it will print live stepcount and total stepcount
			data_atoms = [] # here we are creating empty list for host atom which later will contain atomic coordinates
			#data_atoms_2 = [] # here we are creating empty list for non host atom which later will contain atomic coordinates
			for i,atom in enumerate(self.coordinates[step]):
						data_atoms.append(atom)
			data_atoms = np.array(data_atoms) # just appending in numpy format


			for i, atom1 in enumerate(data_atoms):
				for j in range(self.resolution):
					r1 = ( j + 0.010 ) * dr # this will give the lower radius of the shell
					r2 = r1 + dr # this will give the higher radius of the shell
					v1 = volume_sphere(r1) # this will find the volume from the defined function above.
					v2 = volume_sphere(r2) # same as above
					volumes[j] += v2 - v1  # will find the volume of the shell and will store in a list
					# thus this list will have all the differences of volume of succeding radius.

				for atom2 in data_atoms[i:]: #  this will look for other atoms of same species
					dist = distance(atom1, atom2) # this will call the distance function and compute distance
					index = int(dist / dr) # this can give value either more than self.resolution (if other atom seen is
					# place more than rcutoff) or less than self.resolution( if other atom seen inside rcutoff.)
					if 0 < index < self.resolution: # this will categories only inside atoms of rcutoff.
						self.g_of_r[index] += 2.0 # this will count the pair as two similar atoms

		for i, value in enumerate(self.g_of_r):
			self.g_of_r[i] = value * compute_volume_per_TeO2 / volumes[i]
			# the first line loops over all the elements of g(r) . There are resolution numbers of
			# elements for both g(r) and volumes. Also both values were allowed to run for all the timesteps.
			# thus can be divided. For ex. volume[1] is carrying the difference in vol between (v2 -v1) * total no. of
			# steps. So here we have normalised here volume of the np.erical shells and volume per atom
			# calculated earlier.
		for ik in range(len(self.radii)):
			fw.write("\n  {}          {}  ".format(self.radii[ik],self.g_of_r[ik]))

	def plot(self, filename=""):
		if not self.g_of_r.any():
			print('compute the radial distribution function first\n')
			return
		plt.title('Radial Distribution Function of TeO$_{2}$ glass')
		plt.xlabel('r (A$^0$)')
		plt.ylabel('g$_{ab}$(r)')
		plt.plot(self.radii, self.g_of_r,label = 'Total RDF')
		plt.legend()
		plt.savefig('g_of_r_total.png', dpi=300, bbox='tight', format='png')
		plt.clf()








fo1 = input("Enter the name of trajectory file (.xyz) format : ")
skip_ts = int(input("Enter the integer value to skip steps : "))


print()


fw = open("g_of_r-tot.dat",'w')
TeO2 = Trajectory(fo1,skip_ts)
TeO2.compute_radial_distribution()
TeO2.plot()
fw.close()

print("Your files has been created succesfully")


