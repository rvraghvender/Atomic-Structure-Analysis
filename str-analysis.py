#!/usr/bin/python3.7
#program to calculate the partial distribution funtion and plot it                    ##
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
import os
import logging
from datetime import datetime
import glob
import shutil




def ProgressBar(Total, Progress, BarLength=200, ProgressIcon="#", BarIcon="-"):
	try:
		# You can't have a progress bar with zero or negative length.
		if BarLength <1:
			BarLength = 20
		# Use status variable for going to the next line after progress completion.
		Status = ""
		# Calcuting progress between 0 and 1 for percentage.
		Progress = float(Progress) / float(Total)
		# Doing this conditions at final progressing.
		if Progress >= 1.:
			Progress = 1
			Status = "\r\n"    # Going to the next line
		# Calculating how many places should be filled
		Block = int(round(BarLength * Progress))
		# Show this
		Bar = "[{}] {:.0f}% {}".format(ProgressIcon * Block + BarIcon * (BarLength - Block), round(Progress * 100, 0), Status)
		return Bar
	except:
		return "ERROR"


def ShowBar(Bar):
	sys.stdout.write(Bar)
	sys.stdout.flush()

directory = 'g-of-r'
directory_angle = 'angles'


dirpath = os.path.join(cwd,directory)
dirpath_angle = os.path.join(cwd,directory_angle)


logging.basicConfig(filename='file.log',format = '%(message)s',level=logging.DEBUG,filemode='w')
logging.getLogger('matplotlib.font_manager').disabled = True
logging.warning("Process Started at : {}".format (datetime.now().strftime("%d/%m/%Y %H:%M:%S")))
logging.warning("Version : 1.0 ")
logging.warning("System Information : {}   {}".format(sys.version, sys.platform))                              ##
logging.warning("Build Info : {}  {} ".format(platform.machine(),platform.system()))
logging.warning("")


print("This program calculate and plot the Pair distribution function and Bond Angle distribution function")
print()
logging.warning("")
logging.warning("This program calculate and plot the Pair distribution function and Bond Angle distribution function")
logging.warning("")


print("Note : This program is only valid for orthorhombic cell as of now ")
print()


A = float(input("Enter lattice parameter a : ")) #Entering cell parameters A,B,C
B = float(input("Enter lattice parameter b : "))
C = float(input("Enter lattice parameter c : "))

logging.warning("Lattice parameters are : a = {} , b = {} , c = {} ".format(A,B,C) )

print()

def distance(a,b):
	
	dx =  a[0] - b[0]
	abs_x = abs(dx)
	if abs_x > A :
		dx = abs_x % A
	x  = min(abs_x, (A - abs_x))

	dy =  a[1] - b[1]
	abs_y = abs(dy)
	if abs_y > B :
		dy = abs_y % B
	y  = min(abs_y, (B - abs_y))


	dz =  a[2] - b[2]
	abs_z = abs(dz)
	if abs_z > C :
		dz = abs_z % C
	z  = min(abs_z, (C - abs_z))

	return np.sqrt(x**2 + y**2 + z**2)


def displacement(a,b):

	r_vec = np.zeros(3)

	dx =  (a[0] - b[0])
	abs_x = abs(dx)
	if abs_x > A :
		dx = abs_x % A
	if abs_x > A/2 and dx > 0:
		x = min(abs_x,A - abs_x)
		x = -x
	elif abs_x > A/2 and dx < 0:
		x = min(abs_x,A- abs_x)
		x = x
	else:
		x = dx
	r_vec[0] = x


	dy =  (a[1] - b[1])
	abs_y = abs(dy)
	if abs_y > B :
		dy = abs_y % B
	if abs_y > B/2 and dy > 0:
		y = min(abs_y,A - abs_y)
		y = -y
	elif abs_y > B/2 and dy < 0:
		y = min(abs_y,B- abs_y)
		y = y
	else:
		y = dy
	r_vec[1] = y


	dz =  (a[2] - b[2])
	abs_z = abs(dz)
	if abs_z > C :
		dz = abs_z % C
	if abs_z > C/2 and dz > 0:
		z = min(abs_z,C - abs_z)
		z = -z
	elif abs_z > C/2 and dz < 0:
		z = min(abs_z,C- abs_z)
		z = z
	else:
		z = dz
	r_vec[2] = z

	return r_vec

def volume_sphere(r):
	volume = 4.0 / 3.0 * np.pi * r**3    #  volume of sphere

	if r > min(A,B,C) /2.0 :             # this will remove the spherical cap if we go more than defined.
		volume -= 6 * np.pi * (( 2 * r +  (min(A,B,C) / 2.0) ) * ( r - ( min(A,B,C) / 2.0 ))**2)/ 3 # formula from wiki
	return volume


class Trajectory:
	def __init__(self,filename, skip, resolution = 500):
		self.filename = filename
		self.skip = skip
		self.resolution = resolution
		self.parse_input()
		logging.warning("Resolution : {}".format(self.resolution))

	def parse_input(self):  # for getting coordinates from the trajectory file
		global Atom_type
		with open(self.filename, 'r') as fo: #opening file
			data = fo.readlines()   # reading each line

		self.n_atoms = int(data[0].split()[0]) # reading first line as number of atoms which is expected to be integer
		self.n_steps_total = int(len(data) / (self.n_atoms +2))  # we then find the total timestep by dividing the total
		# no. of lines in trajectory files by n_atoms (no. of atoms) +2 header (2 is for additional two lines (1. no.
		# of atoms 2. comment or blank line))
		logging.warning("Total numbers of atoms found : {}".format(self.n_atoms))
		logging.warning("Total numbers of configuration found : {}".format(self.n_steps_total))


		self.n_steps = self.n_steps_total // self.skip # this allows to check on every interval timestep as atoms does
		#not move appreciably in each step. Of course setting skip = 1, will give 100% accuracy.
		logging.warning("Number of configuration after skipping for which it is analysed : {}". format(self.n_steps))



		self.atom_list = [line.split()[0] for line in data[2:self.n_atoms + 2]]  # here i am trying to capture each atoms
		# i will start reading from 3 rows and will go upto total no. of atoms + 2 that's why data[2:self.n_atoms + 2]]
		# also line.np.it()[0] will read the first element of the list which is indeed atom name

		Atom_type = list(dict.fromkeys(self.atom_list)) ## types of atom



		self.coordinates = np.zeros((self.n_steps, self.n_atoms, 3)) # this will create n_steps 2D matrix with each
		# (n_atoms rows) and (3 columns) all filled with zeros.
		# n_atoms rows can be regarded as total no. of atoms and 3 columns will capture x,y,z coordinates details.
		# similarty n_step matrices will act as timestep .



		for step in range(self.n_steps): # here i am creating a loop with nsteps
			coords = np.zeros((self.n_atoms, 3)) # for each step i am creating zero matrix with n_atoms rows and 3 colu.
			i = step * self.skip * (self.n_atoms + 2)  # now this i will iterate each atom index in increasing timestep.
			# and thus this i represents move forward in timestep starting from 0.

			for j, line in enumerate(data[i + 2 : i + self.n_atoms + 2]): # here i am enumerating in each current
				# timestep "i". Thus this loop scans whole atoms in a single timestep.

				coords[j, :] = [float(value) for value in line.split()[1:]]# here data is a list where number of
				# element is equals to no. of lines in trajectory file. Also, line is temporary variable name
				# assigned to each element of data. And this line is alias for atom name, x ,y and z of each atom.
				#Therefore the list has started from 1 i.e. [1:] so that it ignores the atom name.



			self.coordinates[step] = coords  # this will assign the current value of configuration of timestep
			# to coordinates list. In this way we have n_step matrices with atomic coordinates filled up.


	def compute_volume_per_TeO2(self): # for finding the volume per teO2 molecule. Also we will use this number to
		# normalize our RDF, so that it get saturates to certain value at large distances.
		global atom_name_2
		global atom_name_1
		global compute_volume_per_TeO2

		n_desired_atom = 0 # finding no. of desired atoms in a cell
		for step in range(self.n_steps): # here i am creating a loop with nsteps
			for i, atom in enumerate(self.coordinates[step]):   # this loop is for each timestep's coordinates
				if self.atom_list[i] == atom_name_2: #"or self.atom_list[i] = "O"": will look for other desired atoms.
					n_desired_atom += 1 # after getting any desired atom it will add one to it
		volume = A * B * C  # volume of the cell

		if atom_name_1 == atom_name_2:
			remove = 1   ## this will be useful for subtracting out the no. of host atoms in normalisation
		else:  # for instance if i am looking for Te Te RDF then it will subtract nsteps * one host atoms
			remove = 0 # and if we are interested in Te O then nothing will be subtracted out
		
		average_n_desired_atom = (n_desired_atom - self.n_steps * remove) / self.n_steps   # the above loop will count
		# no. of desired atoms in
		# single timestep as well as total no. of  timestep. Thus IF no. of atoms remains same (NVT,NPT or NVE) it
		# will be no. of atoms times nstep. which we have to normalise thats why divide by nsteps (which logically
		# doesn't make sense, but just to make it general.) May be in future if i simulate a non constant no. of atoms
		# simulation. Let's see. ('.')

		compute_volume_per_TeO2 = volume / average_n_desired_atom

	# -----------------------------------------------------------------------------------------------------------------
	# uffff! Now having sorted out the preliminaries. Let's dive into real calculation part... tk tk tk

	def compute_radial_distribution(self):
		global atom_name_1
		global atom_name_2
		r_cutoff = 9.995
		logging.warning("Cut-off radius is : {}".format(r_cutoff))


		dr = r_cutoff / self.resolution # giving the resolution means how small (thick) i need my bin size.
		self.radii = np.linspace(0.010, r_cutoff, self.resolution)  #this will create an 1-D array with with
		# uniform spacing. Starting from 0.0 go upto rcutoff with self.resolution times of equal spacing.
		# so number of elements will be equal to resolution given.
		logging.warning("Width of infintesimal radius is : {}".format(dr))


		volumes = np.zeros(self.resolution) # this will create an 1-D zero array with resolution no. of elements.
		self.g_of_r = np.zeros(self.resolution) # this will create an 1-D zero array with resolution no. of elements.


		for step in range(self.n_steps):
			print(" Total Steps : {} ".format(self.n_steps), end='\r') #it will print live stepcount and total stepcount
			progressBar = "\rProgress: " + ProgressBar(self.n_steps-1, step, 100)
			ShowBar(progressBar)
			data_atoms_1 = [] # here we are creating empty list for host atom which later will contain atomic coordinates
			data_atoms_2 = []# here we are creating empty list for non host atom which later will contain atomic coordinates
			for i, atom in enumerate(self.coordinates[step]):  # here i am looping as i, also looking the value at i.
				if self.atom_list[i] == atom_name_1 : # if atom_list finds the given atom it will ,//
						data_atoms_1.append(atom)  # append its coordinates in data_atoms list.
				if self.atom_list[i]==  atom_name_2 : # this code will look for other non host atoms
						data_atoms_2.append(atom)
			data_atoms_1 = np.array(data_atoms_1) # just appending in numpy format
			data_atoms_2 = np.array(data_atoms_2) # appending the coordinates of non host atoms

			for i, atom1 in enumerate(data_atoms_1):
				for j in range(self.resolution):
					r1 = ( j + 0.010 )  * dr # this will give the lower radius of the shell
					r2 = r1 + dr # this will give the higher radius of the shell
					v1 = volume_sphere(r1) # this will find the volume from the defined function above.
					v2 = volume_sphere(r2) # same as above
					volumes[j] += v2 - v1 # will find the volume of the shell and will store in a list
					# thus this list will have all the differences of volume of succeding radius.


				for k, atom2 in enumerate(data_atoms_2): #  this will look for other atoms of same np.cies
					dist = distance(atom1, atom2) # this will call the distance function and compute distance
					index = int(dist / dr) # this can give value either more than self.resolution (if other atom seen is
					# place more than rcutoff) or less than self.resolution( if other atom seen inside rcutoff.)

					if 0 < index < self.resolution: # this will categories only inside atoms of rcutoff.
						self.g_of_r[index] += 1.0 # this will count the pair as two similar atoms

		for i, value in enumerate(self.g_of_r):
			self.g_of_r[i] = value * compute_volume_per_TeO2 / volumes[i]# the first line loops over all the
			# elements of g(r) . There are resolution numbers of
			# elements for both g(r) and volumes. Also both values were allowed to run for all the timesteps.
			# thus can be divided. For ex. volume[1] is carrying the difference in vol between (v2 -v1) * total no. of
			# steps. So here we have normalised here volume of the np.erical shells and volume per atom
			# calculated earlier.

		for ik in range(len(self.radii)):
					fw.write("\n  {:1.3f}          {:.10f}  ".format(self.radii[ik],self.g_of_r[ik]))

	def plot(self, filename=""):
		if not self.g_of_r.any():
			print('compute the radial distribution function first\n')
			return
		plt.title('Radial Distribution Function of TeO$_{2}$ glass')
		plt.xlabel('r (A$^0$)')
		plt.ylabel('g$_{ab}$(r)')
		plt.plot(self.radii, self.g_of_r,label = '{}-{}'.format(atom_name_1,atom_name_2))
		plt.legend()
		plt.savefig('g_of_r_{}-{}.png'.format(atom_name_1,atom_name_2), dpi=300, bbox='tight', format='png')
		plt.clf()


	def compute_bond_angle_distribution_funtion(self):
		global atom_name_1
		global atom_name_2
		global atom_name_3
		self.angle = np.linspace(0.0, 200, 200, dtype=int)
		self.bdf = np.zeros(200)
		for step in range(self.n_steps):
			angle_data_atom_1 = []
			angle_data_atom_2 = []
			angle_data_atom_3 = []
			print(" Total Steps : {} ".format(self.n_steps),
				  end='\r')  # it will print live stepcount and total stepcount
			progressBar = "\rProgress: " + ProgressBar(self.n_steps -1 , step, 100)
			ShowBar(progressBar)
			for i, atom in enumerate(self.coordinates[step]):
				if self.atom_list[i] == atom_name_1 :
						angle_data_atom_1.append(atom)
				if self.atom_list[i]==  atom_name_2 :
						angle_data_atom_2.append(atom)
				if self.atom_list[i]==  atom_name_3 :
						angle_data_atom_3.append(atom)
			angle_data_atom_1 = np.array(angle_data_atom_1)
			angle_data_atom_2 = np.array(angle_data_atom_2)
			angle_data_atom_3 = np.array(angle_data_atom_3)

			for i, atom1 in enumerate(angle_data_atom_1):
				for j, atom2 in enumerate(angle_data_atom_2):
					rij = distance(atom1,atom2)
					if 0.0 < rij <= bond_len_cut_pair12:
						for k, atom3 in enumerate(angle_data_atom_3):
							rjk = distance(atom2,atom3)
							if 0.0 < rjk <= bond_len_cut_pair23:
								r12 = displacement(atom2,atom1)
								r32 = displacement(atom2,atom3)
								#print()
								#print(atom2,atom3)
								#print(r12,r32)
								dotprod = np.dot(r12,r32)
								costheta = dotprod / (rij * rjk)
								#print(dotprod, rij, rjk)
								radian_angle = np.arccos(costheta.round(decimals=3))
								degree_angle = np.rad2deg(radian_angle)
								angle_index = int(degree_angle)
								if angle_index != 0:
									self.bdf[angle_index] += 1.0

		for i, value in enumerate(self.bdf):
			self.bdf[i] = value  / ( self.n_steps )
		for ik in range(len(self.angle)):
			fw_angle.write("\n  {}          {}  ".format(self.angle[ik],self.bdf[ik]))

	def bond_plot(self, filename =""):
		if not self.bdf.any():
			print('compute the bond-angle distribution function first\n')
			return
		plt.title('Bond-angle Distribution Function of TeO$_{2}$ glass')
		plt.xlabel('Theta (degrees) ')
		plt.ylabel(' f(theta) ')
		plt.plot(self.angle, self.bdf,label = '{}-{}-{}'.format(atom_name_1,atom_name_2,atom_name_3))
		plt.legend()
		plt.savefig('bdf_{}-{}-{}.png'.format(atom_name_1,atom_name_2,atom_name_3), dpi=300, bbox='tight', format='png')
		plt.clf()



fo1 = input("Enter the name of trajectory file (.xyz) format : ")
skip_ts = int(input("Enter the integer value to skip steps : "))

logging.warning("Name of the trajectory file has been set to : {} ".format(fo1))
logging.warning("Skip (Jump) steps set : {}".format(skip_ts))

print()


logging.warning("")
logging.warning("")


rdf = input("Do you want to calculate 'Radial Distribution Function'  (yes/no) : ")
print()
if rdf == 'yes':
	logging.warning("User opted  to calculate RDF analysis. ")
	TeO2 = Trajectory(fo1,skip_ts)
	logging.warning("Types of atoms are : {}".format(Atom_type))
	logging.warning("")
	for ik, value in enumerate(Atom_type):
		#print(ik,Atom_type[ik],"df")
		for jk, value2 in enumerate(Atom_type[ik:]):
			atom_name_2 = value2
			atom_name_1 = value
			print()
			print("RDF between",atom_name_1,"and",atom_name_2)
			TeO2.compute_volume_per_TeO2()
			fw = open("g_of_r-{}-{}.dat".format(atom_name_1, atom_name_2), 'w')
			TeO2.compute_radial_distribution()
			TeO2.plot()
			print()
			fw.close()
			logging.warning("")
	logging.warning("RDF for each configuration has been created successfully.")
	if os.path.exists(dirpath) and os.path.isdir(dirpath):
		shutil.rmtree(dirpath)
	os.mkdir(dirpath)
	source_path_list = glob.glob(cwd + '/g_of_r*')
	for source_path in source_path_list:
		shutil.move(source_path, dirpath)


logging.warning("")
logging.warning("")

print()

angle_prop = input("Do you want to calculate 'Bond angle distribution'  (yes/no) : ")
print()
if angle_prop == 'yes':
	logging.warning("User opted to calculate Bond Angle analysis")
	while True:
		print("Enter the configuration you are interested in, but be careful about the names of element (Not case sensitive) ")
		print()
		atom_name_1 = input("Enter the name of first atom : ")
		atom_name_2 = input("Enter the name of second (central) atom : ")
		atom_name_3 = input("Enter the name of third atom : ")
		fw_angle = open("bdf_{}-{}-{}.dat".format(atom_name_1, atom_name_2,atom_name_3), 'w')
		bond_len_cut_pair12 = float(input("Enter the bond cut-off length of Atom 1 - {} and Atom 2 - {} : ".format(atom_name_1, atom_name_2)))
		bond_len_cut_pair23 = float(input("Enter the bond cut-off length of Atom 2 - {} and Atom 3 - {} : ".format(atom_name_2, atom_name_3)))
		logging.warning("Cut-off for Caculating Bond Angle Distribution between atoms pair Atom 1 - {} and Atom 2 - {} is {} ".format(atom_name_1, atom_name_2,bond_len_cut_pair12))
		logging.warning(
			"Cut-off for Caculating Bond Angle Distribution between atoms pair Atom 2 - {} and Atom 3 - {} is {} ".format(
				atom_name_2, atom_name_3, bond_len_cut_pair23))

		TeO2 = Trajectory(fo1, skip_ts)
		TeO2.compute_bond_angle_distribution_funtion()
		TeO2.bond_plot()
		print()
		fw_angle.close()
		logging.warning("Bond Angle analysis for configuration - {}-{}-{} has been done. ".format(atom_name_1, atom_name_2,atom_name_3))

		what_s_next = input("Are you interested in finding RDF for other configurations (yes/no) : ")
		print()
		if what_s_next == 'yes':
			continue
		elif what_s_next == 'no':
			break
		else:
			break
	if os.path.exists(dirpath_angle) and os.path.isdir(dirpath_angle):
		shutil.rmtree(dirpath_angle)
	os.mkdir(dirpath_angle)
	source_path_list_angle = glob.glob(cwd + '/bdf*')
	for source_path in source_path_list_angle:
		shutil.move(source_path, dirpath_angle)

logging.warning("")
logging.warning("-------------------------------------------------------------------------------")

logging.warning("Process Completed at : {}".format (datetime.now().strftime("%d/%m/%Y %H:%M:%S")))

print()
print("Your files has been created succesfully")


