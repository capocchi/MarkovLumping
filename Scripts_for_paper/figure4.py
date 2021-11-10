#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import matplotlib.pyplot as plt

### for plot
plot = False

if plot:
	rates = {"Reduction Rate compared to k=n-1":{}, "Reduction Rate compared to all partitions":{}}

def run_command(cmd,matrix):
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout_iterator = iter(process.stdout.readline, b"")
	for line in stdout_iterator:
		print(line)

		if plot:
			txt = str(line)
			if "Reduction Rate compared to k=n-1" in txt:
				val = float(txt.split(':')[1].split('%')[0].replace(" ",""))
				rates["Reduction Rate compared to k=n-1"][matrix].append(round(val,2))
			elif "Reduction Rate compared to all partitions" in txt:
				val = float(txt.split(':')[1].split('%')[0].replace(" ",""))
				rates["Reduction Rate compared to all partitions"][matrix].append(round(val,2))

### filename of the output file that contain the stdout
try:
    fn = sys.argv[1]
except:
    fn = None
else:	
	if fn:
		OUTPUT_FILENAME = fn

		if os.path.isfile(OUTPUT_FILENAME):
			os.remove(OUTPUT_FILENAME)

		class Unbuffered:

		   def __init__(self, stream):

		       self.stream = stream

		   def write(self, data):

		       self.stream.write(data)
		       self.stream.flush()
		       file.write(data)    # Write the data of stdout here to a text file as well

		file = open(OUTPUT_FILENAME,"w")
		sys.stdout=Unbuffered(sys.stdout)

if __name__ == "__main__":

	### for random distribution	
	for dist in ('uniform', 'beta', 'binomial', 'rayleigh', 'weibull'):

		### matrix dim
		for i in (3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100):
			if plot:
				X = []
				Y = []

			#for m in ['0.1','1','5', '10', '50', '100', '500', '1000', '5000', '10000']:
			
				#matrix = f"{i}x{i}_{m}"

			if dist == "uniform":
				matrix = f"{i}x{i}_0.1"
			elif dist == "beta":
				matrix = f"{i}x{i}_1_3"
			elif dist=="binomial":
				matrix = f"{i}x{i}_5_0.5"
			elif dist=="rayleigh":
				matrix = f"{i}x{i}_10"
			elif dist=="weibull":
				matrix = f"{i}x{i}_5"
			else:
				sys.exxit()

			print(f"\n---------------------------- Matrix {matrix} - {dist}")
			
			if plot:
				rates["Reduction Rate compared to k=n-1"][matrix] = []
				rates["Reduction Rate compared to all partitions"][matrix] = []

			pool = "python main_mftp.py "+os.path.join("..","Matrix","bench","diag_max",dist,f"{matrix}.dat")
			run_command(pool, matrix)

			if plot:
				X.append(f"{i}x{i}")
				Y.append(rates["Reduction Rate compared to all partitions"][matrix])

		if plot:
			### comment to display the right curves
			plt.plot(X, Y)
			#plt.plot(X, Y2, label = "compared to all partitions")
			
			# show a legend on the plot
			plt.legend()
			
			plt.title(f"Matrix {i}x{i}")

			plt.grid()
		
			# function to show the plot
			plt.show()
