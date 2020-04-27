# AsterQuant
Multi aster system quantitative analysis using voronoi tessellation

'''
Neha Khetan, 
AsterQuant: 2017 - 
PYTHON SCRIPT FOR VORONOI TESSELLATIONS OF THE ASTERS in the EXPERIMENTS
'''


## This folder contains 3 scripts:
Main file: AsterQuant_V4.py
           and
           plot_voronoi.py
           plot_Simvoronoi.py
-> All the three files should be in the same path




## -------------------------------
To run
>> python AsterQuant_V4.py


## -------------------------------
Version and Libraries
1. Python version 3.7.7
2. numpy, scipy, 
3. shapely, pylab , matplotlib



>> Folder structure:
   a. expt b. sim  // These contain input files and images
   b. output_fig and output_files // to store the outputs from a run
   



## --------------      FOR DEMO PURPOSE: Example files and settings have been chosen as-------------------

a) Demo of Experimental data:

   Option 1: 

	exptdata	         = 1    	# For Experiment =1 ; for simulation = 0
	selectImage              = 1            # Selection of aster centroids: From Image = 1 ; From files = 0
	px2micron                = 0.19         

	
   Option 2:
 	exptdata	         = 1    	# For Experiment =1 ; for simulation = 0
	selectImage              = 0            # Selection of aster centroids: From Image = 1 ; From files = 0
	px2micron                = 0.19         




In option 1: " selectImage  = 1"
-> An image will open. 
-> Select the aster centroids by clicking on it
-> Once all the centroids are selected. Press ENTER
-> Close the Image window.
-> The selected coordinates will be saved in the source folder
-> The code will do the tessellation and files, figures saved in the respective folders: output_files and output_fig

In option 2:  " selectImage  = 0"
-> Ensure that the coordinates of aster centroids are saved in the source folder for the image that you wish to tessellate 
-> You have to do nothing - let the code tessellate - check the output files in the respective folders.



b)Option 3: Demo of simulation data:
	# ------------------ Simulation demo
	exptdata	         = 0    	# For Experiment =1 ; for simulation = 0
	selectImage              = 1        # Selection of aster centroids: From Image = 1 ; From files = 0
	px2micron                = 120/580.0              # if Sim
	'''


# ----------------------------------------------------------------------------------------------------------------------------------



Currently,: px2micron parameter needs to be changed by the user for each image and filename

