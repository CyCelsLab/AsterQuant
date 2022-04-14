#####################################################################################
#   Copyright (C) 2020 Neha Khetan @ khetanneha@gmail.com
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
####################################################################################


import sys
import os
from random import uniform
import numpy as np
import cv2
from scipy.spatial import Voronoi, Delaunay , voronoi_plot_2d
from scipy.spatial import ConvexHull
from shapely.geometry import MultiPoint, Point, Polygon
from shapely.ops import polygonize
from matplotlib.ticker import MaxNLocator
from pylab import plot, ginput, show, axis
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sn
sn.set_style('white')
sn.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2.5})
plt.rcParams['figure.figsize'] = [12, 6]


def plot_AllQdata( polygonSides ,  AllnndVal , AllLengths , polygonArea  , AllEutac , scale , Opath , Opath_fig, outfilename , nR ):

	if scale =='pixel':
		unt  = ' (px)'
		unt2 = ' (px'+ "$^{2}$"+ ')' 
	elif scale =='micron':
		unt  = " ($\mu$m)"
		unt2 = " $\mu$m" + "$^{2}$"+ ')' 
	else:
		print("Check the scale | accepted input is as pixel or micron")


	# Polygon HISTOGRAM ------------------------
	if np.min(polygonSides) < 3:
		minPS = np.min(polygonSides)
	else:
		minPS = 3

	if np.max(polygonSides) > 8:
		maxPS = np.max(polygonSides)
	else:
		maxPS = 8

	nbin =  np.arange( minPS , maxPS+1  , 1 )
	f10 =   plt.figure(1)
	plt.subplot(2,3 ,1 ),sn.histplot(data = polygonSides ,  stat = "probability" , discrete = True  )
	plt.xlabel('Polygonality' , fontsize = 10)
	plt.ylabel('Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	plt.xticks( nbin )
	h1 = plt.legend()
	#plt.show()



	# Cal. Voronoi Entropy n Irregularity: REQD: polygonSides, nbin
	f,b = np.histogram( polygonSides , bins= nbin )
	tmpv = []
	for k in (0 , len(f)-1):
		if f[k] != 0:
			pv = f[k]/sum(f)
			tmpv.append( pv*np.log( pv ) )
	vorEnt =  round( -sum( tmpv ) , 3 ) 
	polygonSides = np.array(polygonSides, dtype = "float" )
	IRR =  round( ( np.std(polygonSides)/np.mean(polygonSides) ) , 3 )





	# All NND
	allnnd = []
	meanNe = []
	for k in AllnndVal:
		allnnd.extend( k )
		meanNe.append( float(np.min(k)) )
	#f4 = plt.figure()
	plt.subplot(2,3 ,2 ),
	sn.histplot( data = allnnd , stat = "probability" )
	plt.xlabel('NND'+ unt , fontsize = 10)
	plt.ylabel('')#'Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	#plt.xticks( np.arange( np.min(allnnd) , np.max(allnnd) , 2 ))
	#plt.legend()
	#plt.show()


	# For circularity: REQD: polygonArea and AllLengths
	circ = []
	regL = []
	al   = []  
	ar   = []
	allen = []
	py    = []
	for pl2 in range(0, len( AllLengths[0] )-1):
		if ( AllLengths[0][pl2] ):
			circ.append( (round( (4*np.pi*polygonArea[0][pl2])/(sum( AllLengths[0][pl2]) )**2 , 3 )) )
			regL.append( np.std(AllLengths[0][pl2])/np.mean(AllLengths[0][pl2]) )
			#al.append( ( polygonArea[0][pl2] , np.std(AllLengths[0][pl2])/np.mean(AllLengths[0][pl2])  ))
			al.append( ( polygonArea[0][pl2] , np.mean(AllLengths[0][pl2])  ))
			ar.append( polygonArea[0][pl2] )
			allen.extend( AllLengths[0][pl2] ) 
			py.append( len( AllLengths[0][pl2] ))

	plt.subplot(2,3 ,3 ),
	sn.histplot( data = ar , stat = "probability" )
	plt.xlabel('Area' + unt2 , fontsize = 10)
	plt.ylabel(' ') #Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)

			
	#f10 = plt.figure()
	plt.subplot(2,3 ,4 ),
	sn.histplot(data = allen ,  stat = "probability" , discrete = False  )
	#sn.boxplot(data = regL )
	#plt.xlim( [ 0 , 1])
	plt.xlabel('Length' + unt , fontsize = 10)
	plt.ylabel('Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	#plt.show()


	plt.subplot(2,3 ,5 ),
	sn.histplot(data = circ ,  stat = "probability" , discrete = False  )
	plt.xlim( [ 0 , 1])
	plt.xlabel('Circularity' , fontsize = 10)
	plt.ylabel(' ') #Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	#plt.show()

	# Eutacticity
	plt.subplot(2,3 , 6 ),
	sn.histplot(data = AllEutac ,  stat = "probability" , discrete = False  )
	plt.xlabel('Eutacticity' , fontsize = 10)
	plt.ylabel(' ') #Frequency' , fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	#plt.xticks( [0,0.5, 1] )
	#plt.show()
	plt.savefig( Opath_fig + 'VorStats_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nR ) + '.pdf'  , dpi=300 )


	f12 = plt.figure()
	plt.subplot( 2 ,3 , 1)
	al = np.array( al )
	sn.scatterplot(x = al[:,0] , y = al[:,1] )
	plt.xlabel('Area' + unt2, fontsize = 10) 
	plt.ylabel('<L>' + unt , fontsize = 10) 
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	#plt.show()

	plt.subplot( 2 ,3 , 2)
	al = np.array( al )
	py = np.array( py )
	sn.scatterplot(x = py , y = al[:,1] )
	plt.xlabel('polygonality', fontsize = 10) 
	plt.ylabel('<L>' + unt, fontsize = 10) 
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	plt.xticks( nbin )

	plt.subplot( 2 ,3 , 3)
	sn.scatterplot(x = py , y = al[:,0] )
	plt.xlabel('polygonality', fontsize = 10) 
	plt.ylabel('Area' + unt2 , fontsize = 10) 
	plt.yticks(fontsize = 10)
	plt.xticks(fontsize = 10)
	plt.xticks( nbin )
	
	plt.savefig( Opath_fig + 'VorStatsCorrelations_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nR ) + '.pdf',  dpi=300 )
	plt.show()
	#plt.close()

	# ---
	ra  = np.sum( meanNe )/len(AllnndVal)
	re  = 0.5*( 1/ np.sqrt(len(AllnndVal)/np.sum(ar)))
	
	# ---
	'''
	print("\n\n ------------------------ \n\n")
	print("Voronoi Entropy :"   , vorEnt )
	print("Irregularity score :", IRR )
	print("<C.V.> in length :", round( np.mean(regL) , 3 ) )
	print( "Measure of random:" ,   round( ra/re , 3 ) )
    '''
    
	import pandas as pd
	outMeasures  = { 'VorEntro': [vorEnt], 'Irregularity': [IRR] , 'CV_Length': round( np.mean(regL) , 3 ) , 'Randomness': round( ra/re , 3 )}
	df  = pd.DataFrame( data = outMeasures )
	df.to_csv( Opath + 'measures.' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nR ) + '.csv', sep = '\t')

