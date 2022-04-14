from pylab import plot, ginput, show, axis
import sys
import numpy as np
import os
from random import uniform
from scipy.spatial import Voronoi, Delaunay , voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from shapely.geometry import MultiPoint, Point, Polygon
from shapely.ops import polygonize
from scipy.spatial import ConvexHull
from matplotlib.ticker import MaxNLocator
import cv2







#------------------------------------------------
# Area of a triangle
def triarea( a , b , c ):
	triarea = 0.5*abs( np.cross( (b - a ) , ( c - a ) )  )
	return triarea



# calculates the area of voronoi region
def calculate_area( vv , p ):
	dpoints = []
	area = 0
	for v in vv.regions[vv.point_region[p]]:
 		dpoints.append(list(vv.vertices[v])) 		

	tri=Delaunay( np.array(dpoints) )
 	#print tri.simplices
 	
 	# adding up the areas of each triangles
	for simplex in tri.simplices:
		area+=triarea(np.array(dpoints[simplex[0]]),np.array(dpoints[simplex[1]]),np.array(dpoints[simplex[2]]))	
	
	return area

# Returns cartesian coordinates from polar
def pointOnCirl( radius  , theta  , w , h , OffsetAU ):
	radius = radius + OffsetAU
	theta =theta*np.pi/180.0
	return radius*np.cos( theta ) + w/2  , radius*np.sin( theta ) + h/2 



# return true if data lies inside the circle, false if outside
def checkCircle( dp , radius , w , h ):
	return( (( dp[0] - (w/2.0))**2 + (dp[1] - (h/2.0))**2 ) < radius**2 )

# return sq. dist	
def get_sdistance(P1, P2):
	square_distance =  (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1])
	return square_distance

# merge the coordinates based on the distance b/w the.
# Not needed for expt data at the moment 
# needed for simulation data - as 2 asters could be over-lapping 
def merge_coordinates(datapoints, dist_cutoff , OffsetAU ):	
	exclude_datapoints = [];
	distance_metrix = [];

	for i in range(0,len(datapoints)-1):
		for j in range(i+1,len(datapoints)):
			p1 = datapoints[i]
			p2 = datapoints[j]

			sdist = get_sdistance(p1,p2)
			distance_metrix.append(np.sqrt(sdist))

			if sdist <= dist_cutoff*dist_cutoff:
				exclude_datapoints.append(j)

	#plt.hist(distance_metrix )
	#plt.show()

	merged_datapoints=[]
	for i in range(0,len(datapoints)):
		if i not in exclude_datapoints:
			merged_datapoints.append(datapoints[i]);
	return merged_datapoints


