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





def plotvoronoiOutput( data_points , CellRad,  nRuns ,   outfilename , OPATH_fig   , bp4viz , w , h , OffsetAU , img ):
	#CellRad  = max(h/2.0 , w/2.0) 
	#CellRad = (( h/2.0 + w/2.0 )/2.0 ) 
	method = 3

	plt.close()
	if method ==3:

		anglesbound = np.linspace( 0 , 360 ,  bp4viz )
		boundary_points = []
	
		for ii in range( 0 , len( anglesbound) ):		
			xval , yval = pointOnCirl(  CellRad , anglesbound[ii] , w , h , OffsetAU )
			#xval , yval = pointOnCirl2(  e_A, e_B , anglesbound[ii] , w , h  )
			boundary_points.append( [ xval , yval ]  )
			
		# ---------------------------------------------------------------------
		points       =  data_points    +  boundary_points 
		vor          =  Voronoi( points )
		#----------------------------

		vorpoints = np.array( vor.points )
		#hull   = ConvexHull( points + inner_bp )
		#points = np.array( points+ inner_bp)
		#hull_xpts =  np.hstack( ( points[hull.vertices,0] , points[hull.vertices[0],0] ) )
		#hull_ypts =  np.hstack( ( points[hull.vertices,1] , points[hull.vertices[0],1] ) )
		

		plt.plot( vorpoints[:,0], vorpoints[:, 1], 'go', ms=4)
		#plt.plot( np.array(inner_bp)[:,0], np.array(inner_bp)[:, 1], 'yo', ms=3)
		plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'bo', ms= 2 )
		#plt.plot( hull_xpts , hull_ypts , 'g' , linewidth = 2 )

		vert_subset = [];

		for vpair in vor.ridge_vertices:
			if vpair[0] >= 0 and vpair[1] >= 0:
				v0 = vor.vertices[vpair[0]]
				v1 = vor.vertices[vpair[1]]

			   	
				chk1 = checkCircle( v0 , CellRad ,w , h , OffsetAU   ) 
				chk2 = checkCircle( v1 , CellRad ,w , h  , OffsetAU ) 
				#chk1 = checkCircle2( v0 , w , h , e_A , e_B )
				#chk2 = checkCircle2( v1 , w , h , e_A , e_B )
				if (chk1 and chk2)==True:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
			   	
				else:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k-.', linewidth = 1 )	


		ax = plt.gca()
		ax.set_aspect('equal', adjustable='box')
		plt.savefig( OPATH_fig + 'display_Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=300 )
		plt.close()




	else:
		vorpoints = np.array( vor.points )
		plt.plot( vorpoints[:,0], vorpoints[:, 1], 'go', ms=4)

		innerbp  = np.array( inner_bp )
		outerbp  = np.array( outer_bp )
		plt.plot( innerbp[:,0], innerbp[:, 1], 'yo', ms=3)
		#plt.plot( outerbp[:,0], outerbp[:, 1], 'ro', ms=3)

		vert_subset = [];

		for vpair in vor.ridge_vertices:
			if vpair[0] >= 0 and vpair[1] >= 0:
				v0 = vor.vertices[vpair[0]]
				v1 = vor.vertices[vpair[1]]
			   	######################################
				chk1 = checkCircle( v0 , CellRad , w , h , OffsetAU ) 
				chk2 = checkCircle( v1 , CellRad , w , h , OffsetAU ) 


			   
				if (chk1 and chk2)==True:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
				   	
				else:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k-.', linewidth = 1 )
					
			   
		plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'bo', ms= 2 )
		
		# for the boundary
		hull = ConvexHull( data_points + inner_bp )
		points = np.array( data_points+ inner_bp)
		hull_xpts =  np.hstack( ( points[hull.vertices,0] , points[hull.vertices[0],0] ) )
		hull_ypts =  np.hstack( ( points[hull.vertices,1] , points[hull.vertices[0],1] ) )
		#plt.plot(, np.hstack(  (points[hull.vertices,1], points[hull.vertices[0],1])    ), 'g' , linewidth = 2 )
		
		plt.plot( hull_xpts , hull_ypts , 'g' , linewidth = 1 )


		ax = plt.gca()
		if ( w == h ):
			ax.set_aspect('equal', adjustable='box')

		plt.savefig( OPATH_fig + 'Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=300 )
		
		plt.clf()
		plt.close()

	



# For Overlay of simulation data:
# Converting the Aster data point from micron scale to pixel scale
def plotvoronoiOverlay( data_points,  CellRad , nRuns ,  outfilename , OPATH_fig ,  nb  , w , h , OffsetAU , data_type , inputimg  ,  px2micron ):

	CellRad  = CellRad*(1/px2micron)
	OffsetAU = OffsetAU*(1/px2micron)
	data_points       = np.array( data_points )
	#print( data_points.shape )
	data_points       = data_points*(1/px2micron)
	data_points[:,0]  = data_points[:,0] + w/2.0
	data_points[:,1]  = data_points[:,1] + h/2.0
	#print("hello!.......", data_points )
	#print("bye.......", data_points )
	anglesbound = np.linspace( 0 , 360 , nb )
	boundary_points = []
	#print( type( data_points))

	for ii in range( 0 , len( anglesbound) ):		
		xval , yval = pointOnCirl(  CellRad , anglesbound[ii] , w , h , OffsetAU )
		#xval , yval = pointOnCirl2(  e_A, e_B , anglesbound[ii] , w , h  )
		boundary_points.append( [ xval , yval ]  )
			
	# ---------------------------------------------------------------------
	points       =  np.vstack( (data_points ,  boundary_points) )
	#print( type( points))
	vor          =  Voronoi( points )
	#----------------------------

	vorpoints   = np.array( vor.points )
	datapoints  = np.array( data_points )
	hull        = ConvexHull( points  )
	#points      = np.array( points + boundary_points )
	#hull_xpts   =  np.hstack( ( points[hull.vertices,0] , points[hull.vertices[0],0] ) )
	#hull_ypts   =  np.hstack( ( points[hull.vertices,1] , points[hull.vertices[0],1] ) )
	

	# -----------------------------------------------------------------------------------------------------------------
	#                 Plotting
	# ------------------------------------------------------------------------------------------------------------------
	fig= plt.figure(1)
	inputimg = inputimg[ ::-1, :]
	plt.imshow( inputimg , origin = 'lower' )
	#datapoints = np.flip( datapoints, axis = 0)
	#datapoints = np.flip( datapoints, axis = 1)
	plt.plot( datapoints[:,0], datapoints[:, 1], 'yo', ms=4)
	#plt.plot( np.array(boundary_points)[:,0], np.array(boundary_points)[:, 1], 'ro', ms=3)
	plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'ro', ms= 3 )
	#plt.plot( hull_xpts , hull_ypts , 'g' , linewidth = 2 )

	vert_subset = [];
	for vpair in vor.ridge_vertices:
		if vpair[0] >= 0 and vpair[1] >= 0:
			v0 = vor.vertices[vpair[0]]
			v1 = vor.vertices[vpair[1]]

			   	
			chk1 = checkCircle( v0 , CellRad ,w , h , OffsetAU ) 
			chk2 = checkCircle( v1 , CellRad ,w , h  , OffsetAU ) 
			#chk1 = checkCircle2( v0 , w , h , e_A , e_B)
			#chk2 = checkCircle2( v1 , w , h , e_A , e_B)
			if (chk1 and chk2)==True:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
			   	
			else:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k-.', linewidth = 1 )	



	ax = plt.gca()
	#ax.set_aspect('equal', adjustable='box')
	#ax.set_facecolor('#210000ff')

	plt.savefig( OPATH_fig + 'Overlays_Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=300 )
	plt.close()
	plt.clf()



def plotvoronoiOverlayRaw( data_points, vor ,  CellRad , nRuns ,  inner_bp , outfilename , OPATH_fig ,  w , h , OffsetAU ,  inputimg   ):


	plt.close()
	vorpoints = np.array( vor.points )
	plt.plot( vorpoints[:,0], vorpoints[:, 1], 'go', ms=4)
	plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'bo', ms= 2 )
	
	vert_subset = [];

	for vpair in vor.ridge_vertices:
		if vpair[0] >= 0 and vpair[1] >= 0:
			v0 = vor.vertices[vpair[0]]
			v1 = vor.vertices[vpair[1]]			   	
			chk1 = checkCircle( v0 , CellRad ,w , h , OffsetAU   ) 
			chk2 = checkCircle( v1 , CellRad ,w , h  , OffsetAU ) 
			if (chk1 and chk2)==True:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
			   	
			else:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k-.', linewidth = 1 )	


	ax = plt.gca()
	ax.set_aspect('equal', adjustable='box')
	

	plt.savefig( OPATH_fig + 'RawOverlays_Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=300 )
	plt.close()



def checkCircle( dp , radius , w , h , OffsetAU   ):
	return( (( dp[0] - (w/2.0))**2 + (dp[1] - (h/2.0))**2 ) < (radius + OffsetAU )**2 )


# Returns cartesian coordinates from polar
def pointOnCirl( radius  , theta  , w , h , OffsetAU ):
	radius = radius + OffsetAU
	theta =theta*np.pi/180.0
	return radius*np.cos( theta ) + w/2  , radius*np.sin( theta ) + h/2 

def pointOnCirl2( A , B , theta , w , h ):
	theta =theta*np.pi/180.0
	return A*np.cos( theta ) + w/2  , B*np.sin( theta ) + h/2 


def checkCircle2( dp , w ,h  , A , B ) :
	return(  ( ( dp[0] -  w/2.0 )**2)/(A**2) +  (( dp[1] - h/2.0)**2)/(B**2)  < 1 ) 
