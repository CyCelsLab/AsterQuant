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





def plotvoronoiOutput( data_points, vor , CellRad, inner_bp, outer_bp , nRuns , outfilename , OPATH_fig , method , nb  , w , h , OffsetAU , data_type , inputimg ):


	if data_type == 'expt':
		CellRad = (( h/2.0 + w/2.0 )/2.0 ) 

	#overlay_plot = 1
	#print "hi", data_type, CellRad


	if method ==3:

		anglesbound = np.linspace( 0 , 360 , nb )
		boundary_points = []
	
		for ii in range( 0 , len( anglesbound) ):		
			xval , yval = pointOnCirl(  CellRad , anglesbound[ii] , w , h , OffsetAU )
			boundary_points.append( [ xval , yval ]  )
			
		# ---------------------------------------------------------------------
		points       =  data_points    +  boundary_points 
		vor          =  Voronoi( points )
		#----------------------------

		vorpoints = np.array( vor.points )
		hull   = ConvexHull( points + inner_bp )
		points = np.array( points+ inner_bp)
		hull_xpts =  np.hstack( ( points[hull.vertices,0] , points[hull.vertices[0],0] ) )
		hull_ypts =  np.hstack( ( points[hull.vertices,1] , points[hull.vertices[0],1] ) )
		

		plt.plot( vorpoints[:,0], vorpoints[:, 1], 'go', ms=4)
		plt.plot( np.array(inner_bp)[:,0], np.array(inner_bp)[:, 1], 'yo', ms=3)
		plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'bo', ms= 2 )
		plt.plot( hull_xpts , hull_ypts , 'g' , linewidth = 2 )

		vert_subset = [];

		for vpair in vor.ridge_vertices:
			if vpair[0] >= 0 and vpair[1] >= 0:
				v0 = vor.vertices[vpair[0]]
				v1 = vor.vertices[vpair[1]]

			   	
				chk1 = checkCircle( v0 , CellRad ,w , h ) 
				chk2 = checkCircle( v1 , CellRad ,w , h  ) 
				if (chk1 and chk2)==True:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
			   	
				else:
					plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k-.', linewidth = 1 )	


		ax = plt.gca()
		ax.set_aspect('equal', adjustable='box')
		plt.savefig( OPATH_fig + 'display_Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=600 )


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
				chk1 = checkCircle( v0 , CellRad , w , h ) 
				chk2 = checkCircle( v1 , CellRad , w , h ) 


			   
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

		plt.savefig( OPATH_fig + 'Voronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=600 )

	



def plotvoronoiOverlay( data_points, vor , CellRad, inner_bp, outer_bp , nRuns , outfilename , OPATH_fig , method , nb  , w , h , OffsetAU , data_type , inputimg ):
	'''
	plt.clf()
	f2 = plt.figure(2) 
	inputimg = inputimg[::-1]         # inversions required if - image & coordinates need to be overlaid
	plt.imshow( inputimg )
	vorpoints = np.array( vor.points )
	plt.plot( vorpoints[:,0], vorpoints[:, 1], 'yo', ms=2)
	innerbp  = np.array( inner_bp )
	
	vert_subset = [];

	for vpair in vor.ridge_vertices:
		if vpair[0] >= 0 and vpair[1] >= 0:
			v0 = vor.vertices[vpair[0]]
			v1 = vor.vertices[vpair[1]]
			######################################
			chk1 = checkCircle( v0 , CellRad , w , h ) 
			chk2 = checkCircle( v1 , CellRad , w , h ) 			   
			if (chk1 and chk2)==True:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'w', linewidth = 1 )
				   	
			else:
				plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'w-.', linewidth = 1 )
					
	

	plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'co', ms= 1 )
	# for the boundary
	hull = ConvexHull( data_points + inner_bp )
	points = np.array( data_points+ inner_bp)
	hull_xpts =  np.hstack( ( points[hull.vertices,0] , points[hull.vertices[0],0] ) )
	hull_ypts =  np.hstack( ( points[hull.vertices,1] , points[hull.vertices[0],1] ) )
	#plt.plot(, np.hstack(  (points[hull.vertices,1], points[hull.vertices[0],1])    ), 'g' , linewidth = 2 )
	plt.plot( hull_xpts , hull_ypts , 'w-.' , linewidth = 1 )
	plt.savefig( OPATH_fig + 'OverlayVoronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=600 )
	'''
	inputimg = inputimg[::-1]         # inversions required if - image & coordinates need to be overlaid
	plt.imshow( inputimg )
	data_points = np.array( data_points )
	plt.plot( data_points[:,0], data_points[:, 1], 'bo', ms=4)

	vorpoints = np.array( vor.points )
	plt.plot( vorpoints[:,0], vorpoints[:, 1], 'yo', ms=2)
	innerbp  = np.array( inner_bp )

	for vpair in vor.ridge_vertices:
	    if vpair[0] >= 0 and vpair[1] >= 0:
	        v0 = vor.vertices[vpair[0]]
	        v1 = vor.vertices[vpair[1]]
	        # Draw a line from v0 to v1.
	        plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth = 1 )
	plt.plot(vor.vertices[:,0], vor.vertices[:, 1], 'co', ms= 3 )

	ax = plt.gca()
	if ( w == h ):
		ax.set_aspect('equal', adjustable='box')
	plt.savefig( OPATH_fig + 'OverlayVoronoi_' + str( "%s" % outfilename )+ '_nRun_' + str( "%s" %  nRuns ) + '.tiff', bbox_inches='tight' , dpi=600 )





def checkCircle( dp , radius , w , h  ):
	return( (( dp[0] - (w/2.0))**2 + (dp[1] - (h/2.0))**2 ) < radius**2 )


# Returns cartesian coordinates from polar
def pointOnCirl( radius  , theta  , w , h , OffsetAU ):
	radius = radius + OffsetAU
	theta =theta*np.pi/180.0
	return radius*np.cos( theta ) + w/2  , radius*np.sin( theta ) + h/2 