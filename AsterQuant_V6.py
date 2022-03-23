'''
Neha Khetan, 
AsterQuant: 2017 - 2022
@ Athale Lab, IISER Pune

Platform to measure aster statistics  i.e. packing and spatial analysis
Developed for Multi - aster project
Can be extended to multiple systems
'''


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
import cv2
import plot_voronoi as plotdata
import plot_Simvoronoi as plotsimdata


import compute_functions


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
def pointOnCirl( radius  , theta  , w , h ):
	radius = radius + OffsetAU
	theta =theta*np.pi/180.0
	return radius*np.cos( theta ) + w/2  , radius*np.sin( theta ) + h/2 



# return true if data lies inside the circle, false if outside
def checkCircle( dp , radius ):
	return( (( dp[0] - (w/2.0))**2 + (dp[1] - (h/2.0))**2 ) < radius**2 )

# return sq. dist	
def get_sdistance(P1, P2):
	square_distance =  (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1])
	return square_distance

# merge the coordinates based on the distance b/w the.
# Not needed for expt data at the moment 
# needed for simulation data - as 2 asters could be over-lapping 
def merge_coordinates(datapoints, dist_cutoff):	
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



#############################################################################################################
# 					                       Specify PATHs, filenames, parameters , criterion
#############################################################################################################


IPATH_Image         = './expt/' #'./Output_DATA/redo_thickLines/' #../rsk10uM_2021/' #NH4Cl/new/' # './../rsk-dose/1uM/' #'./10uM_RSk/' #'./image/'
IPATH                   = './expt/'
filenameImage      = 'cCelegans2018x' # 'ORIT_Celegans2018x' #'ORIT_Celegans2018x'                    #'MAX_Projet_.lif - Series023grt1_cf_0.1851'

OPATH_fig		= './output_fig/' # 
OPATH_files		= './output_files/' # 
outfilename                     = 'test' #'1 microM inh RSK.lif - Image011_CF_0.1804'    # name by which the o/p files should be prefixed or suffixed 

px2micron                       = 1 # 120.0/150.0 #0.284   # 1 px corresponds to um?# work in progress to read metadata here n infer as in matlab - for now - its an input param for user

OffsetAU	                  =  0      # not used now, exists due to historical reasons - pls ignore
typeimage                       ='expt'   # work in progress, - pls ignore for now

vizpt                                = 1  # bool for visual representation of voronoi plot. # 0: minimal bp represented 1: additional points represented but not considered for stats
NumberBoundaryPointsViz = 25
overlay_plot   		     = 1
cuttoff_merge_expt           = 3
exptdata  			     = 0

selectImage 		            = 0  # Boolean
FnameAsterCoor               =  'corr_T_Celegans2018x_scen.out'          # DEFAULT FILEname: 'selected_points.out'

readBoundPoints               = 1 # 1 for reading from file-modified for ;  0 for constructing the boundary points
FnamBoundPoints              = 'corr_T_Celegans2018x_1px_roi.txt'
extn                                 = '.png'




#############################################################################################################
# 					                                            Read in the image & extract aster centroids 
img = cv2.imread( IPATH_Image + filenameImage + extn , 0 ) 
#img  = cv2.imread( IPATH_Image + 'Overlays_Voronoi_1 microM inh RSK.lif - Image009_CF_0.1875_nRun_1' + '.pngf' , 0 )
h = img.shape[0]
w = img.shape[1]
#print( h , w )
# mean half dimension along length n width
CellRad  =  (( h/2.0 + w/2.0 )/2.0 )  
NumberBoundaryPoints  = 35






# Choose the aster centers
if selectImage:
	print("Select points or to read from a file set selectImage to 0 and provide the filename\n" )
	f1 = plt.figure(1), plt.imshow( img , origin = 'lower' )

	# aster centroids interactively chosen by the user.
	# After choosing the points - Press enter and close the image window for the code to continue
	data_points = ginput( -1, show_clicks = True , timeout = 0 ) 
	plt.show()
	plt.close()

	wlines=""
	for dp in data_points:
		wlines+=str(dp[0])+'\t'+str(dp[1])+'\n'
	outf = open( IPATH_Image + 'selected_points.out','w')
	outf.writelines(wlines)
	outf.close()


else:
	
	inf = open( IPATH_Image + FnameAsterCoor , 'r' ) #
	lines = inf.readlines()
	inf.close()

	data_points = []
	for l in lines:
		toks = list(map(float,l.strip().split('\t') ) )
		data_points.append(tuple(toks))
	#print(data_points)

'''
#img = img[::-1]  # inversions required if - image & coordinates need to be overlaid
f1 = plt.figure(1), plt.imshow( img )

# aster centroids interactively chosen by the user.
# After choosing the points - Press enter and close the image window for the code to continue
data_points = ginput( -1, show_clicks = True , timeout = 0 ) 
plt.show()
plt.close()
'''




# not needed in the Current experimental data set as the centroids are manually choosen 
# needed for simulation data - as the coordinates of two asters can be same
data_points =  merge_coordinates( data_points , cuttoff_merge_expt )

print("#################################################")




if readBoundPoints:

	''' 
	6-7-19: 
	Reading in boundary points for Yash data on C elegang
	'''
	boundary_points = []
	f0    =  open( IPATH + FnamBoundPoints  , 'r' )    
	infile0  =  f0.readlines()[0: ] 
	
	xb = []
	yb = []
	for k in range( len( infile0 ) ):
		val = infile0[k].strip().split("	") 
		xb.append( float( val[0])  )
		yb.append( float( val[1]) )  
	boundary_points = zip( xb, yb )
	#print("read bps",  list(boundary_points) )
	boundpt = np.array( list(boundary_points ))

	fig1 = plt.figure()
	plt.plot( boundpt[:,0], boundpt[:,1] , 'ko' , ms=4 )
	plt.savefig( OPATH_fig + 'BoundaryPoints_' + '.pdf', bbox_inches='tight' , dpi=600 )
	method = 4
	
	


else:
	print("Optimizing the number of boundary points required")

	#############################################################################################################
	# 					                        Optimizing the number of boundary points 
	#############################################################################################################


	method = 3   
	# 2: Considering expanison of the offset region as long as the packing fraction of all the voronois are within an error
	#    the manual approach is automated. points are initialized between inner and outer circles , and that choice of distance is determined arbitarily.
	#    while the number of points is optimized. IGNORE THIS FOR NOW. HISTORY
	# 3: Just ensures that voronoi vertices are not inside the cell periphery. The infinit points are accounted in the polygonal calculation
	#    by truncating them at the cell boundary (not implemented while plotting instead for visual representation a higher number of points
	#    are initialized so that infinity points are brought close to the cell periphery)
	# ignore the nomencalture -- historical reasons
	# work under progress for getting rid of history and have a working script ( to do )

	 
	if method == 2 :
		OuterRad = CellRad 
		max_nbp = 0; 
		max_counter = 0;
		while max_nbp <= 1:
			
			
			print("Generating boundary points for outer radius:" , OuterRad)
			for nbp in range(  2  , NumberBoundaryPoints ):
				#for nbp in [NumberBoundaryPoints]:
				#print "in loop:", nbp
				#  Generating boundary points on a circle , if required, 
				anglesbound = np.linspace( 0 , 360 , nbp )

				boundary_points = []
				outer_bp        = []
				for ii in range( 0 , len( anglesbound) ):		
					xval , yval = pointOnCirl(  CellRad , anglesbound[ii] , w , h  )
					boundary_points.append( [ xval , yval ]  )
					xvalx , yvaly = pointOnCirl(  OuterRad , anglesbound[ii] , w , h  )
					outer_bp.append( [ xvalx , yvaly ]  )

				# ---------------------------------------------------------------------
				counter      =  0
				points       =  data_points    +  boundary_points 
				vortmp       =  Voronoi( points )
				#----------------------------
				for k in vortmp.vertices:
					chk1 = checkCircle( k , CellRad ) 
					chk2 = checkCircle( k , OuterRad )

					#print chk1,chk2;
					#if (((chk1 and chk2)==False) and (((chk1 | chk2)==True)) ):

					if ( (chk1 ==False) and (chk2 ==True) ):
						counter+=1;
				
				if ( (counter> max_counter) and (counter - nbp== 0 )):
					max_counter = counter
					max_nbp     = nbp
				#print "!",  max_nbp 


			OuterRad = OuterRad + 1
			#print max_nbp



	else:

		# method 3
		OuterRad = CellRad
		max_nbp = 0; max_counter = 0
		print("Generating boundary points for outer radius:" , OuterRad)
		for nbp in range(  1  , NumberBoundaryPoints ):
		#for nbp in [NumberBoundaryPoints]:
			#print "in loop:", nbp
			#  Generating boundary points on a circle , if required, 
			anglesbound = np.linspace( 0 , 360 , nbp )

			# Generating "nbp" boundary points 
			boundary_points = []
			outer_bp        = []
			for ii in range( 0 , len( anglesbound) ):		
				xval , yval = pointOnCirl(  CellRad , anglesbound[ii] , w , h  )
				boundary_points.append( [ xval , yval ]  )
				xvalx , yvaly = pointOnCirl(  OuterRad , anglesbound[ii] , w , h  )
				outer_bp.append( [ xvalx , yvaly ]  )

			# ---------------------------------------------------------------------
			counter      =  0
			points       =  data_points    +  boundary_points 
			vortmp       =  Voronoi( points )  # Voronoi tessellation done here
			#----------------------------
			
			
			# this checks how many polygon vertices lie outside the cirlce
			for k in vortmp.vertices:
				chk1 = checkCircle( k , CellRad ) 
				
				if ( (chk1 ==False) ):
					counter+=1;
			
			# iteratively, checks whether the points that lie outside is larger than the earlier or not.
			# Idea is to select for the Minimal number of input boundary point for which the maximal number of polygon 
			# vertices are ON the circumference or outside the cell - so that aster centered cells account for the boundary
			if ( (counter > max_counter) ):
				#print (counter - nbp ) , counter
				max_counter = counter
				max_nbp     = nbp
				
			print( "!",  max_nbp , nbp , counter)


	# -------------------------------------------------------------	
	print("No. of points at the boundary:" ,  max_nbp )
	print("-----------------------------------------")
	print("Tesselations with optimized boundary points.")
	print(" Multiple runs for randomized points")




#############################################################################################################
# 					        Statistics on tesellating with the optimized number of boundary points
#############################################################################################################

print("Tessellating and statistics ....")
nRuns = 1
polygonArea = []
polygonSides = []
AllLengths   = []
AllEutac     = []

# Multiple runs with the same number of boundary points which are equidistant to each other - 
# however the starting point is randomized - just to ensure the result for the peripheral asters is 
# not due to the position of the boundary points that are introduced 


for nRuns in range( 0, nRuns ):
	tmp =[]
	tmpArea =[]
	output =[]
	tot_area     = 0
	LengthPolygon = []
	
	
	
	
	if readBoundPoints :
		inner_bp = boundpt
		method = 4
		
		
	else:

		#  Generating boundary points on a circle 
		angles = np.linspace( 0 , 360 , max_nbp  )
		angle_diff = round(angles[2] - angles[1])
		new_angles = angles + angle_diff*np.random.uniform( 0 , 1 )
		angles = new_angles

			
		outer_bp = [] # not reqd here. exists cuz of historical reason
		inner_bp = []
		for ii in range( 0 , len( angles) ):		
			xval , yval = pointOnCirl(  CellRad , angles[ii] , w , h  )
			inner_bp.append( [ xval , yval ]  )
			xvalx , yvaly = pointOnCirl(  OuterRad , angles[ii] , w , h  )
			outer_bp.append( [ xvalx , yvaly ]  )



	# ---------------------------------------------------------------------
	

	points       =  data_points    +  inner_bp.tolist()
	vor           =  Voronoi( points )

	# Iterating over the input points
	for i , p in enumerate( vor.points ):
		#print 'i:=' +str(i) + ' \t indx of vor region:' +  str(  vor.point_region[i] )
		# i is an index of the input point
		# j = vor.point_region[i] gives the index of the region created by ith input point
		# vor.regions[j ] gives the jth region among all the regions. 
		# ----------    if a region is unbound it will have -1 in the list among other elements. 
				
		region = vor.regions[ vor.point_region[i] ]
		area_region = 0
		

		polyLength = []
		# For area calculation: 		
		if -1 not in region:              # if it has infinite bound, region will have a '-1' , Ignore those
			# calculate lengths
			LenCounter = 0
			polygonverticesIndx = vor.regions[vor.point_region[i]]				
			# increment the polygonality if at the boundary
					

			flagLen = 0 # All points inside 
			for k in vor.vertices[polygonverticesIndx]:
				chk1 = checkCircle( k , CellRad ) 
				
				if ( (chk1 ==False) ):
					LenCounter=1
					flagLen = 1
					break

			if (flagLen == 0):
						
				# Lengths of polygon sides
				vertexReg = vor.vertices[region]
				for ivR in range(-1,len(vertexReg)-1):
					ivR1 = vertexReg[ivR]
					ivR2 = vertexReg[ivR + 1 ]
					sqdist =get_sdistance( ivR1 , ivR2 )
					polyLength.append( round(np.sqrt(sqdist)*px2micron ,3 ) )

				# Area of each polygon
				parea =  calculate_area( vor , i )
				tot_area+= parea*px2micron*px2micron			 # calculate total area 
				area_region+=parea*px2micron*px2micron           # calculate area of each region
				#print  'ind pt:', i , 'point:', p , 'N poly:' , len(vor.regions[vor.point_region[i]]) , 'Area of reg:', round( area_region , 2 )
					



			output.append( [ p , len(polygonverticesIndx)+LenCounter , round( area_region , 3 ) ] )
			tmpArea.append( round( area_region , 3 ) )
			tmp.append( ( len(polygonverticesIndx)+LenCounter) )
			LengthPolygon.append( polyLength  )


	# Storing the values of the voronoi areas

	# ------- Output the coords , Polygon , Area of each region
	npol_area = []
	AreaCell = 0
	for j , val in enumerate( output ):	
		npol_area.append( [  val[0][0] , val[0][1] , val[1] , val[2] ])
		AreaCell+=val[2]


	# ===============================================================================
	#  Calculate the nearest neighbour distances
	# for each voronoi point find its region then find all its neighboring regions
	# =================================================================================

	near_neigh = {}
	for p1 in range( 0 , len( vor.points) ):
		near_neigh[p1] = []
		for p2 in range( 0 , len( vor.points)):

			if p1 == p2 :
				continue

			r1 = vor.regions[ vor.point_region[ p1 ] ]
			r2 = vor.regions[ vor.point_region[ p2 ] ]    #vor.regions[ p2 ]
				
			# consider only bounded regions, ignore the unbounded ones
			if ( -1 not in r1 ) & ( -1 not in r2 ): 
					 
				commonEdge = list( set(r1) & set(r2)) 
				if len( commonEdge ) == 0:
					continue

				near_neigh[p1].append(p2)
				#print p1 , p2 ,"\n" , "........" , r1 ,  r2
				#print "Inter:" , "...." , list( set(r1) & set(r2))
				# print "whr am i???"

	near_neigh_distance = {}
	for p1 , v in near_neigh.items():

		if len(v) == 0:
			continue

		#print p1, len(v)
		near_neigh_distance[p1] = [[],0]
		coor1 =  vor.points[ p1 ] 
		for p2 in v:
			coor2=vor.points[ p2 ]
			dist = np.sqrt(  ( coor1[0] - coor2[0] )**2 + ( coor1[1] - coor2[1] )**2 )
			near_neigh_distance[p1][0].append( round( dist , 2 )*px2micron )
			
		mean= round( np.mean( np.array(near_neigh_distance[p1][0]) ) , 2 )
		near_neigh_distance[p1][1]=mean
		#print near_neigh_distance[p1]


	near_neigh_distance = {}
	near_neighborPoints = {}
	for p1 , v in near_neigh.items():
		if len(v) == 0:
			continue

		#print p1, len(v)
		near_neigh_distance[p1] = [[],0]
		near_neighborPoints[p1] = [] 
		coor1 =  vor.points[ p1 ] 
		for p2 in v:
			coor2=vor.points[ p2 ]
			dist = np.sqrt(  ( coor1[0] - coor2[0] )**2 + ( coor1[1] - coor2[1] )**2 )
			near_neigh_distance[p1][0].append( round( dist , 2 )*px2micron )
			near_neighborPoints[p1].append( p2 )
					
		mean= round( np.mean( np.array(near_neigh_distance[p1][0]) ) , 2 )
		near_neigh_distance[p1][1]=mean
		#print near_neigh_distance[p1]


	## Eutactic asters
	for p1 ,n1 in near_neighborPoints.items():
		Mtmat = []
		for p2 in n1:
			coord1 = vor.points[p1]
			coord2 = vor.points[p2]
			Mtmat.append(coord2 - coord1)
				
		Mtmat = np.array( Mtmat)
		Mtmat_t = np.transpose(Mtmat)
		MatS   =  Mtmat @ Mtmat_t
		MatT   = np.trace( MatS)
		eutac  = MatT/(( np.sqrt( np.trace( MatS @ MatS)) * np.sqrt(2)))
		AllEutac.append( np.round(eutac,2))



















	# ==============================================================================
	#                                         PLOTTING
	# print 'hi', CellRad
	# ==============================================================================

	if vizpt:
		bp4viz = NumberBoundaryPointsViz
	else:
		bp4viz = max_nbp

	
	#plotdata.plotvoronoiOutput( data_points , vor , CellRad, inner_bp, outer_bp , nRuns , outfilename , OPATH_fig  , method  , bp4viz , w , h , OffsetAU , typeimage, img )
	print( "Generating the plots....")
	plotdata.plotvoronoiOutput( data_points , vor , CellRad, nRuns , outfilename , OPATH_fig  , method  , bp4viz , w , h , OffsetAU , typeimage, img , boundpt  )
	inputimg = mpimg.imread( IPATH_Image + filenameImage + extn ) #'.tif')
	plotdata.plotvoronoiOverlay( vor , data_points , CellRad , 1 ,  outfilename , OPATH_fig  , bp4viz , w , h , OffsetAU , typeimage, inputimg , method, boundpt ) 	

	if method == 4:
		plotdata.plotvoronoiRaw( vor ,  CellRad , nRuns , boundpt , outfilename , OPATH_fig ,  w , h , OffsetAU ,  inputimg , method  )


	# =====================================================================================
	#
	#   						WRITE INTO THE FILES
	# --------------------------------------------------------------
	

	#np.savetxt(  OPATH_files + 'Phallusia_' +  str( "%s" % outfilename ) + '_nRun_'+ str( "%s" %  nRuns ) + '.out' ,  npol_area  , fmt='%.6f'  ,  delimiter='\t') 
		
	# Near neigbor distances: in um
	outfileNND = 'All_NND' + str( "%s" % outfilename ) + '_nRun_'+ str( "%s" %  nRuns ) + '.out'
	fo = open( OPATH_files +  outfileNND , "w")
	for k, v in near_neigh_distance.items():
		fo.write('\t'.join(map(str, v[0])) + "\n" )
	fo.close()



	outfileEutac = 'AllEutac' + str( "%s" % outfilename ) + '_nRun_'+ str( "%s" %  nRuns ) + '.out'
	#print( AllEutac )
	np.savetxt(  OPATH_files + outfileEutac + '.out' ,  np.array(AllEutac)  , fmt='%.2f'  ,  delimiter='\t') 			


	#print "packing fraction"
	totalCellArea = (CellRad*CellRad*22.0)/7.0
	#print  nRuns , round( ( (AreaCell)/( totalCellArea) )*100 , 2 )

	polygonSides.append( tmp )
	polygonArea.append( tmpArea )
	AllLengths.append( (LengthPolygon) )


	
#print (polygonSides)


# ===================================================================================================================================================================
print("Writing into the files......")
# each row corresponds to one run while asters are in columns

wlines=''
for v , va , vs in zip(AllLengths, polygonArea, polygonSides ):
	for v1 , v1a , vs1 in zip(v,va , vs ) :
		# Polygonality | Area | Lengths of each side
		wlines+= str( vs1 ) + '\t' +  str( v1a ) +  '\t' + '\t'.join(map(str,v1))+'\n'

outf=open(OPATH_files + 'polygonAreaLengths' +  str( "%s" % outfilename ) + '.out','w')
outf.writelines(wlines)
outf.close()




