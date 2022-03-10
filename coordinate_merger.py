import numpy as np

def get_sdistance(P1, P2):
	square_distance =  (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1])
	return square_distance


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

	import matplotlib.pyplot as plt
	#plt.hist(distance_metrix )
	#plt.show()

	merged_datapoints=[]
	for i in range(0,len(datapoints)):
		if i not in exclude_datapoints:
			merged_datapoints.append(datapoints[i]);
	return merged_datapoints


