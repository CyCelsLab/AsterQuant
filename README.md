# tesselPacking
The code that relates to work described in Khetan et al. (2021) and a current manuscript in preparation. The code takes in data points and a boundary and outputs tesselations with packing distributions (histograms of polygon frequency).


REFERENCE:
1) Khetan N, Pruliere G, Hebras C, Chenevert J and Athale CA. (2021) Self-organized optimal packing of kinesin-5-driven microtubule asters scales with cell size. J Cell Sci. 134(10):jcs257543.


-----------------------------------------------------------------------------------------------------
Following are the step-wise instructions to execute the program.

**Version and Libraries dependencies** 
Please install the following:
1. Python version 3.7.7
2. Numpy 
3. Scipy
4. shapely 
5. pylab  
6. matplotlib

**Usuage for master branch: post-jcs-2021**
Execute the following in the shell command line
> python AsterQuant_V6.py

**Output files** generated in the following folders:

<output_fig>
         
         This contains: 
         1. display_Voronoi  >> Data points and Voronoi tessellated image
         2. Overlays_Voronoi >> Overlays above with the raw image if required

<output_files>

         Main file-names   
         1. AllEutactest.out       >> euctactic measure 
         2. All_NND.out            >> near neighbor distance (NND) for each Voronoi cell
         3. polygonAreaLength.out  >> polygonality, area and lengths of/from each polygons
      


-----------------------------------------------------------------------------------------------------
The work pursued 2017- 2021 shared in 3 branches as explained below:

1. Branch (master) :post-jcs-2021. Updates on the many different measures for analysis and from different systems,  post JCS paper in 2021.

2. Branch: v3_jcs. Parts of the work published in Khetan N, Pruliere G, Hebras C, Chenevert J and Athale CA. (2021) Self-organized optimal packing of kinesin- 5-driven microtubule asters scales with cell size. J Cell Sci. 134(10):jcs257543

3. Branch: AsterQuant. Development until 2021, JCS 
------------------------------------------------------------------------------------------------------
Please refer to the file: README_AsterQuant.md in this depository for further details.
