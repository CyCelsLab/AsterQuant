AsterQuant
-----------------------------------------------------------------------------------------------------

AsterQuant: tool for Quantitative analysis of Multi aster systems using computational geometric methods to measure packing statistics and spatial analysis

Developed for Multi - aster project and has been extended to multiple systems in experiments

Primarily Developed by Neha Khetan, <khetanneha@gmail.com>, <neha.khetan@students.iiserpune.ac.in> for the work published as 
_Khetan N, Pruliere G, Hebras C, Chenevert J and Athale CA. (2021) Self-organized optimal packing of kinesin- 5-driven microtubule asters scales with cell size. J Cell Sci. 134(10):jcs257543_

Further improvements, additional methods, measures and analysis added for the developement of the tool.

-----------------------------------------------------------------------------------------------------
The work pursued 2017- 2021 shared in 3 branches 

1. Branch (master) :post-jcs-2021. Updates on the many different measures for analysis and from different systems,  post JCS paper in 2021.

2. Branch: v3_jcs. Parts of the work published in Khetan N, Pruliere G, Hebras C, Chenevert J and Athale CA. (2021) Self-organized optimal packing of kinesin- 5-driven microtubule asters scales with cell size. J Cell Sci. 134(10):jcs257543

3. Branch: AsterQuant. Development until 2021, JCS 


-----------------------------------------------------------------------------------------------------
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
      


**Scripts**
1. AsterQuant_V6.py
2. compute_functions.py
3. plot_Simvoronoi.py
4. plot_voronoi.py
5. coordinate_merger.py


Folder structure
1. expt                         // These contain input files and images
2. sim                          // These contain input files and images
3. output_fig and               // to store the outputs from a run - IMAGES
4. output_files                 // to store the outputs from a run - MEASURES and ANALYSED DATA
   ( Representative output files and images from phallusia, nematodes and simulations included )


