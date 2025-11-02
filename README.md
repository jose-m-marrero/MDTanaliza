# MDTanaliza

MDTanaliza is an open-source application focused on raster data analysis. It was initially designed for natural hazard assessment, though it can be applied to a broader range of uses. In this new version, no GUI is available, so it must be operated from the command line.

The implemented algorithms, internally referred to as analysis strategies, are as follows:

- 1: Modify DEM
- 2: Calculate Sinks
- 3: Calculate Aspect (2 methods)
- 4: Calculate Slope (9 methods)
- 5: Gravitational Flows (4 methods)
- 6: Topohazard
- 7: Multiflow (2 methods)
- 8: Hazard (not yet implemented)
- 9: Interpolation (IDW only)

For the first five strategies, the available reference is:
Marrero, J.; Vasconez, F.; Espín, P.; Ortiz, R.; Yepes, H.; García, A.; Mothes, P. & Estrella, C. (2019). MDTanaliza: understanding digital elevation models when facing gravity-driven flows in a hazard assessment context. Earth Science Informatics, 27, 317-328. https://doi.org/10.1007/s41324-018-0230-y

For the Topohazard strategy, the reference is:
Marrero, J.M; Vasconez, F.J.; Espín-Beón, P.; Sierra, D.; Yepes, H.A. & Mothes, P. (2025). Topohazard, a novel tool for enhancing gravity-driven flows hazard maps: application to Cotopaxi volcano, Ecuador. Natural Hazards, https://doi.org/10.1007/s11069-025-07702-5.

To execute MDTanaliza, use the following command:
bash
$ ./MDTanaliza /path/configfile.cfg
To obtain a clean config file
bash
$ ./MDTanaliza N (where N is the number of the strategy)

To compile the MDTanaliza.c program, refer to the comments in the MDTanaliza.c source file.

# Dependencies
MDTanaliza uses a combination of standard libraries along with some specific libraries developed by the author (located in the lib folder).
Algorithms and analysis strategies are organized in the include folder. 

# Current Version
3.0.2-2025-09-20
Some bugs affecting strategies 5 and 6 in the algorithms have been fixed. Strategy 5 has been improved and optimized. The libraries have also been reviewed.

# Historical versions

*******************************************************
Version 3.0.0-2025-02-19
This new version has been completely reprogrammed, separating the analysis algorithms from the main application into individual files to simplify programming and maintenance. Additionally, new analysis strategies have been incorporated.

*******************************************************
Old version has been discontinued
*******************************************************
Version 2.1-2018-06-04

--GUI--
Improvement, New variables has been added

--MORPHOMETRIC SECTION--

SLOPE-ASPECT
Improvement, add histogram in absolute and relative frequency to facilitate the comparison between different DEMs.

SLOPE-GRADIENT
Improvement, add histogram in frequency and percentage to be compared with other results.

--GRAVITY FLOW PATH SECTION--

Improvement, all flow paths can be generated in a unique global raster or separated in individual rasters using W. Raster Mode parameter (0 join; 1 individual)
Improvement, the output filename has been change including some relevant parameters.
Improvement, a new resume file shows the computing time, distance reached, number of cells affected an the occupied area for each flow and final global raster, the latter if W. Raster Mode is set to 0.

SINGLE FLOW PATH LHM OR SSM
Bug fixed, fill function used incorrect variable to start filling the cell. 
Improvement, a new option have been added to control the interaction or not of different paths during the flow path calculation. This approach uses a control raster, which is reseted or not.
Improvement, the distance to the starting cell is now calculated using coordinates.

DRUNK SAILER MONTECARLO TYPE FLOW PATH
Improvement, a fill function has been added to the algorithm. 
Improvement, the distance to the starting cell is now calculated using coordinates.
Improvement, new variable has been added to control and reduce the repetitions in selected cells.


*******************************************************
Version 2.1-2018-04-06

--DEM PARAMETERS SECTION--

Improvement, DEM Z Max and Min values have been changed to float.

--GRAVITY FLOW PATH SECTION--

SINGLE FLOW PATH LHM OR SSM
Bug fixed, stop flow path calculation if cell is a sink and fill increase parameter is set to 0 
Bug fixed, stop flow path calculation if cell is null

DRUNK SAILER MONTECARLO TYPE FLOW PATH
Bug fixes, stop flow path calculation if cell is null

MULTITRAYECTORY FLOW PATH
Bug fixed, avoid the use of fill increase option.
Bug fixed, Avoid null values in flow path calculation


******************************************************
Version 2.1-2017-10-06

First public version

# Funds
This work has been ostensibly supported in its initial stages by the Secretaría de Educación Superior de Ciencia, Tecnología e Innovación (SENESCYT) of the Government of Ecuador under the PROMETEO Programme (PROMETEO-CEB-004-2015 and PROMETEO-CEB-009-2016)
