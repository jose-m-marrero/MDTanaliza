# MDTanaliza

MDTanaliza is an open source to calculate gravitational flow paths, modify a Digital Elevation Model and calculate aspect and slope using different approaches.

MDTanaliza can be run both as a graphical user interface and as a command line tool:

In command line
$ ./MDTanaliza /path/configfile.cfg
GUI
$ ./MDTa_interface.py

To Compile the MDTanilza.c program:
gcc MDTanaliza.c -o MDTanaliza -lm

A detailed user manual in PDF is available from the releases page.


# Dependencies
The GUI uses Python 2.7 and Tkinter 

# Current Version
2.1-2018-04-06

# Historical versions

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
