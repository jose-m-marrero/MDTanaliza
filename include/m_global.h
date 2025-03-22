/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza  
* Module:  Global_flow.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2022-01-18
* Last update: 2025-01-14
* 
* Description:
* Global variables used by modules of MDTanaliza
* *********************************************************************/
#ifndef _GLOB_F
#define _GLOB_F

/*!
* DEFINES
***********************************************************************/ 
#define PI 3.14159265
#define MAX_PATH     260                                                /*!< String dimension value for some name paths */
#define NUMARR         5                                                /*!< total matrix (raster) used simultaneously by algorithms */ 
#define FIFOSIZE  100000                                                /*!< FIFO size in multiflow model */

#define COLXYZ         4                                                /*!< Expected columns in file init point, format id,x,y,z for DEM modifications */
#define COLVENT        3                                                /*!< Expected columns in file init point, format id,x,y for grav. flow paths */
#define COLINTER       5                                                /*!< Expected columns in file init point, format id,code,x,y,z for interpolation */
#define COLBED         8                                                /*!< Expected columns in file zero level point */

#define TOTXYZ  10000000                                                /*!< total init xyz point allowed */
#define TOTZLP    500000                                                /*!< total Zero Point Level allowed in each Zero Level Trajectory */
#define TOTZLTR      100                                                /*!< total Zero Point Level Trajectories allowed */
#define NUMCROSS   10000                                                /*!< total Cross section allowed */
#define INITVAL  -999999

#endif /* _GLOB_F */


