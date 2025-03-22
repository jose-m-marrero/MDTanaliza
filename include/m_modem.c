/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza 
* Module:  m_modem.sh
* 
* Autores: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2019-10-01
* Last update: 2022-01-19
* 
* Description:
* 
* Modify DEM values following 3 options:
* Only the cell pointed by the new value is changed
* The cell and its 8 surrounded cells, using same values
* (1) The cell and its 8 surrounded cells, using a random value
* (2) the new cell value 
* tip_change is 1, 2 or other value if no 8 cells will change
* 
* *********************************************************************/
#ifndef _MOD_RAS
#define _MOD_RAS

/*!
* MDTanaliza's general modules 
***********************************************************************/ 
/*! Define global variables */
#include "../include/m_global.h"
/*! Define structures */
#include "../include/struc_mdtanaliza.h"
/*!
* General libraries
***********************************************************************/ 
/*! structure to save array's head data from *.grd input raster format*/
#include "struct_head-array.h"
/*! Reading and writing raster functions */
#include "u_arrays.h"
/*! Some geometric calculation and data comparison functions */
#include "u_calculus.h"

/*----------------------------------------------------------------------
 |  Function definition
 ---------------------------------------------------------------------*/
float get_random(float max_interval);


//**********************************************************************
/*! Calc Functions */
//**********************************************************************

/*! Change z values in a input DEM given a set of coordinates
 * In all conditions, same cell is changed by the new Z value
 * id chan_mode != 1 or 2, neighboring cells won't be altered 
 * if chan_mode 1, the eight neighboring cells are changed by same value
 * if chan_mode 2, the eight neighboring cells are changed by random value (0 and the new cell value)*/
double **calc_modem(double **raster,  struct HeadR s_arrayh[], struct coordCSV s_coord[], int tot_xyz, int chan_mode)
{
int i, idxcol, idyrow;
double tx, ty, tz;
	
	
	printf("Modofying Z values from DEM\n");
	for(i=0; i<tot_xyz; i++)
	{
		tx = s_coord[i].ptxcoor;
		ty = s_coord[i].ptycoor;
		tz = s_coord[i].ptzcoor;
		/*!< Check if coordinate is in raster */
		if(tx > s_arrayh[0].hxlo && tx < s_arrayh[0].hxhi && ty > s_arrayh[0].hylo && ty < s_arrayh[0].hyhi) 
		{
			idxcol = calc_Rindex(tx, s_arrayh[0].hxlo, s_arrayh[0].hresx);             /*!< Get Column index */
			idyrow = calc_Rindex(ty, s_arrayh[0].hylo, s_arrayh[0].hresy);             /*!< Get row index */
			raster[idyrow][idxcol] = tz;                                               /*!< Add new z value to the raster */
			/*!< Add same z value to the 8 neighboring cells */
			if (chan_mode == 1)
			{
				//Update perpendicular cells
				raster[idyrow-1][idxcol]   = tz;   //h - 2
				raster[idyrow+1][idxcol]   = tz;   //b - 8
				raster[idyrow][idxcol+1]   = tz;   //f - 4
				raster[idyrow][idxcol-1]   = tz;   //d - 6
				//Update oblique cells
				raster[idyrow+1][idxcol+1] = tz; //c - 7
				raster[idyrow-1][idxcol+1] = tz; //i - 1
				raster[idyrow-1][idxcol-1] = tz; //g - 3
				raster[idyrow+1][idxcol-1] = tz; //a - 9
			}
			/*!< Add random z value to the 8 neighboring cells */
			if (chan_mode == 2)
			{
				//Update perpendicular cells
				raster[idyrow-1][idxcol]   = get_random(tz);   //h - 2
				raster[idyrow+1][idxcol]   = get_random(tz);   //b - 8
				raster[idyrow][idxcol+1]   = get_random(tz);   //f - 4
				raster[idyrow][idxcol-1]   = get_random(tz);   //d - 6
				//Update oblique cells
				raster[idyrow+1][idxcol+1] = get_random(tz); //c - 7
				raster[idyrow-1][idxcol+1] = get_random(tz); //i - 1
				raster[idyrow-1][idxcol-1] = get_random(tz); //g - 3
				raster[idyrow+1][idxcol-1] = get_random(tz); //a - 9
			}
		}
	}
	return raster;
}





#endif /* _MOD_RAS */
