/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza 
* Module:  m_fixsink.h
* 
* Autores: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2018-04-06
* Last update: 2022-09-30
* 
* Description:
* 
* Detect sinks in a DEM. Create an output xy.csv file with each sink
* (1) sink_mode: only check for sink cells
* (2) sink_mode: modify z values if sinks are isolated
* 
* 
* *********************************************************************/
#ifndef _SIN_RAS
#define _SIN_RAS

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
/*! 
* Standard libraries
***********************************************************************/
#include <stdio.h>

/*----------------------------------------------------------------------
 |  Function definition
 ---------------------------------------------------------------------*/


//**********************************************************************
/*! Auxiliar Functions */
//**********************************************************************

int calc_higlow(double h0, double *cellsval);

int calc_higlow16(double **raster, struct HeadR s_arrayh[], int i, int j)
{
double *cellsval, h0;
int k, ci, rj, l, result;

	l=0;
	for(k=0;k<8;k++)
	{
		if (k == 0)                         /*!< if z-diff in zero is higher than 0 */
		{
			ci = i--;                       /*!< Change column index */
			rj = j;                         /*!< Change column index */
		}
		if (k == 1)
		{
			ci = i++;                                
			rj = j;                                  
		}
		if (k == 2)
		{
			ci = i;                                  
			rj = j++;                                
		}
		if (k == 3)
		{
			ci = i;                                  
			rj = j--;                                
		}
		if (k == 4)
		{
			ci = i++;                                
			rj = j++;                                
		}
		if (k == 5)
		{
			ci = i--;                                
			rj = j++;                                
		}
		if (k == 6)
		{
			ci = i--;                                
			rj = j--;                                
		}
		if (k == 7)
		{
			ci = i++;                                
			rj = j--;                                
		}
		//-----------
		if((ci>3 && ci<s_arrayh[0].hn_fy-3 ) && ( rj>3 && rj<s_arrayh[0].hn_cx-3 ))
		{
			h0 = raster[ci][rj];
			cellsval = search_celproxF(raster, ci, rj);
			result   = calc_higlow(h0, cellsval);
		}
		if (result == 8) l++;
	}
	return l;
}

int calc_higlow(double h0, double *cellsval)
{
int i, l, n;
	l = 0;
	n = 0;
	for(i=0;i<8;i++)
	{
		if(cellsval[i] >  h0)l++; //count cells with higher z value
		if(cellsval[i] == h0)n++; //count cells with equal z value 
	}
	return l;
}

//**********************************************************************
/*! Calc Functions */
//**********************************************************************
/*! Search and detect sinks in DEM 
 * (1) sink_mode: only check for sink cells
 * (2) sink_mode: modify z values if sinks are isolated
 * */
double **fix_sinks(double **raster, struct HeadR s_arrayh[], char *outfile, float dem_nulval, int sink_mode)
{
FILE *file;	
int i, j, k, m, n, o, q;
int tisloate;
double h0,  diff;
double txcoor, tycoor;
double dzdx, dzdy, aspect, cell, celda;
float nrun;
int result, result_16, highdif;
double *cellsval, celldif[8];
double c1,c2,c3,c4,c6,c7,c8,c9; //c5 is h0
	
	//srand ( time(NULL) );
	//if(metsink == 1)  printf("%s\n", wrst(22));
	//if(metsink == 2)  printf("%s\n", wrst(23));
	m = n = o = 0;                                                              /*!< count surface depressions */
	//sprintf(outfile, "%sSink_pt.xyz", dir_out);                        /*!< Assign output filename */
    if((file = fopen(outfile,"wt"))== NULL)
    {
        printf("-------ERROR open file--------\n");
        printf("-----------ERROR--------------\n");
        printf("-----------ERROR--------------\n");
        printf("-----------ERROR--------------\n");
        printf("-----------ERROR--------------\n");
        printf("-----------ERROR--------------\n");
        exit(0);
    }
    else
    {
		/**
		* Print the head file 
		*/
		fprintf(file,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
				 "id", "xcoor","ycoor", "zcoorold", "zcoornew", "diffz", "mode", "iso", "tipdif", "ptim1j","pti1j","ptij1","ptijm1","pti1j1","ptim1j1","ptim1jm1","pti1jm1"); 
				 
		for(i=0;i<s_arrayh[0].hn_fy;i++) //column
		{
			for(j=0;j<s_arrayh[0].hn_cx;j++) //row
			{
				highdif=0;
				if((i>3 && i<s_arrayh[0].hn_fy-3 ) && ( j>3 && j<s_arrayh[0].hn_cx-3 ))          /*!< if array index are in working area */
				{
					h0  = raster[i][j];                                   /*!< get z coordinate value */
					if(h0 != dem_nulval)                                   /*!< if z coordinate value is not null */
					{
						cellsval = search_celproxF(raster, i, j);
						for (k=0;k<8;k++)
						{
							celldif[k] = fabsf(h0 - cellsval[k]);
							if (celldif[k] > (2 * s_arrayh[0].hresx)) highdif++;
						} 
						result = calc_higlow(h0, cellsval);
						//-----
						q = 0;                                          /*!< modification mode used */
						if(result == 8)                                 /*!< if center cell is a sink */
						{	
							//--
							if(sink_mode == 1)                            /*!< if algorithm detection is activated */
							{
								diff     = -999;
								q        = -999;
								tisloate = -999;
							}	
							//--
							if(sink_mode == 2)                            /*!< if algorithm correction is activated */
							{
								/**
								* Check if there is not another sink in the 16 neighboring cells 
								*/
								result_16 = calc_higlow16(raster, s_arrayh, i, j);
								if (result_16  >  0) tisloate = 0;              /*!< sink cell not isolate */
								if (result_16  == 0) tisloate = 1;              /*!< isolate sink cell */
								if (tisloate == 1)
								{
									c6 = cellsval[0];   //d - 6
									c4 = cellsval[1];
									c8 = cellsval[2];
									c2 = cellsval[3];   //h - 2
									//oblicuas
									c7 = cellsval[4];
									c9 = cellsval[5];
									c3 = cellsval[6];
									c1 = cellsval[1];
									/**
									* Using the SSM to get the spill cell 
									* http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q900000023000000
									* Burrough and McDonell (1998)
									* dzdx = ((c + 2f + i) - (a + 2d + g)) / 8
									* dzdy = ((g + 2h + i) - (a + 2b + c)) / 8	
									*/	
									dzdx = ((c7 + (2*c4) + c1) - (c9 + (2*c6) + c3)) / 8;	
									dzdy = ((c3 + (2*c2) + c1) - (c9 + (2*c8) + c7)) / 8;  	
									aspect = 180/PI * atan2 (dzdy, - dzdx);   /*!< calc s-aspect in degrees */
									/**
									* Change degree values starting from North as 0 
									*/
									if (aspect >= 0 && aspect <= 90)  cell = 90 - aspect; 
									if (aspect > 90 && aspect <= 180) cell = 360 - (aspect - 90);
									if (aspect >= -180 && aspect < 0) cell =  (aspect*-1) + 90;
									/**
									* Getting the spill cell z value 
									*/
									if (cell >=     0 && cell <   22.5) celda = cellsval[2];
									if (cell >=  22.5 && cell <   67.5) celda = cellsval[4];
									if (cell >=  67.5 && cell <  112.5) celda = cellsval[1];
									if (cell >= 112.5 && cell <  157.5) celda = cellsval[7];
									if (cell >= 157.5 && cell <  202.5) celda = cellsval[3];
									if (cell >= 202.5 && cell <  247.5) celda = cellsval[6];
									if (cell >= 247.5 && cell <  292.5) celda = cellsval[0];
									if (cell >= 292.5 && cell <  337.5) celda = cellsval[5];
									if (cell >= 337.5 && cell <= 360.0) celda = cellsval[2];
									/**
									* if neighboring cell is lower than the center cell 
									*/
									if(h0 > celda)
									{
										raster[i][j] = celda + 0.1;   /*!< exchange z value + 0.1 */
										diff = raster[i][j] - h0;
										q = 1;
									}
									/**
									* if neighboring cell is equal than the center cell 
									*/
									if(celda == h0)
									{
										nrun = (float) (rand()/ (float)RAND_MAX); 
										raster[i][j] = h0 + nrun;     /*!< use a random value */
										diff = nrun;
										q = 2;
									}	
									/**
									* Check if the new center cell produce new surface depressions around (16 neighboring cells) 
									*/
									result_16 = calc_higlow16(raster, s_arrayh, i, j);	
									if(result_16 > 0)                           /*!< if so */
									{
										raster[i][j] = h0;            /*!< do not change the center z value */
										diff = 0;
										q = 0;
									}
									/**
									* For debugging
									* printf("%i %i  -N %i Q %i-  %lf  %lf  %lf\n", i,  j, nalt, q, h0, topoless[i][j], diff);
									*/	
								}	
								if (tisloate == 0)                      /*!< if sink cell is not isolate */
								{
									diff = 0;                           /*!< do not change the z value */
									q    = 0;
								}
							}
							/**
							* Write surface depressions in the output file 
							*/
							txcoor = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, j);                 /*!< get coordinate from array index */
							tycoor = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, i);
							fprintf(file,"%i %lf %lf %lf %lf %lf %i %i %i %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n",
								m,
								txcoor,
								tycoor,
								h0,                                     /*!< original z value */
								raster[i][j],                           /*!< new z value */
								diff,                                   /*!< difference between z values */
								q,                                      /*!< type method used to change x value */
								tisloate,
								highdif,                                /*!< if altitud differences are large or small */
								cellsval[0],                            /*!< i-- j   ptim1j */
								cellsval[1],                            /*!< i++ j   pti1j */
								cellsval[2],                            /*!< i j++   ptij1 */
								cellsval[3],                            /*!< i j--   ptijm1 */
								cellsval[4],                            /*!< i++ j++ pti1j1 */
								cellsval[5],                            /*!< i-- j++ ptim1j1 */
								cellsval[6],                            /*!< i-- j-- ptim1jm1 */
								cellsval[7]                             /*!< i++ j-- pti1jm1 */
								);
							m++;
							
							if (highdif == 0) n++;
							if (highdif > 0) o++;
						}
					}
				}
			}
		}
	}									
	printf("Total surface depressions found %i\n", m);
	printf("Total surface depressions with large differences %i\n", o);
	printf("Total surface depressions with small differences %i\n", n);
	return raster;	
}


#endif /* _SIN_RAS */
