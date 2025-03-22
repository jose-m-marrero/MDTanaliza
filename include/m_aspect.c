/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza  
* Module:  m_aspect.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2018-04-06
* Last update: 2022-09-30
* 
* Description:
* 
* Calculation of aspect from a DEM
* (1) aspec_mode: output only in classes using LHM method
* (2) aspec_mode: output classes and degrees using SSM method
* 
* *********************************************************************/
#ifndef _ASP_RAS
#define _ASP_RAS

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
/*! Handle file names and paths, and other string functions */
#include "u_strings.h"
/*! 
* Standard libraries
***********************************************************************/
#include <string.h>

/*! Almacena los datos de la cabecera de arrays-raster grd */
//struct CabeceraR s_arrayh[NUMARR];


//**********************************************************************
/*! Calc Functions */
//**********************************************************************

/*! SLOPE-ASPECT CALCULATION (write on direct and rast2) 
 * In mode 1 aspect is only calculated in classes using LHM method
 * IN mode 2 aspect is also calculated in degrees using SSM method*/
void calc_saspect(const char *dir_out, const char *namfile, double **raster, struct HeadR s_arrayh[], float dem_nulval, int aspec_mode)
{
FILE *file;
char subfix[50];	
char *file_outrast, *file_histo;
int i, j, k, n;
int celda, cuenta[10];
double h0, alt, superf, aspect, cell, aspor[10];
double dzdx, dzdy;
double *cellsval, celldif[8], **rast_aspg, **rast_aspc; 
double c1,c2,c3,c4,c6,c7,c8,c9;
		
	if(aspec_mode == 2) rast_aspg = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	rast_aspc = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	
	printf("Calculing slope-Aspect method %i\n", aspec_mode);
	for(k=0;k<10;k++)
	{
		cuenta[k] = 0;                                    /*!< Reset values, count cells by s-aspect class */
		aspor[k]  = 0; 
	}
	
	for(i=0;i<s_arrayh[0].hn_fy;i++) //column
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) //row
		{
			/*!< if array index are inside working area  */
			if((i>3 && i<s_arrayh[0].hn_fy-3 ) && ( j>3 && j<s_arrayh[0].hn_cx-3 ))              
			{
				alt = 0;
				n   = 0;                                                /*!< Count if cell has slope-aspect value  */
				h0  = raster[i][j];                                     /*!< Get center cell z value  */
				if(h0 != dem_nulval)                                    /*!< if z coordinate value is not null */     
				{
					cellsval = busca_celproxF(raster, i, j);
					//getmovingcell(i, j, 0);                    			/*!< get 3x3 moving cell z values */
					for (k=0;k<8;k++)
					{
						celldif[k] = h0 - cellsval[k];
					} 
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
					* Using the LHM 
					*/
					if(aspec_mode == 1)                                     
					{
						for(k=0;k<8;k++)
						{
							if(celldif[k] > 0)                          /*!< if z difference is higher than 0  */ 
							{
								if(celldif[k] > alt)                    /*!< if z difference is the higher value */ 
								{
									alt   = celldif[k];                 /*!< max heigh */
									celda = k;                          /*!< Get index */ 
									n++;
								}
							}
						}
					}
					/**
					* Using the SSM 
					* http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q900000023000000
					* Burrough and McDonell (1998)
					* dzdx = ((c + 2f + i) - (a + 2d + g)) / 8
					* dzdy = ((g + 2h + i) - (a + 2b + c)) / 8	
					*/	
					if(aspec_mode == 2)                                    
					{	
						dzdx = ((c7 + (2*c4) + c1) - (c9 + (2*c6) + c3)) / 8;	
						dzdy = ((c3 + (2*c2) + c1) - (c9 + (2*c8) + c7)) / 8;  	
						aspect = 180/PI * atan2 (dzdy, - dzdx); 		/*!< calc s-aspect in degrees  0 located in East*/
						/**
						* Change degree values starting from North as 0 
						*/
						if (aspect >= 0 && aspect <= 90)  cell = 90 - aspect; 
						if (aspect > 90 && aspect <= 180) cell = 360 - (aspect - 90);
						if (aspect >= -180 && aspect < 0) cell =  (aspect*-1) + 90;
						//----
						rast_aspg[i][j] = cell;                                      /*!< Save s-aspect in degrees */
						//----
						if (cell >=     0 && cell <   22.5) celda = 2;
						if (cell >=  22.5 && cell <   67.5) celda = 4;
						if (cell >=  67.5 && cell <  112.5) celda = 1;
						if (cell >= 112.5 && cell <  157.5) celda = 7;
						if (cell >= 157.5 && cell <  202.5) celda = 3;
						if (cell >= 202.5 && cell <  247.5) celda = 6;
						if (cell >= 247.5 && cell <  292.5) celda = 0;
						if (cell >= 292.5 && cell <  337.5) celda = 5;
						if (cell >= 337.5 && cell <= 360.0) celda = 2;
						n++;	
					}					
					if(n == 0)                                          /*!< if cell is flat  */ 
					{
						rast_aspc[i][j] = 255;                          /*!< Save s-aspect in class */
						cuenta[8]++; 
					}	
					if(n > 0) 											/*!< if cell is not flat or null  */
					{
						switch(celda) 									/*!< get the spill neighboring value cell  */
						{
						case 0: 
							rast_aspc[i][j] = 32; 						/*!< Save s-aspect in class */
							cuenta[5]++;                                /*!< Count total cells with class 128 */
							break;
						case 1:
							rast_aspc[i][j] = 2; 
							cuenta[1]++;
							break;
						case 2:
							rast_aspc[i][j] = 128; 
							cuenta[7]++;
							break;
						case 3:
							rast_aspc[i][j] = 8;
							cuenta[3]++;
							break;
						case 4:
							rast_aspc[i][j] = 1;
							cuenta[0]++;
							break;
						case 5:
							rast_aspc[i][j] = 64;
							cuenta[6]++;
							break;
						case 6:
							rast_aspc[i][j] = 16; 
							cuenta[4]++;
							break;
						case 7:
							rast_aspc[i][j] = 4; 
							cuenta[2]++;
							break;
						}	
					}			
				}
				if(h0 == dem_nulval)                                    /*!< if cell is null  */
				{
					rast_aspg[i][j] = dem_nulval;                       /*!< The null valued is transferred  */
					rast_aspc[i][j] = dem_nulval;	 
					cuenta[9]++;
				}	
			}
		}			
	}
	/**
	* Write rasters
	*/
	s_arrayh[1].hn_cx      = s_arrayh[0].hn_cx;
	s_arrayh[1].hn_fy      = s_arrayh[0].hn_fy;
	s_arrayh[1].hxlo       = s_arrayh[0].hxlo;
	s_arrayh[1].hxhi       = s_arrayh[0].hxhi;
	s_arrayh[1].hylo       = s_arrayh[0].hylo;
	s_arrayh[1].hyhi       = s_arrayh[0].hyhi;
	s_arrayh[1].hzlo       = 0;
	s_arrayh[1].hzhi       = 255;                                       /*!< Assign z max value for raster head data */
	s_arrayh[1].hresx      = s_arrayh[0].hresx ;
	s_arrayh[1].hresy      = s_arrayh[0].hresy;
	
	strcpy(subfix,"_aspect-class.grd");
	file_outrast = get_pathnam(dir_out, namfile, subfix, 1);
	write_grdrasterF(file_outrast, rast_aspc, s_arrayh, 1, 1);          /*!< Write S-aspect class */
	
	if(aspec_mode == 2)
	{
		strcpy(subfix,"_aspect-degree.grd");
		file_outrast = get_pathnam(dir_out, namfile, subfix, 1);
		s_arrayh[1].hzhi       = 360;                                   /*!< Assign z max value for raster head data */
		write_grdrasterF(file_outrast, rast_aspg, s_arrayh, 1, 1);      /*!< S-aspect degree raster is only written in mode 2  */
	
	}
	int val[10] = {1,2,4,8,16,32,64,128,255,-9999};                     /*!< array with class values  */	
	/**
	* Write histograms
	*/
	printf("\nWrite histograms - slope-aspect %i\n", aspec_mode);
	i = 0;
	strcpy(subfix,"_aspect-histogram.csv");
	file_histo = get_pathnam(dir_out, namfile, subfix, 1);
	
    if((file = fopen(file_histo,"w"))== NULL)
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
	    fprintf(file,"%s\n",
	         "id freq perc surf"); //primera linea
        for(i=0;i<10;i++)
        {
			/**
			* Calc normalized values
			*/
			if(cuenta[i] == 0)aspor[i]=0;
			if(cuenta[i]  > 0)aspor[i]= (double)cuenta[i] / (double)s_arrayh[0].htotdemcel;
			superf = (cuenta[i] * s_arrayh[0].hresx * s_arrayh[0].hresy)/1000000;					/*!< Calc total surface by s-aspect class in km2 */
			printf("Total cells %i, surface %.2lf km2 with s-aspect %i\n", cuenta[i], superf, val[i]);
			/**
			* Write final values
			*/
			fprintf(file,"%i %i %lf %lf\n",
				i,
				cuenta[i],
				aspor[i],
				superf
				);
		}	
	} 
	fclose(file);
	free_RasterF("Rast aspect class", rast_aspc, s_arrayh[1].hn_fy);
	if(aspec_mode == 2) free_RasterF("Rast aspect degrees", rast_aspg, s_arrayh[1].hn_fy);

}		


#endif /* _ASP_RAS */
