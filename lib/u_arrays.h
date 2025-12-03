/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Module: u_arrays.h
* 
* Author: Jose M. Marrero, Ramon Ortiz.
* 	Code improved following Internet examples
* 
* Version: 0.1.2
* Creation Date: 2019-11-01
* Last update: 2025-12-03
* 
* Description: 
* 
* Get index and coordinates in matrix
* Find neighboring cell values:
* 	By getting proximal values: 8 or 24 neighboring cells 
* 	By cardinal directions
* Smooth raster
* Several methods to create an array (simple - 1 column or double - n column / n rows:
* 	By its dimension (columns and rows)
* 	By centroid and radius
* Dynamical reading of *.grd raster format (Binary or ASCII)
* Dynamical writing of *.grd raster format (Binary or ASCII)
* 
* *********************************************************************/

#ifndef _LEEARR_H
#define _LEEARR_H



/*!
* General libraries
***********************************************************************/  
/*! Some geometric calculation and data comparison functions */
#include "u_calculus.h"
/*! structure to save array's head data from *.grd input raster format*/
#include "struct_head-array.h"
/*! 
* Standard libraries
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*!
* DEFINES
***********************************************************************/ 
#define PI      3.14159265

/**
* DECLARA 
***********************************************************************/
/**
* Declara estructuras
***********************************************************************/


/*----------------------------------------------------------------------
 |  Function definition
 ---------------------------------------------------------------------*/
 
/*! Calcula distancia entre dos puntos, modulo u_calculos.h */
double calc_dist(double x1, double x2, double y1, double y2);
//Local
void Lee_Cabecera(struct HeadR *Cab, short int dimension[2], double coord[6], const char *raster_typetxt);
int *search_celproxI(int **raster, int idfy, int jdcx);
double *search_celproxF(double **raster, int idfy, int jdcx);
double* Crea_1DFarray(int arraySize);
double** Crea_2DFarray(int arraySize_fY, int arraySize_cX);

//**********************************************************************
/*! Get index or coordinates in raster-array */
//**********************************************************************

/*! Get cell coordinate from index
 * Parameters: cell coordinate, corner coordinate of raster (in 0) and resolution*/
int calc_Rindex(double coor, double coor0, double resol)
{
int index;
double dif, dem;

	dif = coor - coor0; 
	index = (int)(dif / resol);
	dem = (index * resol) + coor0; 
	if ((coor - dem) > (resol/2))index++;
	//printf("TEST: dif %lf res %lf, index %i, dem %lf\n", dif, resol, index, dem);
	return index;
}		

/*! Get cell coordinate from index
 * Parameters: cell coordinate, corner coordinate of raster (in 0) and resolution*/
double get_coor(double coor0, int resol, int index)
{
double coord;
	coord = coor0 + (index * resol);
	return coord;
}


/*! Get cell value from coordinate in integer raster
 * Parameters: cell coordinates, raster and header, raster number in struct and null value*/
double get_rastvalueI(double xcoor, double ycoor, int **raster, struct HeadR s_arrayh[], int num_rast, float nulval)
{
int idxfy, idxcx;
int rast_val;
	rast_val = 0;
	idxfy = calc_Rindex(ycoor, s_arrayh[num_rast].hylo, s_arrayh[num_rast].hresy);
	idxcx = calc_Rindex(xcoor, s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hresx);
	if(idxfy > 0 && idxfy < s_arrayh[num_rast].hn_fy && idxcx > 0 && idxcx < s_arrayh[num_rast].hn_cx) rast_val = raster[idxfy][idxcx];
	else rast_val = nulval;
	return rast_val;
}

/*! Get cell value from coordinate in double-float raster
 * Parameters: cell coordinates, raster and header, raster number in struct and null value*/
double get_rastvalueF(double xcoor, double ycoor, double **raster, struct HeadR s_arrayh[], int num_rast, float nulval)
{
int idxfy, idxcx;
double rast_val;
	rast_val = 0;
	idxfy = calc_Rindex(ycoor, s_arrayh[num_rast].hylo, s_arrayh[num_rast].hresy);
	idxcx = calc_Rindex(xcoor, s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hresx);
	if(idxfy > 0 && idxfy < s_arrayh[num_rast].hn_fy && idxcx > 0 && idxcx < s_arrayh[num_rast].hn_cx) rast_val = raster[idxfy][idxcx];
	else rast_val = nulval;
	return rast_val;
}


//**********************************************************************
/*! Search proximal cells in raster-array */
//**********************************************************************

/*! Get value from 8 neighboring cells when raster is integer
 * Important: avoid memory overflow if cells are close to the boundaries  */
int *search_celproxI(int **raster, int idfy, int jdcx)
{
static int celdaI[8];	
	/**
	* Array index named used and 3x3 moving window and equivalent orientation
	* || a b c  || 9  8  7  ||  i+1 j-1 : i+1 j :  i+1 j+1 ||  t5  t2 t4  ||  64  128  1
	* || d e f  || 6  5  4  ||  i  j-1  :  i j  :  i j+1   ||  t0  t  t1  ||  32  255  2
	* || g h i  || 3  2  1  ||  i-1 j-1 : i-1 j :  i-1 j+1 ||  t6  t3 t7  ||  16    8  4
	*/
	/**
	* Get z values from center and 8 neighbouring cells
	*/	
	//buscamos en las 4 celdas proximas perpendiculares
	celdaI[0] = raster[idfy-1][jdcx];   //h - 2
	celdaI[1] = raster[idfy+1][jdcx];   //b - 8
	celdaI[2] = raster[idfy][jdcx+1];   //f - 4
	celdaI[3] = raster[idfy][jdcx-1];   //d - 6
	//buscamos en las 4 celdas proximas oblicuas
	celdaI[4] = raster[idfy+1][jdcx+1]; //c - 7
	celdaI[5] = raster[idfy-1][jdcx+1]; //i - 1
	celdaI[6] = raster[idfy-1][jdcx-1]; //g - 3
	celdaI[7] = raster[idfy+1][jdcx-1]; //a - 9
	
	return (celdaI);
}

/*! Get value from 8 neighboring cells when raster is float
 * Important: avoid memory overflow if cells are close to the boundaries  */
double *search_celproxF(double **raster, int idfy, int jdcx)
{
static double celdaF[8];		
	/**
	* Array index named used and 3x3 moving window and equivalent orientation
	* || a b c  || 9  8  7  ||  i+1 j-1 : i+1 j :  i+1 j+1 ||  t5  t2 t4  ||  64  128  1
	* || d e f  || 6  5  4  ||  i  j-1  :  i j  :  i j+1   ||  t0  t  t1  ||  32  255  2
	* || g h i  || 3  2  1  ||  i-1 j-1 : i-1 j :  i-1 j+1 ||  t6  t3 t7  ||  16    8  4
	*/
	//buscamos en las 4 celdas proximas perpendiculares
	celdaF[0] = raster[idfy-1][jdcx];   //h - 2
	celdaF[1] = raster[idfy+1][jdcx];   //b - 8
	celdaF[2] = raster[idfy][jdcx+1];   //f - 4
	celdaF[3] = raster[idfy][jdcx-1];   //d - 6
	//buscamos en las 4 celdas proximas oblicuas
	celdaF[4] = raster[idfy+1][jdcx+1]; //c - 7
	celdaF[5] = raster[idfy-1][jdcx+1]; //i - 1
	celdaF[6] = raster[idfy-1][jdcx-1]; //g - 3
	celdaF[7] = raster[idfy+1][jdcx-1]; //a - 9
	
	
	return (celdaF);
}	

/*! Get value from 24 neighboring cells when raster is float
* Important: avoid memory overflow if cells are close to the boundaries  */
double *search_celproxF24(double **raster, int idfy, int jdcx)
{
static double celdaF[24];		
	/**
	* Array index named used and 3x3 moving window and equivalent orientation
	* || j k l m n || 10 11 12 13 14 ||i+2 j-2 : i+2 j-1 : i+2 j :  i+2 j+1 : i+2 j+2
	* || o a b c p || 15  9  8  7 16 ||i+1 j-2 : i+1 j-1 : i+1 j :  i+1 j+1 : i+1 j+2    ||  t5  t2 t4  ||  64  128  1
	* || q d e f r || 17  6  5  4 18 ||i   j-2 : i   j-1 : i j   :  i j+1   : i   j+2    ||  t0  t  t1  ||  32  255  2
	* || s g h i t || 19  3  2  1 20 ||i-1 j-2 : i-1 j-1 : i-1 j :  i-1 j+1 : i-1 j+2    ||  t6  t3 t7  ||  16    8  4
	* || u v w x y || 21 22 23 24 25 ||i-2 j-2 : i-2 j-1 : i-2 j :  i-2 j+1 : i-2 j+2
	*/
	//Proximal perpendicular cells
	celdaF[0] = raster[idfy-1][jdcx];   //h - 2
	celdaF[1] = raster[idfy+1][jdcx];   //b - 8
	celdaF[2] = raster[idfy][jdcx+1];   //f - 4
	celdaF[3] = raster[idfy][jdcx-1];   //d - 6
	//Proximal oblique cells 
	celdaF[4] = raster[idfy+1][jdcx+1]; //c - 7
	celdaF[5] = raster[idfy-1][jdcx+1]; //i - 1
	celdaF[6] = raster[idfy-1][jdcx-1]; //g - 3
	celdaF[7] = raster[idfy+1][jdcx-1]; //a - 9
	//Upper row 10-11-12-13-14
	celdaF[8]  = raster[idfy+2][jdcx-2]; 
	celdaF[9]  = raster[idfy+2][jdcx-1]; 
	celdaF[10] = raster[idfy+2][jdcx];
	celdaF[11] = raster[idfy+2][jdcx+1]; 
	celdaF[12] = raster[idfy+2][jdcx+2];
	//Bottom row 21-22-23-24-25
	celdaF[13]  = raster[idfy-2][jdcx-2]; 
	celdaF[14]  = raster[idfy-2][jdcx-1]; 
	celdaF[15]  = raster[idfy-2][jdcx]; 
	celdaF[16]  = raster[idfy-2][jdcx+1]; 
	celdaF[17]  = raster[idfy-2][jdcx+2];
	//Column left side 15-17-19
	celdaF[18]  = raster[idfy+1][jdcx-2];
	celdaF[19]  = raster[idfy][jdcx-2];
	celdaF[20]  = raster[idfy-1][jdcx-2];
	//Column right side 16-18-20
	celdaF[21]  = raster[idfy+1][jdcx+2];
	celdaF[22]  = raster[idfy][jdcx+2];
	celdaF[23]  = raster[idfy-1][jdcx+2];

	return (celdaF);
}	


//**********************************************************************
/*! Get or change values from proximal or distal cells in raster-array */
//**********************************************************************

/*! Get indexes of cell according to defined search
* if cell distance is equal to 1, only will get the 8 neighboring cells
* if cell distance shows variations, a cardinal search will be used */
int* celval_increase(int direction_searching, int row, int column, const int valfy, const int valcx, int celdist)
{
static int indexes[3];

	indexes[0] = row;
	indexes[1] = column;
	indexes[2] = 0;       //check raster boundaries
	
	if (direction_searching == 1) indexes[0] = indexes[0] + celdist;                /*!< up */
	if (direction_searching == 2) indexes[0] = indexes[0] - celdist;                /*!< down */
	if (direction_searching == 3) indexes[1] = indexes[1] + celdist;                /*!< right */
	if (direction_searching == 4) indexes[1] = indexes[1] - celdist;                /*!< left */
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if (direction_searching == 5)                                                   /*!< right up oblique */                                   
	{
		indexes[0] = indexes[0] + celdist;
		indexes[1] = indexes[1] + celdist;
	}
	if (direction_searching == 6)                                                   /*!< left down oblique */
	{
		indexes[0] = indexes[0] - celdist; 
		indexes[1] = indexes[1] - celdist;  
	}
	if (direction_searching == 7)                                                   /*!< left up oblique */
	{
		indexes[0] = indexes[0] + celdist;
		indexes[1] = indexes[1] - celdist;  
	}
	if (direction_searching == 8)                                                   /*!< right down oblique */
	{
		indexes[0] = indexes[0] - celdist; 
		indexes[1] = indexes[1] + celdist;   
	}
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	/*!< check if cell is in raster */
	if(indexes[0] <= 0 || indexes[0] >= valfy) indexes[2] = 1;
	if (indexes[2] == 0)
	{
		if(indexes[1] <= 0 || indexes[1] >= valcx) indexes[2] = 1;
	}
	//printf("TEST: indexes %i : %i %i %i :: %i %i\n", direction_searching, indexes[0], indexes[1], indexes[2], valfy, valcx);
	return (indexes);
}

/*! Get integer values from 8 neighboring cells, excluding nulls and repeated values
 * Parameters: array dimensions, total number of vertex, index of central value
 * null value and raster number */
int* get_valproxI(int **raster, int idfy, int jdcx, int maxver, int id, int grdnul, int num_rast, int filagen, int columgen)
{
int k, i, exist, count;
int *validI;
static int *stored;

	count=0;	
	stored = Crea_1DIarray(maxver);
	
	if ((jdcx > 2 && jdcx < columgen-2) && (idfy > 2 && idfy < filagen-2))  /*! Check if is raster */
	{
		validI = search_celproxI(raster, idfy, jdcx);                   /*! get values from 8 neighboring cells */
		/*! Get results and store them */
		for (k=0;k<8;k++)                                               
		{
			if (validI[k] != id && validI[k] != grdnul)                 /*! if value is not null and different than cell center value */
			{
				exist=0;
				if (count > 0)                                         /*! exclude first element */
				{
					for (i=0;i<count;i++)
					{
						if (stored[i] == validI[k])                     /*! if value already exist - break */
						{
							exist=1;
							break;
						}
					}
				}
				if (exist==0)                                           /*! if value does not exist, add to the array */
				{
					stored[count] = validI[k]; 
					count++;
				}
			}
		}
	}
	return (stored);
}

/*! Get double-float values from 8 neighboring cells, excluding nulls and repeated values
 * Parameters: array dimensions, total number of vertex, index of central value
 * null value and raster number */
double* get_valproxF(double **raster, int idfy, int jdcx, int maxver, int id, int grdnul, int num_rast, int filagen, int columgen)
{
int k, i, exist, count;
double *validI;
static double *stored;

	count=0;	
	stored = Crea_1DFarray(maxver);
	
	if ((jdcx > 2 && jdcx < columgen-2) && (idfy > 2 && idfy < filagen-2))  /*! Check if is raster */
	{
		validI = search_celproxF(raster, idfy, jdcx);                   /*! get values from 8 neighboring cells */
		/*! Get results and store them */
		for (k=0;k<8;k++)                                               
		{
			if (validI[k] != id && validI[k] != grdnul)                 /*! if value is not null and different than cell center value */
			{
				exist=0;
				if (count > 0)                                         /*! exclude first element */
				{
					for (i=0;i<count;i++)
					{
						if (stored[i] == validI[k])                    /*! if value already exist - break */
						{
							exist=1;
							break;
						}
					}
				}
				if (exist==0)                                           /*! if value does not exist, add to the array */
				{
					stored[count] = validI[k]; 
					count++;
				}
			}
		}
	}
	return (stored);
}

/*! Calculate closed z out of a feature, on the floor/street (for buildings heights in exposure models)
 * Two raster: polygon features digitized using ID as z value IDraster
 *   and Digital Superficial Model DSM (including building, trees, etc.)
 * Search in IDraster following one direction until free space is reached
 * If not open space (other features are too close), try in a different direction
 * When open space is reached (a street), get z value form DSM (on the floor)
 * */
double *calc_dirdem(int **Iraster, int columcell, int rowcell, double **Fraster, struct HeadR s_arrayh[], \
	int num_rast, int num_rast2, int direction_searching, double height_feature, int min_height, int max_hieght, int grdnull, int maxsteps)
{
int id_feature, vcol, vrow, steps;
double new_height, dsm_val, diff_height;
double globx, globy;
int *datosid;
static double datos[2];
	steps = 0;          /*!< count distance in cell numbers */
	new_height = 0;     /*!< nearest z value out of the feature */
	/*!< search until maximum searching distance is reached or new_height = 0 */
	do
	{
		if ((columcell > 0 && columcell < s_arrayh[num_rast].hn_cx) && (rowcell > 0 && rowcell < s_arrayh[num_rast].hn_fy))
		{
			id_feature  = Iraster[rowcell][columcell];                  /*!< get feature ID from IDraster */
			if (id_feature == grdnull )                                 /*!< If no ID is found -> open space */
			{
				globx = get_coor(s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hresx, columcell);      /*!< get feature's x-coordinate of in IDraster */
				globy = get_coor(s_arrayh[num_rast].hylo, s_arrayh[num_rast].hresy, rowcell); 
				
				vcol = calc_Rindex(globx, s_arrayh[num_rast2].hxlo, s_arrayh[num_rast2].hresx);       /*!< Get feature column-index in DSM */
				vrow = calc_Rindex(globy, s_arrayh[num_rast2].hylo, s_arrayh[num_rast2].hresy);
				if ((vcol > 0 && vcol < s_arrayh[num_rast2].hn_cx) && (vrow > 0 && vrow < s_arrayh[num_rast2].hn_fy)) /*!< Check if feature is in raster DSM */
				{
					
					dsm_val = Fraster[vrow][vcol];                         /*!< Get z value from DSM */
					if (dsm_val != grdnull || dsm_val > 0)                 /*!< If z value is ok */
					{
						diff_height = abs(dsm_val - height_feature);                      /*!< Calculate height difference (street - roof), its feature's height */
						if (diff_height >= min_height && diff_height <= max_hieght)       /*!< if difference is less than .. or higher than ... ok */
						{
							new_height = dsm_val;
						}
					}
				}
			}
		}
		/**
		* Array index named used and 3x3 moving window and equivalent orientation
		* || a b c  || 9  8  7  ||  i+1 j-1 : i+1 j :  i+1 j+1 ||  t5  t2 t4  ||  64  128  1
		* || d e f  || 6  5  4  ||  i  j-1  :  i j  :  i j+1   ||  t0  t  t1  ||  32  255  2
		* || g h i  || 3  2  1  ||  i-1 j-1 : i-1 j :  i-1 j+1 ||  t6  t3 t7  ||  16    8  4
		*/
		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		/*!< Get a new proximal cell */
		datosid = celval_increase(direction_searching, rowcell, columcell, s_arrayh[num_rast].hn_fy, s_arrayh[num_rast].hn_cx, 1);
		if (datosid[2] == 1) break;  /*!< Out of DSM - break */
		else
		{
			/*!< Reset index and count cell */
			rowcell   = datosid[0];
			columcell = datosid[1];
			steps++;
		}
	}while(steps < maxsteps && new_height == 0);				 
	
	datos[0] = new_height;
	datos[1] = steps;
	return (datos);
}



//**********************************************************************
/*! Change raster-array values */
//**********************************************************************

/*! Smooth output raster as alternative version - mostly for interpolation 
 * mode_soft 1: only 8 neighboring cells
 * mode_soft 2: 24 proximal cells
 * Detect first abrupt differences between data values in raster
 *     Percentages represent the difference between the values of the central cell and its neighboring cells
 *     Smooth is apply if max (r_maxpercent) and min (r_minpercent) percentages are exceeded a certain number of times (r_totoutcells)
 * New smooth value is calculated as the average of all non null proximal cell values, excluding the center cell value
 * Recommended parameters: small percentage differences produce higher changes
 * */
double **smooth_rasterF(double **raster, struct HeadR s_arrayh[], int r_indx, float raster_null, int s_modesmooth, float s_maxpercent, float s_minpercent, int s_totoutcels)
{
int i, j, k;
int tot_cells, cells_mod, count_valid, proxcell_outlim, borderlim;
double *prox_cellval, **rastersoft, max_proxcel, min_proxcel, sum_cellsval, percent;	
	printf("Smoothing raster data in mode %i\n", s_modesmooth);
	rastersoft = Crea_2DFarray(s_arrayh[r_indx].hn_fy, s_arrayh[r_indx].hn_cx);
	
	if (s_modesmooth == 1) 
	{
		borderlim = 2; /*! cells out of calculation to be close to the boundary */
		tot_cells = 8; /*! number of proximal cells for searching */
	}
	if (s_modesmooth == 2) 
	{
		borderlim = 4;
		tot_cells = 24;
	}
	cells_mod=0;
	for(i=0;i<s_arrayh[r_indx].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[r_indx].hn_cx;j++) 
		{
			/*! if cell out of raster boundaries */
			if ((i > borderlim && i < s_arrayh[r_indx].hn_fy -borderlim) && (j > borderlim && j < s_arrayh[r_indx].hn_fy-borderlim))
			{
				/*! if value not null */
				if (raster[i][j] != raster_null)
				{
					count_valid = proxcell_outlim = sum_cellsval = max_proxcel = 0;
					min_proxcel = 100000000000;

					if (s_modesmooth == 1) prox_cellval = search_celproxF(raster, i, j);   /*! Get 8 cell values */
					if (s_modesmooth == 2) prox_cellval = search_celproxF24(raster, i, j); /*! Get 24 cell values */

					/*! check values from proximal cells */
					for(k=0;k<tot_cells;k++) 
					{
						//printf("z next %lf\n", prox_cellval[k]);
						/*! if proximal cell value is not null */
						if (prox_cellval[k] != raster_null)
						{
							/*! calc percentage between cell and proximal cell values */
							percent = (prox_cellval[k] * 100) / raster[i][j];
							
							/*! if percentage out of min and max limits - get maximum or minimum value */
							if (percent > s_maxpercent || percent < s_minpercent)
							{
								if (prox_cellval[k] > max_proxcel) max_proxcel = prox_cellval[k];      /*! get proximal max z value */
								if (prox_cellval[k] < min_proxcel) min_proxcel = prox_cellval[k];      /*! get proximal min z value */
								proxcell_outlim++;                                                     /*! total number of proximal cells out of limits */
							}
							sum_cellsval+= prox_cellval[k];   /*! Sum all non null proximal cell values */
							count_valid++;                    /*! Count all non null proximal cell values */
						}
					}
					//printf("z %lf - sum_cellsval %lf - valid %i :: pe %lf - %i\n", raster[i][j], sum_cellsval, count_valid, percent, tot_cells);
					
					/*! if proximal cell out of limits is higher than a predefined value  */
					if (proxcell_outlim > s_totoutcels)
					{
						rastersoft[i][j] = sum_cellsval / count_valid;             /*! Average calculation between all non null proximal cells */
						//rastersoft[i][j] = max_proxcel + min_proxcel / 2;       /*! Average from maximum and minimum - it does not work well */
						cells_mod++;
					}
					else rastersoft[i][j] = raster[i][j];               /*! if condition is not match, value in raster remains */
				}
				else rastersoft[i][j] = raster[i][j];
			}
			else rastersoft[i][j] = raster[i][j];
		}
	}
	
	printf("Total cell modified: %i\n",cells_mod); 
	
	return rastersoft;
}


//**********************************************************************
/*! Create raster-array */
//**********************************************************************

/*! Create a double dynamic double-float array by its dimension (col and rows) */
double** Crea_2DFarray(int arraySize_fY, int arraySize_cX)
{
int i;	
double** theArray;
	theArray = (double**) malloc(arraySize_fY*sizeof(double*));
	for(i=0; i<arraySize_fY;i++)
	{
		theArray[i] = (double*) malloc(arraySize_cX*sizeof(double));
		if(theArray[i]== NULL)printf("error 2\n");
	}	
	return theArray;
} 

/*! Create a one column dynamic double-float array by its dimension (col) */
double* Crea_1DFarray(int arraySize) 
{
double* db_array;
    db_array = (double*) malloc(arraySize*sizeof(double));
    if (db_array == NULL)
    {       
        printf("Unable to allocate memory, exiting.\n");
        free(db_array);
        exit(0);
    }
    return db_array;
}

/*! Create a double dynamic integer array by its dimension (col and rows) */
int** Crea_2DIarray(int arraySize_fY, int arraySize_cX) 
{
int i;	
int** theArray;
	theArray = (int**) malloc(arraySize_fY*sizeof(int*));
	for(i=0; i<arraySize_fY;i++)
	{
		theArray[i] = (int*) malloc(arraySize_cX*sizeof(int));
		if(theArray[i]== NULL)printf("error 2\n");
	}	
	return theArray;
} 

/*! Create a one column dynamic integer array by its dimension (col) */
int* Crea_1DIarray(int arraySize) 
{
int* int_array;
    int_array = (int*) malloc(arraySize*sizeof(int));

    if (int_array == NULL)
    {       
        printf("Unable to allocate memory, exiting.\n");
        free(int_array);
        exit(0);
    }
    return int_array;
}

/*! Create a double dynamic integer array by centroid and radius 
 * Add equal value to all cells */
int **crea_rasterI(double xcoor, double ycoor, int cell_value, int resenx, int reseny, double radio, struct HeadR s_arrayh[], int num_rast)
{
int i, j, n_cx, n_fy;
double resx, resy;
int **raster;

	s_arrayh[num_rast].hxlo   = xcoor - radio;   /*! minimum x-coordinate */
	s_arrayh[num_rast].hxhi   = xcoor + radio;   /*! maximum x-coordinate */
	s_arrayh[num_rast].hylo   = ycoor - radio;   /*! minimum y-coordinate */
	s_arrayh[num_rast].hyhi   = ycoor + radio;   /*! maximum y-coordinate */
	s_arrayh[num_rast].hzlo   = 0;               /*! minimum z-coordinate */
	s_arrayh[num_rast].hzhi   = cell_value;      /*! maximum z-coordinate */
	resx  = resenx;                              /*! Raster x-resolution */
	resy  = reseny;                              /*! Raster y-resolution */
	n_cx  = (s_arrayh[num_rast].hxlo + (2* radio) - s_arrayh[num_rast].hxlo) / resx;   /*! get number of columns in x */
	n_fy  = (s_arrayh[num_rast].hylo + (2* radio) - s_arrayh[num_rast].hylo) / resy;   /*! get number of rows in y */
	s_arrayh[num_rast].hn_cx   = n_cx;           /*! add number of columns in x */
	s_arrayh[num_rast].hn_fy   = n_fy;           /*! add number of rows in y */
	s_arrayh[num_rast].harea  = (n_cx * resx) * (n_fy * resy);            /*! calc area */
	s_arrayh[num_rast].hperim = (n_cx * resx) * 2 + (n_fy * resy) * 2;    /*! calc perimeter */
	s_arrayh[num_rast].hinvresx   = 1/resx;                               /*! calc inverse of x-resolution */
	s_arrayh[num_rast].hinvresy   = 1/resy;                               /*! calc inverse of y-resolution */
	s_arrayh[num_rast].hresxx     = resx*resx;                            /*! calc double of x-resolution */
	s_arrayh[num_rast].hresyy     = resy*resy;                            /*! calc double of y-resolution */
	s_arrayh[num_rast].htotdemcel = s_arrayh[num_rast].hn_cx * s_arrayh[num_rast].hn_fy; /*! calc total number of cells */
			
    printf("Coordinates xmin = %f -- xmax = %f\n",  s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hxhi);
	printf("Coordinates ymin = %f -- ymax = %f\n",  s_arrayh[num_rast].hylo, s_arrayh[num_rast].hyhi);
    printf("Values zlo = %10.4f -- zhi = %10.4f\n", s_arrayh[num_rast].hzlo, s_arrayh[num_rast].hzhi);
    printf("Raster x-resolution resx = %f\n",       s_arrayh[num_rast].hresx);
    printf("Raster y-resolution resy = %f\n",       s_arrayh[num_rast].hresy);
    printf("Raster inverse resolution invresx = %f and invresy = %f\n", s_arrayh[num_rast].hinvresx, s_arrayh[num_rast].hinvresy);
    printf("Raster double resolution resx2 = %f and resy2 = %f\n",   s_arrayh[num_rast].hresxx, s_arrayh[num_rast].hresyy);
    printf("Perimeter %lf\n",                      s_arrayh[num_rast].hperim);
    printf("Area %lf\n",                           s_arrayh[num_rast].harea);
    
	if (s_arrayh[num_rast].hn_cx > 22000)
	{
		printf("Attention:, raster needs too much memory\n");
		printf("No more than %i columns should be reached\n", n_cx);
		printf("Resize raster with lower resolution or smaller radius\n");
		exit(0);
	}
	raster = Crea_2DIarray(n_fy, n_cx);
	for(i=0;i<n_fy;i++) 
	{
		for(j=0;j<n_cx;j++)raster[i][j] = cell_value; //add equal value to the cell
	}	
	return raster;
}

/*! Create a double dynamic double-float array by centroid and radius 
* Add equal value to all cells */
double **crea_rasterF(double xcoor, double ycoor, int cell_value, int resenx, int reseny, double radio, struct HeadR s_arrayh[], int num_rast)
{
int i, j, n_cx, n_fy;
double resx, resy;
double **raster;

	s_arrayh[num_rast].hxlo   = xcoor - radio;
	s_arrayh[num_rast].hxhi   = xcoor + radio; 
	s_arrayh[num_rast].hylo   = ycoor - radio;
	s_arrayh[num_rast].hyhi   = ycoor + radio;
	s_arrayh[num_rast].hzlo   = 0;
	s_arrayh[num_rast].hzhi   = cell_value;
	resx  = resenx;
	resy  = reseny;
	n_cx  = (s_arrayh[num_rast].hxlo + (2* radio) - s_arrayh[num_rast].hxlo) / resx;
	n_fy  = (s_arrayh[num_rast].hylo + (2* radio) - s_arrayh[num_rast].hylo) / resy;
	s_arrayh[num_rast].hn_cx   = n_cx;
	s_arrayh[num_rast].hn_fy   = n_fy;
	s_arrayh[num_rast].harea  = (n_cx * resx) * (n_fy * resy);
	s_arrayh[num_rast].hperim = (n_cx * resx) * 2 + (n_fy * resy) * 2;
	s_arrayh[num_rast].hinvresx   = 1/resx;
	s_arrayh[num_rast].hinvresy   = 1/resy;
	s_arrayh[num_rast].hresxx     = resx*resx;
	s_arrayh[num_rast].hresyy     = resy*resy;
	s_arrayh[num_rast].htotdemcel = s_arrayh[num_rast].hn_cx * s_arrayh[num_rast].hn_fy;
			
	printf("Coordinates xmin = %f -- xmax = %f\n",   s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hxhi);
	printf("Coordinates ymin = %f -- ymax = %f\n",   s_arrayh[num_rast].hylo, s_arrayh[num_rast].hyhi);
    printf("Values min zlo = %10.4f -- max zhi = %10.4f\n", s_arrayh[num_rast].hzlo, s_arrayh[num_rast].hzhi);
    printf("Raster x-resolution resx = %f\n",        s_arrayh[num_rast].hresx);
    printf("Raster y-resolution resy = %f\n",        s_arrayh[num_rast].hresy);
    printf("Raster 1/2 resolution invresx = %f and invresy = %f\n", s_arrayh[num_rast].hinvresx, s_arrayh[num_rast].hinvresy);
    printf("Raster double resolution resx2 = %f and resy2 = %f\n",  s_arrayh[num_rast].hresxx, s_arrayh[num_rast].hresyy);
    printf("Perimeter %lf\n",  s_arrayh[num_rast].hperim);
    printf("Area %lf\n",  s_arrayh[num_rast].harea);
    
	if (s_arrayh[num_rast].hn_cx > 22000)
	{
		printf("Attention:, raster needs too much memory\n");
		printf("No more than %i columns should be reached\n", n_cx);
		printf("Resize raster with lower resolution or smaller radius\n");
		exit(0);
	}
	raster = Crea_2DFarray(n_fy, n_cx);
	for(i=0;i<n_fy;i++) 
	{
		for(j=0;j<n_cx;j++)raster[i][j] = cell_value; //add equal value to the cell
	}	
	return raster;
}


//**********************************************************************
/*! Read raster-array */
//**********************************************************************

/*! Read and store header *grd ASCII raster
 *  num_rast define position of raster in struct 
 * that number must be retained to access data later*/
void read_headerA(struct HeadR s_arrayh[], int num_rast, FILE* infile)
{
char raster_typetxt[8];
int n_cx, n_fy;
double xlo, xhi, ylo, yhi, zlo, zhi;
double xres, yres;
	/*! read first 5 header lines of .grd ASCII raster */
	fscanf(infile,"%s",      raster_typetxt);    /*! Raster type in text */
    fscanf(infile,"%i %i",   &n_cx, &n_fy);
	fscanf(infile,"%lf %lf", &xlo, &xhi);
	fscanf(infile,"%lf %lf", &ylo, &yhi);
	fscanf(infile,"%lf %lf", &zlo, &zhi);
	/*! Add and calculate values */
	xres = (xhi-xlo)/(double)(n_cx-1);
	yres = (yhi-ylo)/(double)(n_fy-1);
	s_arrayh[num_rast].hn_cx      = n_cx;         /*! add number of columns in x */
	s_arrayh[num_rast].hn_fy      = n_fy;         /*! add number of rows in y */
	s_arrayh[num_rast].hxlo       = xlo;          /*! minimum x-coordinate */
	s_arrayh[num_rast].hxhi       = xhi;          /*! maximum x-coordinate */
	s_arrayh[num_rast].hylo       = ylo;          /*! minimum y-coordinate */
	s_arrayh[num_rast].hyhi       = yhi;          /*! maximum y-coordinate */
	s_arrayh[num_rast].hzlo       = zlo;          /*! minimum z-coordinate */
	s_arrayh[num_rast].hzhi       = zhi;          /*! maximum z-coordinate */
	s_arrayh[num_rast].hresx      = xres;         /*! Raster x-resolution */
	s_arrayh[num_rast].hresy      = yres;         /*! Raster y-resolution */
	s_arrayh[num_rast].hinvresx   = 1/xres;       /*! calc inverse of x-resolution */
	s_arrayh[num_rast].hinvresy   = 1/yres;       /*! calc inverse of y-resolution */
	s_arrayh[num_rast].hresxx     = xres*xres;    /*! calc double of x-resolution */
	s_arrayh[num_rast].hresyy     = yres*yres;    /*! calc double of y-resolution */
	s_arrayh[num_rast].htotdemcel = n_cx * n_fy;  /*! calc total number of cells */
	s_arrayh[num_rast].hperim = (n_cx * xres) * 2 + (n_fy * yres) * 2;  /*! calc perimeter */
	s_arrayh[num_rast].harea  = (n_cx * xres) * (n_fy * yres);          /*! calc area */
	
	printf(".grd Raster header:\n");
	printf("raster type: %s\n", raster_typetxt);
    printf("colx  = %5i -- rowy = %5i\n",         s_arrayh[num_rast].hn_cx, s_arrayh[num_rast].hn_fy);
	printf("xmin = %f -- xmax = %f\n",            s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hxhi);
	printf("ymin = %f -- ymax = %f\n",            s_arrayh[num_rast].hylo, s_arrayh[num_rast].hyhi);
	printf("zlo = %10.4f -- zhi = %10.4f\n",      s_arrayh[num_rast].hzlo, s_arrayh[num_rast].hzhi);
	printf("Raster x-resolution resx = %f\n",     s_arrayh[num_rast].hresx);
	printf("Raster y-resolution resy = %f\n",     s_arrayh[num_rast].hresy);
	printf("Raster inverse resolution invresx = %f and invresy = %f\n", s_arrayh[num_rast].hinvresx, s_arrayh[num_rast].hinvresy);
	printf("Raster double resolution resx2 = %f and resy2 = %f\n",     s_arrayh[num_rast].hresxx, s_arrayh[num_rast].hresyy);
	printf("Perimeter %lf\n",  s_arrayh[num_rast].hperim);
	printf("Area %lf\n",  s_arrayh[num_rast].harea);
	if (s_arrayh[num_rast].hxlo < 0 ||  s_arrayh[num_rast].hylo < 0)
	{
		printf("ATENTION: xlow or ylow corner coordinates are below 0, check raster dimension or geographic proyection\n");
		exit(0);
	}
}

/*! Read and store header Binary *grd data */ 
void read_headerB(struct HeadR s_arrayh[], int num_rast, FILE* infile)
{
char raster_typetxt[8];
short int dimension[4];
double datocoor[8];
double xres, yres;
	fread(&raster_typetxt,4,1,infile);
	fread(&dimension,2*sizeof(short int),1, infile);
	fread(&datocoor,6*sizeof(double),1, infile);
	raster_typetxt[4] = 0;
	xres = (datocoor[1]-datocoor[0])/(double)(dimension[0]-1);
	yres = (datocoor[3]-datocoor[2])/(double)(dimension[1]-1);
	s_arrayh[num_rast].hn_cx      = dimension[0];
	s_arrayh[num_rast].hn_fy      = dimension[1];
	s_arrayh[num_rast].hxlo       = datocoor[0];
	s_arrayh[num_rast].hxhi       = datocoor[1];
	s_arrayh[num_rast].hylo       = datocoor[2];
	s_arrayh[num_rast].hyhi       = datocoor[3];
	s_arrayh[num_rast].hzlo       = datocoor[4];
	s_arrayh[num_rast].hzhi       = datocoor[5];
	s_arrayh[num_rast].hresx      = xres;
	s_arrayh[num_rast].hresy      = yres;
	s_arrayh[num_rast].hinvresx   = 1/xres;
	s_arrayh[num_rast].hinvresy   = 1/yres;
	s_arrayh[num_rast].hresxx     = xres*xres;
	s_arrayh[num_rast].hresyy     = yres*yres;
	s_arrayh[num_rast].htotdemcel = dimension[0] * dimension[1];
	s_arrayh[num_rast].hperim = (dimension[0] * xres) * 2 + (dimension[1] * yres) * 2;
	s_arrayh[num_rast].harea  = (dimension[0] * xres) * (dimension[1] * yres);
	
	printf(".grd Raster header:\n");
	printf("raster type: %s\n", raster_typetxt);
	printf("colx  = %5i -- rowy = %5i\n",         s_arrayh[num_rast].hn_cx, s_arrayh[num_rast].hn_fy);
	printf("xmin = %f -- xmax = %f\n",            s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hxhi);
	printf("ymin = %f -- ymax = %f\n",            s_arrayh[num_rast].hylo, s_arrayh[num_rast].hyhi);
	printf("zlo = %10.4f -- zhi = %10.4f\n",      s_arrayh[num_rast].hzlo, s_arrayh[num_rast].hzhi);
	printf("Raster x-resolution resx = %f\n",       s_arrayh[num_rast].hresx);
	printf("Raster y-resolution = %f\n",       s_arrayh[num_rast].hresy);
	printf("Raster inverse resolution invresx = %f and invresy = %f\n", s_arrayh[num_rast].hinvresx, s_arrayh[num_rast].hinvresy);
	printf("Raster double resolution resx2 = %f and resy2 = %f\n",     s_arrayh[num_rast].hresxx, s_arrayh[num_rast].hresyy);
	printf("Perimeter %lf\n",  s_arrayh[num_rast].hperim);
	printf("Area %lf\n",  s_arrayh[num_rast].harea);
	
	if (s_arrayh[num_rast].hxlo < 0 ||  s_arrayh[num_rast].hylo < 0)
	{
		printf("ATENTION: xlow or ylow corner coordinates are below 0, check raster dimension or geographic proyection\n");
		exit(0);
	}
}

/*! Read ASCII or Binary integer *grd raster 
 * ntype_raster: Raster type must be defined (1-binary 2-ASCII)
 * Other inputs parameters: raster name, null value and maximum and minimum valid data */
int **read_grdrasterI(const char *name_grd, float grdnul, float grdmin, float grdmax, struct HeadR s_arrayh[], int ntype_raster, int num_rast)
{
FILE *in;
float datofl[4], dato;
int datfin;
int i, j, contnul, contgod;
int nyfilas, nxcolum;
double resolx, resoly;
int **raster;

	printf("\n***Lectura archivo grd, %s ***\n", name_grd);
	printf("Tipo raster %i\n", ntype_raster);
    if((in=fopen(name_grd,"rb"))==NULL)
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
		* Reading golden software format head in binary
		* En la cabecera las "x" representan las columnas y las "y" las filas
		* Se define siempre array[fila-y][columna-x]
		*/
		if (ntype_raster == 1)
		{
			/**
			* Envia datos leidos para que sean almacenados en la estructura HeadR
			* Se utiliza un puntero para operar con la estructura en la funcion destino
			* En esta se opera normal
			*/
			read_headerB(s_arrayh, num_rast, in);
		}
		if (ntype_raster == 2)
		{
			read_headerA(s_arrayh, num_rast, in);
		}
		/**
		* Creando arrays y asignando valores
		*/
		printf("Creando array\n");
		contnul=0;
		contgod=0;
		nyfilas = s_arrayh[num_rast].hn_fy;
		nxcolum = s_arrayh[num_rast].hn_cx;
		resolx  = s_arrayh[num_rast].hresx;
		resoly  = s_arrayh[num_rast].hresy;
		printf("datos devueltos %i %i %lf %lf\n", nyfilas, nxcolum, resolx, resoly);
		printf("Valid data interval: %f - %f, both included\n", grdmin, grdmax);
		raster = Crea_2DIarray(nyfilas, nxcolum);
		
		for(i=0;i<nyfilas;i++) //filas
		{
			for(j=0;j<nxcolum;j++) //columnas
			{           
				if (ntype_raster == 1) 
				{
					fread(&datofl,sizeof(float),1,in);                      /*!< Lee *.grd en formato binario */
					datfin = (int)datofl[0];
				}
				
				if (ntype_raster == 2) 
				{
					fscanf(in,"%f", &dato);                                /*!< Lee *.grd en formato ASCII */
					datfin = dato;
				}
				if((datfin > grdmax) || (datfin < grdmin))                 /*!< Evalua si el valor esta dentro de los limites definidos */
				{
					datfin = grdnul;                                        /*!< Se asigna valor nulo si esta fuera de los limites */
					contnul+=1;
				}
				else  contgod +=1;                                          /*!< count celdas validas */
				
				raster[i][j]  = datfin;                                /*!< Almacena valor en raster-array */ 			
			}                                          
		}    
		fclose(in);
		/**
		* Print results
		*/
		printf("---------------------------------\n\n");
		printf("Tot cells ok leidas %i\n", contgod);
		printf("Area de trabajo %lf\n", contgod * resolx * resoly);
		printf("Tot cells nulls in DEM = %i\n", contnul);
		printf("Area nula %lf\n", contnul * resolx * resoly);
		printf("end read DEM file\n");  
		printf("---------------------------------\n\n");
	}
    return raster;
}


/*! Read ASCII or Binary double-float *grd raster 
 * ntype_raster: Raster type must be defined (1-binary 2-ASCII)
 * Other inputs parameters: raster name, null value and maximum and minimum valid data */
double **read_grdrasterF(const char *name_grd, float grdnul, float grdmin, float grdmax, struct HeadR s_arrayh[], int ntype_raster, int num_rast)
{
FILE *in;
float datofl[4];
double dato, datfin;
int i, j, contnul, contgod;
int nyfilas, nxcolum;
double resolx, resoly;
//struct HeadR array_headF;
double **raster;

	printf("\n***Lectura archivo grd, %s ***\n", name_grd);
	printf("Tipo raster %i\n", ntype_raster);
	if((in=fopen(name_grd,"rb"))==NULL)
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
		* Reading golden software format head in binary
		* En la cabecera las "x" representan las columnas y las "y" las filas
		* Se define siempre array[fila-y][columna-x]
		*/
		if (ntype_raster == 1)
		{
			/**
			* Envia datos leidos para que sean almacenados en la estructura HeadR
			* Se utiliza un puntero para operar con la estructura en la funcion destino
			* En esta se opera normal
			*/
			read_headerB(s_arrayh, num_rast, in);
		}
		if (ntype_raster == 2)
		{
			read_headerA(s_arrayh, num_rast, in);
		}
		/**
		* Creando arrays y asignando valores
		*/
		printf("Creando array\n");
		contnul=0;
		contgod=0;
		nyfilas = s_arrayh[num_rast].hn_fy;
		nxcolum = s_arrayh[num_rast].hn_cx;
		resolx  = s_arrayh[num_rast].hresx;
		resoly  = s_arrayh[num_rast].hresy;
		printf("datos devueltos %i %i %lf %lf\n", nyfilas, nxcolum, resolx, resoly);
		printf("Valid data interval: %f - %f, both included\n", grdmin, grdmax);
		raster = Crea_2DFarray(nyfilas, nxcolum);
		
		for(i=0;i<nyfilas;i++) //filas
		{
			for(j=0;j<nxcolum;j++) //columnas
			{           
				if (ntype_raster == 1) 
				{
					fread(&datofl,sizeof(float),1,in);                      /*!< Lee *.grd en formato binario */
					datfin = datofl[0];
				}
				if (ntype_raster == 2) 
				{
					fscanf(in,"%lf", &dato);                                /*!< Lee *.grd en formato ASCII */
					datfin = dato;
				}
				if((datfin > grdmax) || (datfin < grdmin))                 /*!< Evalua si el valor esta dentro de los limites definidos */
				{
					datfin = grdnul;                                        /*!< Se asigna valor nulo si esta fuera de los limites */
					contnul+=1;
				}
				else  contgod +=1;                                          /*!< count celdas validas */
				raster[i][j]  = datfin;                                     /*!< Almacena valor en raster-array */ 			
			}                                          
		}    
		fclose(in);
		printf("datos temp entrada: %f %f %f %i %i cont: %i %i Dim %lf %lf\n", grdnul, grdmin, grdmax, ntype_raster, num_rast, contnul, contgod, resoly, resolx);
		/**
		* Print results
		*/
		printf("---------------------------------\n\n");
		printf("Tot cells ok leidas %i\n", contgod);
		printf("Area de trabajo %lf\n", contgod * resolx * resoly);
		printf("Tot cells nulls in DEM = %i\n", contnul);
		printf("Area nula %lf\n", contnul * resolx * resoly);
		printf("end read DEM file\n");  
		printf("---------------------------------\n\n");
	}
	return raster;
}


/*! Read ASCII or Binary double-float *grd raster in clip mode
 * ntype_raster: Raster type must be defined (1-binary 2-ASCII)
 * Other inputs parameters: raster name, null value and maximum and minimum valid data
 * clips requires corner coordinates  */
double **read_grdrasterF_clip(const char *name_grd, float grdnul, float grdmin, float grdmax, \
	struct HeadR s_arrayh[], int ntype_raster, int num_rast, \
	double clip_xlow, double clip_xhig, double clip_ylow, double clip_yhig)
{
FILE *in;
float datofl[4];
double dato, datfin;
int i, j, contnul, contgod;
//int nyfilas, nxcolum;
double resolx, resoly;
//struct HeadR array_headF;
double **raster;
//clip
int big_ncx, big_nfy;
double big_xlow, big_xhig, big_ylow, big_yhig;
int col_min_clip, col_max_clip, row_min_clip, row_max_clip, clip_ncx, clip_nfy;
	
	printf("\n***Lectura archivo grd and CLIP, %s ***\n", name_grd);
	printf("Tipo raster %i\n", ntype_raster);
	if((in=fopen(name_grd,"rb"))==NULL)
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
		* reading input raster header
		*/
		if (ntype_raster == 1)
		{
			read_headerB(s_arrayh, num_rast, in);
		}
		if (ntype_raster == 2)
		{
			read_headerA(s_arrayh, num_rast, in);
		}
		/**
		* Calculating new clip-raster dimension
		*/
		//resolution remain the same
		resolx  = s_arrayh[num_rast].hresx;
		resoly  = s_arrayh[num_rast].hresy;
		printf("Transfering big arraday header \n");
		big_ncx  = s_arrayh[num_rast].hn_cx;         /*! add number of columns in x */
		big_nfy  = s_arrayh[num_rast].hn_fy;         /*! add number of rows in y */
		big_xlow = s_arrayh[num_rast].hxlo;          /*! minimum x-coordinate */
		big_xhig = s_arrayh[num_rast].hxhi;          /*! maximum x-coordinate */
		big_ylow = s_arrayh[num_rast].hylo;          /*! minimum y-coordinate */
		big_yhig = s_arrayh[num_rast].hyhi;          /*! maximum y-coordinate */
		printf("Calculating new clip header \n");
		//**************************************************************
		// check if xy clip-raster limits exceed big-raster dimension */
		if (clip_xlow < big_xlow) clip_xlow = big_xlow;
		if (clip_xhig > big_xhig) clip_xhig = big_xhig;
		if (clip_ylow < big_ylow) clip_ylow = big_ylow;
		if (clip_yhig > big_yhig) clip_yhig = big_yhig;
		// calc clip-raster indexes
		col_min_clip = (int)((clip_xlow - big_xlow) / resolx);
		col_max_clip = (int)((clip_xhig - big_xlow) / resolx);
		row_min_clip = (int)((clip_ylow - big_ylow) / resoly);
		row_max_clip = (int)((clip_yhig - big_ylow) / resoly);
		// check if clip-raste indexes exceed big-raster dimension
		if (col_min_clip < 0) col_min_clip = 0;
		if (row_min_clip < 0) row_min_clip = 0;
		if (col_max_clip >= big_ncx) col_max_clip = big_ncx - 1;
		if (row_max_clip >= big_nfy) row_max_clip = big_nfy - 1;
		// calc total columns and rows of clip-raster
		clip_ncx = col_max_clip - col_min_clip + 1;
		clip_nfy = row_max_clip - row_min_clip + 1;
		// Update clip-raster header
		s_arrayh[num_rast].hn_cx = clip_ncx;
		s_arrayh[num_rast].hn_fy = clip_nfy;
		s_arrayh[num_rast].hxlo  = big_xlow + col_min_clip * resolx;
		s_arrayh[num_rast].hxhi  = s_arrayh[num_rast].hxlo + (clip_ncx - 1) * resolx;
		s_arrayh[num_rast].hylo  = big_ylow + row_min_clip * resoly;
		s_arrayh[num_rast].hyhi  = s_arrayh[num_rast].hylo + (clip_nfy - 1) * resoly;
		//**************************************************************
		printf("New clip header \n");
		printf("array dimension and resolution %i %i %lf %lf\n", s_arrayh[num_rast].hn_fy, s_arrayh[num_rast].hn_cx, resolx, resoly);
		printf("Valid data interval: %f - %f, both included\n", grdmin, grdmax);
		printf("Big-raster dimension: ncx %i nfy %i Clip-raster dimesinon: ncx %i nfy %i\n", big_ncx, big_nfy, clip_ncx, clip_nfy);
		printf("Clip limits: col %i %i row  %i %i\n", col_min_clip, row_max_clip, row_min_clip, row_max_clip);
		//**************************************************************
		printf("Creando clip array \n");
		raster = Crea_2DFarray(clip_nfy, clip_ncx);
		//**************************************************************
		/**
		* adding valid values to the new clip-raster
		*/
		contnul=0;
		contgod=0;
		//read big-raster
		for(i=0;i<big_nfy;i++) //rows of big array
		{
			for(j=0;j<big_ncx;j++) //columns of big array
			{
				//read each data
				if (ntype_raster == 1) 
				{
					if (fread(&datofl[0], sizeof(float), 1, in) != 1)   /*!< reading *.grd in binary format */
					{
						printf("Error reading data\n");
						exit(0);
					}    
					dato = (double)datofl[0];
				}
				if (ntype_raster == 2) 
				{
					if (fscanf(in, "%lf", &dato) != 1)                  /*!< reading *.grd in ASCII format */
					{
						printf("Error reading ASCII data\n");   
						exit(0);
					}
					//fscanf(in,"%lf", &dato);
				}
				datfin = dato;
				
				//check if data is in limits
				if((datfin > grdmax) || (datfin < grdmin))              /*!< check if cell value is within limits */
				{
					datfin = grdnul;                                    /*!< add null value if not */
					contnul++;
				}
				else  contgod ++;
				//when data are saved in clipped array
				
				/*! if cell is within clip-raster limits add value */
				if ((i >= row_min_clip && i <= row_max_clip) && (j >= col_min_clip && j < col_max_clip))
				{
					int i_clip = i - row_min_clip;
					int j_clip = j - col_min_clip; 
					raster[i_clip][j_clip]  = datfin;
					//printf("dat %lf : big %i %i : clip i %i j %i\n", datfin, i, j, i_clip, j_clip);
				}
			}
		}
		fclose(in);
		printf("Input raster parameters: %f %f %f %i %i cont: %i %i Dim %lf %lf\n", grdnul, grdmin, grdmax, ntype_raster, num_rast, contnul, contgod, resoly, resolx);
		
		/**
		* Print header again
		*/
		printf("New .grd Clip Raster header:\n");
		printf("raster type: %i\n", ntype_raster);
		printf("colx  = %5i -- rowy = %5i\n",         s_arrayh[num_rast].hn_cx, s_arrayh[num_rast].hn_fy);
		printf("xmin = %f -- xmax = %f\n",            s_arrayh[num_rast].hxlo, s_arrayh[num_rast].hxhi);
		printf("ymin = %f -- ymax = %f\n",            s_arrayh[num_rast].hylo, s_arrayh[num_rast].hyhi);
		printf("zlo = %10.4f -- zhi = %10.4f\n",      s_arrayh[num_rast].hzlo, s_arrayh[num_rast].hzhi);
		printf("Raster x-resolution resx = %f\n",     s_arrayh[num_rast].hresx);
		printf("Raster y-resolution resy = %f\n",     s_arrayh[num_rast].hresy);
		/**
		* Print results
		*/
		printf("---------------------------------\n\n");
		printf("Tot cells ok leidas %i\n", contgod);
		printf("Area de trabajo %lf\n", contgod * resolx * resoly);
		printf("Tot cells nulls in DEM = %i\n", contnul);
		printf("Area nula %lf\n", contnul * resolx * resoly);
		printf("end read DEM file\n");  
		printf("---------------------------------\n\n");
	}
	return raster;
}



//**********************************************************************
/*! Escribe raster-array */
//**********************************************************************

/*! Escribe grd raster 
 * Requiere el numero de raster a utilizar (1-3), el formato (1-entero 2-double) y el tipo (1-binario 2-ASCII) */
int write_grdrasterI(const char *name_grd, int **raster, struct HeadR s_arrayh[], int num_rast, int ntype_raster)
{
char   buffer[255], tipofile[10];
short int buff_int[4];
double buff_double[32];
int *buff_datint;
int i, j, ncxfin, nfyfin;
FILE *out;
    
	printf("Escribe raster generado %s\n", name_grd);		
	/**
	* Captura valores para cabecera
	*/
	if (ntype_raster == 1)
	{
		sprintf(buffer,"DSBB");
		sprintf(tipofile, "wb");            /*!< Tipo de apertura de archivo, escritura en binario */
	}
	if (ntype_raster == 2)
	{
		sprintf(buffer,"ASBB"); 
		sprintf(tipofile, "wt");            /*!< Tipo de apertura de archivo, escritura en ascii */
	}
	buff_int[0]    = s_arrayh[num_rast].hn_cx;
	buff_int[1]    = s_arrayh[num_rast].hn_fy;
	buff_double[0] = s_arrayh[num_rast].hxlo;
	buff_double[1] = s_arrayh[num_rast].hxhi;
	buff_double[2] = s_arrayh[num_rast].hylo;
	buff_double[3] = s_arrayh[num_rast].hyhi;
	buff_double[4] = s_arrayh[num_rast].hzlo;
	buff_double[5] = s_arrayh[num_rast].hzhi;
	ncxfin = s_arrayh[num_rast].hn_cx;
	nfyfin = s_arrayh[num_rast].hn_fy;
	
	/**
	* Crea y abre archivo de salida 
	*/
	if((out=fopen(name_grd,tipofile))==NULL)
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
		printf("%s\n", buffer);
		printf("%i %i\n", buff_int[0], buff_int[1]);
		printf("%lf %lf\n", buff_double[0], buff_double[1]);
		printf("%lf %lf\n", buff_double[2], buff_double[3]);
		printf("%lf %lf\n", buff_double[4], buff_double[5]);
		if (ntype_raster == 1)
		{
			/*!> Crea buffer dinamico para almacenar las filas de datos */
			buff_datint = (int*)malloc(ncxfin*sizeof(int));
			/**
			* Escribe datos de cabecera
			*/
			fwrite(buffer, 4, 1, out);
			fwrite(buff_int,    sizeof(short int)*2, 1, out);
			fwrite(buff_double, sizeof(double)*6, 1, out);
			/**
			* Lee y escribe linea a linea 
			*/
			for (i=0;i<nfyfin;i++)
			{
				for(j=0;j<ncxfin;j++)
				{
					buff_datint[j] = (int)raster[i][j];       /*!< Si es entero, almacena toda la linea en un buffer */	
					/*
					//TEST**********************************************
					if (j <  ncxfin)
					{	
						if (buff_datint[j] > 0) printf("%i", buff_datint[j]); //TEST
					}
					else 
					{
						if (buff_datint[j] > 0) printf("%i\n", buff_datint[j]);
					}
					*/
				}
				fwrite(buff_datint, sizeof(int)*ncxfin, 1, out);       /*!< Si es entero, escribe el buffer en el archivo de salida */
			}
		}	
		if (ntype_raster == 2)
		{
			/**
			* Escribe datos de cabecera
			*/			
			fprintf(out,"%s\n", buffer);
			fprintf(out,"%i %i\n", buff_int[0], buff_int[1]);
			fprintf(out,"%lf %lf\n", buff_double[0], buff_double[1]);
			fprintf(out,"%lf %lf\n", buff_double[2], buff_double[3]);
			fprintf(out,"%lf %lf\n", buff_double[4], buff_double[5]);
			for (i=0;i<nfyfin;i++)
			{
				for(j=0;j<ncxfin;j++)
				{
					fprintf(out,"%i ",raster[i][j]);
				}
				fprintf(out,"\n");
			}
		}			 
	}
	fclose(out);
	printf("Finaliza escritura raster generado\n");
	return 1;
}

/*! Escribe grd raster 
 * Requiere el numero de raster a utilizar (1-3), el formato (1-entero 2-double) y el tipo (1-binario 2-ASCII) */
int write_grdrasterF(const char *name_grd, double **raster, struct HeadR s_arrayh[], int num_rast, int ntype_raster)
{
char   buffer[255], tipofile[10];
short int buff_int[4];
double buff_double[32];
float  *buff_datdob;
int i, j, ncxfin, nfyfin;
FILE *out;
    
	printf("Write raster %s\n", name_grd);		
	/**
	* Captura valores para cabecera
	*/
	if (ntype_raster == 1)
	{
		sprintf(buffer,"DSBB");
		sprintf(tipofile, "wb");            /*!< Tipo de apertura de archivo, escritura en binario */
	}
	if (ntype_raster == 2)
	{
		sprintf(buffer,"ASBB"); 
		sprintf(tipofile, "w");            /*!< Tipo de apertura de archivo, escritura en ascii */
	}
	buff_int[0]    = s_arrayh[num_rast].hn_cx;
	buff_int[1]    = s_arrayh[num_rast].hn_fy;
	buff_double[0] = s_arrayh[num_rast].hxlo;
	buff_double[1] = s_arrayh[num_rast].hxhi;
	buff_double[2] = s_arrayh[num_rast].hylo;
	buff_double[3] = s_arrayh[num_rast].hyhi;
	buff_double[4] = s_arrayh[num_rast].hzlo;
	buff_double[5] = s_arrayh[num_rast].hzhi;
	ncxfin = s_arrayh[num_rast].hn_cx;
	nfyfin = s_arrayh[num_rast].hn_fy;
	
	/**
	* Crea y abre archivo de salida 
	*/
	if((out=fopen(name_grd,tipofile))==NULL)
	{
		printf("-------ERROR open file--------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("ATENCION: LA RUTA DE SALIDA DEL ARCHIVO NO ES CORRECTA\n");
		exit(0);
	}
	else
	{
		printf("%s\n", buffer);
		printf("%i %i\n", buff_int[0], buff_int[1]);
		printf("%lf %lf\n", buff_double[0], buff_double[1]);
		printf("%lf %lf\n", buff_double[2], buff_double[3]);
		printf("%lf %lf\n", buff_double[4], buff_double[5]);
		if (ntype_raster == 1)
		{
			/*!> Crea buffer dinamico para almacenar las filas de datos */
			buff_datdob = (float *)malloc(sizeof(float)*ncxfin);
			/**
			* Escribe datos de cabecera
			*/
			fwrite(buffer, 4, 1, out);
			fwrite(buff_int,    sizeof(short int)*2, 1, out);
			fwrite(buff_double, sizeof(double)*6, 1, out);
			/**
			* Lee y escribe linea a linea 
			*/
			for (i=0;i<nfyfin;i++)
			{
				for(j=0;j<ncxfin;j++)
				{
					buff_datdob[j] = (float)raster[i][j];                  /*!< Si es double, almacena toda la linea en un buffer */
					/*
					//TEST**********************************************
					if (j <  ncxfin)
					{	
						if (buff_datdob[j] > 0) printf("%lf", buff_datdob[j]); 
					}
					else 
					{
						if (buff_datdob[j] > 0) printf("%lf\n", buff_datdob[j]);
					}
					*/
				}

				fwrite(buff_datdob, sizeof(float)*ncxfin, 1, out);     /*!< Si es double, escribe el buffer en el archivo de salida */
				//printf("\n");
			}
		}
		if (ntype_raster == 2)
		{
			/**
			* Escribe datos de cabecera
			*/
			fprintf(out,"%s\n", buffer);
			fprintf(out,"%i %i\n", buff_int[0], buff_int[1]);
			fprintf(out,"%lf %lf\n", buff_double[0], buff_double[1]);
			fprintf(out,"%lf %lf\n", buff_double[2], buff_double[3]);
			fprintf(out,"%lf %lf\n", buff_double[4], buff_double[5]);
			for (i=0;i<nfyfin;i++)
			{
				for(j=0;j<ncxfin;j++)
				{
					fprintf(out,"%lf ",raster[i][j]);
				}
				fprintf(out,"\n");
			}
		}			 
	}
	fclose(out);
	printf("Finaliza escritura raster generado\n");
	return 1;
}


//**********************************************************************
/*! Libera raster-array */
//**********************************************************************

/*! Libera memoria del raster 2D entero */
void free_RasterI(const char *donde, int **raster, int nfilas)
{
int i;
	printf("Liberando raster en %s\n", donde);
	for(i=0;i<nfilas;i++)
	{
		free(raster[i]);
	}
	free(raster);
	raster = NULL;
	printf("raster liberado\n");
}

/*! Libera memoria del raster 2D doble */
void free_RasterF(const char *donde, double **raster, int nfilas)
{
int i;
	printf("Liberando raster en %s\n", donde);
	for(i=0;i<nfilas;i++)
	{
		free(raster[i]);
	}
	free(raster);
	raster = NULL;
	printf("raster liberado\n");
}


/*! Libera cabecera raster */
void free_Headraster(struct HeadR s_arrayh[], int num_rast)
{
	s_arrayh[num_rast].hn_cx      = 0;
	s_arrayh[num_rast].hn_fy      = 0;
	s_arrayh[num_rast].hxlo       = 0;
	s_arrayh[num_rast].hxhi       = 0;
	s_arrayh[num_rast].hylo       = 0;
	s_arrayh[num_rast].hyhi       = 0;
	s_arrayh[num_rast].hzlo       = 0;
	s_arrayh[num_rast].hzhi       = 0;
	s_arrayh[num_rast].hresx      = 0;
	s_arrayh[num_rast].hresy      = 0;
	s_arrayh[num_rast].hinvresx   = 0;
	s_arrayh[num_rast].hinvresy   = 0;
	s_arrayh[num_rast].hresxx     = 0;
	s_arrayh[num_rast].hresyy     = 0;
	s_arrayh[num_rast].htotdemcel = 0;
	s_arrayh[num_rast].hperim     = 0;
	s_arrayh[num_rast].harea      = 0;
}

#endif /* _LEEARR_H */




/*!
 * UPDATEs
 * 2025/12/03: Update clip function based on IA recommendation to solve raster dimension
 * 2025/08/11: New read raster function for clip option
 */
