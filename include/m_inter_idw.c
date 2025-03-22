/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza 
* Module:  m_inter_idw.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2023-01-31
* Last update: 2023-01-31
* 
* Description:
* 
* Interpolate points in a raster using Inverse Distance Weighting (IDW)
* parameter p to control the weight of known pt
* parameter distance to control searching radius and total numbers considered
* in each assessment. The higher the softer the result.
* If needed, data could be differentiated by seven types:
* Type 1: point over continuous isopleth/lines (continuous)
* Type 2: point over discontinous isopleth/lines less confident (discontinuous)
* Type 3: point on site - a measure taken in specific locations (onsite)
* Type 4: point over interpreted-added isopleth/line to complete the previous (interpret)
* Type 5: point over boundary line (boundaryline)
* Type 6: discretize boundary point (boundarypoint)
* Type 7: removable point - config file allows to select set of data * (rm-xxxx)
* *********************************************************************/
#ifndef _INT_IDW
#define _INT_IDW

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


/*----------------------------------------------------------------------
 |  Section Global variables
 ---------------------------------------------------------------------*/
/*!< The raster is defined in main module */
extern double **rast_inter;

/*!< By closet points, not ready yet */
void calc_idwbypts(struct HeadR s_arrayh[], float dem_nulval, coordCSV s_coord[], int n_xyz, int idw_searpts, float idw_power)
{
int i,j,k,l,p;
int pti, ptj, cont, comp, id, nfirst, com_id, cont_cells;
double ptx, pty, ptz, cellx, celly, dist, com_dist;
double inptz, sum_sup, sum_inf;
double **newxyz;
	
	cont_cells=0;
	p=100;
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) 
		{
			if (rast_inter[i][j] != dem_nulval)
			{
				cellx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, j);  /*!< getting the coordinate from raster index */
				celly = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, i);
				
				for(k=0;k<n_xyz;k++)                                        /*!< for each known point */
				{
					ptx = s_coord[k].ptxcoor;
					pty = s_coord[k].ptycoor;
					ptz = s_coord[k].ptzcoor;
					s_coord[k].ptdist  = 0;
					s_coord[k].ptjerar = 0;
					/*!< calc index in raster */
					ptj = calc_Rindex(ptx, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
					pti = calc_Rindex(pty, s_arrayh[0].hylo, s_arrayh[0].hresy);
					/*!< if pt is in cell get value */
					if (ptj == j && pti == i)
					{
						rast_inter[i][j] = ptz;
						break;
					}
					else
					{
						/*!< calc dist */
						dist = calc_dist(cellx, ptx, celly, pty);
						s_coord[k].ptdist = dist;
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						/*!< compare  distance */
						for(l=0;l<k;l++)                                    /*!< only for those pt already checked */
						{
							cont = 0;
							if (cont == 0)                                  /*!< only for those pt already checked */
							{
								com_id   = s_coord[l].ptjerar;
								com_dist = s_coord[l].ptdist;
								comp = compara_val(dist, com_dist);
								
								if (l == 0)
								{
									if (comp == 1) 
									s_coord[l].ptjerar = com_id;
									s_coord[l].ptjerar++;
									cont = 1;
								}
								if (l > 0)
								{
									if (comp == 3) 
									{
										s_coord[l].ptjerar = com_id;
										s_coord[l].ptjerar++;
										cont = 1;
									}
								}
							}
							else
							{
								if (s_coord[l].ptjerar > com_id) s_coord[l].ptjerar++;
							}
						}
					}
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				/*!< sort the xyz by dist */
				newxyz = Crea_2DFarray(n_xyz, 2);
				for(k=0;k<n_xyz;k++)
				{
					ptz  = s_coord[k].ptzcoor;
					dist = s_coord[k].ptdist;
					id   = s_coord[k].ptjerar;
					if (id < n_xyz)
					{
						newxyz[id][0] = ptz;
						newxyz[id][1] = dist;
					}
					else
					{
						printf("Error, Id %i is higher than total expected pt %i\n", id, n_xyz);
						exit(0);
					}
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				nfirst = idw_searpts / 2;
				sum_sup = sum_inf = 0;
				cont = 0;
				for(k=0;k<n_xyz;k++)
				{
					ptz  = s_coord[k].ptzcoor;
					dist = s_coord[k].ptdist;
					if (k == 0) inptz = ptz;
					if (k < nfirst)
					{
						sum_sup =  ptz / pow(dist, idw_power);
						sum_inf =  1 / pow(dist, idw_power);
					}
					if (k >= nfirst)
					{
						if (ptz != inptz)
						{
							sum_sup =  ptz / pow(dist, idw_power);
							sum_inf =  1 / pow(dist, idw_power);
							cont++;
						}
						if (cont==nfirst)break;
					}
				}
				rast_inter[i][j] = sum_sup / sum_inf;
				
				cont_cells++;
			}
			if (cont_cells == p)
			{
				printf("Tot cells process %i, total cells %i (including nulls)\n", cont_cells, s_arrayh[0].htotdemcel);
				p = p*2;
			}
		}
	}
}


/*!< Only by distance, min searching pts, the value remain */
void calc_idwbydist(struct HeadR s_arrayh[], float dem_nulval, coordCSV s_coord[], int n_xyz, float idw_distan, int idw_searpts, float idw_power)
{
int i,j,k,p;
int cont, cont_cells;
double ptx, pty, ptz, cellx, celly, dist;
double sum_sup, sum_inf, sum_control;
	
	cont_cells=0;
	p=100;
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) 
		{
			if (rast_inter[i][j] != dem_nulval)
			{
				cellx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, j);  /*!< getting the coordinate from raster index */
				celly = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, i);                           /*!< Init search distance */
				sum_control = sum_sup = sum_inf = cont = 0;
				for(k=0;k<n_xyz;k++)                                /*!< for each known point */
				{
					ptx  = s_coord[k].ptxcoor;
					pty  = s_coord[k].ptycoor;
					ptz  = s_coord[k].ptzcoor;
					//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					dist = calc_dist(cellx, ptx, celly, pty);       /*!< calc distance cell - point */
					if (dist <= idw_distan)                     /*!< If distance is in radious search */
					{
						sum_control += ptz;
						sum_sup     +=  ptz / pow(dist, idw_power);     /*!< sum zval / dist^ for each searched point  */
						sum_inf     +=  1 / pow(dist, idw_power);       /*!< sum 1 / dist^ for each searched point  */
						cont++;
					}
					//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					if (idw_searpts > 0)                          /*!< End loop if found points matchs */
					{
						if (cont == idw_searpts)break;
					}
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				if (cont >0) 
				{
					if (sum_control == 0)rast_inter[i][j] = 0;          /*!< if all z search points is 0, then idw = 0 */
					else rast_inter[i][j] = sum_sup / sum_inf;          /*!< Calc interpolated value */
				}
				else rast_inter[i][j] = 0;
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				cont_cells++;
			}
			if (cont_cells == p)
			{
				printf("Tot cells process %i, total cells %i (including nulls)\n", cont_cells, s_arrayh[0].htotdemcel);
				p = p*2;
			}
		}
	}						
}


/*!< By distance, increase searching distance if not enough points are found */
void calc_idwbydistup(struct HeadR s_arrayh[], float dem_nulval, coordCSV s_coord[], int n_xyz, float idw_distan, int idw_searpts, float idw_power, int idw_bytype)
{
int i,j,k,p;
int cont, cont_cells, type;
double ptx, pty, ptz, cellx, celly, dist;
double sum_sup, sum_inf, sum_control, idw_distancalc;
	
	cont_cells=0;
	p=100;
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) 
		{
			if (rast_inter[i][j] != dem_nulval)
			{
				cellx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, j);  /*!< getting the coordinate from raster index */
				celly = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, i);
				idw_distancalc  = idw_distan;                           /*!< Init search distance */
				do
				{
					sum_control = sum_sup = sum_inf = cont = 0;
					for(k=0;k<n_xyz;k++)                                /*!< for each known point */
					{
						ptx  = s_coord[k].ptxcoor;
						pty  = s_coord[k].ptycoor;
						ptz  = s_coord[k].ptzcoor;
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						if (idw_bytype > 0)
						{
							type = s_coord[i].ptype;
							if (idw_bytype == 1)
							{
								if (type == 4)idw_power = idw_power / 2;
							}
						}
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						dist = calc_dist(cellx, ptx, celly, pty);       /*!< calc distance cell - point */
						if (dist <= idw_distancalc)                     /*!< If distance is in radious search */
						{
							sum_control += ptz;
							sum_sup     +=  ptz / pow(dist, idw_power);     /*!< sum zval / dist^ for each searched point  */
							sum_inf     +=  1 / pow(dist, idw_power);       /*!< sum 1 / dist^ for each searched point  */
							cont++;
						}
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						if (idw_searpts > 0)                          /*!< End loop if found points matchs */
						{
							if (cont == idw_searpts)break;
						}
					}
					if (idw_searpts > 0)                                /*!< check if loop end with less expected points */
					{
						if (cont < idw_searpts)cont= 0;                 /*!< if found points < min point, try again with higher radioush */
					}
					idw_distancalc = idw_distancalc + idw_distancalc;   /*!< If no point, increase radious search */
					
				}while(cont<1);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				if (cont >0) 
				{
					if (sum_control == 0)rast_inter[i][j] = 0;          /*!< if all z search points is 0, then idw = 0 */
					else rast_inter[i][j] = sum_sup / sum_inf;          /*!< Calc interpolated value */
				}
				else rast_inter[i][j] = 0;
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				cont_cells++;
			}
			if (cont_cells == p)
			{
				printf("Tot cells process %i, total cells %i (including nulls)\n", cont_cells, s_arrayh[0].htotdemcel);
				p = p*2;
			}
		}
	}						
}



#endif /* _INT_IDW */
