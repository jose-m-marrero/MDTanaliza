/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> and Hugo Yepes <
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza
* Module:  m_topohazard.h
* 
* Authors: Jose M. Marrero 
*          Hugo Yepes
* 
* Version: 0.1.1
* Creation Date: 2018-04-06
* Last update: 2025-02-19
* 
* Description:
* 
* Calc differential of topography between a given point and all cells in a raster following:
* direc_calc 1: horizontal
* direc_calc 2: vertical
* direc_calc 3: oblique
* from_center 1: to get differences for Zero Level Trajectory to the maximum extension
* If maximum distance is 0, all cells in DEM will be processed from ZLP. If >0 calculation will stop when distance in reached 
* from_limask 1: Mask is needed. Two output raster in mask, and from mask limits
* Topohazard calculate first in one direction and then the opposite
* Important: if mask is used, the amount of RAM will increase by three times
* 
* *********************************************************************/
#ifndef _TOPO_RAS
#define _TOPO_RAS

/*!
* MDTanaliza's general modules 
***********************************************************************/ 
/*! Define global variables */
#include "../include/m_global.h"
/*! Define structures */
#include "../include/struc_mdtanaliza.h"
/*! Common calculation functions */
#include "../include/m_comuntools.h"
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

/*----------------------------------------------------------------------
 |  Section Global variables
 ---------------------------------------------------------------------*/

/*!< Topohazard needs a raster-base to save results
	 * This raster is kept until last ZLP has been processed
	 * The raster is raster every new trajectory 
	 * The raster are defined in main module */
extern double **raster_topohaz; 
extern double **inhaz_topo, **outhaz_topo;

//**********************************************************************
/*! Auxiliar Functions */
//**********************************************************************


/*! Check if cell was previously procesed */
int check_ifriberbed(int nrio, int idxfy, int idxcx, struct Zeropoint s_bedriver[], int n_bedxyz)
{
int i;
int idx_fy, idx_cx, n_rio;
int if_center;
	if_center = 0;
	for(i=0;i<n_bedxyz;i++)
	{
		idx_fy = s_bedriver[i].bidxfy;
		idx_cx = s_bedriver[i].bidxcx;
		n_rio  = s_bedriver[i].brio;
		if (n_rio == nrio)
		{
			if (idxfy == idx_fy && idxcx == idx_cx)if_center = 1;
		}
	}
	return if_center;
}

/*! Check if cell is in mask */
int check_ifinmask(int **mask, struct HeadR s_arrayh[], int idx_fy, int idx_cx)
{
int idxm_fy, idxm_cx, mask_bed;
double txm, tym;
//double txm2, tym2;
	mask_bed = 0;
	/*!< get coordinates from topo DEM raster*/
	txm = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_cx);
	tym = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
	/*!< When mask and DEM have different size, coordinate must be checked first to avoid invalid idex, 
	 * Check is in mask raster */
	if(tym > s_arrayh[1].hylo && tym <= s_arrayh[1].hyhi && txm > s_arrayh[1].hxlo && txm <= s_arrayh[1].hxhi)
	{
		/*!< get index from coordinates over mask raster*/
		idxm_fy = abs(calc_Rindex(tym, s_arrayh[1].hylo, s_arrayh[1].hresy));
		idxm_cx = abs(calc_Rindex(txm, s_arrayh[1].hxlo, s_arrayh[1].hresx));
		//txm2 = get_coor(s_arrayh[1].hxlo, s_arrayh[1].hresx, idxm_cx);
		//tym2 = get_coor(s_arrayh[1].hylo, s_arrayh[1].hresy, idxm_fy);
		
		//printf("mask %i %i :: %i %i\n", idxm_fy, idxm_cx, s_arrayh[1].hn_fy, s_arrayh[1].hn_cx);
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< chek if index are inside mask raster */
		if (idxm_fy > 0 && idxm_fy < s_arrayh[1].hn_fy && idxm_cx > 0 && idxm_cx < s_arrayh[1].hn_cx)
		{
			mask_bed  = mask[idxm_fy][idxm_cx];
			//printf("mask val %i \n", mask_bed);
		}
		else mask_bed = -5555; /*!< cell is out of mask raster */
	}
	else mask_bed = -5555; /*!< cell is out of mask raster */
	
	//printf("%i ::: xy %lf %lf :: idx %i %i : oriidx %i %i :: orixy %lf %lf\n", mask_bed, txm2, tym2, idxm_fy, idxm_cx, idx_fy, idx_cx, txm, tym);
	
	return mask_bed;
}

/*! Calculation of differential of topography */
double calc_diff(double **raster, double z_cell, double z_zlp, float dem_nulval, int idx_fy, int idx_cx, int typcal)
{
double diff_topo;
	
	diff_topo=0;
	if(z_cell != dem_nulval)                                  
	{
		diff_topo = z_cell - z_zlp;                                       /*!< Calc z diff between zlp and next cells */
	}
	else                                                                    /*!< If z value is null */
	{
		diff_topo = get_newzval(raster, dem_nulval, idx_fy, idx_cx, z_zlp, typcal); /*!< Get z value from 8 cells if they are not null */
	}
	return diff_topo;
}

//**********************************************************************
/*! Tophazard Functions */
//**********************************************************************

/*! Calculate horizontal (East-West), vertical (North-South) or oblique direction
 * direc_calc control vertical or horizontal, no efect on oblique
 * direc_count control positive or negative cell count, from ZLP to the edge of the raster or defined distance
 * */
void central_topohazard(double **raster, struct HeadR s_arrayh[], float dem_nulval, \
		int idxfy, int idxcx, double txori, double tyori, double tzori, struct Zeropoint s_bedriver[], int n_bedxyz, \
		int nrio, double raddx, double raddy, int direc_calc, int direc_count, float topmax_dist)
{
int n, idx_fy, idx_cx;
int count_loops, out_riberbed, cont_err, cont_fix, in_demlimit;
double z_cell, diff_topo;
double p4x, p4y;
int tot_idxcx, tot_idxfy;
float raster_discre, crossdist;
	n = count_loops = cont_fix = cont_err = 0;
	raster_discre = s_arrayh[0].hresx / 2;
	tot_idxcx = s_arrayh[0].hn_cx - 2;
	tot_idxfy = s_arrayh[0].hn_fy - 2;
	idx_fy = idxfy;
	idx_cx = idxcx;
	raster_topohaz[idx_fy][idx_cx] = 0;
	in_demlimit = 1;
	/**
	* Moving from center to left/right, up/dow or oblique
	*/
	do
	{
		if (direc_calc == 1)                                            /*!< horizontal (East-West) direction */
		{
			if (direc_count == 1)idx_cx--;                              /*!< Reducing row value to move from center to the left (West) */   //<----------------------
			if (direc_count == 2)idx_cx++;                              /*!< Reducing row value to move from center to the right (East) */
			p4x = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_cx);
			p4y = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
		}
		if (direc_calc == 2)                                            /*!< Vertical (North-South) direction */
		{
			if (direc_count == 1)idx_fy--;                              /*!< Reducing row value to move from center to botton (South) */
			if (direc_count == 2)idx_fy++;                              /*!< Reducing row value to move from center to top (North) */
			p4x = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_cx);
			p4y = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
		}
		if (direc_calc == 3)
		{
			if(direc_count == 1)n++;                                        /*!< if positive count */
			if(direc_count == 2)n--;
			p4x  = txori + (raddx * raster_discre * n);  /*!< Get oblique coordinates */
			p4y  = tyori + (raddy * raster_discre * n);
			idx_cx = calc_Rindex(p4x, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
			idx_fy = calc_Rindex(p4y, s_arrayh[0].hylo, s_arrayh[0].hresy);
			//printf("%lf %lf -- %i %i\n", p4x, p4y, idx_cx, idx_fy);
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< Check if index if they are in the raster */
		if(idx_fy > 2 && idx_fy < tot_idxfy -2  && idx_cx > 2 && idx_cx < tot_idxcx - 2)     /*!< far from raster limits, just in case is needed  to find 8 cells */
		{
			in_demlimit = 1;   
			//printf("in_limit %i\n", in_limits);
			//printf("%i %i : %i %i\n", idx_fy, tot_idxfy, idx_cx, tot_idxcx);
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			/*!< If we are in a 0 cell value over the bedriver */
			if (count_loops < 10)
			{
				/*!< Check if new cell belong to the riberbed */
				out_riberbed = check_ifriberbed(nrio, idx_fy, idx_cx, s_bedriver, n_bedxyz);
			}
			else out_riberbed = 0;
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			/*!< if we are far from bedriver */
			if (out_riberbed == 0)
			{
				z_cell  = raster[idx_fy][idx_cx];                                          /*!< Get new topo value */
				diff_topo = calc_diff(raster, z_cell, tzori, dem_nulval, idx_fy, idx_cx, 2);  /*!< calc diff topo */
				raster_topohaz[idx_fy][idx_cx] = diff_topo;
			}
			else raster_topohaz[idx_fy][idx_cx] = 0;
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (topmax_dist > 0 && count_loops > 1) 
			{
				crossdist = calc_dist(p4x, txori, p4y, tyori);
				if (crossdist > topmax_dist) break;                         /*!< break if distance > than max distance */
			}
		}
		else 
		{
			in_demlimit = 0; 
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		count_loops++;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	}while(in_demlimit > 0);                                               /*!< While in raster  */
	printf("Loops resolved with nloops %i and count direction %i\n", count_loops, direc_count);
	//printf("Cells fixed with 8 neighbour cells %i and no fixed %i\n", cont_fix, cont_err);
	printf("-----------------------------------------------------------\n");
}


/*! Calculate horizontal (East-West) or vertical (North-South) direction
 * direc_calc control vertical or horizontal 
 * direc_count control positive or negative cell count, 
 * IN from ZLP to the edge of the mask
 * OUT from mask's limits to the edge of the raster or defined distance
 * */
void mask_topohazard(double **raster, int **mask, struct HeadR s_arrayh[], float dem_nulval, int mas_nulval, \
		int idxfy, int idxcx, double txori, double tyori, double tzori, struct Zeropoint s_zerolevpt[], int n_bedxyz, \
		int nrio, double raddx, double raddy, int direc_calc, int direc_count, float topmax_dist)
{
int n, idx_fy, idx_cx;
int count_loops, out_riberbed, count_out, cont_infix, cont_inerr, cont_outfix, cont_outerr, in_demlimit;
int inmask_count, outmask_count, outlimask_count;
double z_cell, diff_topo;
int tot_idxfy, tot_idxcx;
float raster_discre, crossdist;
double p4x, p4y;
int mask_bed;
double tzori_lim, txori_lim, tyori_lim;

	n = count_loops = cont_infix = cont_outfix =  cont_inerr = cont_outerr = 0;
	raster_discre = s_arrayh[0].hresx / 2;
	tot_idxcx = s_arrayh[0].hn_cx - 2;
	tot_idxfy = s_arrayh[0].hn_fy - 2;
	idx_fy = idxfy;
	idx_cx = idxcx;
	raster_topohaz[idx_fy][idx_cx] = 0;
	inhaz_topo[idx_fy][idx_cx]     = 0;                                     /*!< Add diff topo to hazin */
	outhaz_topo[idx_fy][idx_cx]    = 0;
	in_demlimit = 1;
	count_out   = 0;
	tzori_lim   = -9999;
	inmask_count = outmask_count = outlimask_count = 0;
	//---
	inhaz_topo[idx_fy][idx_cx]  = 0;
	//outhaz_topo[idx_fy][idx_cx] = 0;
	/**
	* Moving from center to left/right or up/dow
	*/
	do
	{
		if (direc_calc == 1)                                            /*!< horizontal (East-West) direction */
		{
			if (direc_count == 1)idx_cx--;                              /*!< Reducing row value to move from center to the left (West) */   //<----------------------
			if (direc_count == 2)idx_cx++;                              /*!< Reducing row value to move from center to the right (East) */
			p4x = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_cx);
			p4y = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
		}
		if (direc_calc == 2)                                            /*!< Vertical (North-South) direction */
		{
			if (direc_count == 1)idx_fy--;                              /*!< Reducing row value to move from center to botton (South) */
			if (direc_count == 2)idx_fy++;                              /*!< Reducing row value to move from center to top (North) */
			p4x = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_cx);
			p4y = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);;
		}
		if (direc_calc == 3)
		{
			if(direc_count == 1)n++;                                        /*!< if positive count */
			if(direc_count == 2)n--;
			p4x  = txori + (raddx * raster_discre * n);  /*!< Get oblique coordinates */
			p4y  = tyori + (raddy * raster_discre * n);
			idx_cx = calc_Rindex(p4x, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
			idx_fy = calc_Rindex(p4y, s_arrayh[0].hylo, s_arrayh[0].hresy);
		}	
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< if we are still in dem raster */
		if(idx_fy > 2 && idx_fy < tot_idxfy && idx_cx > 2 && idx_cx < tot_idxcx)
		{
			//printf("%i %i :: %i %i\n", idx_fy, tot_idxfy, idx_cx, tot_idxcx);
			in_demlimit = 1;                                            
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			/*!< chek if in mask */
			mask_bed  = check_ifinmask(mask, s_arrayh, idx_fy, idx_cx);
			//if (count_loops < 50)printf("valmask %i  --  masknul %i\n", mask_bed, mas_nulval);
			/*!< if inside mask raster */
			if  (mask_bed != -5555)                                 
			{
				//printf("L %i valmask %i  %lf \n", count_loops, mask_bed, z_cell);
				/*!< Is inside mask*/
				if  (mask_bed != mas_nulval)    
				{
					//printf("TEST: calc inside with mask no null %i\n", mask_bed);
					/*!< If near original pt */
					if (count_loops < 10)
					{
						/*!< Check if new cell belong to the riberbed */
						out_riberbed = check_ifriberbed(nrio, idx_fy, idx_cx, s_zerolevpt, n_bedxyz);
					}
					else out_riberbed = 0;
					//printf("out_riberbed %i\n", out_riberbed);
					//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					/*!< if we are far from center */
					if (out_riberbed == 0)
					{
						z_cell  = raster[idx_fy][idx_cx];                                          /*!< Get new topo value */
						diff_topo = calc_diff(raster, z_cell, tzori, dem_nulval, idx_fy, idx_cx, 2);  /*!< calc diff topo */
						inhaz_topo[idx_fy][idx_cx]  = diff_topo;                                     /*!< Add diff topo to hazin */
						outhaz_topo[idx_fy][idx_cx] = dem_nulval;                                    /*!< While in mask,  outmask is null */
						inmask_count++;
					}
					else 
					{
						inhaz_topo[idx_fy][idx_cx]  = 0;
					}
					
				}					
				/*!< Is out mask but inside raster mask*/
				if  (mask_bed == mas_nulval)   
				{
					//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					z_cell  = raster[idx_fy][idx_cx];                                             /*!< Get new topo value */
					if(z_cell != dem_nulval && tzori_lim == -9999)                                /*!< Get z in mask limit */
					{
						/*!< get init values of mask limits */
						tzori_lim = z_cell;
						txori_lim = p4x;
						tyori_lim = p4y;
					}
					
					diff_topo = calc_diff(raster, z_cell, tzori_lim, dem_nulval, idx_fy, idx_cx, 2);  /*!< calc diff topo */
					outhaz_topo[idx_fy][idx_cx] = diff_topo;                                         /*!< Add diff topo to hazout */
					inhaz_topo[idx_fy][idx_cx]  = dem_nulval;                                        /*!< While out mask,  inmask is null */
					outmask_count++;
				}
			}
			/*!< if cell is out of mask raster but it was inside before */
			if  (mask_bed == -5555 && tzori_lim  > 0) 
			{
				/*!< Continuous out of mask until end of DEM or max distand are reached */
				z_cell  = raster[idx_fy][idx_cx]; 
				diff_topo = calc_diff(raster, z_cell, tzori_lim, dem_nulval, idx_fy, idx_cx, 2);
				outhaz_topo[idx_fy][idx_cx] = diff_topo;  
				inhaz_topo[idx_fy][idx_cx]  = dem_nulval; 
				outlimask_count++;
			}
			/*!< if cell is out of mask raster and new mask zlim is not available, 
			 * when no mask exist in some points */
			if  (mask_bed == -5555 && tzori_lim  == -9999)
			{
				z_cell  = raster[idx_fy][idx_cx];
				if(z_cell != dem_nulval && tzori_lim == -9999)                                      /*!< If z value is not null */
				{
					/*!< get init values of mask limits */
					tzori_lim = z_cell;
					txori_lim = p4x;
					tyori_lim = p4y;
				}
				/*!< Continuous out of mask until end of DEM or max distand are reached */
				diff_topo = calc_diff(raster, z_cell, tzori, dem_nulval, idx_fy, idx_cx, 2);
				outhaz_topo[idx_fy][idx_cx] = diff_topo;  
				inhaz_topo[idx_fy][idx_cx]  = dem_nulval;
			}
				//printf("loop %i count = %i in %lf  --  out %lf :: %lf %lf %lf\n", count_loops, count_out, inhaz_topo[idx_fy][idx_cx], outhaz_topo[idx_fy][idx_cx],  diff_topo, z_cell, tzori_lim);			
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (topmax_dist > 0 && count_out > 1) 
			{
				crossdist = calc_dist(p4x, txori_lim, p4y, tyori_lim);
				if (crossdist > topmax_dist) break;                         /*!< break if distance from mask lim > than max distance */
			}
		}
		else  in_demlimit = 0;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		count_loops++;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	}while(in_demlimit > 0);                                               /*!< While in raster  */
	printf("MASK Loops resolved with nloops %i and count direction %i\n", count_loops, direc_count);
	
	printf("Cells In mask %i\n", inmask_count);
	printf("Cells out mask %i\n", outmask_count);
	printf("Cells out raster mask %i\n", outlimask_count);
	
	
	//printf("In mask, cells fixed with 8 neighbour cells %i and no fixed %i\n", cont_infix, cont_inerr);
	//printf("Out mask, cells fixed with 8 neighbour cells %i and no fixed %i\n", cont_outfix, cont_outerr);
	printf("-----------------------------------------------------------\n");
}

/*! Differential of topography by nearest point
 * Calculation direction number 4
 * Not well implement yet, working in progress
 */
void nearest_topohazard(double **raster, struct HeadR s_arrayh[], float dem_nulval, struct Zeropoint s_zerolevpt[], int n_zlps, int nrio)
{
int i, j, k;
int brio; 
//idxfy, idxcx;
double x_zlp, y_zlp, z_zlp;
double x_cell, y_cell, z_cell;
double min_distance, dif_distance, endz_zlp;
	/**
	* Start topohazard calculation
	*/	
	for(i=0;i<s_arrayh[0].hn_fy;i++) //column
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) //row
		{
			/*!< if array index are in working area */
			if((i>1 && i<s_arrayh[0].hn_fy-1 ) && ( j>1 && j<s_arrayh[0].hn_cx-1 ))
			{
				
				z_cell  = raster[i][j];
				/*!< get coordinates using from index */
				x_cell = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, j);                       /*!< Calcula coordenada a partir de indice */
				y_cell = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, i);
				
				min_distance = 1000000000.0;
				dif_distance = 0;
				for(k=0;k<n_zlps;k++)
				{
					/*!< getting ZLP values */
					brio  = s_zerolevpt[k].brio;
					x_zlp = s_zerolevpt[k].bxcoor;                                   
					y_zlp = s_zerolevpt[k].bycoor;
					z_zlp = s_zerolevpt[k].bzcoor;
					//idxfy = s_zerolevpt[k].bidxfy;
					//idxcx =	s_zerolevpt[k].bidxcx;
					
					if (z_zlp != dem_nulval && brio  == nrio)
					{
						dif_distance = calc_dist(x_cell, x_zlp, y_cell, y_zlp);
						/*!< while distance is lower every time, continuous
						 * if distance increase, stop calculation */
						//if (k > 0 && dif_distance > min_distance) break;
						
						if (dif_distance < min_distance)
						{
							endz_zlp = z_zlp;
							min_distance = dif_distance;
						}
					}
				}
				
				/*!< add differential */
				raster_topohaz[i][j] = z_cell - endz_zlp;
				
				//printf("TEST: exit en k %i de n_zlps %i en %i %i totales %i %i con zdfi %lf de %lf\n", k, n_zlps, i, j, s_arrayh[0].hn_fy, s_arrayh[0].hn_cx, raster_topohaz[i][j], z_cell);
			}
		}
		if (i == (int)(s_arrayh[0].hn_fy * 5) / 100) printf("5%% columns evaluated\n");
		if (i == (int)(s_arrayh[0].hn_fy * 10) / 100) printf("10%% columns evaluated\n");
		if (i == (int)(s_arrayh[0].hn_fy * 25) / 100) printf("25%% columns evaluated\n");
		if (i == (int)(s_arrayh[0].hn_fy * 50) / 100) printf("50%% columns evaluated\n");
		if (i == (int)(s_arrayh[0].hn_fy * 75) / 100) printf("75%% columns evaluated\n");
		if (i == (int)(s_arrayh[0].hn_fy * 90) / 100) printf("90%% columns evaluated\n");
	}
}


#endif /* _TOPO_RAS */
