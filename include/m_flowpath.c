/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza 
* Module:  m_flowpath.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.11
* Creation Date: 2018-04-06
* Last update: 2025-12-03
* 
* Description:
* 
* Calc gravitational flows for morphology DEM analysis
* flow_mode 1: calculate single path flow using LHM 
* flow_mode 2: calculate single paht flow using SSM
* flow_mode 3: Drunk sailor multi path model
* flow_mode 4: Multipath, find all possible paths
* 
* Force interaction, in flow_mode 1 and 2
* for_inter 0: Find the most probable path
* for_inter 1: See if cell was selected in previous path and search for alternative path 
* 
* Flow_mode 1-2, In GIS open output csv, select a path by number, create a polyline path (point to path), 
* smooth path (simplify), get points (points along geometry), add attributes, export as csv
* 
* IMPORTANT: This algorithm triples RAM usage
* *********************************************************************/
#ifndef _FLO_RAS
#define _FLO_RAS

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
#include <time.h>

/*! Almacena los datos de la cabecera de arrays-raster grd */
struct Zeropoint s_zerolevpt[TOTZLP];

/*----------------------------------------------------------------------
 |  Section Global variables
 ---------------------------------------------------------------------*/
/*!< Flow calculation needs rasters to add new data and remember the result after finish
	 * the path from each init pt. We need ex tern variables to control the whole cycle
	 * We need almost three times memory to carry out this process, 
	 * so check the RAM availability */
extern double **raster_path; 
extern int **raster_control;
int FIFO[FIFOSIZE][2], ncell_infifo;

//**********************************************************************
/*! Auxiliar Functions */
//**********************************************************************

/*! DEPRESSION FILL FUNCTION FOR SINGLE PATH MODELS */
float *calc_levelfill(double nwval[8], int idx_fy, int idx_xc, struct HeadR s_arrayh[], int flow_mode)
{
int l, j, n, celda;	
double alt;
static float values[4];

	for(l=0;l<4;l++)values[l] = 0;
	
	n = j = 0;
	alt  = 0;
	for(l=0;l<8;l++)
	{
		/**
		* Check surface type
		*/
		if(nwval[l]==0)	j++;                                            /*!< if cell used or flats - count */
		if(nwval[l] <0)	n++; 	                                        /*!< if Surface depression - count */
		/**
		* Get new xy using LHM
		*/
		if(flow_mode == 1) 
		{
			if(nwval[l]>0)                                              /*!< if z diff is higher than 0 : cell is up */
			{
				values[0] += nwval[l];                                  /*!< Sum all z-diff higher than 0 */
				/**
				* Get max z-diff value
				*/
				if(nwval[l] > alt)                                        
				{
					alt = nwval[l];                                     /*!< Set max value */
					celda = l;                                          /*!< Set cell number */
				}	
			}
		}
		/**
		* Get new xy using SSM
		*/
		if(flow_mode == 2)
		{
			if(nwval[l]>0)                                              /*!< if z diff is higher than 0 */
			{
				values[0] += nwval[l];                                  /*!< Sum all z-diff higher than 0 */
				/**
				* Depending on neighboring cell position
				*/
				if(l <  4)nwval[l] = nwval[l] / (s_arrayh[0].hresx / 3);  /*!< Get slope-gradient */
				if(l >= 4)nwval[l] = nwval[l] / (s_arrayh[0].hresx / 2);
				/**
				* Get max s-gradient value
				*/
				if(nwval[l] > alt) 
				{
					alt = nwval[l];                                     /*!< Set max value */
					celda = l;                                          /*!< Set cell number */
				}	
			}
		}	
	}
	/**
	* According to surface type, set sum1 var
	*/
	if(n+j == 8) values[0] = -3;                                        /*!< Mix surface flat-sink */
	if(j == 8)   values[0] = -1;                                        /*!< Mix surface all_cell_used-flat */
	if(n == 8)   values[0] = -2;                                        /*!< Sink */
	if(values[0] == 0) values[0] = -4;                                  /*!< When no available cells to go, all have been used before */
	/**
	* Set new cell array index
	*/
	switch(celda) 
		{
		case 0:
			idx_fy--;
			break;
		case 1:
			idx_fy++;
			break;
		case 2:
			idx_xc++;
			break;
		case 3:
			idx_xc--;
			break;
		//--------	
		case 4:
			idx_fy++;
			idx_xc++;
			break;
		case 5:
			idx_fy--;
			idx_xc++;
			break;
		case 6:
			idx_fy--;
			idx_xc--;
			break;
		case 7:
			idx_fy++;
			idx_xc--;
			break;
		}
		/**
		* Set global variables
		*/
		values[1] = idx_fy;
		values[2] = idx_xc;
		/**
		* Calc distance between center cell and selected cell
		*/
		if(celda <  4) values[3] = s_arrayh[0].hresx;
		if(celda >= 4) values[3] = s_arrayh[0].hresx + (s_arrayh[0].hresx / 2);
	
	
	/**
	* This function returns
	* In value[0] sum of all z-diff higher than 0
	* In value[1] new file index
	* In value[2] new column index
	* In value[3] Distance to the selected cell  
	*/
	return (values);
}


/*! DEPRESSION FILL FUNCTION FOR RAMDOM PATH MODEL */
float calc_levelfill2(double *nwdiftopo, int *filein, int *coluin)
{
int i, j, l;
int f1, f2, c1, c2, buq;	
float sum1;
	
	buq = 0;
	for(i=0;i<8;i++)
	{
		/**
		*Check last selected cells, if they are the same, end filling depression
		*/
		f1 = filein[i];
		c1 = coluin[i];
		for(j=0;j<8;j++)
		{   
			f2 = filein[j];
			c2 = coluin[j];
			if(i != j)
			{
				if(f1 == f2 && c1 == c2) buq++;
			}
		}
	}
	if (buq < 4)         /*!< if repited cells are less than  */
	{
		sum1 = 0;
		for(l=0;l<8;l++)
		{
			if(nwdiftopo[l]>0)                                                 /*!< if z diff is higher than 0 */
			{
				sum1 += nwdiftopo[l];                                          /*!< Sum all z-diff higher than 0 */	
			}	
		}
	}
	if(sum1 < 0) sum1 = -2;         /*!< if still sink */
	if(buq  > 4) sum1 = -1;	        /*!< if selected cells higher than .. finish filling depression */
	
	return sum1;
}

/*! Return new altitude value after adding some critical height */
float increas_critheight(float cri_heig, float heig_add, double **raster, int idx_fy, int idx_xc)
{
float new_height, new_altitud;
	/**
	* Change Critical Height if necessary
	*/
	//printf("%i %i\n", idx_fy, idx_xc);
    new_height = cri_heig + heig_add;
	new_altitud = new_height+raster[idx_fy][idx_xc];  
	return new_altitud;   
}

int previous_cell(int **ocupcell, int idx_fy, int idx_xc, int pasos)
{
int c, r, cr;
int prefil, precol;
	prefil = ocupcell[pasos][0];                     /*!< previous selected cell's index row  */
	precol = ocupcell[pasos][1];                     /*!< previous selected cell's index col */
	c = 0;
	r = 0;
	if(prefil == idx_fy-1) c = 1;
	if(prefil == idx_fy+1) c = 2;
	if(prefil == idx_fy)   c = 3;
	if(precol == idx_xc-1) r = 10;
	if(precol == idx_xc+1) r = 20;
	if(precol == idx_xc)   r = 30;
	cr = c + r;
	return cr;
}

double *nw_difftopo(float dem_nulval, double *topoval, float h0_medval)
{
int l;
static double diftopo[8];
	
	for(l=0;l<8;l++)
	{
		if (topoval[l] != dem_nulval)                          /*!< if z is a right value */
		{
			diftopo[l] = h0_medval - topoval[l];
		}
	}
	return (diftopo);
}

/*! FIFO INIT FUNCTION */
void fifo_ini(void)
{
int i;
    /**
	* FIFO initialitation. 
	*/
    for(i=0;i<FIFOSIZE;i++) 
    {
        FIFO[i][0] = -1;                                                /*!< -1 means empty cell */
        FIFO[i][1] = -1;
    }
}
//----------------------------------------------------------------------
/*! FIFO SUBTRACT FUNCTION */
void fifo_move(void)
{
int i;
    i = 0;
    /**
	* Substract values from FIFO and recal its index. 
	*/
    do
    {
        FIFO[i][0] = FIFO[i+1][0];                                      /*!< Move current value to next position */
        FIFO[i][1] = FIFO[i+1][1];
        i++;
    }while(FIFO[i][0] != -1);                                           /*!< While current value is not -1 */
    FIFO[FIFOSIZE-1][0] = -1;
    FIFO[FIFOSIZE-1][1] = -1;
    ncell_infifo = i-1;                                                 /*!< 1 must be subtracted to avoid the insertion of -1 values in the middle */
}
//----------------------------------------------------------------------
/*! FIFO LOAD FUNCTION */
void fifo_load(int indyf, int indxc)
{
    /**
	* Add values to FIFO 
	*/
    if(ncell_infifo < FIFOSIZE)
    {
        FIFO[ncell_infifo][0] = indyf;
        FIFO[ncell_infifo][1] = indxc;
        ncell_infifo ++;
    }
}
//----------------------------------------------------------------------




//**********************************************************************
/*! Calc Functions */
//**********************************************************************

/*! Single flow path, calculate the most direct path
 * flow_mode 1 use LHM method to calculate new path
 * flow mode 2 use SSM method to calculate new path
 * for_inter 0 paths are calculated as individual
 * for_inter 1 previous path are considered to calculate the new one */
void calc_singflow(int nflow, const char *dir_out, const char *namfile, const char *namxyven, double **raster, struct HeadR s_arrayh[], \
	float dem_nulval, int indyf, int indxc, int ptid, double txori, double tyori, double tzori, int flow_mode, float max_dist, \
	float cri_heig, float heig_add, int num_repit, int for_inter, int dist_type)
{	
int i, j, k, n_gflowpath, n, o, maxbdpt, ext, out;
unsigned char done;
double ll, tx, ty, tz;

char *subfix;
char *file_outpath;
int ras_cx, ras_fy;
int idx_fy, idx_xc, status, nwsearch;
float *values, distpt, hl2;
float newheight;
double *topoval, nwval[8], areaflow;
double dist_2d, diff_he, dist_3d;
int *flowval;

	done = 0;                                                           /*!< Control the loop */
	o    = 0;                                                           /*!< Count sinks crossed */
	n_gflowpath = 0;                                                    /*!< Count selected cells by gravitational flow path */
	nwsearch   = 0;                                                         
	hl2  = heig_add;
	maxbdpt = (TOTZLP * 80) / 100;                                    /*!< Decrease the max value of points available for each flow path to avoid memory overflow */
	idx_fy= indyf;                                                      /*!< Init column index */
	idx_xc= indxc;                                                      /*!< Init row index */
	raster_control[indyf][indxc] = 1;                                   /*!< Starting point value : only when force intercation 1*/
	raster_path[indyf][indxc]    = 1;                                   /*!< Starting point value */
	ext=0;                                                              /*!< Control exit mode type */
	out=0;
	//-----------------
	
	ras_cx = s_arrayh[0].hn_cx;
	ras_fy = s_arrayh[0].hn_fy;
	printf("Iniciate loop for flow path %i %i\n", indyf, indxc);
	do
	{
		n = 0;                                                          /*!< Count how many times the z is increased when sink */
		for(k=0;k<8;k++)nwval[k]=0;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//printf("height, %lf\n", raster[idx_fy][idx_xc]);
		newheight = increas_critheight(cri_heig, hl2, raster, idx_fy, idx_xc);  /*!< Calc new altitud in the new position increasing its value */
		topoval   = search_celproxF(raster, idx_fy, idx_xc);               /*!< Get 8 cels value in DEM */  
		flowval   = search_celproxI(raster_control, idx_fy, idx_xc);       /*!< Get 8 cels value in Flow-control rater */
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		for (k=0;k<8;k++)
		{
			if (flowval[k] == 0 && topoval[k] != dem_nulval) nwval[k] = newheight - topoval[k];  /*!< Get z-diff if celval is > 0 */
			else  nwval[k] = 0;
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		values = calc_levelfill(nwval, idx_fy, idx_xc, s_arrayh, flow_mode);  /*!< Get a new cell in the flow path */
		status = values[0];
		idx_fy = (int)values[1];
		idx_xc = (int)values[2];
		distpt = values[3];
		nwsearch++;
		//printf("%i : %i %i : %f\n", status, idx_fy, idx_xc, distpt);
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if(idx_fy<5||idx_fy>ras_fy-5||idx_xc<5||idx_xc>ras_cx-5||raster[idx_fy][idx_xc]== 0||raster[idx_fy][idx_xc]== dem_nulval)done = 0;    /*!< if the new cell index is out of limits - stop loop */
		else
		{
			if(status == -1)                                            /*!< All 8 cells used or they have same altitud */	
			{
				if(for_inter == 0)                                      /*!< if no force interaction */
				{
					status = -2;
					if(out > num_repit)                                 /*!< out higher than total number iteration allowed */
					{
						ext = 6;
						out = 0;
						break; 
					}
					out++;
				}
				if(for_inter == 1)	                                    /*!< if force interaction */
				{
					ext = 1;
					break; 
				}	                                       				/*!< if all neighboring cell were used before - stop */	
			}	
			if(status == -2 || status == -3)                            /*!< if the cell is a sink */
			{

				if(heig_add > 0)										/*!< if fill increase is higher than 0 */
				{
					hl2 += heig_add;                                    /*!< increase the critical height value using Filling Increase var */
					n++;
				}
				else break;		
			}
			if(status == -4) done = 0;                                  /*!< exit, no available cells to go */
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			if(status > 0)                                                   /*!< if cell is not a sink */
			{
				tx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_xc);  /*!< getting the coordinate from DEM index */
				ty = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
				tz = raster[idx_fy][idx_xc];                                 /*!< getting original z value */                            
				//------------------------
				if (tz == dem_nulval) 
				{
					ext = 3;
					break;								/*!< Avoid z null values */
				}	
				//------------------------
				/**
				* Saving data in caucerio struc 
				*/
				s_zerolevpt[n_gflowpath].bjerar  = 1;
				s_zerolevpt[n_gflowpath].brio    = nflow;
				s_zerolevpt[n_gflowpath].btramo  = ptid;
				s_zerolevpt[n_gflowpath].bidpt   = n_gflowpath; 
				if(n_gflowpath == 0)s_zerolevpt[n_gflowpath].bdist = 0;
				if(n_gflowpath >  0)s_zerolevpt[n_gflowpath].bdist = s_zerolevpt[n_gflowpath-1].bdist + distpt;
				s_zerolevpt[n_gflowpath].bxcoor  = tx;
				s_zerolevpt[n_gflowpath].bycoor  = ty;
				s_zerolevpt[n_gflowpath].bzcoor  = tz;
				s_zerolevpt[n_gflowpath].bdx     = distpt;
				s_zerolevpt[n_gflowpath].bdt     = 0;
				s_zerolevpt[n_gflowpath].bquality = n;
				raster_path[idx_fy][idx_xc]      +=1;                                /*!< total cells selected */
				raster_control[idx_fy][idx_xc]   = 1;                                /*!< Control raster - this cell was selected */
				//ll = sqrt((tx-txori)*(tx-txori)+(ty-tyori)*(ty-tyori));              /*!< Calc distance from cell to init cell */
				
				dist_2d = calc_dist(tx, txori, ty, tyori);                /*!< Calc 2d distance to init */ 
				diff_he = fabs(tzori - tz);                                /*!< Calc height difference to init */
				dist_3d = calc_hipo(dist_2d, diff_he);                    /*!< Calc 3d distance to init */ 
				
				if (dist_type == 1) ll = dist_2d;
				if (dist_type == 2) ll = dist_3d;
				if (dist_type == 3) ll = n_gflowpath * s_arrayh[0].hresx;  /*!< Calc cell distance based on number of cells covered */
				
				//ll = s_zerolevpt[n_gflowpath].bdist;
				//--------
				if(ll > max_dist)
				{
					ext  = 4;
					done = 0;                                           /*!< If distance is higher than Max. Dist - stop loop */
					printf("TEST: ending by maximum distance %lfi of %lf\n", ll, max_dist);
				}
				else done = 1;  
				if(n>0)o++;                                             /*!< if sink, add count */
				n_gflowpath++;                                          /*!< Count a new cell in grav. flow path */
				n = 0; 													/*!< Reset Fill increase count */
				hl2 = 0;                                                /*!< Reset Fill increase value */
			}
			if(n_gflowpath == maxbdpt)                                  /*!< if total cells allowed per grav. flow path are reached */
			{
				ext  = 5;
				done = 0; 
				printf("TEST: ending by maximum cells reached %i of %i\n", n_gflowpath, maxbdpt);
				break;                                                     
			}
		}
			
		//if (nflow == 5090)printf("TEST:status %i out %i ext %i max %i, dist %lf done %i hei %lf::%lf search %i\n", status, out, ext, m, ll, done, tz, hl2, nwsearch);
	}while(done>0);                                                     /*!< Keep the loop on while done > 0 */
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	
	printf("nfile %i exit type %i hl2 %f\n", nflow, ext, hl2);
	//resflow = m;
	if(n_gflowpath>0)                                                   /*!< If the flow path exist */
	{
		//write the first flow path
		subfix = (char *)malloc(2+2+strlen(namxyven)+20); 
		sprintf(subfix,"_pathflow_mod-%i-%i_%s.csv", flow_mode, for_inter, namxyven);
		if (nflow == 0)
		{
			file_outpath = get_pathnam(dir_out, namfile, subfix, 1);
			write_flowpath(file_outpath, "w", nflow, s_zerolevpt, n_gflowpath);
		}
		//add the next flow paths
		if (nflow > 0)
		{
			file_outpath = get_pathnam(dir_out, namfile, subfix, 1);
			write_flowpath(file_outpath, "a", nflow, s_zerolevpt, n_gflowpath);
		}
		areaflow = n_gflowpath * s_arrayh[0].hresx * s_arrayh[0].hresy;
		printf("Exit by %i\n", ext);
		printf("Maximum distance reached %.2lf\n", ll);
		printf("Total points/cells per flow path calculated %i\n", n_gflowpath);
		printf("Flow area %lf\n", areaflow);
		printf("Total sinks per flow path found %i\n", o);
		printf("Total searchs %i\n", nwsearch);
		printf("----------------------------------------------\n");
		printf("----------------------------------------------\n");
		printf("----------------------------------------------\n");
		free(subfix);
	}
	/*!< If force interacction is equal to 0 means, no interaction will be taken
	 * each flow-paht will follow the main route without taken into accout other simulated paths
	 * To do that, we must to reset the raster_control */
	if(for_inter == 0)
	{
		for(i=0;i<s_arrayh[0].hn_fy;i++) 
		{
			for(j=0;j<s_arrayh[0].hn_cx;j++)raster_control[i][j] = 0;
		}
		
	}
}

/*! Drunk salier montecalor type flow path
 * Based on Felpeto et al., 2001 (https://dx.doi.org/10.1023/A%3A1011112330766)
 * and Marrero et al., 2013 (https://dx.doi.org/10.1007/s11069-013-0672-4)
 * flow_mode 3
 * for_inter 0 Control path is not reset after finishing all iteration of each pt vent
 * for_inter 1 Control path is reset after finishing all iteration of each pt vent*/
void calc_drunksailflow(double **raster, struct HeadR s_arrayh[], \
	float dem_nulval, int indyf, int indxc, double txori, double tyori, double tzori, float max_dist, int num_itera, \
	int num_repit, float radiu_val, float cri_heig, float heig_add, int for_inter, int dist_type)
{
int i, l, m, q, o;
int cr, cn;
int limext; 
int n_gflowpath, limpass, rep;
double s[8];
double pr;
unsigned char done;
double ll;

int okval, n_attem, nitera_done, ocu_rows;
int filein[8], coluin[8];
int step_reach, maxd_reach, nowayout, badcell;
int idx_fy, idx_xc, ras_cx, ras_fy, **ocupcell;
float h0_medval, search_fradius, search_xradius;
float min_indyf, max_indyf, min_indxc, max_indxc;
double sum, sum1;
double tx, ty, h0, *topoval, diftopo[8], sumtopo, *nwdiftopo;
double dist_2d, diff_he, dist_3d;
	/**
	* Init var
	*/
	//iniepoch = clock();                                                 /*!< Get init time */
	//if(nfile == 0) firsepoch = clock();
    for(i=0;i<8;i++)s[i]= 0; 
    cn         = 0;                          
    o          = 0;                                                     /*!< Count total number of sinks in flow path */
    m          = 0;                                                     /*!< Count total number of cells in flow path */
    n_attem    = 1;                                                     /*!< Control the numbers of flow calculation attempts in complicated cells */
    rep        = 0;                                                     /*!< Control noway-out */
    nitera_done= 0;                                                     /*!< Count total number of iterations */
    nowayout   = 0;                                                     /*!< Count cells selected high numbers of times consecutively */
    badcell    = 0;														/*!< Count cells with null value or 0 in altitude */
    step_reach = 0;														/*!< Count number of times where the total number of steps by flow simulation is reached */
    maxd_reach = 0;														/*!< Count number of times where the maximum distance is reached */
	limext     = 0;														/*!< Count number of times where the raster's limits are reached */
	limpass    = (max_dist / s_arrayh[0].hresx) * num_itera;            /*!< Total cells permitted per flow path calculation */
	ocu_rows   = limpass;
	ocupcell   = Crea_2DIarray(limpass, 2);                             /*!< Array to save indexs of selected cells */
	ras_cx = s_arrayh[0].hn_cx;
	ras_fy = s_arrayh[0].hn_fy;
	
	search_fradius = radiu_val /  s_arrayh[0].hresy;                    /*!< Convert search radius in row cells */
	search_xradius = radiu_val /  s_arrayh[0].hresx;                    /*!< Convert search radius in column cells */
	min_indyf = indyf - search_fradius;                                 /*!< Minumun decrease row cells */
	max_indyf = indyf + search_fradius;                                 /*!< Maximun increase row cells */
	min_indxc = indxc - search_xradius;                                 /*!< Minumun decrease column cells */
	max_indxc = indxc + search_xradius;                                 /*!< Maximun increase column cells */
    /**
	* First Loop - Control iterations and init zpoint
	*/
	do
	{
		n_gflowpath = 0;                                                 /*!< Count cells selected per path */
		rep     = 0;
		done    = 0;                                                    /*!< Control the loop - reset to 0 */
		/**
		* Change original xy index using a random function to find a new vent in each iteration
		*/
		idx_fy = (int)random_float(min_indyf, max_indyf);               /*!< New index row based on a minimun and maximum distance */
		idx_xc = (int)random_float(min_indxc, max_indxc);               /*!< New index column based on a minimun and maximum distance */
		//k = (int)get_random(999.0);                                     // index position of numer of cells in which the original index will be moved
		//idx_fy = indyf + gauss_belt[k];                                 //moving original index
		//k = (int)get_random(999.0);
		//idx_xc = indxc + gauss_belt[k];
		/**
		* Check if indexs are in DEM
		*/
		//printf("TEST: fy %i xc %i :: %lf %lf\n", idx_fy, idx_xc, min_indyf, max_indyf);
		if(idx_fy<5||idx_fy>ras_fy-5||idx_xc<5||idx_xc>ras_cx-5)done=0;
		else
		{
			/**
			* Get coordinates of original xy
			*/ 
			tx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_xc);  /*!< getting the coordinate from DEM index */
			ty = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
			/**
			* Second Loop - Get the flow path in each iteration
			*/
			do
			{
				/**
				* get cell value
				*/
				h0 = raster[idx_fy][idx_xc];
				if(h0 == 0 || h0 == dem_nulval)
				{
					badcell++;
					done = 1;  /*!< if z == 0 if sink s[7]=0 or s[7]=hl and sum=0 then done=1 out loop 2 */
				}
				else
				{
					/*!< check if we are in the raster */
					if(idx_fy<5||idx_fy>ras_fy-5||idx_xc<5||idx_xc>ras_cx-5)done=1;
					else
					{
						/**
						* find 8 cells
						*/
						topoval = search_celproxF(raster, idx_fy, idx_xc);       /*!< get 3x3 moving cell z values */
						sum     = 0;                                            /*!< reset sum var */
						sumtopo = h0;
						okval   = 1;                                   
						for(l=0;l<8;l++)
						{
							if (topoval[l] != dem_nulval)                          /*!< if z is a right value */
							{
								sumtopo   += topoval[l];                            /*!< sum all z values, including center cell */
								diftopo[l] = h0 - topoval[l];                       /*!< get height difference */
								if(diftopo[l]>0)sum += diftopo[l];                  /*!< Sum all z-diff higher than 0 */
								s[l]       =sum;                                    /*!< store sum var in s[] array for each element */ 
								okval++;                                            /*!< Num the z values ok */
							}  
							else 
							{
								diftopo[l]=0; 
								s[l] = 0; 
							}                                                                              
						}
						/**
						* Calc z-diff media value
						*/
						h0_medval = sumtopo / okval;                               /*!< Medium z value from all valid cells */
						/**
						* Check loop conditions - sink - z value - sum
						*/
						if(sum==0)                                                 /*!< central cell is a sink */
						{
							o++;
							if(heig_add > 0)									   /*!< if fill increase is higher than 0 */
							{
								if(n_gflowpath > 1)                                 /*!< if sink is located after two cells in the new path */
								{
									cr = previous_cell(ocupcell, idx_fy, idx_xc, n_gflowpath-2);    /*!< calc direcction of previous cell */
									if(cr > 0 && cr != 33)                                         /*!< if coming from antoher cell and it is not the same one */
									{
										nwdiftopo = nw_difftopo(dem_nulval, topoval, h0_medval);   /*!< Calc z-diff using the z medium value */
										raster[idx_fy][idx_xc] = h0_medval;                        /*!< Change the z value with the new medium value */
										
										sum1 = calc_levelfill2(nwdiftopo, filein, coluin);
										if(sum1 > 0)
										{
											sum=0;
											for(l=0;l<8;l++)
											{
												if(nwdiftopo[l]>0)
												{
													sum += nwdiftopo[l];                      /*!< Sum all z-diff higher than 0 */
													s[l]=sum;                                 /*!< store sum var in s[] array for each element */
												}
												else s[l]=0;
											}
										}
										if(sum1 == -1)
										{
											done = 2;
											break;
										}
									}
									else done = 2;				
								}
								else done = 2;
							}	
							else done = 2;                                      /*!< if sink can not be filled */
						}	
					}	
				}
				if(done==0)
				{
					/**
					* select the new cell from the 8 neighbours cells
					*/
					pr = s[7]*((double)(rand()%1000)/1001.0);               /*!< val from 0 to sum max */
					q = 7;                                                  /*!< starting value for q var */
					if(pr<s[6])q=6;                                         /*!< if pr random var is less than s[x] */
					if(pr<s[5])q=5;
					if(pr<s[4])q=4;
					if(pr<s[3])q=3;
					if(pr<s[2])q=2;
					if(pr<s[1])q=1;
					if(pr<s[0])q=0;
					/**
					* Set new cell array index
					*/
					switch(q) 
					{
					case 0:
						idx_fy--;
						break;
					case 1:
						idx_fy++;
						break;
					case 2:
						idx_xc++;
						break;
					case 3:
						idx_xc--;
						break;
					case 4:
						idx_fy++;
						idx_xc++;
						break;
					case 5:
						idx_fy--;
						idx_xc++;
						break;
					case 6:
						idx_fy--;
						idx_xc--;
						break;
					case 7:
						idx_fy++;
						idx_xc--;
						break;
					}
					/**
					* Check loop conditions to see if new cell is inside array
					* borders are not considered
					*/
					if((idx_fy<2 || idx_fy>ras_fy-3) || (idx_xc<2 || idx_xc>ras_cx-3))
					{
						limext++;
						done = 3;     /*!< if outside then done= 3 out loop 2 */
					}	
					else
					{
						/*!< Count new cells in m
						if(for_inter == 0 && heig_add  > 0) 
						{
							if(raster_control[idx_fy][idx_xc]==0)m++;
						}	
						else
						{
							if (nitera_done == 0)
							{
								if(raster_control[idx_fy][idx_xc]==0)m++;
							}	
							if (nitera_done > 0)
							{
								if(raster_path[idx_fy][idx_xc]==0)m++;
							}
						}
						*/									
						raster_path[idx_fy][idx_xc]   += 1;							/*!< sum 1 to the cell */
						raster_control[idx_fy][idx_xc] = 1;                                    /*!< Control raster */
						//rast3[i][j] ++;                                     /*!< sum 1 to the cell */
						//if(maximun < raster_path[i][j])maximun = raster_path[i][j];     /*!< Get the maximum value */
						ocupcell[n_gflowpath][0] = idx_fy;                                /*!< Save array index from selected cells */
						ocupcell[n_gflowpath][1] = idx_xc;
						filein[cn] = idx_fy;
						coluin[cn] = idx_xc;
						if(cn == 8) cn = 0;                                 /*!< Save array index of last eight cells -control filling depression */
						if(cn  < 8) cn++;
						/**
						* Get coordinates of original xy
						*/ 
						tx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, idx_xc);     /*!< getting the coordinate from DEM index */
						ty = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, idx_fy);
						
						
						dist_2d = calc_dist(tx, txori, ty, tyori);                /*!< Calc 2d distance to init */ 
						diff_he = fabs(tzori - tzori);                                /*!< Calc height difference to init */
						dist_3d = calc_hipo(dist_2d, diff_he);                    /*!< Calc 3d distance to init */ 
						
						if (dist_type == 1) ll = dist_2d;
						if (dist_type == 2) ll = dist_3d;
						if (dist_type == 3) ll = n_gflowpath * s_arrayh[0].hresx;  /*!< Calc cell distance based on number of cells covered */
						

						ll = sqrt((tx-txori)*(tx-txori)+(ty-tyori)*(ty-tyori));         /*!< Calc distance from init */
						if(ll>max_dist)                                     /*!< Check if distance is higher than permitted */
						{
							maxd_reach++;
							done = 5;                                       /*!< if yes then done=5 - exit second loop */
						}	
						if(n_gflowpath > limpass)	                            /*!< Check if number of cells is higher than permitted */
						{
							step_reach++;
							done = 6;                                       /*!< if yes then done=6 - exit second loop */
						}
						if(n_gflowpath > 2) 
						{
							if(ocupcell[n_gflowpath-2][0] == idx_fy && ocupcell[n_gflowpath-2][1] == idx_fy) rep++;   /*!< Check if cell has been selected recently */
							
							if(rep > num_repit)                                  /*!< if cell has been selected recently more than N times */
							{ 
								nowayout++;
								done = 4;                                   /*!< if yes then done=4 - exit second loop */
							}	
						}	
						n_gflowpath++;	                                        /*!< Count selected cells */
					}	
				}
			}
			while(!done);                                                   /*!< While done = 0 keep loop 2 */
			/**
			* check if remain in the first Loop
			*/
			if(done == 3 || done == 4 || done == 5) done=0;
			if(done == 6)
			{
				limpass = limpass / 2;                                      /*!< Reduce the number of allowed steps */
				done=0;	
			}	
			if(done == 1 || done == 2)                                      /*!< if done=1-2 */ 
			{
				if(n_attem < 100)done=0;
				n_attem++;   
			}	
			nitera_done++;                                              /*!< Add new iteration */ 
			//printf("TEST: nitera %i done %i: %i %i\n",  nitera_done, done, idx_fy, idx_xc);
		} 
		/**
		 * Done = 0 Continuous the calculation
		 * Done = 1 if z value is null or 0
		 * Done = 2 if sink is not solved
		 * Done = 3 if out of raster limits
		 * Done = 4 If num_repit is exceeded
		 * Done = 5 if max_dist is exceeded
		 * Done = 6 if num pasos is exceeded
		*/                                                       
	}while(!done && nitera_done < num_itera);                                    /*!< While done = 0 and cont1 lower than total iterations then keep loop 1 */

	/**
	* Write and Reset raster
	*/                                         
    printf("Total cells selected per flow path %i\n", m);
	printf("Total sinks found %i\n", o);
	printf("Total interactions %i\n", num_itera);
	printf("Total number of times that max distance was reached %i\n", maxd_reach);
	printf("Total number of times that max step was reached %i\n", step_reach);
	printf("Total number of times without exit %i\n", nowayout);
	printf("Total null cells %i\n", badcell);
    //-----
    free_RasterI("Array ocupcell", ocupcell, ocu_rows);
    /*!< If force interacction is equal to 0 means, no interaction will be taken
	 * each flow-paht will follow the main route without taken into accout other simulated paths
	 * To do that, we must to reset the raster_control */
}



/*! MULTITRAYECTORY FLOW PATH (4) 
 * Based on Felpeto et al., 2001 (https://dx.doi.org/10.1023/A%3A1011112330766)
 * and Marrero et al., 2013 (https://dx.doi.org/10.1007/s11069-013-0672-4)
 * */
void calc_multyflow(double **raster, struct HeadR s_arrayh[], \
	float dem_nulval, int indyf, int indxc, double txori, double tyori, double tzori, float max_dist, float rest_heig, int dist_type)
{
int getval;	
double ll;	

int l, idx_fy, idx_xc, okdif, nwfy, nwxc;
int ras_fy, ras_cx, n_gflowpath, cont_sinks, cont_loops;
float sumdif, por_diftopo;
double *topoval, diftopo[8], h0;
double dist_2d, diff_he, dist_3d;
double tx, ty;
	
	//iniepoch = clock();                                                 /*!< Get init time */
	//if(nfile == 0) firsepoch = clock();
	n_gflowpath   = 0;                                                  /*!< Count total number of selected cell per flow path */
	cont_loops   = 0;
	ncell_infifo = 0;
	cont_sinks   = 0;													/*!< Count sinks */ 
	/**
	* Get init array index and z value
	*/ 
	idx_fy= indyf;                                                      /*!< Init column index */
	idx_xc= indxc;                                                      /*!< Init row index */    
	ras_cx = s_arrayh[0].hn_cx;
	ras_fy = s_arrayh[0].hn_fy;                                                    
	fifo_ini();                                                          /*!< Init FIFO */
	/**
	* Set init z point to 1
	*/                                                                  
	raster_control[indyf][indxc] = 1;                                   /*!< Store value */   
	raster_path[indyf][indxc]    = 1;                                   /*!< Store value */
	fifo_load(indyf, indxc);                                            /*!< Sent first data to the FIFO */
	n_gflowpath++;
	/**
	* Start Loop
	*/
	do
	{
		if (cont_loops > 0)
		{
			idx_fy = FIFO[0][0];                                        /*!< Get from FIFO col value */
			idx_xc = FIFO[0][1];                                        /*!< Get from FIFO row value */
			if(idx_fy == -1 && idx_xc == -1)break;                      /*!< if FIFO is empty */
		}
		cont_loops++;
		h0     = raster[idx_fy][idx_xc];                                /*!< Get z value from DEM raster */
		getval = raster_path[idx_fy][idx_xc];                           /*!< Get rast path value */
		if (h0 != dem_nulval)						                    /*!< If z val is not null */
		{
			/**
			* Check if cell is inside array
			*/
			if((idx_fy < 2 || idx_fy >= ras_fy-2) || (idx_xc < 2 || idx_xc >= ras_cx-2))break; 		/*!< if not then stop loop*/
			else
			{
				okdif   = 0;
				sumdif  = 0.0;
				topoval = search_celproxF(raster, idx_fy, idx_xc);
				for(l=0;l<8;l++)
				{
					if (topoval[l] != dem_nulval)                       /*!< if z is a right value */
					{
						diftopo[l] = h0 - topoval[l];                   /*!< get height difference */
						if (diftopo[l] > 0) 
						{
							sumdif += diftopo[l];
							okdif++;
						}
					}
					else
					{
						diftopo[l] = 0;
					}
				}
				if (okdif > 0)
				{
					for(l=0;l<8;l++)
					{
						if (diftopo[l] > 0) 
						{
							if (l == 0)                                 /*!< if z-diff in cero is higher than 0 */
							{
								nwfy = idx_fy--;                        /*!< Change column index */
								nwxc = idx_xc;
							}
							if (l == 1)
							{
								nwfy = idx_fy++;   
								nwxc = idx_xc;                                                              
							}
							if (l == 2)
							{
								nwfy = idx_fy;
								nwxc = idx_xc++;                                
							}
							if (l == 3)
							{                                  
								nwfy = idx_fy;
								nwxc = idx_xc--;                                
							}
							if (l == 4)
							{
								nwfy = idx_fy++;                                
								nwxc = idx_xc++;                                
							}
							if (l == 5)
							{
								nwfy = idx_fy--;                                
								nwxc = idx_xc++;                                
							}
							if (l == 6)
							{
								nwfy = idx_fy--;                                
								nwxc = idx_xc--;                                
							}
							if (l == 7)
							{
								nwfy = idx_fy++;                                
								nwxc = idx_xc--;                                
							}
							por_diftopo = (diftopo[l] * 100) / sumdif;	           /*!< Calc percentage of total sum */
							if( raster_control[nwfy][nwxc] == 0 )                  /*!< if it is a new cell in the flow path */    
							{
								if(por_diftopo > rest_heig)                        /*!< if percentage is higher than Restric. Multiflow var then ... */       
								{  
									raster_control[nwfy][nwxc] = 1;                /*!< Set new value in rast1 */ 
									raster_path[nwfy][nwxc]    = getval+1;         /*!< Set new value in rast3 */ 
									fifo_load(nwfy,nwxc);                          /*!< Load array index in FIFO */ 
									n_gflowpath++;	                               /*!< Count selected cell */  
									
									//printf("loop %i ncell %i l %i :: %i %i : %i %i :: r %lf :: p %f\n", cont_loops, n_gflowpath, l, idx_fy, idx_xc, nwfy, nwxc, raster_path[nwfy][nwxc], por_diftopo);                                                
								}	
							}
							/**
							* Check distance to init cell
							*/
							tx = get_coor(s_arrayh[0].hxlo, s_arrayh[0].hresx, nwxc);  /*!< getting the coordinate from DEM index */
							ty = get_coor(s_arrayh[0].hylo, s_arrayh[0].hresy, nwfy);
							
							dist_2d = calc_dist(tx, txori, ty, tyori);                /*!< Calc 2d distance to init */ 
							diff_he = fabs(tzori - h0);                                /*!< Calc height difference to init */
							dist_3d = calc_hipo(dist_2d, diff_he);                    /*!< Calc 3d distance to init */ 
							
							if (dist_type == 1) ll = dist_2d;
							if (dist_type == 2) ll = dist_3d;
							if (dist_type == 3) ll = n_gflowpath * s_arrayh[0].hresx;  /*!< Calc cell distance based on number of cells covered */

							if(ll>max_dist)break;                                     /*!< Check if the distance is higher than Max. Dist. */	
						}
					}	
					fifo_move();	                                    /*!< subtract cell from FIFO */	
				}
				if (okdif == 0)
				{			
					cont_sinks++;
					fifo_move();                                        /*!< if not, subtract cell from FIFO */                    		
				}
			}
		}
		else	fifo_move();	   	                                    /*!< subtract cell from FIFO */
	}while(ncell_infifo > 1 || idx_fy != -1 || idx_xc != -1);		    /*!< While FIFO has cells */		
			
	printf("Maximun distance reached %.2lf m\n", ll);
	printf("Total points/cells selected %i\n", n_gflowpath);
	printf("Total sinks found %i\n", cont_sinks);
	printf("--------------------------------\n\n");


}



#endif /* _FLO_RAS */

/*!
 * UPDATEs
 * 2025/12/03: changed abs by fabs function to avoid compilation issues
 */
