/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza
* Module:  m_ultiflow.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2018-04-06
* Last update: 2022-09-30
* 
* Description:
* Calculate cross section of a river valley morphology or Schilling equation implementation 
 * multi_mode 1: cross section 
 * multi_mode 2: Schilling equations 
* 
* *********************************************************************/
#ifndef _MUL_RAS
#define _MUL_RAS

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
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------
 |  Section Global variables
 ---------------------------------------------------------------------*/

struct  CrossSecc s_crosec[NUMCROSS];
extern double **rast_multi;


/*! Calc crosection in one direction
 * This function must be call twice */
int calc_cross(const char *file_cross, double **raster, struct HeadR s_arrayh[], float dem_nulval, \
		double vnewx, double vnewy, int direc_count, float cros_maxwith, float cros_height, \
		int njerar, int nrio, int ntramo, float ndist, double txori, double tyori, double tzori, int idxfy, int idxcx, float tdx)
{
FILE *file;
int n, nwpt;
int idx_cx, idx_fy, tot_idxcx, tot_idxfy, in_demlimit;
double p4x, p4y, p4z, diff_topo;
float raster_discre, crossdist;
	n = nwpt      = 0;
	raster_discre = s_arrayh[0].hresx / 2;
	in_demlimit   = 1;
	tot_idxcx     = s_arrayh[0].hn_cx - 2;
	tot_idxfy     = s_arrayh[0].hn_fy - 2;
	crossdist     = tdx;
	diff_topo     = 0;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if((file = fopen(file_cross, "a"))== NULL)
	{
		printf("-------ERROR open file--------\n");
		exit(0);
	}
	else
	{
		fprintf(file,"%i %i %i %.2f %lf %lf %lf %i %i %.2f %.2f\n",
				
				njerar,
				nrio,
				ntramo,
				ndist,
				txori,
				tyori,
				tzori,
				idxfy,
				idxcx,
				crossdist,
				diff_topo			
				);	
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		do
		{
			if(direc_count == 1)n++;                                        /*!< if positive count */
			if(direc_count == 2)n--;
			p4x  = txori + (vnewx * raster_discre * n);  /*!< Get oblique coordinates */
			p4y  = tyori + (vnewy * raster_discre * n);
			idx_cx = calc_Rindex(p4x, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
			idx_fy = calc_Rindex(p4y, s_arrayh[0].hylo, s_arrayh[0].hresy);
			crossdist = calc_dist(p4x, txori, p4y, tyori);
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			/*!< Check if index if they are in the raster */
			if(idx_fy > 2 && idx_fy < tot_idxfy && idx_cx > 2 && idx_cx < tot_idxcx) 
			{
				in_demlimit = 1;    
				//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				p4z  = raster[idx_fy][idx_cx]; //get z value
				if (p4z != dem_nulval)diff_topo = (float)p4z - tzori;
				else                  diff_topo = 0.0;
				/*!< write data in output file */
				fprintf(file,"%i %i %i %.2f %lf %lf %lf %i %i %.2f %.2f\n",
					njerar,
					nrio,
					ntramo,
					ndist,
					p4x,
					p4y,
					p4z,
					idx_fy,
					idx_cx,
					crossdist,
					diff_topo			
					);	
				nwpt++;
				//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				if (crossdist > cros_maxwith) break;                        /*!< break if cross distance > max distance */
				//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				if (cros_height >0)                                     /*!< if crossection height is activated */
				{
					diff_topo = p4z - (tzori + cros_height);            /*!< check if the new cell es below or above to tzori+added heigh */
					if (diff_topo > 0) break;
				}
				//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			}
			else in_demlimit = 0; 

		}while(in_demlimit > 0); //calcula mientras estemos dentro de la distancia maxima
	}
	fclose(file);
	return nwpt;
}


/*! Get crossection from specific zero level trajectory */
void calc_crossfile(const char *file_cross, double **raster, struct HeadR s_arrayh[], float dem_nulval, \
	struct Zeropoint s_zerolevpt[], int n_bedxyz, int idx_pt, float cros_maxwith, float cros_height, \
	int njerar, int nrio, int ntramo, float ndist, double txori, double tyori, double tzori, int idxfy, int idxcx, float tdx, int tend, int triodir)
{
int nextpt, prevpt,negpt,postpt;
int increas;
double pretx,prety,nextx,nexty;
double *xy_fromvec, vnewx, vnewy;

	negpt=postpt=0;
	//points located in the middel----------------------------------
	if(tend == 0)
	{
		//0 top/botton area
		if (s_arrayh[0].hresx <= 5)                            increas = 2;
		if (s_arrayh[0].hresx >  5 && s_arrayh[0].hresx <= 10) increas = 3;
		if (s_arrayh[0].hresx > 10 && s_arrayh[0].hresx <= 20) increas = 4;
		if (s_arrayh[0].hresx > 20)                            increas = 5;
		if (idx_pt < increas || idx_pt > (n_bedxyz-increas))   increas = 1;
		
		if(triodir == 1) 
		{
			nextpt = idx_pt - increas;
			prevpt = idx_pt + increas;
		}	
		else
		{
			nextpt = idx_pt + increas;
			prevpt = idx_pt - increas;
		}
		pretx = s_zerolevpt[prevpt].bxcoor;
		prety = s_zerolevpt[prevpt].bycoor;
		nextx = s_zerolevpt[nextpt].bxcoor;
		nexty = s_zerolevpt[nextpt].bycoor;
		//get perpendicular vector
		xy_fromvec = calc_vector(txori, tyori, nextx, pretx, nexty, prety, 0, 0, 2);
		vnewx = xy_fromvec[0];
		vnewy = xy_fromvec[1];
		/*!< Negative search */
		negpt = calc_cross(file_cross, raster, s_arrayh, dem_nulval, \
		vnewx, vnewy, 1, cros_maxwith, cros_height, \
		njerar, nrio, ntramo, ndist, txori, tyori, tzori, idxfy, idxcx, tdx);
		/*!< Positive search */
		postpt = calc_cross(file_cross, raster, s_arrayh, dem_nulval, \
		vnewx, vnewy, 2, cros_maxwith, cros_height, \
		njerar, nrio, ntramo, ndist, txori, tyori, tzori, idxfy, idxcx, tdx);
	}
	printf("End cross section calculation in rio %i and pt %i\n", nrio, idx_pt);
	printf("total cross section points added %i, neg %i post %i\n", negpt+postpt, negpt, postpt);
	printf("---------------------------------\n");			
		
}	


/*! Auxiliary function of schilling implementation
 */
double loop_laharz(int isfirst, const char *file_cross, double **raster, struct HeadR s_arrayh[], float dem_nulval, struct Zeropoint s_zerolevpt[], \
	float cros_heincre, double area_section, float maning_coef, double vnewx, double vnewy, \
	int idx_pt, int idpt, int nrio, double txori, double tyori, double tzori, double slope_por, float tdx)
{
FILE *file;
float incre, discre;
double nwtzori;
int i, n, p, k, stopn, stopp, tot_pt, end_both, end_mainloop, n_loops;
int idx_cx, idx_fy, tot_idxcx, tot_idxfy;
double p4x, p4y, p4z, diff_topo, depth_topo, sum_areasec;
float crossdist;
double wetted_perimeter, crosec_area, hidrau_radius, speed_flow;
double planim_area, caudal_flow, arriv_time;
	
	tot_idxcx     = s_arrayh[0].hn_cx - 2;
	tot_idxfy     = s_arrayh[0].hn_fy - 2;
	if(tdx == 1)  discre = tdx;
	else          discre = tdx/2; //discrtize value
	incre         = 0;
	crossdist     = 0;
	diff_topo     = 0;
	planim_area   = 0;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if((file = fopen(file_cross, "a"))== NULL)
	{
		printf("-------ERROR open file--------\n");
		exit(0);
	}
	else
	{
		fprintf(file,"%i %i %lf %lf %lf %.2f %.2f\n",
				idpt,
				nrio,
				txori,
				tyori,
				tzori,
				crossdist,
				diff_topo			
				);
	
		/*!< Try calc section until the area is reached */
		end_mainloop=0;
		do
		{
			incre      += cros_heincre;      /*!< crossection height increase rate */ 
			nwtzori     = tzori + incre;     /*!< height from the tzori to the top of the flow */
			sum_areasec = incre * discre;    /*!< Area of first cell (alt*side) */
			stopn = stopp = tot_pt = k = n_loops = end_both = 0;
			n = p = 0;
			wetted_perimeter = 0;
			/*!< Calc section several times until reach the adecuate area */
			do
			{								
				if (stopn == 0)
				{
					p4x  = txori + (vnewx * discre * n);  /*!< Get oblique coordinates */
					p4y  = tyori + (vnewy * discre * n);
					idx_cx = calc_Rindex(p4x, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
					idx_fy = calc_Rindex(p4y, s_arrayh[0].hylo, s_arrayh[0].hresy);
					crossdist = calc_dist(p4x, txori, p4y, tyori);                                                                /*!< Distance between previous and current crossection point*/
					if(idx_fy > 2 && idx_fy < tot_idxfy && idx_cx > 2 && idx_cx < tot_idxcx) 
					{
						depth_topo=-9999;
						p4z  = raster[idx_fy][idx_cx]; //get z value
						if (p4z != dem_nulval)
						{
							diff_topo = (float)p4z - tzori;
							depth_topo = (float)p4z - nwtzori;
						}
						else                  
						{
							/*!< there could be a fail in the DEM values */		
							diff_topo  = get_newzval(raster, dem_nulval, idx_fy, idx_cx, tzori, 2);			
							depth_topo = get_newzval(raster, dem_nulval, idx_fy, idx_cx, nwtzori, 2);                                /*!< if positivo, celda above height*/
						}
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
						//printf("%lf %lf:: %lf %lf %lf :: %f %f\n", diff_topo, depth_topo, p4z, tzori, nwtzori, incre, cros_heincre);
						if (depth_topo != dem_nulval && depth_topo < 0)
						{
							sum_areasec += fabs(depth_topo * discre);
							s_crosec[tot_pt].cxcoor = p4x;
							s_crosec[tot_pt].cycoor = p4y;
							s_crosec[tot_pt].czcoor = p4z;
							s_crosec[tot_pt].ccrdis = crossdist;
							s_crosec[tot_pt].cdifto = diff_topo;
							s_crosec[tot_pt].cdepth = depth_topo;
							s_crosec[tot_pt].cidxfy = idx_fy;
							s_crosec[tot_pt].cidxcx = idx_cx;
							tot_pt++;
							wetted_perimeter += calc_hipo(crossdist, fabs(diff_topo));
							n--;
							//printf("n %i, and sum %lf ar of %lf :: k %i n %i :: %i :: %lf\n", n, sum_areasec, area_section, k, stopn, n_loops, incre);
						}
						else stopn = 1;
					}
				}
				//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				if (stopp == 0)
				{
					p4x  = txori + (vnewx * discre * p);  /*!< Get oblique coordinates */
					p4y  = tyori + (vnewy * discre * p);
					idx_cx = calc_Rindex(p4x, s_arrayh[0].hxlo, s_arrayh[0].hresx);                                               /*!< Get array index of new point*/
					idx_fy = calc_Rindex(p4y, s_arrayh[0].hylo, s_arrayh[0].hresy);
					crossdist = calc_dist(p4x, txori, p4y, tyori);
					if(idx_fy > 2 && idx_fy < tot_idxfy && idx_cx > 2 && idx_cx < tot_idxcx) 
					{
						depth_topo=-9999;
						p4z  = raster[idx_fy][idx_cx]; //get z value
						if (p4z != dem_nulval)
						{
							diff_topo = (float)p4z - tzori;
							depth_topo = (float)p4z - nwtzori;
						}
						else                  
						{
							/*!< there could be a fail in the DEM values */		
							diff_topo  = get_newzval(raster, dem_nulval, idx_fy, idx_cx, tzori, 2);			
							depth_topo = get_newzval(raster, dem_nulval, idx_fy, idx_cx, nwtzori, 2);                                /*!< if positivo, celda above height*/
						}
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
						//printf("%lf %lf:: %lf %lf %lf :: %f %f\n", diff_topo, depth_topo, p4z, tzori, nwtzori, incre, cros_heincre);
						if (depth_topo != dem_nulval && depth_topo < 0)
						{
							sum_areasec += fabs(depth_topo * discre);
							s_crosec[tot_pt].cxcoor = p4x;
							s_crosec[tot_pt].cycoor = p4y;
							s_crosec[tot_pt].czcoor = p4z;
							s_crosec[tot_pt].ccrdis = crossdist;
							s_crosec[tot_pt].cdifto = diff_topo;
							s_crosec[tot_pt].cdepth = depth_topo;
							s_crosec[tot_pt].cidxfy = idx_fy;
							s_crosec[tot_pt].cidxcx = idx_cx;
							tot_pt++;
							wetted_perimeter += calc_hipo(crossdist, fabs(diff_topo));
							p++;
							//printf("p %i and sum %lf ar of %lf :: k %i p %i :: %i :: %lf\n", p, sum_areasec, area_section, k, stopp, n_loops, incre);
						}
						else stopp = 1;
					}
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				if(stopp == 1 && stopn == 1)                                /*!< if out of raster or null values are found y both side before reach the correct area, start againg with a new height */
				{
					end_both = 1;                                           /*!< To keep searching or end the two loops when area is reached */
					break;                        
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				if (tot_pt > NUMCROSS)                                      /*!< Check total amount of crossec pt allowed */
				{
					printf("Total cross pt %i higher that allowed %i, reset define var\n", tot_pt, NUMCROSS);
					exit(0);
				}
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			}while(sum_areasec <= area_section);
			//printf("-------- pt %i n %i, p %i and sum %lf ar of %lf :: k %i n %i p %i end %i :: %i :: %lf\n", idx_pt, n, p, sum_areasec, area_section, k, stopn, stopp, end_both, n_loops, incre);
			
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			if (n_loops > 1000) break;                                 /*!< if too much loops are used to solve the area, stop the loop */
			n_loops++;
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			if (end_both == 1 && sum_areasec  < area_section) end_mainloop=0;  /*!< if both sides are close but area_section was not reached, continous search */
			if (end_both == 1 && sum_areasec >= area_section) end_mainloop=1;  /*!< if both sides are close and area_section is reached, stop searching */
			if (end_both == 0 && sum_areasec >= area_section) end_mainloop=1;  /*!< if both sides are open and area_section was reached, stop searching */
			if (end_both == 0 && sum_areasec  < area_section) end_mainloop=0;  /*!< if both sides are open and area_section was not reached, continous searching */                                             
			
		}while(end_mainloop != 1);	
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< write parameter in riverbed pt */
		if (isfirst == 1)                                               /*!< Only first time crossection is calculated */
		{
			speed_flow = 0;
			wetted_perimeter = wetted_perimeter + (tot_pt * discre);
			crosec_area      = sum_areasec;
			//planim_area      = (double)tot_pt * (2.0 * discre);
			hidrau_radius    = crosec_area / wetted_perimeter;
			//-----
			if (slope_por == 0) slope_por = 1;
			speed_flow       = (1 / maning_coef) * pow(hidrau_radius, 2.0/3.0) * pow((slope_por/100), 0.5);          /*!< Slope in tantos por 1: formula de maning https://es.wikipedia.org/wiki/F%C3%B3rmula_de_Manning */
			arriv_time       = tdx / speed_flow;                                                                     /*!< Time = space (m) / velocity (m/s) */
			//caudal_flow      = (1 / maning_coef) * (pow(crosec_area, (1/5)) * pow(wetted_perimeter, (1/2)));         /*!< Q = Volume (m3) / Time (seg)  Q = Area Tras-Sec (m2) * Velocidad (m/s)*/
			caudal_flow      = crosec_area * speed_flow;                                                             /*!< Q = Area Tras-Sec (m2) * Velocidad (m/s)*/
			//printf("---TEST mancoes %lf rad %lf slop %lf  = spee %lf\n", maning_coef, hidrau_radius, slope_por, speed_flow);
			s_zerolevpt[idx_pt].bcrossperim  = wetted_perimeter;       /*!< cross permimeter in m */
			s_zerolevpt[idx_pt].bcrossarea   = crosec_area;            /*!< cross area in m2 */
			s_zerolevpt[idx_pt].bareaplanim  = planim_area;            /*!< cumulated planimetric area in m2 */
			s_zerolevpt[idx_pt].bcrosshidrad = hidrau_radius;          /*!< radio hidraulico */
			s_zerolevpt[idx_pt].bmaninspeed  = speed_flow;             /*!< manning speed in m/s */
			s_zerolevpt[idx_pt].bttime       = arriv_time;             /*!< cumulated arrival time in seg */
			s_zerolevpt[idx_pt].bcaudal      = caudal_flow;            /*!< caudal in m3/seg */
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< write crosection */
		for(i=0;i<tot_pt;i++)
		{
			/*!< Write structure */ 
			fprintf(file,"%i %i %lf %lf %lf %.2f %.2f %.2f\n",
					idpt,
					nrio,
					s_crosec[i].cxcoor,
					s_crosec[i].cycoor,
					s_crosec[i].czcoor,
					s_crosec[i].ccrdis,
					s_crosec[i].cdifto,
					s_crosec[i].cdepth			
					);
			/*!< Add to raster */
			idx_fy = s_crosec[i].cidxfy;
			idx_cx = s_crosec[i].cidxcx;
			if (rast_multi[idx_fy][idx_cx] == -9999)planim_area += 	s_arrayh[0].hresx * s_arrayh[0].hresy;	 /*!< Only cell selected the first time */
			rast_multi[idx_fy][idx_cx] = 1;
			/*!< Reset structure */ 	
			s_crosec[i].cxcoor = 0;
			s_crosec[i].cycoor = 0;
			s_crosec[i].czcoor = 0;
			s_crosec[i].ccrdis = 0;
			s_crosec[i].cdifto = 0;
			s_crosec[i].cdepth = 0;
		}
		fclose(file);
	}
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	return planim_area;
}


/*! 
 * Based on Schilling 1998 equations (http://pubs.usgs.gov/of/1998/0638/report.pdf)
 * Similar results than LaharZ with a simplified implementation    
 * hazard zones US Geological Survey; Information Services [distributor],
 * US Geological Survey; Information Services 1998
 * Crosection area A = 0.05V^(2 / 3) 
 * Total planimetric area B = 200V^(2 / 3) 
 * V = volume
 * DFLOWZ
 * Matteo Berti, Alessandro Simoni. 2007. Prediction of debris flow inundation
 * areas using empirical mobility relationships. Geomorphology 90, 144 – 161
 * doi:10.1016/j.geomorph.2007.01.014
 * Crosection area A = 0.08V^(2 / 3) 
 * Total planimetric area B = 17V^(2 / 3) 
 * V = volume
 * 
 * Calc lahar in a specific zero level trajectory */ 
void calc_laharz(const char *file_cross, double **raster, struct HeadR s_arrayh[], float dem_nulval, \
	struct Zeropoint s_zerolevpt[], int n_bedxyz, int simple_flow, float cros_maxwith, float cros_heincre, double area_section, float maning_coef, double a0, double a1, double a2, \
	int idx_pt, int idpt, int nrio, double txori, double tyori, double tzori, double slope_por, int idxfy, int idxcx, float tdx, int tend, int triodir)
{
int increas, nextpt, prevpt;
double pretx, prety, nextx, nexty;
double *xy_fromvec, vnewx, vnewy;
double planim_area;

//Calc lahar values
//arealah  = 0.05 * pow(volumen, 2/3); //en metros cuadrados
//planilah = 200 * pow(volumen, 2/3);  //en metros cuadrados
	if (s_arrayh[0].hresx <= 5)                            increas = 2;
	if (s_arrayh[0].hresx >  5 && s_arrayh[0].hresx <= 10) increas = 3;
	if (s_arrayh[0].hresx > 10 && s_arrayh[0].hresx <= 20) increas = 4;
	if (s_arrayh[0].hresx > 20)                            increas = 5;
	if (idx_pt < increas || idx_pt > (n_bedxyz-increas))   increas = 1;
	
	/*! First approach, perpendicular to the flow direction */ 
	if(triodir == 1) 
	{
		nextpt = idx_pt - increas;
		prevpt = idx_pt + increas;
	}	
	else
	{
		nextpt = idx_pt + increas;
		prevpt = idx_pt - increas;
	}
	pretx = s_zerolevpt[prevpt].bxcoor;
	prety = s_zerolevpt[prevpt].bycoor;
	nextx = s_zerolevpt[nextpt].bxcoor;
	nexty = s_zerolevpt[nextpt].bycoor;
	//calculamos vector perpendicular
	xy_fromvec = calc_vector(txori, tyori, nextx, pretx, nexty, prety, 0, 0, 2);
	vnewx = xy_fromvec[0];
	vnewy = xy_fromvec[1];
	
	planim_area = loop_laharz(1, file_cross, raster, s_arrayh, dem_nulval, s_zerolevpt, \
	cros_heincre, area_section, maning_coef, vnewx, vnewy, \
	idx_pt, idpt, nrio, txori, tyori, tzori, slope_por, tdx);
	/*! Recalc crossections in oblicual direction to fill afected area */ 
	if (simple_flow != 1)
	{
		/*! Second approach, perpendicular to the flow direction */ 
		pretx = s_zerolevpt[prevpt].bxcoor - increas;
		prety = s_zerolevpt[prevpt].bycoor - increas;
		nextx = s_zerolevpt[nextpt].bxcoor + increas;
		nexty = s_zerolevpt[nextpt].bycoor + increas;
		xy_fromvec = calc_vector(txori, tyori, nextx, pretx, nexty, prety, 0, 0, 2);
		vnewx = xy_fromvec[0];
		vnewy = xy_fromvec[1];
		
		planim_area += loop_laharz(0, file_cross, raster, s_arrayh, dem_nulval, s_zerolevpt, \
		cros_heincre, area_section, maning_coef, vnewx, vnewy, \
		idx_pt, idpt, nrio, txori, tyori, tzori, slope_por, tdx);
		
		/*! Third approach, perpendicular to the flow direction */ 
		pretx = s_zerolevpt[prevpt].bxcoor + increas;
		prety = s_zerolevpt[prevpt].bycoor + increas;
		nextx = s_zerolevpt[nextpt].bxcoor - increas;
		nexty = s_zerolevpt[nextpt].bycoor - increas;
		xy_fromvec = calc_vector(txori, tyori, nextx, pretx, nexty, prety, 0, 0, 2);
		vnewx = xy_fromvec[0];
		vnewy = xy_fromvec[1];
		
		planim_area += loop_laharz(0, file_cross, raster, s_arrayh, dem_nulval, s_zerolevpt, \
		cros_heincre, area_section, maning_coef, vnewx, vnewy, \
		idx_pt, idpt, nrio, txori, tyori, tzori, slope_por, tdx);
	}
	/*! Add planimetric area */
	s_zerolevpt[idx_pt].bareaplanim  = planim_area;     //based on total number of new cells ocupated by the flow
}

#endif /* _MUL_RAS */
