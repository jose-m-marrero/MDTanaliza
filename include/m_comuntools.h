/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza  
* Module:  m_commontools.sh
* 
* Autores: Jose M. Marrero 
*          Hugo Yepes
* 
* Version: 0.1.11
* Creation Date: 2018-04-06
* Last update: 2025-12-03
* 
* Description:
* 
* Common auxiliary functions of MDTanaliza
* 
* *********************************************************************/
#ifndef _COM_RAS
#define _COM_RAS


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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>


/*----------------------------------------------------------------------
 |  Function definition
 ---------------------------------------------------------------------*/
double get_newzval(double **raster, float dem_nulval, int idx_fy, int idx_cx, double tzori, int typcal);


/*----------------------------------------------------------------------
 |  MDTanaliza's basic auxiliary functions
 ---------------------------------------------------------------------*/

/*! Function to check if a string is a valid integer -> provide by ChatGPT */
int is_integer(const char *str) {
    if (str[0] == '-' || str[0] == '+') str++; // Handle negative/positive sign

    if (*str == '\0') return 0; // Empty string after sign is invalid

    while (*str) {
        if (!isdigit(*str)) return 0; // If any character is not a digit, return false
        str++;
    }
    return 1; // Valid integer
}


/*! Random variate generator */
double box_muller(double m, double s)	/* normal random variate generator */
{				                        /* mean m, standard deviation s */
double x1, x2, w, y1;   
	do
	{
		x1 = 2.0*(double)rand()/(double)RAND_MAX-1.0;
		x2 = 2.0*(double)rand()/(double)RAND_MAX-1.0;
		w  = x1*x1 + x2*x2;
	}while (w >= 1.0);
	w  = sqrt((-2.0*log(w))/w);
	y1 = x1*w;
    return(m + y1*s);
}

/*! Random Gauss function */
int *calc_gaussb(struct HeadR s_arrayh[])
{
int i, dlbox;
static int gaussb[1000];
    dlbox = 100 * 2  * s_arrayh[0].hresx;
    for(i=0;i<1000;i++)gaussb[i] = (int)(box_muller(0.0,dlbox)/s_arrayh[0].hresx);
    return (gaussb);
}

/*! Random Gauss function with radius */
int *calc_gauss(struct HeadR s_arrayh[], float radiu_val)
{
int i, dlbox;
static int gauss[1000];  //mil posibles opciones
    dlbox = radiu_val * s_arrayh[0].hresx;  //distancia maxima en celdas
    for(i=0;i<1000;i++)gauss[i] = (int)(box_muller(0.0,dlbox)/s_arrayh[0].hresx);  //numero de celdas a recorrer 
	return (gauss);
}

/*! Auxiliary function to fix Z values in zero level trajectory */
double get_newzval(double **raster, float dem_nulval, int idx_fy, int idx_cx, double tzori, int typcal)
{
int l;
int cont_notnull, cont_topo;
double *topoval, sumcell, topo_cells[8], diff_topo, max_val;

	cont_notnull = cont_topo = sumcell = 0;
	diff_topo = -9999;
	max_val = 0;
	topoval = search_celproxF(raster, idx_fy, idx_cx);    /*!< get z valuels  from 8 neighbour cells */
	
	for(l=0;l<8;l++)
	{
		if (topoval[l] > max_val) max_val=topoval[l];    /*!< get max z value */
	}
	for(l=0;l<8;l++)
	{
		if (topoval[l] != dem_nulval)                    /*!< if z valuels is not null */
		{
			if (fabs(topoval[l]-max_val) < 100)           /*!< avoid error values with anormal zvalues */
			{
				topo_cells[cont_topo] = topoval[l];                   
				sumcell += topoval[l]; 
				cont_topo++;                             /*!< Count ok */
			}  
			else cont_notnull++;                  
		}
		else cont_notnull++;                             /*!< Count ok */
	}
	if (cont_notnull <= 3)                               /*!< if nulls are less than */
	{
		if (typcal == 1)                                 /*!< interpolated from media value */
		{
			diff_topo = sumcell / cont_topo;
		}
		if (typcal == 2)                                 /*!< diff topo from random selection */
		{
			int idx = get_randomI(0, cont_topo);         /*!< select random value */
			diff_topo = topo_cells[idx] - tzori;
		}
		
	}
	return diff_topo;
}

/*----------------------------------------------------------------------
 |  MDTanaliza's complex auxiliary functions
 ---------------------------------------------------------------------*/

/*! Alternative initial mode of MDTanaliza 
 * Create a default configuration input file
 * Set strategy number as input MDTanaliza's argument
 */
void print_config_file(int num_strategy)
{
FILE *file;
char namefile[55];
	sprintf(namefile,"default_global_configfile_strategy_%i.cfg", num_strategy);
	printf("\n Default configuration file %s\n", namefile);
	
	printf("Please, adapt parameter values according to your needs\n");

	if((file = fopen(namefile, "w"))== NULL)
	{
		printf("-------ERROR open file--------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		printf("-----------ERROR--------------\n");
		exit(0);
	}
	
	fprintf(file,"%s\n", "#*****INI********");
	fprintf(file,"%s%i\n", "STRATEGY:",num_strategy);
	fprintf(file,"%s\n", "SHORTNAMES:1");
	fprintf(file,"%s\n", "#*****DEM********");
	fprintf(file,"%s\n", "DEMFILE:/pathwithouthomeuser/DEMfilename.grd");
	fprintf(file,"%s\n", "DEM_TYPRAS:1");
	fprintf(file,"%s\n", "DEM_NULVAL:-9999");
	fprintf(file,"%s\n", "DEM_MINVAL:0");
	fprintf(file,"%s\n", "DEM_MAXVAL:2000");
	fprintf(file,"%s\n", "DEM_OUTFIL:0");
	fprintf(file,"%s\n", "DEM_CLIPED:0");
	fprintf(file,"%s\n", "DEM_NWXLOW:0.0");
	fprintf(file,"%s\n", "DEM_NWXHIG:0.0");
	fprintf(file,"%s\n", "DEM_NWYLOW:0.0");
	fprintf(file,"%s\n", "DEM_NWYHIG:0.0");
	
	fprintf(file,"%s\n", "#*****MASK********");
	fprintf(file,"%s\n", "MAS_IFUSED:0");
	fprintf(file,"%s\n", "MAS_TYPRAS:1");
	fprintf(file,"%s\n", "MAS_NULVAL:-9999");
	fprintf(file,"%s\n", "MAS_MINVAL:0");
	fprintf(file,"%s\n", "MAS_MAXVAL:5");
	fprintf(file,"%s\n", "MASFILE:/pathwithouthomeuser/MASKfilename.grd");
	
	if (num_strategy == 0 || num_strategy == 1)
	{
		fprintf(file,"%s\n", "#*XYZDEMOD********");
		fprintf(file,"%s\n", "TYPCHANG:1");
		fprintf(file,"%s\n", "#Format 1: id x y z");
		fprintf(file,"%s\n", "XYZFILE:/pathwithouthomeuser/coordfilename.csv");
	}
	if (num_strategy == 0 || num_strategy == 2)
	{
		fprintf(file,"%s\n", "#*****SINK********");
		fprintf(file,"%s\n", "MODESINK:1");
	}
	
	if (num_strategy == 0 || num_strategy == 3)
	{
		fprintf(file,"%s\n", "#*****ASPEC********");
		fprintf(file,"%s\n", "MODEASP:2");
	}
	if (num_strategy == 0 || num_strategy == 4)
	{
		fprintf(file,"%s\n", "#*****SLOPE********");
		fprintf(file,"%s\n", "MODESLOP:1");
	}
	
	
	if (num_strategy == 0 || num_strategy == 5)
	{
		fprintf(file,"%s\n", "#*****FLOWPATH********");
		fprintf(file,"%s\n", "MODEFLOW:2");
		fprintf(file,"%s\n", "DISTMAX:10000.0");
		fprintf(file,"%s\n", "CRITHEIG:2.0");
		fprintf(file,"%s\n", "HEIGINCR:0.25");
		fprintf(file,"%s\n", "ALLOWREP:5000");
		fprintf(file,"%s\n", "FORCINTE:0");
		fprintf(file,"%s\n", "NUMITER:100");
		fprintf(file,"%s\n", "RADIUS:100");
		fprintf(file,"%s\n", "RESTRICH:1.0");
		fprintf(file,"%s\n", "MODEDIST:1");
		
		fprintf(file,"%s\n", "#**ZEROINITPOINT********");
		fprintf(file,"%s\n", "#Format 2: id x y");
		fprintf(file,"%s\n", "ZEROINITPT:/pathwithouthomeuser/INITfilename.csv");
	}
	
	if (num_strategy == 0 || num_strategy == 6)
	{
		fprintf(file,"%s\n", "#****---------------TOPOHAZARD------------------------*****");
		fprintf(file,"%s\n", "DIRECTOPO:3");
		fprintf(file,"%s\n", "FROMCENTER:1");
		fprintf(file,"%s\n", "FROMLIMASK:1");
		fprintf(file,"%s\n", "TOPDISTMAX:0");
		
		fprintf(file,"%s\n", "#**ZEROLEVELPOINT********");
		fprintf(file,"%s\n", "ZCORRECTION:1");
		fprintf(file,"%s\n", "ZDIST:100");
		fprintf(file,"%s\n", "#Format ZLP: jerar traject line ptdits x y discrex discret");
		fprintf(file,"%s\n", "ZEROLEVPT:/pathwithouthomeuser/ZLPfilename.csv");
	}
	
	if (num_strategy == 0 || num_strategy == 7)
	{
		fprintf(file,"%s\n", "#****---------------MULTIFLOW------------------------*****");
		fprintf(file,"%s\n", "MODEMULTI:2");
		fprintf(file,"%s\n", "CROSMAXWITH:100");
		fprintf(file,"%s\n", "CROSHEIGHT:0");
		fprintf(file,"%s\n", "CROSJUMPT:3");
		fprintf(file,"%s\n", "CROSHEINCRE:0.2");
		fprintf(file,"%s\n", "LAHARVOLUM:30000.0");
		fprintf(file,"%s\n", "SIMPLEDIREC:1");
		fprintf(file,"%s\n", "FLOWMAXDIS:10000.0");
		fprintf(file,"%s\n", "MANINGCOEF:0.35");
		fprintf(file,"%s\n", "SEARCHRADIUS:4");
		
		fprintf(file,"%s\n", "#**ZEROLEVELPOINT********");
		fprintf(file,"%s\n", "ZCORRECTION:1");
		fprintf(file,"%s\n", "ZDIST:100");
		fprintf(file,"%s\n", "#Format ZLP: jerar traject line ptdits x y discrex discret");
		fprintf(file,"%s\n", "ZEROLEVPT:/pathwithouthomeuser/ZLPfilename.csv");
	}
	if (num_strategy == 0 || num_strategy == 8)
	{
		fprintf(file,"%s\n", "#****---------------VOLHAZ------------------------*****");
		fprintf(file,"%s\n", "MODVOLHAZ:1");
		fprintf(file,"%s\n", "OLEALTCOL:1");
		fprintf(file,"%s\n", "OLEANGLE:1");
		
		fprintf(file,"%s\n", "#**ZEROINITPOINT********");
		fprintf(file,"%s\n", "ZEROINITPT:/pathwithouthomeuser/INITfilename.csv");
	}
	if (num_strategy == 0 || num_strategy == 9)
	{
		fprintf(file,"%s\n", "#****---------------INTER------------------------*****");
		fprintf(file,"%s\n", "MODEINTER:1");
		fprintf(file,"%s\n", "IDWAPROACH:2");
		fprintf(file,"%s\n", "IDWSEARCHPT:0");
		fprintf(file,"%s\n", "IDWDIST:4000.0");
		fprintf(file,"%s\n", "IDWPOWER:4");
		fprintf(file,"%s\n", "IDWBYTYPE:0");
		fprintf(file,"%s\n", "IDWCONTIYES:1");
		fprintf(file,"%s\n", "IDWDISCOYES:1");
		fprintf(file,"%s\n", "IDWONSITEYES:1");
		fprintf(file,"%s\n", "IDWINTERYES:1");
		fprintf(file,"%s\n", "IDWBOUNLYES:1");
		fprintf(file,"%s\n", "IDWBOUNPYES:1");
		
		fprintf(file,"%s\n", "#**ZEROINITPOINT********");
		fprintf(file,"%s\n", "#Format 3: id type x y z");
		fprintf(file,"%s\n", "ZEROINITPT:/pathwithouthomeuser/INITfilename.csv");
	}
	
	if (num_strategy == 1 || num_strategy == 3 || num_strategy == 4 || num_strategy == 9)
	{
		fprintf(file,"%s\n", "#****---------SMOOTHING-----------------*****");
		fprintf(file,"%s\n", "SMOOTHMODE:2");
		fprintf(file,"%s\n", "SMOOTHPMIN:90");
		fprintf(file,"%s\n", "SMOOTHPMAX:110");
		fprintf(file,"%s\n", "SMOOTHTCEL:6");
	}
	fprintf(file,"%s\n", "#*****OUT********");
	fprintf(file,"%s\n", "DIR_OUTGEN:/pathwithouthomeuser/");
	fprintf(file,"%s\n", "#****ENDFILE****");

}

/**
* Calculating the strainght line by least squares
*/
double *get_linebyleastsquares(struct HeadR s_arrayh[], struct Zeropoint s_zerolevpt[], int n_xyz, double sumx, double sumy, double sumxy, double sumxx)
{
double medx, medy, slope1, slope, slope2;
int nextpt, prevpt;
double alfrad, alfang;
double txmed, tymed;
double brect;
double nextx, pretx, nexty, prety;
static double addcoor[4], *addfromvec;
//double angmax, angmin, angmax2, angmin2;
 /**
	* Calculating the strainght line by least squares
	*/
    medx    = sumx / n_xyz;                 //half distance in x
    medy    = sumy / n_xyz;                 //half distance in y
    slope1  = sumxy - ((sumx * sumy) /  n_xyz);
    slope2  = sumxx - (pow(sumx, 2) / n_xyz);
    slope   = slope1 / slope2;
    brect   = medy - (slope * medx);         //straight line b parameter
    alfrad  = atan (slope);
	alfang  = alfrad * 180 / PI;
	if(alfang < 0)alfang += 180;
	//angmax  = (alfang - 90) + 15;
	//angmin  = (alfang - 90) - 15;
	//angmax2 = ((alfang - 90) + 180) + 15;
	//angmin2 = ((alfang - 90) + 180) - 15;
	txmed   = s_zerolevpt[n_xyz/2].bxcoor;    //x-coordinate of the midpoint of the zero level trajectory
	tymed   = (slope * txmed) + brect;        //y-coordinate of the midpoint of the zero level trajectory
	nextpt  = (n_xyz/2)-10;
	prevpt  = (n_xyz/2)+10;
	//nextx = s_zerolevpt[nextpt].bxcoor;
	//pretx = s_zerolevpt[prevpt].bxcoor;
	//nexty = s_zerolevpt[nextpt].bycoor;
	//prety = s_zerolevpt[prevpt].bycoor;
	/*!< Calc next point in the same straight line y=ax+b */
	nextx = medx - (s_arrayh[0].hresx * 4);
	nexty = slope * nextx + brect;                  
	/*!< Calc previous point in the same straight line y=ax+b */
	pretx = medx + (s_arrayh[0].hresx * 4);
	prety = slope * pretx + brect;
	
	printf("%lf %lf\n", s_zerolevpt[nextpt+1].bxcoor, s_zerolevpt[prevpt-1].bycoor);
	printf("para vec %lf %lf %lf %lf %lf %lf %lf %lf :: %i %i\n", txmed, tymed, nextx, pretx, nexty, prety, slope, brect, nextpt, prevpt);
	addfromvec = calc_vector(medx, medy, nextx, pretx, nexty, prety, slope, brect, 2);
	
	addcoor[0] = addfromvec[0];
	addcoor[1] = addfromvec[1];
	addcoor[2] = slope;
	addcoor[3] = brect;
	
	printf("sumx %lf - sumy %lf - medx %lf - medy %lf\n", sumx, sumy, medx, medy);
	printf("slope %lf - b %lf - sumxy %lf - sumxx %lf\n", slope, brect, sumxy, sumxx);
	printf("alfrad %lf - alfang %lf\n", alfrad, alfang);
	printf("addx %lf - addy %lf\n", addcoor[0], addcoor[1]);
	printf("------------\n");
	return (addcoor);
}

/*! Fix Z values of zero level points - Strategies Topohazard and Multiflow  */ 
double **fix_zval(double **raster, struct HeadR s_arrayh[], float dem_nulval, struct Zeropoint s_zerolevpt[], int n_zlps, float xyz_zdist)
{
int i;
//int j, was_fix, npost fward;
int idxfy, idxcx;
int typfix, detc, fix_8, tot_prob, fix_pp;
//float tdisc;
double tdist, tz, new_zcoor, diff_coor;
//double tx, ty; 
	
	printf("Fixing possible zval error in cero-level trajectory\n");
	
	fix_8 = tot_prob = fix_pp = 0;
	/*!< total Zero Level Points */
	for(i=0;i<n_zlps;i++)                                             
	{
		//tx    = s_zerolevpt[i].bxcoor;                                   
		//ty    = s_zerolevpt[i].bycoor;
		tz    = s_zerolevpt[i].bzcoor;  /*!< ZLP's Z value */
		idxfy = s_zerolevpt[i].bidxfy;
		idxcx =	s_zerolevpt[i].bidxcx;
		tdist = s_zerolevpt[i].bdist;   /*!< distance between ZLP */
		//tdisc = s_zerolevpt[i].bdx;
		typfix = detc = 0;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
		/*!< DETECT PROBLEM */
		/*!< int the first 50 m */
		if (tdist <= 50)
		{
			if (tz == dem_nulval)
			{
				if(idxfy > 2 && idxfy < s_arrayh[0].hn_fy-2 && idxcx > 2 && idxcx < s_arrayh[0].hn_cx-2) typfix = 1;
				detc=1;
			}	
		}
		/*!< after the first 50 m */
		if (tdist > 50)
		{
			/*!< Check difference with pre */
			diff_coor = fabs(tz - s_zerolevpt[i-1].bzcoor);
			if (diff_coor > xyz_zdist)
			{
				detc=1;
				/*!< Check difference between pre and post if not null */
				if (s_zerolevpt[i+1].bzcoor != dem_nulval && s_zerolevpt[i-1].bzcoor != dem_nulval)
				{
					diff_coor = fabs(s_zerolevpt[i+1].bzcoor - s_zerolevpt[i-1].bzcoor);
					if (diff_coor < xyz_zdist)  typfix = 2;
					else 
					{
						/*!< otherwise use 8 cells */
						if(idxfy > 2 && idxfy < s_arrayh[0].hn_fy-2 && idxcx > 2 && idxcx < s_arrayh[0].hn_cx-2) typfix = 1;
					}
				}
				else 
				{
					/*!< otherwise use 8 cells */
					if(idxfy > 2 && idxfy < s_arrayh[0].hn_fy-2 && idxcx > 2 && idxcx < s_arrayh[0].hn_cx-2) typfix = 1;
				}
			}
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< SOLVE PROBLEM */
		/*!< using 8 cells */
		if (typfix == 1)
		{
			new_zcoor = get_newzval(raster, dem_nulval, idxfy, idxcx, 0, 1);
			if (new_zcoor != dem_nulval)
			{
				s_zerolevpt[i].bzcoor = new_zcoor;
				raster[idxfy][idxcx] = new_zcoor;
				s_zerolevpt[i].bzfix  = typfix;
				fix_8++;
			}
		}
		/*!< using pre and post pt zvalues */
		if (typfix == 2)
		{
			new_zcoor            = (s_zerolevpt[i+1].bzcoor + s_zerolevpt[i-1].bzcoor) / 2;
			s_zerolevpt[i].bzcoor = new_zcoor;
			raster[idxfy][idxcx] = new_zcoor;
			s_zerolevpt[i].bzfix  = typfix;
			fix_pp++;
		}
		if (detc > 0) tot_prob++;
		
		if (typfix == 1 || typfix == 2)
		{
			if (i > 0)
			{
				printf("Zval error detected in %i. Init zval %lf new zval %lf pre %lf :: post %lf, using approach %i\n", i, tz, new_zcoor, s_zerolevpt[i-1].bzcoor, s_zerolevpt[i+1].bzcoor, typfix);
			}
		}
	}
	printf("Total zval error detected  %i :: 8cells %i prepost %i ::: solved %i, Not solv %i\n", tot_prob, fix_8, fix_pp, (fix_8 + fix_pp), (tot_prob-(fix_8 + fix_pp)));
	printf("---------------------------------\n\n");
	return (raster);
}


/*! Buffer calculation - for zero level trajectory connections */ 
void calc_buff(struct Zeropoint s_zerolevpt[], int n_zlps)
{
int i, irio, ifbuff;	
int j, jrio;
float idx, radius_buff;
double ix, iy, jx, jy, disfin;
int ok_conected, cont_conect, cont_nocon;
int nw_idbuf;
	
	printf("Calculating connections in Zero Level Trajectory\n");
	/**
	* Variable Initialization
	*/	
	cont_conect = cont_nocon = 0;           /*!< Count connected and isolated buffers */
	nw_idbuf=1;                             /*!< Add buffer number */
	for(i=0;i<n_zlps;i++)
    {
		/**
		* Get data from zpoint structure in i
		*/
		irio    = s_zerolevpt[i].brio;
		ifbuff  = s_zerolevpt[i].bis_buff;
        ix      = s_zerolevpt[i].bxcoor;
		iy      = s_zerolevpt[i].bycoor;
		idx     = s_zerolevpt[i].bdx;          /*!< discretization value  */
		s_zerolevpt[i].bnum_buf = nw_idbuf;   /*!< Add buff number */
		radius_buff = idx + idx/2;           /*!< Radius of the buffer, slightly longer than the discritization  */
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		/*!< If we are in the end of the bedriver and it is the first time  */
		if (ifbuff == 0)   
		{
			/**
			* find a different pt in the search radius
			*/
			ok_conected = 0;
			for(j=0;j<n_zlps;j++)
			{
				/**
				* Get data from zpoint structure in j
				*/
				jrio    = s_zerolevpt[j].brio;
				jx      = s_zerolevpt[j].bxcoor;
				jy      = s_zerolevpt[j].bycoor;
				/*!< If the stream in j is different and its hierarchy is lower  */
				if(irio != jrio)   
				{
					/*!< Find the nearest z point  */
					disfin = calc_dist(ix, jx, iy, jy);
					if(disfin < radius_buff)                /*!< if inside the buff radius  */
					{
						s_zerolevpt[j].bis_buff = 1;
						s_zerolevpt[j].bnum_buf = nw_idbuf;
						s_zerolevpt[j].btyp_buf = 1;
						ok_conected = 1;						
					}
				}
			}
			/**
			* if buffer is conected with other rivers
			*/
			if (ok_conected == 1)
			{
				s_zerolevpt[i].bis_buff = 1;
				s_zerolevpt[i].bnum_buf = nw_idbuf;
				s_zerolevpt[i].btyp_buf = 1;
				cont_conect++;
			}
			else cont_nocon++;
			nw_idbuf++;
		}
	}
	printf("total number of analyzed buffer = %i of total pt %i\n", nw_idbuf-1, n_zlps);
	printf("Total number of buffer connected = %i\n", cont_conect);
	printf("Total number of pt isolated = %i\n", cont_nocon);
	printf("End buffer calculation\n");
	printf("---------------------------------\n\n");		
}

/*! CALC SLOPE-GRADIENT AND RADIOUS - ZPOINTS */
void calc_bedslope(struct HeadR s_arrayh[], struct Zeropoint s_zerolevpt[], int n_zlps)
{
int k, l, m, errores, pt_out;
int i, irio, nextpt, prevpt, increas, idxfy, idxcx;	
float idist, spa_dis, dist_increas;
double jz, tslope, diffz, dist, dist2;
double tx, ty, tz, px, py, nx, ny, midx, midy;
int tis_end;
	printf("Calculating Slope in Zero Level Trajectory\n");
	/**
	* Variable Initialization
	*/
	k = 0; 				      /*!< count processed zpoints  */
	l = 0;              	  /*!< count errors  */
	m = 0;         			  /*!< count duplicates  */
	pt_out = errores = 0;     /*!< count total errors  */
	for(i=0;i<n_zlps;i++)
	{
		/**
		* Get data from zpoint structure in i
		*/
		irio     = s_zerolevpt[i].brio;
		idist    = s_zerolevpt[i].bdist;       /*!< accumulated distance */
		spa_dis  = s_zerolevpt[i].bdx;         /*!< spatial discretize  */
		tx       = s_zerolevpt[i].bxcoor;
		ty       = s_zerolevpt[i].bycoor;
		tz       = s_zerolevpt[i].bzcoor;
		tis_end  = s_zerolevpt[i].bends;
		idxfy    = s_zerolevpt[i].bidxfy;
		idxcx    = s_zerolevpt[i].bidxcx;
		
		if(idxfy > 0 && idxfy < s_arrayh[0].hn_fy && idxcx > 0 && idxcx < s_arrayh[0].hn_cx) 
		{
			if(tis_end == 0)                         /*!< if middel point  */
			{
				/*!< the higher the resolution the larger the multiplier  */
				if (s_arrayh[0].hresx <= 5)   dist_increas = s_arrayh[0].hresx * 9;   /*!< previous or next point should be a given distance to improve slope calculation  */
				if (s_arrayh[0].hresx  > 5)   dist_increas = s_arrayh[0].hresx * 6;   
				increas = dist_increas / s_arrayh[0].hresx;                           /*!< number of points to reach such distance  */  
				if (i < increas || i > (n_zlps-increas))      increas = 1;            /*!< if pt is close to extremes we reduce the distance to 1 pt to avoid memory overflow  */

				prevpt = i - increas;						/*!< get previous zpoint index in i  */
				nextpt = i + increas;                       /*!< get next zpoint index in i  */
				/**
				* Calc Slope gradient 
				*/
				//printf("TEST: Interior i=%i  nextpt=%i  dist=%f\n",i, nextpt, dist_increas);
				jz                      = s_zerolevpt[nextpt].bzcoor;           /*!< get z value of the nex zpoint in i  */ 
				diffz                   = fabs(tz - jz);                       /*!< Calc z difference  */
				dist                    = calc_hipo(spa_dis, diffz);               /*!< dist based in height differences  */
				if(diffz  > 0)tslope    = atan(diffz/(dist))*180/PI;
				if(diffz == 0)tslope    = 0.00; 
				s_zerolevpt[i].bslope    = tslope;                              /*!< Save slope-grandient in degrees  */
				s_zerolevpt[i].bslopepor = tslope * 100 / 45;                /*!< Saving s-gradient in percentages  */
				k++;
				l++;
				m++;
				/**
				* Calc radious H = dist2 and W3 = dist H/2 + W2/8H
				*/ 
				px    = s_zerolevpt[prevpt].bxcoor;											/*!< Get coordinates  */
				py    = s_zerolevpt[prevpt].bycoor;
				nx    = s_zerolevpt[nextpt].bxcoor;
				ny    = s_zerolevpt[nextpt].bycoor;
				dist  = calc_dist(px, nx, py, ny);                          /*!< Calc distance between next and previous zpts  */
				midx  = (px + nx) / 2;                                      /*!< Calc middle point in the line formed by next and previous zpts  */
				midy  = (py + ny) / 2;
				dist2 = calc_dist(midx, tx, midy, ty);                      /*!< Calc distance between middle and central zpts  */
				if (dist2 < 0.0001) s_zerolevpt[i].bradius = 10000;								            /*!< Set radious Almost a straight line  */
				else                s_zerolevpt[i].bradius = (dist2 / 2) + ((pow(dist, 2)) / (8*dist2));	 	/*!< Calc radious */
				
			}
			if(tis_end == 1)          /*!< if extreme  */
			{
				if(idist == 0)nextpt = i + 1;      /*!< if it is the first zpoint */
				if(idist  > 0)nextpt = i - 1;      /*!< if it is the last zpoint */
				//printf("TEST: Exterior i=%i  nextpt=%i\n",i, nextpt);
				jz  = s_zerolevpt[nextpt].bzcoor;   /*!< get z value of the nex zpoint according to flag  */
				/**
				* Calc Slope gradient 
				*/
				diffz                   = fabs(tz - jz);
				dist                    = calc_hipo(spa_dis, diffz);          /*!< dist based in height differences  */
				if(diffz  > 0)tslope    = atan(diffz/(dist))*180/PI;
				if(diffz == 0)tslope    = 0.00; 
				s_zerolevpt[i].bslope    = tslope;	                        /*!< Save slope-grandient in degrees  */
				s_zerolevpt[i].bslopepor = tslope * 100 / 45;                /*!< Saving s-gradient in percentages  */
				s_zerolevpt[i].bradius   = -9999;
				k++;
				l++;
				m++;
			}
			/**
			* Detect errors	 
			*/																
			if(m == 2)printf("error en %i\n",i);   /*!< the zpoint has been calculated twice  */
			if(l == 0)                             /*!< if the slope-gradient is not calculated  */
			{
				s_zerolevpt[i].bslope  = -1;
				s_zerolevpt[i].bradius = -9999;
				printf("Error slope in trajectory number = %i and ZLP = %i and distance to init pt  %f\n", irio, i, idist);	
			}
			if(m==2 || l == 0) errores++;                 /*!< count total errors  */
			m=0;
			l=0;
		}
		else pt_out++;
			
	}
	printf("total bedriver point procesed %i of %i\n", k, n_zlps);	
	printf("total error detected %i\n", errores);
	printf("Out of raster %i\n", pt_out);
	printf("End slope calculation\n");
	printf("---------------------------------\n\n");			
}	


//**********************************************************************
/*! General writing functions */
//**********************************************************************


/*! Write gravitational flow path LHM or SSM from strategy 5  */
int write_flowpath(const char *namefile, const char *mode, int nfile, struct Zeropoint s_zerolevpt[], int n_gflowpath)
{
FILE *file;
int i;

	printf("\nWrite flow path in xyz file %i\n", nfile);
	i = 0;
	if((file = fopen(namefile, mode))== NULL)
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
		if (nfile == 0)
		{ 
			fprintf(file,"%s\n",
				"njear nrios ntramo ndist xcoor ycoor zcoor dx dt npt quality"); //primera linea
		}
		
		for(i=0;i<n_gflowpath;i++)
		{		
			fprintf(file,"%i %i %i %.2f %lf %lf %lf %.2f %.2f %i %i\n",
				
				s_zerolevpt[i].bjerar,
				s_zerolevpt[i].brio,
				s_zerolevpt[i].btramo,
				s_zerolevpt[i].bdist,
				s_zerolevpt[i].bxcoor,
				s_zerolevpt[i].bycoor,
				s_zerolevpt[i].bzcoor,
				s_zerolevpt[i].bdx,
				s_zerolevpt[i].bdt,
				s_zerolevpt[i].bidpt,
				s_zerolevpt[i].bquality			
				);	
		}
	}
	fclose(file);
	printf("Output file name = %s\n",namefile);
	printf("End writing output file\n");
	printf("---------------------------------\n");
	return 0;
}	

/*! Write topohazard ZLP from strategy 6  */
int write_zltopohaz(const char *namefile, struct Zeropoint s_zerolevpt[], int n_zlps)
{
FILE *file;
int i;

	i = 0;
	if((file = fopen(namefile, "w"))== NULL)
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
			"njear nrios ntramo ndist xcoor ycoor zcoor zfix dx dt npt quality idxfy idxcx idxmfy idxmcx direcc ifend isbuff numbuff typbuf slope bradius"); //primera linea
		for(i=0;i<n_zlps;i++)
		{		
			fprintf(file,"%i %i %i %.2f %lf %lf %lf %i %.2f %.2f %i %i %i %i %i %i %i %i %i %i %i %f %f\n",
				
				s_zerolevpt[i].bjerar,
				s_zerolevpt[i].brio,
				s_zerolevpt[i].btramo,
				s_zerolevpt[i].bdist,
				s_zerolevpt[i].bxcoor,
				s_zerolevpt[i].bycoor,
				s_zerolevpt[i].bzcoor,
				s_zerolevpt[i].bzfix,
				s_zerolevpt[i].bdx,
				s_zerolevpt[i].bdt,
				s_zerolevpt[i].bidpt,
				s_zerolevpt[i].bquality,	
				s_zerolevpt[i].bidxfy,             /*!< row index of DEM */
				s_zerolevpt[i].bidxcx,             /*!< col index of DEM */
				s_zerolevpt[i].bidxmasfy,          /*!< row index of MASK */
				s_zerolevpt[i].bidxmascx,          /*!< col index of MASK */
				s_zerolevpt[i].bstreamdir,         /*!< flow direction of riverbed line */
				s_zerolevpt[i].bends,              /*!< if pt is an ends */  
				s_zerolevpt[i].bis_buff,           /*!< if pt is conected with other bedrives */   
				s_zerolevpt[i].bnum_buf,           /*!< buffer numer */
				s_zerolevpt[i].btyp_buf,           /*!< Type of buffer connection */
				s_zerolevpt[i].bslope,             /*!< Slope-grandient in degrees of each bedriver pts */
				s_zerolevpt[i].bradius            /*!< Radious of each bedriver pts, in Slope-gradient */		
				);	
		}
	}
	fclose(file);
	printf("Output file name = %s\n",namefile);
	printf("End writing output file\n");
	printf("---------------------------------\n");
	return 0;
}	

/*! Write mulfiflow ZLP from strategy 7  */
int write_zltmulti(const char *namefile, struct Zeropoint s_zerolevpt[], int n_zlps, int nriver)
{
FILE *file;
int i;

	i = 0;
	if((file = fopen(namefile, "w"))== NULL)
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
			"njear nrios ntramo ndist xcoor ycoor zcoor dx dt npt quality idxfy idxcx idxmfy idxmcx direcc ifend isbuff numbuff typbuf slopedeg bradius cperim carea chidrad acumpla inplanic flowspeed time timearr caudal"); //primera linea
		for(i=0;i<n_zlps;i++)
		{		
			if (nriver == s_zerolevpt[i].brio)
			{
				fprintf(file,"%i %i %i %.2f %lf %lf %lf %.2f %.2f %i %i %i %i %i %i %i %i %i %i %i %f %f %lf %lf %lf %lf %i %lf %lf %lf %lf\n",
					
					s_zerolevpt[i].bjerar,
					s_zerolevpt[i].brio,
					s_zerolevpt[i].btramo,
					s_zerolevpt[i].bdist,
					s_zerolevpt[i].bxcoor,
					s_zerolevpt[i].bycoor,
					s_zerolevpt[i].bzcoor,
					s_zerolevpt[i].bdx,
					s_zerolevpt[i].bdt,
					s_zerolevpt[i].bidpt,
					s_zerolevpt[i].bquality,	
					s_zerolevpt[i].bidxfy,             /*!< row index of DEM */
					s_zerolevpt[i].bidxcx,             /*!< col index of DEM */
					s_zerolevpt[i].bidxmasfy,          /*!< row index of MASK */
					s_zerolevpt[i].bidxmascx,          /*!< col index of MASK */
					s_zerolevpt[i].bstreamdir,         /*!< flow direction of riverbed line */
					s_zerolevpt[i].bends,              /*!< if pt is an ends */  
					s_zerolevpt[i].bis_buff,           /*!< if pt is conected with other bedrives */   
					s_zerolevpt[i].bnum_buf,           /*!< buffer numer */
					s_zerolevpt[i].btyp_buf,           /*!< Type of buffer connection */
					s_zerolevpt[i].bslope,             /*!< Slope-grandient in degrees of each bedriver pts */
					s_zerolevpt[i].bradius,            /*!< Radious of each bedriver pts, in Slope-gradient */	
					
					s_zerolevpt[i].bcrossperim,         /*!< cross permimeter in m */
					s_zerolevpt[i].bcrossarea,          /*!< cross area in m2 */
					s_zerolevpt[i].bcrosshidrad,        /*!< radio hidraulico */
					s_zerolevpt[i].bareaplanim,         /*!< cumulated planimetric area in m2 */
					s_zerolevpt[i].binplanim,           /*!< is in cumulated planimetric area */
					s_zerolevpt[i].bmaninspeed,         /*!< manning speed in m/s */
					s_zerolevpt[i].bttime,              /*!< displacement time per section in seg */
					s_zerolevpt[i].btimearrive,         /*!< arrival time in seg */
					s_zerolevpt[i].bcaudal              /*!< caudal in m3/seg */
					);	
			}
		}
	}
	fclose(file);
	printf("Output file name = %s\n",namefile);
	printf("End writing output file\n");
	printf("---------------------------------\n");
	return 0;
}	


/*! write output raster from interpolate strategy 9*/
void write_outputinter(double **raster, struct HeadR s_arrayh[], float dem_nulval, int r_index, char *out_namfile, int typmeasure)
{
int i, j, celltot;
double max,sumtot,volum, value;
	/**
	* Detect max value
	*/
	max=volum=celltot=sumtot=0;
	for(i=0;i<s_arrayh[r_index].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[r_index].hn_cx;j++)
		{
			if (raster[i][j] > max) max = raster[i][j];
			
			if (raster[i][j] != dem_nulval)
			{
				sumtot+=raster[i][j];
				if (typmeasure == 1) value = raster[i][j]/1000;  /*! mm */
				if (typmeasure == 2) value = raster[i][j]/100;   /*!cm */
				if (typmeasure == 3) value = raster[i][j]/10;    /*!dc */
				if (typmeasure == 4) value = raster[i][j]*10;    /*!dca */
				if (typmeasure == 5) value = raster[i][j]*100;   /*!hecto */
				if (typmeasure == 6) value = raster[i][j]*1000;  /*!km */
				volum+= value * s_arrayh[r_index].hresx * s_arrayh[r_index].hresx;
				celltot++;
			}
		}
	}
	printf ("Total sum = %lf\n", sumtot);
	printf ("Valid cells = %i\n", celltot);
	printf ("volume in cubic meters = %lf\n", volum);
	/**
	* Change max value in zmax
	*/
	s_arrayh[r_index].hzhi = max;
	/**
	* Write raster path
	*/
	//char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
	write_grdrasterF(out_namfile, raster, s_arrayh, 0, 1);
	//free(file_outpath);
}

#endif /* _COM_RAS */

/*!
 * UPDATEs
 * 2025/12/03: changed abs by fabs function to avoid compilation issues
 */
