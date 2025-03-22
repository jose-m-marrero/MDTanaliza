/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza  
* Module:  struc_mdtanaliza.h
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2022-01-18
* Last update: 2025-01-14
* 
* Description:
* struc variables used by modules of MDTanaliza
* *********************************************************************/


#ifndef _XYZ_F
#define _XYZ_F


/*! For input init points in modify DEM, gravitational flows and interpolation
 * Three possible input formats
 * 1 xyz for DEM modification: id x y z
 * 2 xy for gravitational flow analysis: id x y
 * 3 idtypexy for interpolation: id code x y z
 */ 
typedef struct coordCSV 
{
	int ptidx;            /*!< Internal index */
	int ptid;             /*!< ID of each feature */
	int ptype;
	double ptxcoor;
	double ptycoor;
	double ptzcoor;
	//---
	float  ptidwp;       /*!< IDW' p parameter */
	float  ptidwdis;     /*!< IDW' searching distance parameter */
	//---
	int    ptjerar;
	int    ptrio;
	int    pttramo;
	float  ptdist;
	float  ptdx;
	float  ptdt;
	int    ptzmask;
	double ptlong;
	double ptlat;
	int    ptquality;
	
}coordCSV;


/*! For input zero level trajectories as input file in topohazard and multiflow algorithms
 * Gravitational flow path as output file
 */ 
typedef struct Zeropoint
{  
	int    bidx;               /*!< global index zpoint */
	int    bidpt;              /*!< global id zpoint */
	int    bidptinbed;         /*!< id zpoint in each bedriver */
	int    bjerar;             /*!< bedriver hierarchy */
	int    brio;                /*!< id bedriver */
	int    btramo;             /*!< id link of each bedriver */
	float  bdist;              /*!< Accumulated distance of each zpoint per bedriver */
	double bxcoor;
	double bycoor;
	double bzcoor;
	int    bzfix;              /*!< if zval was fixed, 1 by null, 2 by high */
	int    bzmask;             /*!< Inside mask area */
	double blong;
	double blat;
	float  bdx;                /*!< pt discretize distance */
	float  bdt;                /*!< pt discretize time */
	int    bquality;
	int    bidxfy;             /*!< row index of DEM */
	int    bidxcx;             /*!< col index of DEM */
	int    bidxmasfy;          /*!< row index of MASK */
	int    bidxmascx;          /*!< col index of MASK */
	int    bstreamdir;         /*!< flow direction of riverbed line */
	int    bends;              /*!< if pt is an ends */  
	int    bis_buff;           /*!< if pt is conected with other bedrives */   
	int    bnum_buf;           /*!< buffer numer */
	int    btyp_buf;           /*!< Type of buffer connection */
	float  bslope;             /*!< Slope-grandient in degrees of each bedriver pts */
	float  bslopepor;          /*!< Slope-grandient in % of each bedriver pts */
	float  bradius;            /*!< Radious of each bedriver pts, in Slope-gradient */
	//multiflow
	double bcrossperim;         /*!< cross permimeter in m */
	double bcrossarea;          /*!< cross area in m2 */
	double bcrosshidrad;        /*!< radio hidraulico */
	double bareaplanim;         /*!< cumulated planimetric area in m2 */
	int    binplanim;           /*!< if in cumulated planimetric area */
	double btimearrive;         /*!< arrival time in seg */
	double bttime;
	double bcaudal;             /*!< caudal in m3/seg */
	double bmaninspeed;         /*!< manning speed in m/s */
}Zeropoint;

/*! complementary structure for zero level trajectories */

typedef struct Traject
{
	int ridx;
	int rid;
	int rnio;
	double raddx;
	double raddy;
	double rslope;
	double rbrect;
	int    rtotpt;        /*!< total pt in bedrive */ 
	double rlenght;       /*!< lenght of the bedrive */ 
}Traject;


/*! Create cross sections with hydrographic parameters */

typedef struct  CrossSecc
{
	double cxcoor;
	double cycoor;
	double czcoor;
	double ccrdis;
	double cdifto;
	double cdepth;
	int    cidxfy;
	int    cidxcx;
	//double cttime;

	int    czhaz1;
	double czhafin;
	int  nrowhaz;
	int  ncolhaz;
	
}CrossSecc;


#endif /* _XYZ_F */

