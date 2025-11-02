/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza
* Module:  Main module
* 
* Author: Jose M. Marrero 
*          
* Version: 3.0.2
* Creation Date: 2022-01-18
* Last update: 2025-09-20
* 
* Compilation line: gcc MDTanaliza.c -o MDTanaliza -I/home/jmarrero/GIS/HAZARD_MODELS/SOFTWARE/c_libraries -lm
* Compilation line in repository: gcc MDTanaliza.c -o MDTanaliza -I./lib -lm
* Execution line:  ./MDTanaliza configfile.cfg  or ./MDTanaliza N (number to create configfile according to selected strategy)
* 
* Dependencies:
* General libraries for raster reading and writing (u_arrays.h), 
* some calculations (u_calculos),  handle strings for file names (u_strings) and 
* general purpose (u_generales). All of them in general_lib folder
* 
* Description:
* 
* Main module to manage input data and selected functions 
* Input configuration file is needed or a number to generate such file
* Number is referred to the strategy:
* 
* 1: Modify DEM
* 2: Calculate sinks
* 3: Calculate aspect (2 methods)
* 4: Calculate slope (9 methods)
* 5: Gravitational flows (4 methods)
* 6: Topohazard
* 7: Multiflow (2 methods)
* 8: Hazard, not implemented yet
* 9: Interpolation (IDW only)
* 
* Not implemented yet:
* 
* Multiflow: debris-flow
* hazard: lava flow; energy cone; ash
* 
* IMPORTANT: 
* Check TOTXYZ, TOTZLP, TOTZLTR, NUMARR values in glob_flow.h to avoid memory overflow
* 
* *********************************************************************/

/*----------------------------------------------------------------------
 |  Section Include
 ---------------------------------------------------------------------*/

/*!
* MDTanaliza's general modules 
***********************************************************************/ 
/*! Define global variables */
#include "./include/m_global.h"
/*! Define structures */
#include "./include/struc_mdtanaliza.h"
/*! Common calculation functions */
#include "./include/m_comuntools.h"
/*! 
 * MDTanaliza's basic strategy modules
 * Modules to assess Digital Elevation Model characteristics
***********************************************************************/ 
/*! Modify specific areas of DEM using a xyz file,  strategy module (1) */
#include "./include/m_modem.c"
/*! Fix or detect sink cells,  strategy module (2) */
#include "./include/m_sinks.c"
/*! Aspect calculation in degree or classes,  strategy module (3) */
#include "./include/m_aspect.c"
/*! Slope calculation in degree, classes and percentages,  strategy module (4) */
#include "./include/m_slope.c"

/*!
* MDTanaliza's specific modules
***********************************************************************/ 
/*! Gravity Flow path calculation, strategy module (5) */
#include "./include/m_flowpath.c"
/*! Topohazard,  strategy module (6) */
#include "./include/m_topohazard.c"
/*! Multiflow (Cross-section and Schilling Equations), strategy module (7) */
#include "./include/m_multiflow.c"
/*! Interpolation IDW, strategy module (9) */
#include "./include/m_inter_idw.c"

/*!
* General libraries
***********************************************************************/ 
/*! structure to save array's head data from *.grd input raster format*/
#include "struct_head-array.h"
/*! Some geometric calculation and data comparison functions */
#include "u_calculus.h"
/*! Reading and writing raster functions */
#include "u_arrays.h"
/*! Handle file names and paths, and other string functions */
#include "u_strings.h"
/*! Some general purpose functions */
#include "u_global.h"

/*! 
* Standard libraries
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------
 |  Section Global variables
 ---------------------------------------------------------------------*/

/*!
* Global variables declaration
*/ 
/*!
* Path names
***********************************************************************/ 
int n_strategy;
char cwd[MAX_PATH], *dir_out, *file_outdem;
char *gb_namefile, *namfile, *pathfile, *extfile;
//output
int nsuffix_short;
/*!
* Digital elevation model
***********************************************************************/ 
char *file_dem;
int dem_typras, dem_cliped, dem_outfil;
float dem_nulval, dem_minval, dem_maxval;
double dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig;
double **rast_dem;
/*!
* Mask raster file
***********************************************************************/ 
char *file_mas;
int mas_ifused, mas_typras;
float mas_nulval, mas_minval, mas_maxval;
int  **rast_mask;
/*!
* Modify DEM
***********************************************************************/ 
char *file_xyz, *namxyz;
int n_xyz, chan_mode;
/*!
* Zero Point Level
***********************************************************************/ 
int xyz_zfix;
float xyz_zdist;
/*!
* Digital Elevation Model characteristics
***********************************************************************/ 
int sink_mode, aspec_mode, slope_mode;
/*!
* Flow path or gravitational flows
***********************************************************************/ 
char *file_xyven, *namxyven;
int flow_mode, num_repit, for_inter, num_itera;
float max_dist, cri_heig, incre_heig, radiu_val, rest_heig;
double  **raster_path;
int **raster_control;
/*!
* Topohazard
***********************************************************************/ 
char *file_xybed, *namxybed;
int n_zlts, n_zlps, direc_calc, from_center, from_limask;
float topmax_dist;
double **raster_topohaz, **inhaz_topo, **outhaz_topo; 
/*!
* Multiflow
***********************************************************************/ 
int multi_mode, simple_flow;
float cros_maxwith, cros_jump, maning_coef, search_radi, cros_height;
float maxdist_flow, cros_heincre;
double area_section, plani_lahar, volum_lahar;
double peak_dischar, height_sec; //peak_discharge Qp
double a0, a1, a2;
double **rast_multi;
/*!
* Flow path or gravitational flows and Multiflow
***********************************************************************/ 
int dist_type;
/*!
* Volcanic hazard (not implemented yet)
***********************************************************************/ 
int volhaz_mode;
float ole_altcolu, ole_angle;
double **rast_volhaz;
/*!
* IDW interpolating
***********************************************************************/ 
int inter_mode, idw_approach, idw_searpts, idw_fastexit, idw_bytype;
int idw_contyes, idw_distyes, idw_onsiyes, idw_inteyes, idw_boulyes, idw_boupyes, idw_fillyes; 
float idw_distan, idw_power;
double **rast_inter;
//soft
int s_modesmooth, s_totoutcels;
float s_minpercent, s_maxpercent;


//old
int flowtyp, forcevar, nitera, mod, huso, hemis, ncentros;
float lmax, hl, distran, incre;

/*----------------------------------------------------------------------
 |  Section Struct definition
 ---------------------------------------------------------------------*/

/*!
* In MDTanaliza's libraries
* IMPORTANT: Check TOTXYZ, TOTZLP, TOTZLTR, values in glob_flow.h to avoid memory overflow
* ***********************************************************************/
/*! General input coordinates, in struc_mdtanaliza.h */
struct coordCSV s_coord[TOTXYZ];

/*! Input Zero level point (ZLP) or initial point in struc_mdtanaliza.h
 * Depending on the selected strategy, they can be initial point (DEM modify, flow path)
 * or Zero Level Point (topohazard, multiflow) */
struct Zeropoint s_zerolevpt[TOTZLP];

/*! Total Zero Level Trajectories under analysis in struc_mdtanaliza.h */
struct Traject s_trajectory[TOTZLTR];

/*!
* In General libraries
* IMPORTANT: Check NUMARR, in glob_flow.h to avoid memory overflow
***********************************************************************/
/*! Save array head data from grd format -> in struct_array-cabecera.h*/
struct HeadR s_arrayh[NUMARR];

/*----------------------------------------------------------------------
 |  Function definition
 ---------------------------------------------------------------------*/

void print_config_file(int num_strategy);
int read_xyz(char *namfile, int format_tipe);
double **calc_modem(double **raster,  struct HeadR s_arrayh[], struct coordCSV s_coord[], int tot_xyz, int tip_change);
double **fix_sinks(double **raster, struct HeadR s_arrayh[], char *outfile, float dem_nulval, int sink_mode);
void calc_saspect(const char *dir_out, const char *namfile, double **raster, struct HeadR s_arrayh[], float dem_nulval, int aspec_mode);
void chk_cfgfile(int n_strategy);
void smooth_funtion(int n_strategy, const char *dir_out, const char *namefile, float name_param1, \
	float name_param2, double **raster, struct HeadR s_arrayh[], int r_indx, float dem_nulval, \
	int r_softmod, float r_softmax, float r_softmin, float r_softcel);
void write_output(double **raster, struct HeadR s_arrayh[], float dem_nulval, int r_index, char *suffix, int typmeasure);


/*----------------------------------------------------------------------
 |  Reading input data section
 ---------------------------------------------------------------------*/

/*! Reading config file
 * Depending on the selected strategy, only the require parameters will be used */
int read_cfg(const char *nombre, const char *dirhome)
{
char suffix[50];
FILE *file;
float two_thirds, one_thirds;
	
	n_strategy = INITVAL;
	dem_typras = dem_outfil = dem_cliped = INITVAL;
	dem_nulval = dem_minval = dem_maxval = INITVAL;
	dem_nwxlow = dem_nwxhig = dem_nwylow = dem_nwyhig = INITVAL;
	mas_ifused = mas_typras =  INITVAL;
	mas_nulval = mas_minval = mas_maxval = INITVAL;
	chan_mode  = sink_mode  = aspec_mode = slope_mode = INITVAL;
	//--
	flow_mode = for_inter = num_itera = num_repit = dist_type = INITVAL;
	max_dist = cri_heig = incre_heig = radiu_val = rest_heig = dist_type = INITVAL;
	//--
	direc_calc = from_center = from_limask = INITVAL;
	topmax_dist = INITVAL;
	//--
	inter_mode = idw_approach = idw_searpts = idw_distan = idw_power = idw_bytype = INITVAL;
	idw_contyes = idw_distyes = idw_onsiyes = idw_inteyes = idw_boulyes = idw_boupyes = idw_fillyes = INITVAL;
	//---
	printf("\n***Reading cfg file, %s in folder %s***\n", nombre, dirhome);
	if ((file = fopen(nombre,"rt"))== NULL)
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
		char line[1024] = { 0 };
		char line2[1024] = { 0 };
		int  fin=0;                   /*!< If file is ending */ 
		int  nline=0;                 /*!< Count for lines */ 
		 /*!< feof needs one line more to detect file ending, 
		  * thus a null value is generate at the end */
		while (!feof(file))            
		{
			/*!< size string */ 
			memset(line, 0, 1024);
			memset(line2, 0, 1024);		 
			/*!< read line */        
			fgets(line, 1024, file);
			/*!< get comments and print them in command line */
			if (line[0] == '#')
			{
				printf("%s", line);
				if ((strcmp(line,"#****ENDFILE****\n") == 0)) fin=1;
				continue;
			}
			else
			{
				/*!< Copy line value before any change */
				strcpy(line2,line);
				/*!< strtok function changes original value, so we generate two tokens
				 * First to compare label
				 * Second, to get the parameter value */
				char *token  = strtok(line, ":");
				char *token2 = strtok(line2, ":");

				if (token !=NULL)
				{
					/*!< Get parameter value after colon */
					token2 = strtok(NULL, ":");
					
					if (token2 !=NULL)
					{
						/*!< remove \n and replace with string \0 */
						int size = strlen(token2);
						token2[size-1] = '\0';
						printf("%s: %s\n",token, token2);
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
						/*! Set strategy */
						if ((strcmp(token,"STRATEGY")      == 0))
						{
							n_strategy    = atoi(token2);    /*!< selected analysis strategy  */
							if (n_strategy == 1)
							{
								printf("Selected strategy MODIFY DEM in %i\n", n_strategy);
								/*! change input DEM file by adding a suffix */
								strcpy(suffix,"_modif.grd");
							}
							if (n_strategy == 2) printf("Selected strategy SINK %i\n", n_strategy);   /*!< Available algorithms 1 Detect only 2 Modify */
							if (n_strategy == 3) printf("Selected strategy ASPECT %i\n", n_strategy); /*!< Available algorithms 1 LHM 2 SSM */ 
							if (n_strategy == 4) printf("Selected strategy SLOPE %i\n", n_strategy);  /*!< Available algorithms 1-9 */
							if (n_strategy == 5) printf("Selected strategy GRAVITATIONAL FLOWS %i\n", n_strategy);
							if (n_strategy == 6) printf("Selected strategy TOPOHAZARD %i\n", n_strategy);
							if (n_strategy == 7) printf("Selected strategy MULTIFLOW %i\n", n_strategy);
							if (n_strategy == 8) printf("Selected strategy VOLCANIC HAZARDS %i\n", n_strategy);  //NOT IMPLEMENTED YET
							if (n_strategy == 9) printf("Selected strategy INTERPOLA %i\n", n_strategy);
							//if (n_strategy == 10) printf("Estrategia seleccionada RADIOVISUAL en %i\n", n_strategy);
						}
						if ((strcmp(token,"SHORTNAMES")    == 0)) nsuffix_short  = atoi(token2); /*!< if output short name will be used (1) */ 
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
						//Digital Elevation Model or main input raster
						if ((strcmp(token,"DEMFILE")    == 0))
						{
							file_dem = get_pathnam(dirhome, "none", token2, 2);  /*!< get filepath for DEM file in *.grd format */ 
							
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							char *gb_namefile = (char *)malloc(strlen(token2)+1); 
							if (gb_namefile== NULL)
							{
								printf("Attention, not memory space was generated for %s\n", gb_namefile);
								free(gb_namefile);
								exit(0);
							}
							else strcpy(gb_namefile, get_name(token2));
							printf("---gb_name %s\n",gb_namefile);
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							namfile  = remove_extension(gb_namefile);
							extfile  = get_exten(gb_namefile);							
							printf("---namefile %s\n",namfile);
							printf("---extefile %s\n",extfile);
                           								
							free(gb_namefile);
						}
						if ((strcmp(token,"DEM_TYPRAS")    == 0)) dem_typras     = atoi(token2);        /*!< Set input raster type, 1 bynary and 2 ASCII */
						if ((strcmp(token,"DEM_NULVAL")    == 0)) dem_nulval     = atof(token2);        /*!< Set null value, suggested -9999 */
						if ((strcmp(token,"DEM_MINVAL")    == 0)) dem_minval     = atof(token2);        /*!< Set mininum valid value */
						if ((strcmp(token,"DEM_MAXVAL")    == 0)) dem_maxval     = atof(token2);        /*!< Set maximum valid value */
						if ((strcmp(token,"DEM_OUTFIL")    == 0)) dem_outfil     = atoi(token2);        /*!< If DEM will be written in an output file, for testing */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						if ((strcmp(token,"DEM_CLIPED")    == 0)) dem_cliped     = atoi(token2);        /*!< If DEM will be clipped value 1-0 */
						if ((strcmp(token,"DEM_NWXLOW")    == 0)) dem_nwxlow     = atof(token2);        /*!< X min. coord. bottom-left corner */
						if ((strcmp(token,"DEM_NWXHIG")    == 0)) dem_nwxhig     = atof(token2);        /*!< X max. coord. bottom-right corner */
						if ((strcmp(token,"DEM_NWYLOW")    == 0)) dem_nwylow     = atof(token2);        /*!< Y min. coord. bottom-left corner */
						if ((strcmp(token,"DEM_NWYHIG")    == 0)) dem_nwyhig     = atof(token2);        /*!< Y max. coord. bottom-right corner */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						//Maask raster
						if ((strcmp(token,"MAS_IFUSED")    == 0)) mas_ifused     = atoi(token2);        /*!< Whether mask file is used or not (1) */ 
						if ((strcmp(token,"MAS_TYPRAS")    == 0)) mas_typras     = atoi(token2);        /*!< Set input raster type, 1 bynary and 2 ASCII */
						if ((strcmp(token,"MAS_NULVAL")    == 0)) mas_nulval     = atof(token2);        /*!< Set null value, suggested -9999 */
						if ((strcmp(token,"MAS_MINVAL")    == 0)) mas_minval     = atof(token2);        /*!< Set mininum valid value */
						if ((strcmp(token,"MAS_MAXVAL")    == 0)) mas_maxval     = atof(token2);        /*!< Set maximum valid value */
						if ((strcmp(token,"MASFILE")       == 0)) file_mas       = get_pathnam(dirhome, "none", token2, 2);
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						//Strategy 1
						if ((strcmp(token,"TYPCHANG")      == 0)) chan_mode      = atoi(token2);        /*!< Type mode (1) all 8 cells same value (2) all 8 cell random value */
						if ((strcmp(token,"XYZFILE")       == 0)) 
						{
							file_xyz = get_pathnam(dirhome, "none", token2, 2);
						
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							char *gb_xyz = (char *)malloc(strlen(token2)+1); 
							if (gb_xyz == NULL)
							{
								printf("Attention, not memory space was generated for %s\n", gb_xyz);
								free(gb_xyz);
								exit(0);
							}
							else strcpy(gb_xyz, get_name(token2));
							printf("---gb_name %s\n",gb_xyz);
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							namxyz  = remove_extension(gb_xyz);
							free(gb_xyz);
						}
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Strategies 2, 3 and 4
						if ((strcmp(token,"MODESINK")      == 0)) sink_mode  = atoi(token2);            /*!< Find sink and get and xy sink file. Modes: With no changes in DEM (1) or with changes (2)  */
						if ((strcmp(token,"MODEASP")       == 0)) aspec_mode = atoi(token2);            /*!< Calc aspect. Modes LHM (1) or Using the SSM (2) */
						if ((strcmp(token,"MODESLOP")      == 0)) slope_mode = atoi(token2);            /*!< If > 0 Calc slope. Modes 1 to 9 */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						//Strategy 5 (gravitational flows)
						if ((strcmp(token,"MODEFLOW")      == 0)) flow_mode  = atoi(token2);            /*!< Set strategy - Type of gravitational flow 1 (Single LHM) 2 (Single SSM) 3 (Drunk sailor) 4 (multipath) */
						if ((strcmp(token,"DISTMAX")       == 0)) max_dist   = atof(token2);            /*!< Max. distance in meters of flow path (flow_mode 1-4) */
						if ((strcmp(token,"CRITHEIG")      == 0)) cri_heig   = atof(token2);            /*!< Critical height in meters (flow_mode 1-4) */
						if ((strcmp(token,"HEIGINCR")      == 0)) incre_heig = atof(token2);            /*!< Increase height in meters in each iteration (flow_mode 1-3)*/
						if ((strcmp(token,"FORCINTE")      == 0)) for_inter  = atoi(token2);            /*!< Force interaction between different paths (1 yes - 0 No) (flow_mode 1,2) */
						if ((strcmp(token,"NUMITER")       == 0)) num_itera  = atoi(token2);            /*!< Total number of Interations (flow_mode 3) */
						if ((strcmp(token,"ALLOWREP")      == 0)) num_repit  = atoi(token2);            /*!< Total number of repeticion (flow_mode 3) */
						if ((strcmp(token,"RADIUS")        == 0)) radiu_val  = atof(token2);            /*!< Desplacement in meters of xy original vent-pt (flow_mode 3) */
						if ((strcmp(token,"RESTRICH")      == 0)) rest_heig  = atof(token2);            /*!< Restrict multiplow according to % of heigh difference (flow_mode 4)*/  
						if ((strcmp(token,"MODEDIST")      == 0)) dist_type  = atoi(token2);            /*!< How to calculate maximum runout distance reached by flow path: 1 (in 2dimension) 2 (in 3d) or 3 (by number cells) (flow_mode 1-4) */            
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Strategy 6 (Topohazard)
						if ((strcmp(token,"DIRECTOPO")       == 0)) direc_calc    = atoi(token2);        /*!< Diftop calc direction == 1 horiz - 2 vert - 3 cross - 4 Nearest */
						if ((strcmp(token,"FROMCENTER")      == 0)) from_center   = atoi(token2);        /*!< From riverbed */
						if ((strcmp(token,"FROMLIMASK")      == 0)) from_limask   = atoi(token2);        /*!< From mask limits */
						if ((strcmp(token,"TOPDISTMAX")      == 0)) topmax_dist   = atof(token2);        /*!< Max. perpendicular distance in meters if > 0, otherwise raster limit */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						//Strategy 7 (Multiflow)
						if ((strcmp(token,"MODEMULTI")       == 0)) multi_mode    = atoi(token2);        /*!< 1 Croos section only; 2 LaharZ */
						if ((strcmp(token,"CROSMAXWITH")     == 0)) cros_maxwith  = atof(token2);        /*!< Maximum lenght in meters of cross section to stop calculation, for flat areas */
						if ((strcmp(token,"CROSHEIGHT")      == 0)) cros_height   = atof(token2);        /*!< Maximum height in meters of cross section in multi_mode 1 */
						if ((strcmp(token,"CROSJUMPT")       == 0)) cros_jump     = atof(token2);        /*!< Increase height in meters in each iteration to jump obstacles, only considered if > 0*/ 
						if ((strcmp(token,"CROSHEINCRE")     == 0)) cros_heincre  = atof(token2);        /*!< Height added to the center cell to calculate area and permiter of crosection in flow simulations */
						if ((strcmp(token,"LAHARVOLUM")      == 0)) volum_lahar   = atof(token2);        /*!< lahar volume in cubic meters multi_mode 1-2 */
						if ((strcmp(token,"SIMPLEDIREC")     == 0)) simple_flow   = atoi(token2);        /*!< If 1, only perpendicular. Other than 1 use triple perpendicular direction | / \ */
						//if ((strcmp(token,"CROSMAXWITH")     == 0)) cros_maxwith  = atof(token2);        /*!< Maximum lenght in meters of cross section to stop calculation, for flat areas */
						if ((strcmp(token,"FLOWMAXDIS")      == 0)) maxdist_flow  = atof(token2);        /*!< Maximun leght in meters of flows, only considered if > 0 */
						if ((strcmp(token,"MANINGCOEF")      == 0)) maning_coef   = atof(token2);
						if ((strcmp(token,"SEARCHRADIUS")    == 0)) search_radi   = atof(token2);        /*!< search distance in meters to fill gaps when finish output multiflow raster - NOT IMPLEMENTED YET*/
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Strategy 8 - Not implemented yet
						if ((strcmp(token,"MODVOLHAZ")       == 0)) volhaz_mode   = atoi(token2); 
						if ((strcmp(token,"OLEALTCOL")       == 0)) ole_altcolu   = atof(token2);        /*!< Heigh in meters of column colapse, ~0 for block & ash */
						if ((strcmp(token,"OLEANGLE")        == 0)) ole_angle     = atof(token2);        /*!< Angle in degrees of... */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Strategy 9 (IDW Interpolation)
						if ((strcmp(token,"MODEINTER")       == 0)) inter_mode    = atoi(token2);         /*!< Interpolation method, only IDW (1) has been implemented */
						if ((strcmp(token,"IDWAPROACH")      == 0)) idw_approach  = atoi(token2);         /*!< Search approach 1 by search distance 2 by group of specific pts */
						if ((strcmp(token,"IDWSEARCHPT")     == 0)) idw_searpts   = atoi(token2);         /*!< Number of min pt needed for idw calculation (only if > 0) */
						if ((strcmp(token,"IDWDIST")         == 0)) idw_distan    = atoi(token2);         /*!< Searching distance */
						if ((strcmp(token,"IDWPOWER")        == 0)) idw_power     = atoi(token2);         /*!< P parameter to increase weight of known pts */
						//if ((strcmp(token,"IDWFASTEXIT")     == 0)) idw_fastexit  = atoi(token2);         /*!< When searching distance exit loop if searched pt are reached (1) otherwise, find more pts */
						if ((strcmp(token,"IDWBYTYPE")       == 0)) idw_bytype    = atoi(token2);          /*!< Change power or searching parameter according to pt type (1) */ 
						/*!<
						 * The following parameters were designed to interpolate isoplet maps
						 * Depending on type of input point characteristics, it can be considered
						 * or not in the interpolation process. If no distinction is available, set
						 * all values to 0. You need a column to filter them
						 */
						if ((strcmp(token,"IDWCONTIYES")     == 0)) idw_contyes   = atoi(token2);         /*!< point over continous isopleth (value 1 - yes 1) */ 
						if ((strcmp(token,"IDWDISCOYES")     == 0)) idw_distyes   = atoi(token2);         /*!< point over discontinous isopleth (value 2 - yes 1) */ 
						if ((strcmp(token,"IDWONSITEYES")    == 0)) idw_onsiyes   = atoi(token2);         /*!< point on site (value 3 - yes 1) */  
						if ((strcmp(token,"IDWINTERYES")     == 0)) idw_inteyes   = atoi(token2);         /*!< point over interpreted-added isopleth (value 4 - yes 1) */ 
						if ((strcmp(token,"IDWBOUNLYES")     == 0)) idw_boulyes   = atoi(token2);         /*!< boundary line added (value 5 - yes 1) */ 
						if ((strcmp(token,"IDWBOUNPYES")     == 0)) idw_boupyes   = atoi(token2);         /*!< boundary discretize point (value 6 - yes 1) */
						if ((strcmp(token,"IDWFILLEDYES")    == 0)) idw_fillyes   = atoi(token2);         /*!< Filled discretize point (value 7 - yes 1) */
						
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//common methods and parameters
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						//smoothing is a complemetary method valid for all strategies when output is a raster
						if ((strcmp(token,"SMOOTHMODE")        == 0)) s_modesmooth    = atoi(token2);         /*!< Smooth final raster data if needed !=0. */
						if ((strcmp(token,"SMOOTHPMIN")        == 0)) s_minpercent    = atof(token2);         /*!< Maximun porcentual difference allowed */
						if ((strcmp(token,"SMOOTHPMAX")        == 0)) s_maxpercent    = atof(token2);         /*!< Minimun porcentual difference allowed */
						if ((strcmp(token,"SMOOTHTCEL")        == 0)) s_totoutcels    = atoi(token2);         /*!< Minimun cell allowed overpassing max or min limit */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Initial points for strategies 5 and 8
						if ((strcmp(token,"ZEROINITPT")       == 0)) 
						{
							file_xyven = get_pathnam(dirhome, "none", token2, 2);      /*!< XY init pt to calculate a gravitational flow or other dynamics*/
						
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							char *gb_xyven = (char *)malloc(strlen(token2)+1); 
							if (gb_xyven== NULL)
							{
								printf("Attention, not memory space was generated for %s\n", gb_xyven);
								free(gb_xyven);
								exit(0);
							}
							else strcpy(gb_xyven, get_name(token2));
							printf("---gb_name %s\n",gb_xyven);
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							namxyven  = remove_extension(gb_xyven);
							free(gb_xyven);
						}
						
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//z correction of trajectory - Valid for st 6 and 7
						if ((strcmp(token,"ZCORRECTION")      == 0)) xyz_zfix   = atoi(token2);     /*!< if z value is fixed in both zero level pt input file and DEM (1) or not */
						if ((strcmp(token,"ZDIST")            == 0)) xyz_zdist  = atof(token2);     /*!< Minimum vertical difference between two pts to considere an anomaly in a zero level trajectory */
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
						//Zero level points in trajectories for strategies 6 and 7
						if ((strcmp(token,"ZEROLEVPT")       == 0)) 
						{
							file_xybed  = get_pathnam(dirhome, "none", token2, 2); /*!< zero level point filename*/
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							char *gb_xybed = (char *)malloc(strlen(token2)+1); 
							if (gb_xybed== NULL)
							{
								printf("Attention, not memory space was generated for %s\n", gb_xybed);
								free(gb_xybed);
								exit(0);
							}
							else strcpy(gb_xybed, get_name(token2));
							printf("---gb_name %s\n",gb_xybed);
							//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
							namxybed  = remove_extension(gb_xybed);
							free(gb_xybed);
						}
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						if ((strcmp(token,"DIR_OUTGEN")     == 0))
						{
							dir_out     = get_pathnam(dirhome, "none", token2, 2);                     /*!< Get output directory path */  						
						}
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						nline++;
					}
					else
					{
						if (fin == 0)
						{
							printf("Attention: reading error has been found in cfg file\n");
							printf("Check line %i with text: %s\n", nline, line);
							printf("Expected format should be LABEL:value, without black space");
							exit(0);
						}
					}
				}
			}
		}		
		printf("END PROCESSING cfg file\n\n");
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (n_strategy == 7)
		{
			two_thirds  = 2.0/3.0; 
			one_thirds  = 1.0/3.0; 
			if(multi_mode == 2)    //Lahar - laharz
			{
				area_section = 0.05 * pow(volum_lahar, two_thirds); //in m2
				plani_lahar  = 200 * pow(volum_lahar, two_thirds);  //in m2
				peak_dischar = 0.0135 * pow(volum_lahar, 0.780);    //for granular flows
				height_sec   = 0;
			}
			if(multi_mode == 3) //	Debris flow - DFLOWZ
			{
				area_section = 0.07 * pow(volum_lahar, two_thirds);     //cross section area in m2 - confined
				plani_lahar  = 18 * pow(volum_lahar, two_thirds);       //planimetric area in m2
				height_sec   = 0.06 * pow(volum_lahar, one_thirds);     //cross section heigh in m - not confined
				peak_dischar = 0.0188 * pow(volum_lahar, 0.790);        // for clay-rich flows
			}
			if(multi_mode == 4)
			{
				area_section  = 0; 
				plani_lahar   = 0;
				peak_dischar  = 0;
			}	
						
			if(peak_dischar >= 1000000) // peak discharge m3/s
			{
				a0 =  0.00909467;
				a1 =  0.0000904871;
				a2 =  0.000282666;
			}	
			if(peak_dischar >= 10000 && peak_dischar < 1000000)
			{
				a0 = -0.2710086;
				a1 =  0.0378719;
				a2 =  0.000110375;
			}
			if(peak_dischar >= 1000 && peak_dischar < 10000)
			{	
				a0 =  0.087511;
				a1 = -0.00889418;
				a2 =  0.0015254;
			}	
			if(peak_dischar >= 100 && peak_dischar < 1000)
			{	
				a0 =  0.300674;
				a1 = -0.0179581;
				a2 =  0.002099817;
			}
			if(peak_dischar >= 0 && peak_dischar < 100) //No data for peak discharge < 100 m3/s
			{	
				a0 = 0;
				a1 = 0;
				a2 = 0;
			}
			printf("Simulation data in mode %i: Area-sec %lf m2:: Plani-Area %lf m2:: Height sec %lf m:: peak-disc %lf\n", multi_mode, area_section,plani_lahar, height_sec, peak_dischar);
			printf("----------------------------------------------------\n");
			printf("----------------------------------------------------\n\n");
		}
	}
	fclose(file);
	/*!check input parameters */
	chk_cfgfile(n_strategy);
	return 1;		
}


/*! Checking config file parameters according to the selected strategy */
void chk_cfgfile(int n_strategy)
{
	printf("Checking strategy number\n");
	if (n_strategy < 1 && n_strategy > 9)
	{
		printf("Attention, no strategy has been set: %i\n", n_strategy);
			printf("Value must be between 1 and 9\n");
			printf("Strategy 1: Modify DEM\n");
			printf("Strategy 2: Sink detection or correction\n");   
			printf("Strategy 3: Aspect\n");   
			printf("Strategy 4: Slope\n");  
			printf("Strategy 5: Gravitational flows\n");
			printf("Strategy 6: Topohazard\n");
			printf("Strategy 7: Multiflow\n");
			printf("Strategy 8: Hazards - NOT IMPLEMENTED YET\n");
			printf("Strategy 9: Interpolate IDW\n");
			exit(0);
	}
	//******************************************************************
	//dem
	printf("Checking DEM input data\n");
	if (dem_typras == INITVAL || (dem_typras < 1 &&  dem_typras > 2))
	{
		printf("Attention, no type DEM has been set: %i\n", dem_typras);
		printf("Value must be 1 for binary and 2 for ASCII\n");
		exit(0);
	}
	if (dem_minval == INITVAL || dem_maxval == INITVAL)
	{
		printf("Attention, no minimum or maximum value in DEM has been set: %lf ;; %lf\n", dem_minval, dem_maxval);
		exit(0);
	}
	if (dem_cliped == INITVAL && dem_cliped != 1)
	{
		printf("Attention, clip activation with wrong value: %i\n", dem_cliped);
		printf("Value must be 1\n");
		exit(0);
	}
	if (dem_cliped == 1)
	{
		if (dem_nwxlow == INITVAL || dem_nwxhig == INITVAL)
		{
			printf("Attention, no X coordinate in clip has been set: %lf ;; %lf\n", dem_nwxlow, dem_nwxhig);
			exit(0);
		}
		if (dem_nwylow == INITVAL || dem_nwyhig == INITVAL)
		{
			printf("Attention, no Y coordinate in clip has been set: %lf ;; %lf\n", dem_nwylow, dem_nwyhig);
			exit(0);
		}
	}
	printf("No error were found in DEM input data\n\n");
	//******************************************************************
	//mask
	printf("Checking mask input data\n");
	if (mas_ifused == INITVAL && mas_ifused != 1)
	{
		printf("Attention, mask activation with wrong value: %i\n", mas_ifused);
		printf("Value must be 1\n");
		exit(0);
	}
	if (mas_ifused == 1)
	{
		if (mas_typras == INITVAL || (mas_typras < 1 &&  mas_typras > 2))
		{
			printf("Attention, no type MASK has been set: %i\n", mas_typras);
			printf("Value must be 1 for binary and 2 for SCII\n");
			exit(0);
		}
		if (mas_minval == INITVAL || mas_maxval == INITVAL)
		{
			printf("Attention, no minimum or maximum value in MASK has been set: %lf ;; %lf\n", mas_minval, mas_maxval);
			exit(0);
		}
		printf("No error were found in mask input data\n\n");
	}
	//******************************************************************
	//Modify DEM
	if (n_strategy == 1)
	{
		printf("Checking modify DEM mode (1)\n");
		if (chan_mode == INITVAL || (chan_mode < 1 &&  chan_mode > 2))
		{
			printf("Attention, wrong value in modify DEM strategy: %i\n", chan_mode);
			printf("Value must be 1 (same value) or 2 (random value)\n");
			exit(0);
		}
		printf("No error were found in modify DEM strategy\n\n");
	}
	//Sink
	if (n_strategy == 2)
	{
		printf("Checking sink strategy (2)\n");
		if (sink_mode == INITVAL || (sink_mode < 1 &&  sink_mode > 2))
		{
			printf("Attention, wrong value in sink mode: %i\n", sink_mode);
			printf("Value must be 1 (detection only) or 2 (detection and correction)\n");
			exit(0);
		}
		printf("No error were found in sink strategy\n\n");
	}
	//Aspect
	if (n_strategy == 3)
	{
		printf("Checking Aspect strategy (3)\n");
		if (aspec_mode == INITVAL || (aspec_mode < 1 &&  aspec_mode > 2))
		{
			printf("Attention, wrong value in aspect mode: %i\n", aspec_mode);
			printf("Value must be 1 (LHM) or 2 (SSM)\n");
			exit(0);    
		}
		printf("No error were found in Aspect strategy\n\n");
	}
	//Slope
	if (n_strategy == 4)
	{
		printf("Checking Slope strategy (4)\n");
		if (slope_mode == INITVAL || (slope_mode < 1 &&  slope_mode > 9))
		{
			printf("Attention, wrong value in slope mode: %i\n", slope_mode);
			printf("Value must be between 1-9\n");
			exit(0);    
		}
		printf("No error were found in Slope strategy\n\n");
	}
	//******************************************************************
	//gravitational flows
	if (n_strategy == 5)
	{
		printf("Checking strategy 5 - gravitational flows\n");
		if (flow_mode == INITVAL || (flow_mode < 1 &&  flow_mode > 4))
		{
			printf("Attention, wrong value in gravitational flow mode: %i\n", flow_mode);
			printf("Value must be 1 (Single LHM) 2 (Single SSM) 3 (Drunk sailor) 4 (multipath)\n");
			exit(0);
		}
		if (flow_mode >=1 &&  flow_mode <= 4)
		{
			//all models
			if (max_dist == INITVAL || max_dist < 0)
			{
				printf("Attention, wrong max-distance value or value not set: %lf\n", max_dist);
				printf("Value must be positive\n");
				exit(0);
			}
			if (cri_heig == INITVAL || cri_heig < 0)
			{
				printf("Attention, wrong critical-height value or value not set: %lf\n", cri_heig);
				printf("Value must be positive\n");
				exit(0);
			}
			if (radiu_val == INITVAL || radiu_val < 0)
			{
				printf("Attention, wrong initial location modification value or value not set: %lf\n", radiu_val);
				printf("Value must be 0 (no modification) or positive\n");
				exit(0);
			}
			if (dist_type == INITVAL || (dist_type < 1 &&  dist_type > 2))
			{
				printf("Attention, wrong type of calculation runout distance or value not set: %i\n", dist_type);
				printf("Value must be 1 (in 2d) 2 (in 3d) or 3 (cell number) \n");
				exit(0);
			}
			//specific for models
			if (flow_mode <= 3)
			{
				if (incre_heig == INITVAL || incre_heig < 0)
				{
					printf("Attention, wrong incremental-height value or value not set: %lf\n", incre_heig);
					printf("Value must be 0 or >0\n");
					exit(0);
				}
			}
			if (flow_mode ==1 || flow_mode == 2)
			{
				if (for_inter == INITVAL || (for_inter < 0 &&  for_inter > 1))
				{
					printf("Attention, wrong value in force interaction: %i\n", for_inter);
					printf("Value must be 1\n");
					exit(0);
				}
			}
			if (flow_mode == 3)
			{
				if (num_itera == INITVAL || num_itera < 0)
				{
					printf("Attention, wrong total number of iterations or value not set: %i\n", num_itera);
					printf("Value must be positive\n");
					exit(0);
				}
				if (num_repit == INITVAL || num_repit < 0)
				{
					printf("Attention, wrong total number of repetitions or value not set: %i\n", num_repit);
					printf("Value must be positive\n");
					exit(0);
				}
			}
			if (flow_mode == 4)
			{
				if (rest_heig == INITVAL || rest_heig < 0)
				{
					printf("Attention, wrong height difference or value not set: %lf\n", rest_heig);
					printf("Value must be positive\n");
					exit(0);
				}
			}
		} 
		printf("No error were found in gravitational flows strategy\n\n");
	}
	//******************************************************************
	//topohazard
	if (n_strategy == 6)
	{
		printf("Checking topohazard strategy %i\n", n_strategy);
		if (direc_calc == INITVAL || (direc_calc < 1 &&  direc_calc > 3))
		{
			printf("Attention, wrong value in direction of calculation: %i\n", direc_calc);
			printf("Value must be 1 (horizontal/parallel) ,2 (vertical/meridian) or 3 (oblique)\n");
			exit(0);
		}
		if (from_center == INITVAL || (from_center < 0 &&  from_center > 1))
		{
			printf("Attention, wrong value in from center: %i\n", from_center);
			printf("Value must be 1 or 0\n");
			exit(0);
		}
		if (from_limask == INITVAL && from_limask != 1)
		{
			printf("Attention, wrong value in from mask limit: %i\n", from_limask);
			printf("Value must be 1\n");
			exit(0);
		}
		if (from_limask == 1 && mas_ifused == 0)
		{
			printf("Attention: to use mask limits in topohazard a raster mask file is needed\n");
			printf("Check in configfile FROMLIMASK = %i and MAS_IFUSED = %i\n", from_limask, mas_ifused);
			exit(0);
		}
		if (topmax_dist == INITVAL && topmax_dist < 0)
		{
			printf("Attention, wrong value in maximum lateral distance: %lf\n", topmax_dist);
			printf("Value must be 0 ( con limitation) or > 0\n");
			exit(0);
		}
		printf("No error were found in Topohazard strategy\n\n");
	}
	//******************************************************************
	//multiflow
	if (n_strategy == 7)
	{
		printf("Checking Multiflow strategy (7)\n");
		if (multi_mode == INITVAL || (multi_mode < 1 &&  multi_mode > 2))
		{
			printf("Attention, wrong value in mode: %i\n", multi_mode);
			printf("Value must be 1 (cross section) or 2 (LahaZ type)\n");
			exit(0);
		}
		
		if (cros_maxwith == INITVAL || cros_maxwith < 0)
		{
			printf("Attention, wrong value in maximum cross section distance: %lf\n", cros_maxwith);
			printf("Value must be > 0\n");
			exit(0);
		}
		
		if (volum_lahar == INITVAL || volum_lahar < 0)
		{
			printf("Attention, wrong value in lahar volume: %lf\n", volum_lahar);
			printf("Value must be > 0\n");
			exit(0);
		}
		
		//cross se
		if (multi_mode == 1)
		{
			if (cros_height == INITVAL || cros_height < 0)
			{
				printf("Attention, wrong value in maximum height: %lf\n", cros_height);
				printf("Value must be > 0\n");
				exit(0);
			}
			if (cros_jump == INITVAL || cros_jump < 0)
			{
				printf("Attention, wrong value in increase height: %lf\n", cros_jump);
				printf("Value must be 0 or > 0\n");
				exit(0);
			}
		}
		if (multi_mode == 2)
		{
			if (cros_heincre == INITVAL || cros_heincre < 0) 
			{
				printf("Attention, wrong value in maximum height: %lf\n", cros_heincre);
				printf("Value must be > 0\n");
				exit(0);
			}
			if (maning_coef == INITVAL || maning_coef < 0)
			{
				printf("Attention, wrong value in Maning coefficient: %lf\n", maning_coef);
				printf("Value must be > 0\n");
				exit(0);
			}
			if (simple_flow == INITVAL || simple_flow != 1)
			{
				printf("Attention, wrong value in simple flow: %i\n", simple_flow);
				printf("Value must be 1\n");
				exit(0);
			}
			if (maxdist_flow == INITVAL || maxdist_flow < 0)
			{
				printf("Attention, wrong value in maximum flow distance: %lf\n", maxdist_flow);
				printf("Value must be > 0\n");
				exit(0);
			}
		}
	}
	//******************************************************************
	//not implemented yet
	if (n_strategy == 8)
	{
		printf("Attention, strategy 8 not implemented yet\n");
		printf("Try a different approach\n");
		exit(0);
	}
	//******************************************************************
	//interpolation
	if (n_strategy == 9)
	{
		printf("Checking Interpolation strategy (9)\n");
		if (inter_mode == INITVAL || inter_mode != 1)
		{
			printf("Attention, wrong value in interpolation method: %i\n", inter_mode);
			printf("Value must be 1 for IDW\n");
			exit(0);
		}
		if (idw_approach == INITVAL || (idw_approach < 1 &&  idw_approach > 3))
		{
			printf("Attention, wrong value in mode: %i\n", idw_approach);
			printf("Value must be 1 (search by distance), 2 (by pt fix) or 3 (by pt, varying value)\n");
			exit(0);
		}
		if (idw_searpts == INITVAL || idw_searpts < 0)
		{
			printf("Attention, wrong value in min pt needed: %i\n", idw_searpts);
			printf("Value must be 0  or > 0\n");
			exit(0);
		}
		if (idw_distan == INITVAL || idw_distan < 0)
		{
			printf("Attention, wrong value in searching distance: %lf\n", idw_distan);
			printf("Value must be > 0\n");
			exit(0);
		}
		if (idw_power == INITVAL || idw_power < 0)
		{
			printf("Attention, wrong value in weighted power: %lf\n", idw_power);
			printf("Value must be > 0, recommended around 4\n");
			exit(0);
		}
		if (idw_bytype == INITVAL || (idw_bytype < 0 &&  idw_bytype > 1))
		{
			printf("Attention, wrong value in parameter power shift by pt type: %i\n", idw_bytype);
			printf("Value must be 0 (no selection) or 1 (selected), only implemented for pt 4 in IDW approach 3\n");
			exit(0);
		}
		if (idw_contyes == INITVAL || (idw_contyes < 0 &&  idw_contyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_contyes);
			printf("Value must be 0 (no selection) or 1 (selected), continuous isoplet\n");
			exit(0);
		}
		if (idw_distyes == INITVAL || (idw_distyes < 0 &&  idw_distyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_distyes);
			printf("Value must be 0 (no selection) or 1 (selected), discontinuous isoplet\n");
			exit(0);
		}
		if (idw_onsiyes == INITVAL || (idw_onsiyes < 0 &&  idw_onsiyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_onsiyes);
			printf("Value must be 0 (no selection) or 1 (selected), onsite sample\n");
			exit(0);
		}
		if (idw_inteyes == INITVAL || (idw_inteyes < 0 &&  idw_inteyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_inteyes);
			printf("Value must be 0 (no selection) or 1 (selected), interpreted added isoplet\n");
			exit(0);
		}
		if (idw_inteyes == INITVAL || (idw_inteyes < 0 &&  idw_inteyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_inteyes);
			printf("Value must be 0 (no selection) or 1 (selected), interpreted added isoplet\n");
			exit(0);
		}
		if (idw_boulyes == INITVAL || (idw_boulyes < 0 &&  idw_boulyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_boulyes);
			printf("Value must be 0 (no selection) or 1 (selected), boundary line\n");
			exit(0);
		}
		if (idw_boupyes == INITVAL || (idw_boupyes < 0 &&  idw_boupyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_boupyes);
			printf("Value must be 0 (no selection) or 1 (selected), boundary discrete pt line\n");
			exit(0);
		}
		if (idw_fillyes == INITVAL || (idw_fillyes < 0 &&  idw_fillyes > 1))
		{
			printf("Attention, wrong value in pt type: %i\n", idw_fillyes);
			printf("Value must be 0 (no selection) or 1 (selected), filled discrete pt line\n");
			exit(0);
		}
	}
}

/*! Read initial point xyz data form DEM modification and flow paths
 * format_tipe =1 modify DEM, format id,x,y,z
 * format_tipe =2 init point, format id,x,y 
 * CHECK COLXYZ, COLVENT, COLINTER and TOTXYZ in glob_flow to control format columns and file size*/
int read_xyz(char *namfile, int format_tipe)
{
FILE *file;
char texto[256];
int i, j, ncol, tid, type;
double tx, ty, tz;

	i=j=0;	
	printf("\n***Reading initial xyz file: %s***\n", namfile);
	if ((file = fopen(namfile,"rt"))== NULL)
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
		fgets(texto,256,file); 
		/**
		* Reading xyz for DEM modification id,x,y,z
		*/ 
		if (format_tipe == 1)
		{
			 
			while ((ncol = fscanf(file,"%i %lf %lf %lf", 	
				&tid,
				&tx,
				&ty,
				&tz
				)) == COLXYZ)
			{		
				s_coord[i].ptidx     = i-1;
				s_coord[i].ptid      = tid;
				s_coord[i].ptxcoor   = tx;
				s_coord[i].ptycoor   = ty;
				s_coord[i].ptzcoor   = tz;
				//printf("%i %lf %lf %lf\n", tid, tx, ty, tz);
				i++;
			}
		}
		/**
		* Reading xy for gravitational flow analysis id,x,y
		*/
		if (format_tipe == 2)
		{
			while ((ncol = fscanf(file,"%i %lf %lf", 	
				&tid,
				&tx,
				&ty
				)) == COLVENT)
			{		
				s_coord[i].ptidx     = i-1;
				s_coord[i].ptid      = tid;     
				s_coord[i].ptxcoor   = tx;
				s_coord[i].ptycoor   = ty;
				//printf("%i %lf %lf\n", tid, tx, ty);
				i++;
			}
		}
		/**
		* Reading xy for interpolation id,code,x,y,z
		*/
		if (format_tipe == 3)
		{
			while ((ncol = fscanf(file,"%i %i %lf %lf %lf", 	
				&tid,
				&type,
				&tx,
				&ty,
				&tz
				)) == COLINTER)
			{		
				/*!< remove data by typology if type n is 0 */
				if (type == 1 && idw_contyes == 0)type = 8;
				if (type == 2 && idw_distyes == 0)type = 8;
				if (type == 3 && idw_onsiyes == 0)type = 8;
				if (type == 4 && idw_inteyes == 0)type = 8;
				if (type == 5 && idw_boulyes == 0)type = 8;
				if (type == 6 && idw_boupyes == 0)type = 8;
				if (type == 7 && idw_fillyes == 0)type = 8;
				
				if (type != 8)
				{
					s_coord[i].ptidx     = i-1;
					s_coord[i].ptid      = tid;
					s_coord[i].ptype     = type;
					s_coord[i].ptxcoor   = tx;
					s_coord[i].ptycoor   = ty;
					s_coord[i].ptzcoor   = tz;
					//printf("%i %lf %lf\n", tid, tx, ty);
					i++;
				}
				else j++;
			}
		}
	}
	fclose(file);
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	n_xyz = i;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	
	if (format_tipe == 1) chk_readfiles("xyz mod csv", n_xyz, TOTXYZ, 0, 0, ncol, COLXYZ, 1, 0, 0);
	if (format_tipe == 2) chk_readfiles("vents csv", n_xyz, TOTXYZ, 0, 0, ncol, COLVENT, 1, 0, 0);
	if (format_tipe == 3) 
	{
		printf("Total pt removed from the list %i from %i\n", j, j+i);
		chk_readfiles("inter csv", n_xyz, TOTXYZ, 0, 0, ncol, COLINTER, 1, 0, 0);
	}	
	printf("end read file\n");
	printf("---------------------------------\n\n");
	return 1;
}

/*! Read Zero Level Points from trajectories
 * The input file is obtained from data conversion into GIS (polyline to point) 
 * There could be more than one Zero Level Trajectory, but they will be analyze independently
 * For each pt, the z value, the mask value the direction of the flow and several
 * other variables will be also calculated while reading 
 * CHECK COLBED and TOTZLP in glob_flow to control format columns and file size*/
int read_zerolevpoint(char *namfile, double **rast_dem, struct HeadR s_arrayh[], float dem_nulval)
{
FILE *file;
char texto[256];
int i, j, k, ncol, idxfy, idxcx;
//int tot_null, fix_null, tot_high, fix_high;  //fixing data function removed from here
double tx, ty;
int tjerar, trio, trio2, ttramo, npts, okz, pt_out, pt_in;
float tdist, tdxpt, tdtpt;
int mask_in, n_zerolevtr;
int trioini;
double sumx, sumxx, sumxy, sumy, sumlengh, zcoorini, difzcoor;
double *rectaval;
	//tot_null=fix_null=tot_high=fix_high=0;
	i=j=0;	
	pt_in = pt_out = mask_in = n_zerolevtr = 0;
	sumx = sumxx = sumxy = sumy = 0;
	npts=1;
	printf("\n***Lectura archivo proximidad %s***\n", namfile);
	if ((file = fopen(namfile,"rt"))== NULL)
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
		fgets(texto,256,file); 

		/**
		* Reading xy for topohazar and laharz analysis
		*/      
		while ((ncol = fscanf(file,"%i %i %i %f %lf %lf %f %f", 
				&tjerar,
				&trio,
				&ttramo,
				&tdist,
				&tx,
				&ty, 
				&tdxpt,
				&tdtpt          
				)) == COLBED)
		{
			s_zerolevpt[i].bidx     = i-1;
			s_zerolevpt[i].bidpt    = i;
			s_zerolevpt[i].bjerar   = tjerar;
			s_zerolevpt[i].brio     = trio;
			s_zerolevpt[i].btramo   = ttramo;
			s_zerolevpt[i].bdist    = tdist;
			s_zerolevpt[i].bxcoor   = tx;
			s_zerolevpt[i].bycoor   = ty;
			s_zerolevpt[i].bdx      = tdxpt;
			s_zerolevpt[i].bdt      = tdtpt;
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			idxfy = calc_Rindex(ty, s_arrayh[0].hylo, s_arrayh[0].hresy);
			idxcx = calc_Rindex(tx, s_arrayh[0].hxlo, s_arrayh[0].hresx);
			if(idxfy > 0 && idxfy < s_arrayh[0].hn_fy && idxcx > 0 && idxcx < s_arrayh[0].hn_cx) 
			{
				s_zerolevpt[i].bzcoor  = rast_dem[idxfy][idxcx];
				s_zerolevpt[i].bidxfy  = idxfy;
				s_zerolevpt[i].bidxcx  = idxcx;
				//printf("pt %i in raster: Y %lf %lf %lf %i --- X %lf %lf %lf %i -- Z %lf\n", i, ty, s_arrayh[0].hylo, s_arrayh[0].hresy, idxfy, tx, s_arrayh[0].hxlo, s_arrayh[0].hresx, idxcx, s_zerolevpt[i].bzcoor);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				/**
				* If mask, get index raster from mask
				*/
				if (mas_ifused == 1)
				{
					idxfy = calc_Rindex(tx, s_arrayh[1].hylo, s_arrayh[1].hresy);
					idxcx = calc_Rindex(ty, s_arrayh[1].hxlo, s_arrayh[1].hresx);
					if(idxfy > 0 && idxfy < s_arrayh[1].hn_fy && idxcx > 0 && idxcx < s_arrayh[1].hn_cx) 
					{
						s_zerolevpt[i].bzmask     = rast_mask[idxfy][idxcx];
						s_zerolevpt[i].bidxmasfy  = idxfy;
						s_zerolevpt[i].bidxmascx  = idxcx;
						mask_in++;
					}
				}
				/**
				* Get stream extremes and flow direction
				*/
				if (pt_in == 0)										    /*!< first pt inside raster, avoid mistake when pt outside raster exist */
				{
					trioini  = trio;                                    /*!< Get trajectory number */
					//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					s_trajectory[n_zerolevtr].ridx    = 0;
					s_trajectory[n_zerolevtr].rid     = n_zerolevtr+1;
					s_trajectory[n_zerolevtr].rnio    = trio;
					//printf("river hier %i -- idrio %i\n", s_trajectory[n_zerolevtr].rid,  s_trajectory[n_zerolevtr].rnio );
					//--- hierarchy defined internaly. Both zero-level trajectories and global trajectory need same id 
					//s_zerolevpt[i].bjerar   = n_zerolevtr+1;
					//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					s_zerolevpt[i].bends = 1;  						    /*!< trajectory extremes */	
					n_zerolevtr++;			
				}		
				if (trioini != trio || feof(file))                      /*!< When trajectory change */
				{
					s_trajectory[n_zerolevtr].ridx    = n_zerolevtr;
					s_trajectory[n_zerolevtr].rid     = n_zerolevtr+1;
					s_trajectory[n_zerolevtr].rnio    = trioini;
					//-----
					//s_zerolevpt[i].bjerar   = n_zerolevtr+1;
					/*!< reset values  */
					trioini    = trio;
					n_zerolevtr++;  
					/*!< first step buff calc */
					s_zerolevpt[i-1].bends = 1;                         /*!< Add Buffer value in trajectory extremes */
					s_zerolevpt[i].bends   = 1;
					npts=1;
				}
				s_zerolevpt[i].bidptinbed = npts;                   	/*!< Add recno by ZLP */
				npts++;
				pt_in++;
			}
			else
			{
				//printf("pt %i out of raster: Y %lf %lf %lf %i --- X %lf %lf %lf %i\n", i, ty, s_arrayh[0].hylo, s_arrayh[0].hresy, idxfy, tx, s_arrayh[0].hxlo, s_arrayh[0].hresx, idxcx);
				s_zerolevpt[i].bzcoor  = -9999;
				s_zerolevpt[i].bidxfy  = -9999;
				s_zerolevpt[i].bidxcx  = -9999;
				s_zerolevpt[i].bends   = -9999;
				pt_out++;
			}
			//printf("%i :: %lf %lf :: %i %i  : %lf\n",i, tx, ty, idxfy, idxcx, s_zerolevpt[i].bzcoor);
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			if (i > TOTZLP)	
			{
				printf("-----------ERROR--------------\n");
				exit(0);
			}
			i++;                                                           
		}
	}
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//Add last extreme pt to the list
	s_zerolevpt[i-1].bends = 1;
	//----
	fclose(file);
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	n_zlps = i;                     /*!< total Zero Level Points  */
	n_zlts = n_zerolevtr;           /*!< total Zero Level Trajectories  */
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	if (direc_calc == 3)
	{
		for(i=0;i<n_zlts;i++)
		{
			npts=0;
			trio = s_trajectory[i].rnio;
			sumx =  sumxx = sumxy = sumy = sumlengh = npts = okz = 0;
			for(j=0;j<n_zlps;j++)
			{
				if (okz == 0)
				{
					zcoorini = s_zerolevpt[j].bzcoor;                    /*!< Get z value  */
					okz = 1;
				}
				trio2  = s_zerolevpt[j].brio;
				tx     = s_zerolevpt[j].bxcoor;
				ty     = s_zerolevpt[j].bycoor;
				tdist  = s_zerolevpt[j].bdist;
				if (trio == trio2)
				{
					s_trajectory[i].rtotpt  ++;
					sumx     += tx;
					sumxx    += pow(tx, 2);
					sumxy    += tx * ty;
					sumy     += ty;
					sumlengh += tdist;
					npts     ++;
				}
			}
			s_trajectory[i].rtotpt = npts;
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			rectaval = get_linebyleastsquares(s_arrayh, s_zerolevpt, npts, sumx, sumy, sumxy, sumxx);
			s_trajectory[i].raddx  = rectaval[0];
			s_trajectory[i].raddy  = rectaval[1];
			s_trajectory[i].rslope = rectaval[2];
			s_trajectory[i].rbrect = rectaval[3];
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			difzcoor = zcoorini - s_zerolevpt[j].bzcoor;           /*!< Calculate stream flow with a flag */
			for(k = npts; k > 0; k--) 						       /*!< reverse loop n j */
			{
				if (difzcoor > 0)s_zerolevpt[j-k].bstreamdir = 1; 	/*!< high values upstream direction */
				if (difzcoor < 0)s_zerolevpt[j-k].bstreamdir = 2; 	/*!< low values downstream direction */
			}
			
		}
	}
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	chk_readfiles("Check xy pts csv", n_zlps, TOTXYZ, 0, 0, ncol, COLVENT, 1, 0, 0);
	printf("Total zero-level trajectories %i\n", n_zlts);
	printf("Total ZLP inside of raster %i of %i\n", pt_in, n_zlps);
	printf("Total ZLP outside of raster %i of %i\n", pt_out, n_zlps);
	if (pt_out == n_zlps)
	{
		printf("Attention: all ZLP located outside the raster\n");
		printf("Check coordinates\n");
		exit(0);
	}
	//printf("Total high-diff values fixed %i de %i\n", fix_high, tot_high);
	printf("End read file\n");
	printf("---------------------------------\n\n");
	return 1;
}


/*----------------------------------------------------------------------
 |  Raster analysis strategies
 ---------------------------------------------------------------------*/
 
/*! Modify DEM
 * Change z values in a input DEM given a set of coordinates
 * In all conditions, same cell is changed by the new Z value
 * id chan_mode != 1 or 2, neighboring cells won't be altered 
 * if chan_mode 1, the eight neighboring cells are changed by same value
 * if chan_mode 2, the eight neighboring cells are changed by random value (0 and the new cell value)
 * */
void mod_dem(char *file_xyz, const char *file_demn, const char *dir_out, char *namfile)
{
char suffix[50];
	/*!< Read xyz file to change DEM values */
	read_xyz(file_xyz, 1);
	/*!< Read original DEM file */
	rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	/*!< Modify DEM file */
	rast_dem = calc_modem(rast_dem, s_arrayh, s_coord, n_xyz, chan_mode);
	/*!< Write DEM final file */
	strcpy(suffix,"_modificate.grd");
	char *file_outmod = get_pathnam(dir_out, namfile, suffix, 1);
	write_grdrasterF(file_outmod, rast_dem, s_arrayh, 0, 1);
	/*!< Generate smooth output raster - Optional */
	if (s_modesmooth > 0)
	{
		printf("--------------------------\n\n");
		smooth_funtion(n_strategy, dir_out, file_outmod, 0, 0, rast_dem, s_arrayh, 0, dem_nulval, s_modesmooth, s_maxpercent, s_minpercent, s_totoutcels);
	}
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);
}


/*! Search and detect sinks in DEM 
 * (1) sink_mode: only check for sink cells
 * (2) sink_mode: modify z values if sinks are isolated
 * */
void mod_sinks(const char *file_demn)
{
char suffix[50];
	printf("Calculating sinks\n");
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	
	/*!< Detect and create sink xyz file */
	strcpy(suffix,"_Sink_pt.xyz");
	char *file_outsink = get_pathnam(dir_out, namfile, suffix, 1);
	/*!< Calc fix sink */
	rast_dem = fix_sinks(rast_dem, s_arrayh, file_outsink, dem_nulval, sink_mode);
	if (sink_mode == 2)
	{
		strcpy(suffix,"_Sink.grd");
		char *file_outsink = get_pathnam(dir_out, namfile, suffix, 1);
		write_grdrasterF(file_outsink, rast_dem, s_arrayh, 0, 1);
	}
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);
}


/*! Calc aspect in classes and degrees from original DEM 
 * (1) aspec_mode only in classes
 * (2) aspec_mode classes and degrees
 * */
void mod_aspect(const char *file_demn, const char *dir_out, const char *namfile)
{
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	/*!< Calc aspect */
	calc_saspect(dir_out, namfile, rast_dem, s_arrayh, dem_nulval, aspec_mode);
	/*!< Generate smooth output raster - Optional */
	if (s_modesmooth > 0)
	{
		printf("--------------------------\n\n");
		smooth_funtion(n_strategy, dir_out, namfile, 0, 0, rast_dem, s_arrayh, 0, dem_nulval, s_modesmooth, s_maxpercent, s_minpercent, s_totoutcels);
	}
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);
}

/*! Calc slope in degrees, classes and percentages from original DEM 
 * (1) slope_mode Using the SSM Burrough and McDonell (1998)
 * (2) slope_mode Second-order finite difference 2FD Fleming and Hoffer (1979)
 * (3) slope_mode Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Distance  3FDWRD Unwin (1981)
 * (4) slope_mode Three-order  Finite  Difference,  Linear  regression  plan  3FD Sharpnack et al (1969)
 * (5) slope_mode Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Squared  Distance 3FDWRSD Horn (1981)
 * (6) slope_mode Frame Finite difference FFD Chu and Tsai (1995)
 * (7) slope_mode Maximum Max Travis et al. (1975)
 * (8) slope_mode Simple difference Simple-D Jones (1998)
 * (9) slope_mode Constrained Quadratic Surface Quadsurface (Wood Evans [13])
 * see m_slope to know more about these methods
 * */
void mod_slope(const char *file_demn, const char *dir_out, const char *namfile)
{
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	
	/*!< Calc aspect */
	calc_slope(dir_out, namfile, rast_dem, s_arrayh, dem_nulval, slope_mode);
	/*!< Generate smooth output raster - Optional */
	if (s_modesmooth > 0)
	{
		printf("--------------------------\n\n");
		smooth_funtion(n_strategy, dir_out, namfile, 0, 0, rast_dem, s_arrayh, 0, dem_nulval, s_modesmooth, s_maxpercent, s_minpercent, s_totoutcels);
	}
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);
}


/*! Calc gravitational flows for morphology DEM analysis
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
*/
void mod_gravflows(char *file_xyven, const char *file_demn, const char *dir_out, const char *namfile, int flow_mode, float max_dist)
{
//char suffix[100];
char *suffix;
int k, i, j, ptid, count_ok;
int indxc, indyf, max;
//int *gauss_belt;
double tx, ty, tz;
	count_ok=0;  /*!< count processed init xy */
	/*!< Read xy vent file to calculate gravitational flows */
	read_xyz(file_xyven, 2);
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	/*!< Create working rasters */
	raster_path    = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	raster_control = Crea_2DIarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++)
		{
			raster_path[i][j]    = 0;
			raster_control[i][j] = 0;
		}
	}
	/**
	* For each vent, calc the path using a gravitational flow model 
	*/
	for(k=0;k<n_xyz;k++)                                             /*!< total init z points */
	{
		tx   = s_coord[k].ptxcoor;                                   /*!< getting init z points values */
		ty   = s_coord[k].ptycoor;
		ptid = s_coord[k].ptid;
		//printf("Limits x %.2lf  xmin %.2lf xmax %.2lf, y %.2lf  ymin %.2lf ymax %.2lf\n", xc, xmin, xmax, yc, ymin, ymax);
		if(tx > s_arrayh[0].hxlo && tx < s_arrayh[0].hxhi && ty > s_arrayh[0].hylo && ty < s_arrayh[0].hyhi)        /*!< if the init z point is inside the working area */
		{
			indxc = calc_Rindex(tx, s_arrayh[0].hxlo, s_arrayh[0].hresx);
			indyf = calc_Rindex(ty, s_arrayh[0].hylo, s_arrayh[0].hresy);
			//printf("x %lf y %lf xlo %lf resx %lf\n", tx, ty, s_arrayh[0].hxlo, s_arrayh[0].hresx);
			tz = rast_dem[indyf][indxc];
			if (tz != dem_nulval)
			{
				/**
				* cells located in border won't be considered here
				*/
				if( (indyf > 2 && indyf < s_arrayh[0].hn_fy-2) && (indxc > 2 && indxc < s_arrayh[0].hn_cx-2) ) 
				{
					printf("Init Z point %3i de %i  fily= %i colx= %i  tz=%.2lf\n", k, n_xyz, indyf, indxc, tz);
					
					/**
					* Single flow in modes 1 (LHM) or 2 (SSM) to calculate paht
					* Force interaction return different result
					*/				
					if(flow_mode == 1 || flow_mode == 2)
					{
						printf("Iniciate Single flow path calculation\n");
						calc_singflow(k, dir_out, namfile, namxyven, rast_dem, s_arrayh, \
							dem_nulval, indyf, indxc, ptid, tx, ty, tz,\
							flow_mode, max_dist, cri_heig, incre_heig, num_repit, for_inter, dist_type);
					}
					/**
					* Multiflow is based in the drunk sailor montecarlo model
					* Force interaction is not a big deal yet
					*/
					if(flow_mode == 3) 
					{
						printf("Iniciate Multiflow Drunk Sailor path calculation\n");
						calc_drunksailflow(rast_dem, s_arrayh, dem_nulval, \
							indyf, indxc, tx, ty, tz, max_dist, num_itera, \
							num_repit, radiu_val, cri_heig, incre_heig, for_inter, dist_type);
					}
					/**
					* Multipath, find all possible path using a FIFO
					* Force interaction is not a big deal yet
					*/
					if(flow_mode == 4) 
					{
						calc_multyflow(rast_dem, s_arrayh, dem_nulval, \
						indyf, indxc, tx, ty, tz, max_dist, rest_heig, dist_type);
					}
					count_ok++;
				}
			}
			else printf("ZLP with null z value. Check and coordinates DEM\n");			
		}
		else printf("Point outside the DEM: %lf %lf\n", tx,  ty);
	}
	/**
	* If not init xy processed
	*/
	if (count_ok == 0)
	{
		printf("Attention: No Init xy processed, something went wrong\n");
		exit(0);
	}
	/**
	* Detect max value
	*/
	max=0;
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++)
		{
			if (raster_path[i][j] > max) max = raster_path[i][j];
		}
	}
	/**
	* Normalize
	*/
	for(i=0;i<s_arrayh[0].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++)
		{
			if (raster_path[i][j] != 0) raster_path[i][j] = (double)raster_path[i][j] / (double)max;
		}
	}
	/**
	* Change max value in zmax
	*/
	s_arrayh[0].hzhi = max;
	/**
	* Write dem raster for testing
	*/
	if (dem_outfil == 1)
	{
		suffix = (char *)malloc(strlen("_test.grd")+1);
		sprintf(suffix,"_test.grd");
		char *file_outdem = get_pathnam(dir_out, namfile, suffix, 1);
		write_grdrasterF(file_outdem, rast_dem, s_arrayh, 0, 1);
	}
	/**
	* Write raster path
	*/
	suffix = (char *)malloc(2+2+strlen(namxyven)+20); 
	sprintf(suffix,"_pathflow_mod-%i-%i_%s.grd", flow_mode, for_inter, namxyven);
	char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
	write_grdrasterF(file_outpath, raster_path, s_arrayh, 0, 1);
	/**
	* Free raster memory
	*/
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);
	free_RasterF("Rast flow", raster_path, s_arrayh[0].hn_fy);
	free_RasterI("Rast control", raster_control, s_arrayh[0].hn_fy);
	free(suffix);
}

/*! Calc differential of topography between a given point and all cells in a raster following
* direc_calc 1: horizontal
* direc_calc 2: vertical
* direc_calc 3: oblique
* from_center 1: to get differences for Zero Level Trajectory to the maximum extension
* If maximum distance is 0, all cells in DEM will be processed from ZLP. If >0 calculation will stop when distance in reached 
* from_limask 1: Mask is needed. Two output raster in mask, and from mask limits
* Topohazard calculate first in one direction and then the opposite
* Important: if mask is used, the amount of RAM will increase by three times
 * */
void mod_topohazard(char *file_xybed, const char *file_demn, const char *dir_out, const char *namfile)
{
char *suffix;
int k, i, j, num, count_ok;
int idxfy, idxcx, nrio, brio;
//float dx;
double tx, ty, tz;
double max, min, traddx, traddy;
double min_in,min_out,max_in,max_out;
	count_ok=0;                     /*!< count processed ZLPs */
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);

	num = 1;
	if (mas_ifused == 1)
	{
		rast_mask = read_grdrasterI(file_mas, mas_nulval, mas_minval, mas_maxval, s_arrayh, mas_typras, 1);
		num = 2;
	}
	/*!< Due to max and min values of output raster, new head is needed */
	s_arrayh[num].hn_cx      = s_arrayh[0].hn_cx;
	s_arrayh[num].hn_fy      = s_arrayh[0].hn_fy;
	s_arrayh[num].hxlo       = s_arrayh[0].hxlo;
	s_arrayh[num].hxhi       = s_arrayh[0].hxhi;
	s_arrayh[num].hylo       = s_arrayh[0].hylo;
	s_arrayh[num].hyhi       = s_arrayh[0].hyhi;
	s_arrayh[num].hzlo       = 0;
	s_arrayh[num].hzhi       = 0;                                      
	s_arrayh[num].hresx      = s_arrayh[0].hresx ;
	s_arrayh[num].hresy      = s_arrayh[0].hresy;
	/*!< Read xy vent file to calculate gravitational flows */
	read_zerolevpoint(file_xybed, rast_dem, s_arrayh, dem_nulval);
	/*!< Fix z val from DEM and bedriverpt */
	if (xyz_zfix == 1) rast_dem = fix_zval(rast_dem, s_arrayh, dem_nulval, s_zerolevpt, n_zlps, xyz_zdist);
	/*!< Calc buffer connections in ZLP */
	calc_buff(s_zerolevpt, n_zlps);
	/*!< Calc slope of each ZLP */
	calc_bedslope(s_arrayh, s_zerolevpt, n_zlps);
	/*!< For each ZLP */
	for(k=0;k<n_zlts;k++)                                                /*!< for each ZLT calculate tophazard */
	{
		printf("Creating base raster\n"); 
		raster_topohaz  = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
		
		if (mas_ifused == 1)
		{
			inhaz_topo  = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
			outhaz_topo = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
		}
		for(i=0;i<s_arrayh[0].hn_fy;i++) 
		{
			for(j=0;j<s_arrayh[0].hn_cx;j++) 
			{
				raster_topohaz[i][j] = -9999;                 /*!< all values of the working raster will be initiated to -9999 */
				if (mas_ifused == 1)
				{
					inhaz_topo[i][j]  = -9999;
					outhaz_topo[i][j] = -9999;
				}
			}
		}
		printf("Calculating topohazard in ZLT %i\n", k);
		nrio   = s_trajectory[k].rnio;
		traddx = s_trajectory[k].raddx;
		traddy = s_trajectory[k].raddy;
		/*!< Modes 1, 2 and 3 start calculation from ZLP and corresponding cell
		 * Mode 4 start calculation from any given cell in raster, and check nearest
		 * from all available ZLP in a given ZLT */
		if (direc_calc < 4)
		{
			/*!< For each ZLPs in trajectory */
			for(i=0;i<n_zlps;i++)
			{
				/*!< getting ZLP values */
				brio  = s_zerolevpt[i].brio;
				tx    = s_zerolevpt[i].bxcoor;                                   
				ty    = s_zerolevpt[i].bycoor;
				tz    = s_zerolevpt[i].bzcoor;
				idxfy = s_zerolevpt[i].bidxfy;
				idxcx =	s_zerolevpt[i].bidxcx;
				//dx    = s_zerolevpt[i].bdx;
				//printf("z %lf and rio %i - %i\n", tz, nrio, brio);
				if (tz != dem_nulval && nrio == brio)                               /*!<  if ZLP belongs trajectory */
				{
					if (idxfy > 1 && idxfy < s_arrayh[0].hn_fy-1 && idxcx > 1 && idxcx < s_arrayh[0].hn_cx-2)  /*!<  if in raster */
					{
						if(from_center == 1)
						{
							printf("Calculating topohazard in river %i with direction %i and pt riverbed %i in row %i and col %i\n", k, direc_calc, i, idxfy, idxcx);		
							/*!< Negative search */
							central_topohazard(rast_dem, s_arrayh, dem_nulval, \
								idxfy, idxcx, tx, ty, tz, s_zerolevpt, n_zlps, \
								nrio, traddx, traddy, direc_calc, 1, topmax_dist);
							/*!< Positive search */
							central_topohazard(rast_dem, s_arrayh, dem_nulval, \
								idxfy, idxcx, tx, ty, tz, s_zerolevpt, n_zlps, \
								nrio, traddx, traddy, direc_calc, 2, topmax_dist);
						}
						
						if (mas_ifused == 1)
						{
							printf("Mask---------\n");
							printf("Calculating topohazard with mask in %i with direction %i and pt riverbed %i in row %i and col %i\n", k, direc_calc, i, idxfy, idxcx);	
							/*!< Negative search */
							mask_topohazard(rast_dem, rast_mask, s_arrayh, dem_nulval, mas_nulval, \
								idxfy, idxcx, tx, ty, tz, s_zerolevpt, n_zlps, \
								nrio, traddx, traddy, direc_calc, 1, topmax_dist);
							/*!< Positive search */
							mask_topohazard(rast_dem, rast_mask, s_arrayh, dem_nulval, mas_nulval, \
								idxfy, idxcx, tx, ty, tz, s_zerolevpt, n_zlps, \
								nrio, traddx, traddy, direc_calc, 2, topmax_dist);
						}
						count_ok++;
					}
				}
			}
			/**
			* If not ZLP processed
			*/
			if (count_ok == 0)
			{
				printf("Attention: No ZLP processed, something went wrong\n");
				exit(0);
			}
		}
		if (direc_calc == 4) 
		{
			getime(); 
			nearest_topohazard(rast_dem, s_arrayh, dem_nulval, s_zerolevpt, n_zlps, nrio);
			getime(); 
		}
		//**************************************************************
		//**************************************************************
		/**
		* Detect max value
		*/
		printf("Find max value\n");
		max=min=0;
		min_in=min_out=max_in=max_out;
		for(i=0;i<s_arrayh[0].hn_fy;i++) 
		{
			for(j=0;j<s_arrayh[0].hn_cx;j++)
			{
				if (raster_topohaz[i][j] != dem_nulval &&  raster_topohaz[i][j] < min) 
				{
					min = raster_topohaz[i][j];
				}
				if (raster_topohaz[i][j] > max) 
				{
					max = raster_topohaz[i][j];
				}
				if (mas_ifused == 1)
				{
					if (inhaz_topo[i][j] != dem_nulval &&  inhaz_topo[i][j] < min_in) min_in = inhaz_topo[i][j];
					if (inhaz_topo[i][j] > max_in)                                    max_in = inhaz_topo[i][j];
					if (outhaz_topo[i][j] != dem_nulval &&  outhaz_topo[i][j] < min_out) min_out = outhaz_topo[i][j];
					if (outhaz_topo[i][j] > max_out)                                     max_out = outhaz_topo[i][j];
				}
			}
		}
		s_arrayh[num].hzlo = min;
		s_arrayh[num].hzhi = max;
		/**
		* Write dem raster
		*/
		if (dem_outfil == 1)
		{
			suffix = (char *)malloc(strlen("_test.grd")+1);
			sprintf(suffix,"_test.grd");
			char *file_outdem = get_pathnam(dir_out, namfile, suffix, 1);
			write_grdrasterF(file_outdem, rast_dem, s_arrayh, 0, 1);
		}
		/**
		* Write topohazard raster
		*/
		printf("Write topohazard raster\n");
		suffix = (char *)malloc(2+2+strlen(namxybed)+25); 
		sprintf(suffix,"_topocenter-%i_nrio-%i_%s.grd", direc_calc, nrio, namxybed);
		printf("suffix %s -- %li)\n", suffix, strlen(namxybed));
		char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
		write_grdrasterF(file_outpath, raster_topohaz, s_arrayh, num, 1);
		/**
		* Free topohazard raster
		*/
		free_RasterF("Rast topocenter", raster_topohaz, s_arrayh[0].hn_fy);
		free(suffix);
		
		if (mas_ifused == 1)
		{
			s_arrayh[num].hzlo = min_in;
			s_arrayh[num].hzhi = max_in;
			/**
			* Write topohazard raster
			*/
			printf("Write topohazard int raster\n");
			suffix = (char *)malloc(2+2+strlen(namxybed)+30); 
			sprintf(suffix,"_topomaskcenter-%i_nrio-%i_%s.grd", direc_calc, nrio, namxybed);
			file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
			write_grdrasterF(file_outpath, inhaz_topo, s_arrayh, num, 1);
			/**
			* Free topohazard raster
			*/
			free_RasterF("Rast topomaskcenter", inhaz_topo, s_arrayh[0].hn_fy);
			free(suffix);
			
			s_arrayh[num].hzlo = min_out;
			s_arrayh[num].hzhi = max_out;
			/**
			* Write topohazard raster
			*/
			printf("Write topohazard int raster\n");
			suffix = (char *)malloc(2+2+strlen(namxybed)+30); 
			sprintf(suffix,"_topomasklimit-%i_nrio-%i_%s.grd", direc_calc, nrio, namxybed);
			file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
			write_grdrasterF(file_outpath, outhaz_topo, s_arrayh, num, 1);
			/**
			* Free topohazard raster
			*/
			free_RasterF("Rast topomasklimit", outhaz_topo, s_arrayh[0].hn_fy);
			free(suffix);
		}
	}
	/**
	* Write riverbed xy file
	*/
	char *file_outbed = get_pathnam(dir_out, namxybed, ".csv", 1);
	write_zltopohaz(file_outbed, s_zerolevpt, n_zlps);
	free_RasterF("Rast DEM", rast_dem, s_arrayh[0].hn_fy);		
	
}

/*!Calc cross section or Schilling equation implementation 
 * multi_mode 1: cross section 
 * multi_mode 2: Schilling equations 
 * */
void mod_multiflow(char *file_xybed, const char *file_demn, const char *dir_out, const char *namfile)
{
FILE *file;
char *suffix, *file_cross, *file_bedriver;
int k, i, j, num, jumpt;
int idxfy, idxcx;
int idpt, njerar, nrio, brio, ntramo, tend, triodir;
float ndist, tdx, distance;
double tx, ty, tz, max,min, tslopep;
	
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	num = 1;
	if (mas_ifused == 1)
	{
		rast_mask = read_grdrasterI(file_mas, mas_nulval, mas_minval, mas_maxval, s_arrayh, mas_typras, 2);
		num = 2;
	}
	s_arrayh[num].hn_cx      = s_arrayh[0].hn_cx;
	s_arrayh[num].hn_fy      = s_arrayh[0].hn_fy;
	s_arrayh[num].hxlo       = s_arrayh[0].hxlo;
	s_arrayh[num].hxhi       = s_arrayh[0].hxhi;
	s_arrayh[num].hylo       = s_arrayh[0].hylo;
	s_arrayh[num].hyhi       = s_arrayh[0].hyhi;
	s_arrayh[num].hzlo       = 0;
	s_arrayh[num].hzhi       = 0;                                      
	s_arrayh[num].hresx      = s_arrayh[0].hresx ;
	s_arrayh[num].hresy      = s_arrayh[0].hresy;
	/*!< Read xy vent file to calculate gravitational flows */
	read_zerolevpoint(file_xybed, rast_dem, s_arrayh, dem_nulval);
	/*!< Fix z val from DEM and bedriverpt */
	if (xyz_zfix == 1) rast_dem = fix_zval(rast_dem, s_arrayh, dem_nulval, s_zerolevpt, n_zlps, xyz_zdist);
	/*!< Calc buffer connection of each riverbed pt */
	calc_buff(s_zerolevpt, n_zlps);
	/*!< Calc slope of each riverbed pt */
	calc_bedslope(s_arrayh, s_zerolevpt,n_zlps);
	/*!< For each ZLP */
	for(k=0;k<n_zlts;k++)  
	{
		printf("Creating base raster\n");
		rast_multi  = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
		for(i=0;i<s_arrayh[0].hn_fy;i++) 
		{
			for(j=0;j<s_arrayh[0].hn_cx;j++) rast_multi[i][j] = -9999;
		}
		nrio   = s_trajectory[k].rnio;		
		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		/*!< Create and write first line to save cross section data */
		if (multi_mode == 1)    
		{
			printf("Write cross section file\n");
			suffix = (char *)malloc(2+2+strlen(namxybed)+40); 
			sprintf(suffix,"_multiflow-%i_nrio-%i_crosection_%s.csv", multi_mode, nrio, namxybed);
			file_cross = get_pathnam(dir_out, namfile, suffix, 1);
			if((file = fopen(file_cross, "w"))== NULL)
			{
				printf("-------ERROR open file--------\n");
				exit(0);
			}
			else
			{
				fprintf(file,"%s\n",
				"njear nrios ntramo ndist xcoor ycoor zcoor idx_fy idx_cx horzdist difftopo"); //primera linea
			}
			fclose(file);
			free(suffix);
		}
		if (multi_mode == 2)    
		{
			printf("Write laharZ cross section file\n");
			suffix = (char *)malloc(2+2+strlen(namxybed)+40); 
			sprintf(suffix,"_multiflow-%i-%i_nrio-%i_crosection_%s.csv", multi_mode, (int)volum_lahar, nrio, namxybed);
			file_cross = get_pathnam(dir_out, namfile, suffix, 1);
			if((file = fopen(file_cross, "w"))== NULL)
			{
				printf("-------ERROR open file--------\n");
				exit(0);
			}
			else
			{
				fprintf(file,"%s\n",
				"idpt nrios xcoor ycoor zcoor horzdist difftopo depth"); //primera linea
			}
			fclose(file);
			free(suffix);
		}
		//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		if (cros_jump > 0) jumpt=cros_jump;
		jumpt=1;
		/*!< For each ZLP in trajectory */
		for(i=0;i<n_zlps;i++)                                           
		{
			/*!< getting ZLP values */
			idpt   = s_zerolevpt[i].bidpt;
			njerar = s_zerolevpt[i].bjerar,
			brio   = s_zerolevpt[i].brio,
			ntramo = s_zerolevpt[i].btramo,
			ndist  = s_zerolevpt[i].bdist,
			tx     = s_zerolevpt[i].bxcoor;                                   
			ty     = s_zerolevpt[i].bycoor;
			tz     = s_zerolevpt[i].bzcoor;
			idxfy  = s_zerolevpt[i].bidxfy;
			idxcx  = s_zerolevpt[i].bidxcx;
			tdx    = s_zerolevpt[i].bdx;
			tend   = s_zerolevpt[i].bends;
			triodir= s_zerolevpt[i].bstreamdir;
			tslopep= s_zerolevpt[i].bslopepor;
			/*!<  if ZLP belongs trajectory */
			if (tz != -9999 && nrio == brio)
			{
				if (idxfy > 1 && idxfy < s_arrayh[0].hn_fy-1 && idxcx > 1 && idxcx < s_arrayh[0].hn_cx-2)
				{
					if (i > 0)
					{
						/*!< Calc cross section file */
						if (multi_mode == 1 && jumpt == i)                            
						{
							
							calc_crossfile(file_cross, rast_dem, s_arrayh, dem_nulval, \
								s_zerolevpt, n_zlps, i, cros_maxwith, cros_height, \
								njerar, brio, ntramo, ndist, tx, ty, tz, idxfy, idxcx, tdx, tend, triodir);
							if (cros_jump > 0) jumpt += cros_jump;
							else jumpt=i;
						}
					}
					/*!< Calc Schilling equation*/
					if (multi_mode == 2 || multi_mode == 3)                                
					{
						//distkm = (float)ndist/1000;                     
						//ttime  = a0 + a1*distkm + a2*distkm*distkm;     /*!< Arrival time laharZ */
						//s_zerolevpt[i].bttime = ttime;
						printf("Calculating laharz in river %i and pt riverbed %i in row %i and col %i -- slop %lf\n", k, i, idxfy, idxcx, tslopep);
						calc_laharz(file_cross, rast_dem, s_arrayh, dem_nulval, \
							s_zerolevpt, n_zlps, simple_flow, cros_maxwith, cros_heincre, area_section, maning_coef, a0, a1, a2, \
							i, idpt, brio, tx, ty, tz, tslopep, idxfy, idxcx, tdx, tend, triodir);
							
						/*!< cummulative planimetric area, calculated in two times */
						if (i > 0) 
						{
							s_zerolevpt[i].bareaplanim = s_zerolevpt[i-1].bareaplanim + s_zerolevpt[i].bareaplanim;
							s_zerolevpt[i].btimearrive = s_zerolevpt[i-1].btimearrive + s_zerolevpt[i].bttime;
							
							if (s_zerolevpt[i].bareaplanim < plani_lahar) s_zerolevpt[i].binplanim = 1;
						}
						if (maxdist_flow > 0)
						{
							distance = calc_dist(tx, s_zerolevpt[0].bxcoor, ty, s_zerolevpt[0].bycoor);
							if (distance > maxdist_flow) break;
						}
						
						
					}
					//if (multi_mode == 3)                                /*!< Calc debris-flow */
					
					//if (multi_mode == 4)                                /*!< Calc flow */
					
				}
			}
		}
		if (multi_mode == 2 || multi_mode == 3)
		{
			/**
			* Write bedriver csv
			*/
			suffix = (char *)malloc(2+2+strlen(namxybed)+40); 
			sprintf(suffix,"_multiflow-%i_nrio-%i_%s.csv", multi_mode, nrio, namxybed);
			file_bedriver = get_pathnam(dir_out, namfile, suffix, 1);
			write_zltmulti(file_bedriver, s_zerolevpt, n_zlps, nrio);
			free(suffix);
			/**
			* Detect max value
			*/
			printf("Find max value\n");
			max=min=0;
			for(i=0;i<s_arrayh[0].hn_fy;i++) 
			{
				for(j=0;j<s_arrayh[0].hn_cx;j++)
				{
					if (rast_multi[i][j] != dem_nulval && rast_multi[i][j] < min) 
					{
						min = rast_multi[i][j];
					}
					if (rast_multi[i][j] > max) 
					{
						max = rast_multi[i][j];
					}
				}
			}
			s_arrayh[num].hzlo = min;
			s_arrayh[num].hzhi = max;
			/**
			* Write multiflow raster
			*/
			suffix = (char *)malloc(2+2+strlen(namxybed)+40); 
			printf("Write topohazard raster\n");
			if (multi_mode == 1) sprintf(suffix,"_multiflow-%i_nrio-%i_%s.grd", multi_mode, nrio, namxybed);
			if (multi_mode == 2)sprintf(suffix,"_multiflow-%i-%i_nrio-%i_%s.grd", multi_mode, (int)volum_lahar, nrio, namxybed);
			
			char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
			write_grdrasterF(file_outpath, rast_multi, s_arrayh, num, 1);
			/**
			* Free multiflow raster
			*/
			free_RasterF("Rast topocenter", rast_multi, s_arrayh[0].hn_fy);
			free(suffix);
		}
	}
	
}

//**********************************************************************
/*! Interpolate functions */
//**********************************************************************
/*! Interpolate points in a raster using Inverse Distance Weighting (IDW)
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
 * */
void mod_interpolate(int n_strategy, char *file_xyz, const char *file_demn, const char *dir_out, char *namfile)
{
char *suffix, *file_outpath;
	printf("Reading input filesf\n");
	
	/*!< Read original DEM file */
	if (dem_cliped != 1) rast_dem = read_grdrasterF(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0);
	else rast_dem = read_grdrasterF_clip(file_demn, dem_nulval, dem_minval, dem_maxval, s_arrayh, dem_typras, 0, dem_nwxlow, dem_nwxhig, dem_nwylow, dem_nwyhig);
	/*!< Read original DEM file */
	read_xyz(file_xyz, 3);
	
	printf("Calculating interpolation method %i with distance %f\n", n_strategy, idw_distan);
	
	if (inter_mode == 1)
	{
		if (idw_approach == 1)
		{
			/*!< By closet points, not ready yet */
			calc_idwbypts(s_arrayh, dem_nulval, s_coord, n_xyz, idw_searpts, idw_power);
		}
		if (idw_approach == 2)
		{
			/*!< Only by distance, min searching pts, the value remain */
			calc_idwbydist(s_arrayh, dem_nulval, s_coord, n_xyz, idw_distan, idw_searpts, idw_power);
		}
		if (idw_approach == 3)
		{
			/*!< By distance, increase searching distance if not enough points are found */
			calc_idwbydistup(s_arrayh, dem_nulval, s_coord, n_xyz, idw_distan, idw_searpts, idw_power, idw_bytype);
		}
	}
	/**
	* Write interpola raster
	*/
	if (nsuffix_short == 1)
	{
		suffix = (char *)malloc(strlen("_intidw_.grd")+strlen(namxyz)+10);
		sprintf(suffix,"_intidw_%s.grd", namxyz);
	}
	else
	{
		suffix = (char *)malloc(strlen("_intidw_XXX_XXX_.grd")+strlen(namxyz)+10);
		sprintf(suffix,"_intidw_%.lf_%.lf_%s.grd", idw_power, idw_distan, namxyz);
	}
	
	file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
	write_outputinter(rast_inter,s_arrayh, dem_nulval, 0, suffix, 1);
	//write_grdrasterF(file_outpath, rast_inter, s_arrayh, 0, 1);
	free(suffix);
	free(file_outpath);
	if (s_modesmooth > 0)
	{
		printf("--------------------------\n\n");
		smooth_funtion(n_strategy, dir_out, namxyz, idw_power, idw_distan, rast_inter, s_arrayh, 0, dem_nulval, s_modesmooth, s_maxpercent, s_minpercent, s_totoutcels);
		
		/*
		raster_smooth = soft_rasterF(rast_inter, s_arrayh, 0, dem_nulval, r_softmod, r_softmax, r_softmin, r_softcel);
		printf("End smoothing\n");
		
		if (nsuffix_short == 1)
		{
			suffix = (char *)malloc(strlen(namxyz)+strlen("_intidw-soft_.grd")+10);
			sprintf(suffix,"_intidw-soft_%s.grd", namxyz);
		}
		else
		{
			suffix = (char *)malloc(strlen(namxyz)+strlen("_intidw-soft_XXX_XXX_XX.grd")+10);
			sprintf(suffix,"_intidw-soft_%.lf_%.lf_%i_%s.grd", idw_power, idw_distan, r_softmod, namxyz);
		}
		write_output(raster_smooth, s_arrayh, dem_nulval, 0, suffix, 1);
		*/
	}
	free_RasterF("Rast INTERPOLA", rast_inter, s_arrayh[0].hn_fy);
}


/*! Smooth raster output depending on strategy */
void smooth_funtion(int n_strategy, const char *dir_out, const char *namefile, float name_param1, \
	float name_param2, double **raster, struct HeadR s_arrayh[], int r_indx, float dem_nulval, \
	int mode_soft, float r_maxpercent, float r_minpercent, float r_totoutcells)
{
char *suffix;
double **raster_smooth;
	printf("--------------------------\n\n");
	raster_smooth = smooth_rasterF(raster, s_arrayh, r_indx, dem_nulval, mode_soft, r_maxpercent, r_minpercent, r_totoutcells);
	printf("End smoothing\n");
	
	if (nsuffix_short == 1)
	{
		suffix = (char *)malloc(strlen(namefile)+strlen("strat-X_soft_.grd")+10);
		sprintf(suffix,"strat-%i_soft_%s.grd", n_strategy, namefile);
	}
	else
	{
		suffix = (char *)malloc(strlen(namxyz)+strlen("_soft_XXX_XXX_XX.grd")+10);
		sprintf(suffix,"strat-%i_soft_%.lf_%.lf_%i_%s.grd", n_strategy, name_param1, name_param2, mode_soft, namxyz);
	}
	char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
	write_outputinter(raster_smooth, s_arrayh, dem_nulval, 0, file_outpath, 1);
	free(file_outpath);
	
}




/*
void write_output(double **raster, struct HeadR s_arrayh[], float dem_nulval, int r_index, char *out_namfile, int typmeasure)
{
int i, j, celltot;
double max,sumtot,volum, value;

	max=volum=celltot=sumtot=0;
	for(i=0;i<s_arrayh[r_index].hn_fy;i++) 
	{
		for(j=0;j<s_arrayh[r_index].hn_cx;j++)
		{
			if (raster[i][j] > max) max = raster[i][j];
			
			if (raster[i][j] != dem_nulval)
			{
				sumtot+=raster[i][j];
				if (typmeasure == 1) value = raster[i][j]/1000;  
				if (typmeasure == 2) value = raster[i][j]/100;   
				if (typmeasure == 3) value = raster[i][j]/10;    
				if (typmeasure == 4) value = raster[i][j]*10;    
				if (typmeasure == 5) value = raster[i][j]*100;   
				if (typmeasure == 6) value = raster[i][j]*1000; 
				volum+= value * s_arrayh[r_index].hresx * s_arrayh[r_index].hresx;
				celltot++;
			}
		}
	}
	printf ("suma total = %lf\n", sumtot);
	printf ("celdas validas = %i\n", celltot);
	printf ("volume in cubic meters = %lf\n", volum);

	s_arrayh[r_index].hzhi = max;

	//char *file_outpath = get_pathnam(dir_out, namfile, suffix, 1);
	write_grdrasterF(out_namfile, raster, s_arrayh, 0, 1);
	//free(file_outpath);
}
*/


//**********************************************************************
/*! MAIN */
//**********************************************************************

int main(int argn, char **args)
{
int call, es_med;
char home[MAX_PATH];
	/**
	* Init random feed
	*/
	srand(time(NULL));
	es_med=call=0;
	printf("Initiating MDTanaliza\n");
	getime(); 
	printf("Check working directory\n");
	strcpy(home, getenv("HOME"));
	if (home != NULL) 
	{
		if (getcwd(cwd, sizeof(cwd)) != NULL)                           /*!< get working directory */
		{
			es_med = es_enmedia(cwd);                                   /*!< check  if external disk */
			if (es_med == 1)
			{
				strcpy(home, new_home(cwd));                            /*!< sift string home */
				printf("New home is %s\n", home);
			}
		} 
	}
	printf("Define output outcome based on argument input\n");
	if (is_integer(args[1])) 
	{
		printf("Argument is a number\n");
		int arg_num = atoi(args[1]);
		if (arg_num < 0 || arg_num > 9)
		{
			printf("Error: Argument to generate configuration file must be 1-9\n");
			exit(0);
		}
		else
		{
			if (arg_num == 8)
			{
				printf("Error: strategy 8 has not been implemented yet. Try again, please\n"); 
				exit(0);
			}
			else
			{
				printf("Generating configuration file as \n");
				print_config_file(arg_num);
			}
		}
	}
	else
	{
		printf("Argument is a configuration file\n");
		if (args[1]) 
		{
			if (es_med == 0) printf("Working directory in: %s\n", home);
			if (es_med == 1) printf("Alternative working directory in: %s\n", home);
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
			call = read_cfg(args[1], home);
			printf("Check main cfg principal, if return value different than 1: %i\n\n", call);
			if (call == 1)
			{
				printf("Initiating calculation in %i\n", n_strategy);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				if (n_strategy == 1) mod_dem(file_xyz, file_dem, dir_out, namfile);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//SINK mode
				if (n_strategy == 2) mod_sinks(file_dem);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//ASPECT  mode
				if (n_strategy == 3) mod_aspect(file_dem, dir_out, namfile);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//SLOPE mode
				if (n_strategy == 4) mod_slope(file_dem, dir_out, namfile);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//Gravitational Flow path mode
				if (n_strategy == 5)mod_gravflows(file_xyven, file_dem, dir_out, namfile, flow_mode, max_dist);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//TOPOHAZARD mode
				if (n_strategy == 6)mod_topohazard(file_xybed, file_dem, dir_out, namfile);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//MULTIFLOW mode
				if (n_strategy == 7)mod_multiflow(file_xybed, file_dem, dir_out, namfile);
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//HAZARD mode (not implemented yet)
				//if (n_strategy == 8)
				//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				//INTERPOLATION mode
				if (n_strategy == 9)
				{
					//IDW
					if (inter_mode == 1) mod_interpolate(n_strategy, file_xyz, file_dem, dir_out, namfile);
				}
			}
		}
		else
		{
			printf("ERROR: insert configuration path and filename\n");
			exit(0);
		}
	}
	
}




/*!
 * PENDING
 * New option in strategy 1 using a mask instead xy file only.
 * Normalized yes or no in output rasters
 * 
 * UPDATES
 * 2025/09/20: Fixed gravity flow in mode 3. Error in random function implementation and search vents
 * 2025/09/12: random function (general lib) in gravity flow in mode 3 did not work well.
 * 2025/08/11: Implementation clip raster option in all strategies, but modify DEM
 * 2025/08/10: Implementation mode 4 in Topohazards (nearest). Not finish
 * 	Results are not quite good yet
 */





