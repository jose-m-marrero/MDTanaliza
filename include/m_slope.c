/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza
* Module:  m_slope
* 
* Author: Jose M. Marrero 
*          
* Version: 0.1.1
* Creation Date: 2018-04-06
* Last update: 2022-09-30
* 
* Description:
* Calc slope in degrees, classes and percentages from original DEM 
* (1) slope_mode Using the SSM Burrough and McDonell (1998)
* (2) slope_mode Second-order finite difference 2FD Fleming and Hoffer (1979)
* (3) slope_mode Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Distance  3FDWRD Unwin (1981)
* (4) slope_mode Three-order  Finite  Difference,  Linear  regression  plan  3FD Sharpnack et al (1969)
* (5) slope_mode Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Squared  Distance 3FDWRSD Horn (1981)
* (6) slope_mode Frame Finite difference FFD Chu and Tsai (1995)
* (7) slope_mode Maximum Max Travis et al. (1975)
* (8) slope_mode Simple difference Simple-D Jones (1998)
* (9) slope_mode Constrained Quadratic Surface Quadsurface (Wood Evans [13])
* 
* *********************************************************************/
#ifndef _SLP_RAS
#define _SLP_RAS


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

//**********************************************************************
/*! Calc Functions */
//**********************************************************************

/*! SLOPE-GRADIENT CALCULATION (write on rast1 rast2 rast3) 
 * Mode:
 * 1 Using the SSM  
 * 2 Second-order finite difference 2FD 
 * 3 Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Distance  3FDWRD
 * 4 Three-order  Finite  Difference,  Linear  regression  plan  3FD
 * 5 Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Squared  Distance 3FDWRSD
 * 6 Frame Finite difference FFD
 * 7 Maximum Max
 * 8 Simple difference Simple-D
 * 9 Constrained Quadratic Surface Quadsurface (Wood Evans [13])
 * */
void calc_slope(const char *dir_out, const char *namfile, double **raster, struct HeadR s_arrayh[], float dem_nulval, int slope_mode)
{
FILE *file;
char subfix[50];
char *file_outrast, *file_histo;
int i, j, k, degslo[19];
double topo[8], degpor[19];
double fx, fy, fx2, fy2, rise;
double slprad, slpdeg, superf;
double inival, finx, finz;
double *cellsval;
double c1,c2,c3,c4,c5,c6,c7,c8,c9;
double **rast_degree, **rast_porcen, **rast_class;

	rast_degree = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	rast_porcen = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	rast_class  = Crea_2DFarray(s_arrayh[0].hn_fy, s_arrayh[0].hn_cx);
	
	printf("Calculating slope-gradient method %i\n", slope_mode);
	/**
	* Reset to 0 to calculate the histogram
	*/
	for(i=0;i<19;i++)
	{
		degslo[i] = 0;
		degpor[i] = 0;
	}
	/**
	* Start slope calculation
	*/	
	for(i=0;i<s_arrayh[0].hn_fy;i++) //column
	{
		for(j=0;j<s_arrayh[0].hn_cx;j++) //row
		{
			/*!< if array index are in working area */
			if((i>3 && i<s_arrayh[0].hn_fy-3 ) && ( j>3 && j<s_arrayh[0].hn_cx-3 ))               
			{
				c5  = raster[i][j];                                     /*!< Get center cell z value  */
				if(c5 != dem_nulval)                                    /*!< if z coordinate value is not null */
				{                 
					//getmovingcell(i, j, 0);                                /*!< get 3x3 moving cell z values */
					cellsval = search_celproxF(raster, i, j);
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
					* Using the SSM 
					* http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q900000023000000
					* Burrough and McDonell (1998)
					* dzdx = ((c + 2f + i) - (a + 2d + g)) / 8 * g
					* dzdy = ((g + 2h + i) - (a + 2b + c)) / 8 * g 
					*/
					if(slope_mode == 1)
					{
						fx = ((c7 + (2*c4) + c1) - (c9 + (2*c6) + c3)) / (8 * s_arrayh[0].hresx);
						fy = ((c3 + (2*c2) + c1) - (c9 + (2*c8) + c7)) / (8 * s_arrayh[0].hresy);
					}
					/**
					* Second-order finite difference 2FD 
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Fleming and Hoffer (1979)
					* fx=(z6-z4)/2g  
					* fy=(z8-z2)/2g  
					*/
					if(slope_mode == 2)
					{	
                        fx = (c6 - c4) / (2*s_arrayh[0].hresx);
                        fy = (c8 - c2) / (2*s_arrayh[0].hresy);
					}
					/**
					* Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Distance  3FDWRD
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Unwin (1981)
					* fx=(z3-z1+√2(z6-z4)+z9-z7)/(4+2√2)g   
					* fy=(z8-z2)/2g 
					*/	
					if(slope_mode == 3)
					{	 
                        fx = ((c3 - c1) + (sqrt(2)*(c6-c4)) + (c9 - c7)) / ((4 + (2*sqrt(2)) )* s_arrayh[0].hresx);
                        fy = (c8 - c2) / (2*s_arrayh[0].hresy);  
                        
					}	
					/**
					* Three-order  Finite  Difference,  Linear  regression  plan  3FD
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Sharpnack et al (1969)
					* fx=(z3-z1+z6-z4+z9-z7)/6*g  
					* fy=(z7-z1+z8-z2+z9-z3)/6*g 
					*/		
					if(slope_mode == 4)
					{		
						fx = (((((c3-c1)+c6)-c4)+c9)-c7) / (6*s_arrayh[0].hresx);
						fy = (((((c7-c1)+c8)-c2)+c9)-c3) / (6*s_arrayh[0].hresy);
					}
					/**
					* Three-order  Finite  Difference  Weighted  by  Reciprocal  of  Squared  Distance 3FDWRSD
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Horn (1981)
					* fx=(z3-z1+2 (z6-z4)+z9-z7/8*g  
					* fy=(z7-z1+2(z8-z2)+z9-z3)/8*g 
					*/
					if(slope_mode == 5)
					{ 
						fx = ((c3-c1)+(2*(c6-c4))+(c9-c7)) / (8*s_arrayh[0].hresx);
						fy = ((c7-c1)+(2*(c8-c2))+(c9-c3)) / (8*s_arrayh[0].hresy);	
					}
					/**
					* Frame Finite difference FFD
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Chu and Tsai (1995)
					* fx=(z3-z1+z9-z7)/4*g 
					* fy=(z7-z1+z9-z3)/4*g
					*/
					if(slope_mode == 6)
					{	
						fx=(c3-c1+c9-c7)/(4*s_arrayh[0].hresx);
						fy=(c7-c1+c9-c3)/(4*s_arrayh[0].hresy);
					}
					/**
					* Maximum Max
					* https://babel.hathitrust.org/cgi/pt?id=umn.31951d029867321;view=1up;seq=17
					* Travis et al. (1975)
					* max(abs((z5-z1)/(√2×g)),abs((z5-z2)/g), abs((z5-z3)/(√2×g)),abs((z5-z9)/(√2×g)), 
					* abs((z5-z7)/(√2×g)),abs((z5-z6)/g), abs((z5-z8)/g), abs((z5-z4)/g))
					*/
					if(slope_mode == 7)
					{	
						topo[0] = fabs((c5-c1)/(sqrt(2)*s_arrayh[0].hresx));
						topo[1] = fabs((c5-c2)/s_arrayh[0].hresx);
						topo[2] = fabs((c5-c3)/(sqrt(2)*s_arrayh[0].hresx));
						topo[3] = fabs((c5-c9)/(sqrt(2)*s_arrayh[0].hresx)); 
						topo[4] = fabs((c5-c7)/(sqrt(2)*s_arrayh[0].hresx));
						topo[5] = fabs((c5-c6)/s_arrayh[0].hresx);
						topo[6] = fabs((c5-c8)/s_arrayh[0].hresx); 
						topo[7] = fabs((c5-c4)/s_arrayh[0].hresx);
						inival = 0.0;
						for(k=0;k<8;k++)
						{	
							if(topo[k] > inival)
							{
								if(k==0||k==2||k==3||k==4)finx = sqrt(2)*s_arrayh[0].hresx;
								if(k==1||k==5||k==6||k==7)finx = s_arrayh[0].hresx;
								if(k==0)finz = fabs(c5-c1);
								if(k==1)finz = fabs(c5-c2);
								if(k==2)finz = fabs(c5-c3);
								if(k==3)finz = fabs(c5-c9);
								if(k==4)finz = fabs(c5-c7);
								if(k==5)finz = fabs(c5-c6);
								if(k==6)finz = fabs(c5-c8);
								if(k==7)finz = fabs(c5-c4);
							}	
						}
						fx = finx; //distance
						fy = finz; //dif height
					}
					/**
					* Simple difference Simple-D
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Jones (1998)
					* fx=(z5-z4)/g
					* fy=(z5-z2)/g
					*/	                                    
					if(slope_mode == 8)
					{	
						fx=(c5-c4)/s_arrayh[0].hresx;
						fy=(c5-c2)/s_arrayh[0].hresy;	
					}
					/**
					* Constrained Quadratic Surface Quadsurface (Wood Evans [13])
					* //https://www.witpress.com/Secure/elibrary/papers/RM11/RM11013FU1.pdf
					* Wood (1996)
					* F(x,y)=ax^2 +by^2 +cxy+dx+ey+f
					* δz/δy =2by + cx + e                                
                    * si x,y = 0
                    * δz/δx = d
                    * δz/δy = e 
                    * D = [(c6 + c4) /2 - c5] / resx^2
					* E = [(c8 + c2) /2 - c5] / resy^2
					*/
					if(slope_mode == 9)
					{	
						fx=(((c6 + c4) /2) - c5) / pow(s_arrayh[0].hresx,2);
						fy=(((c8 + c2) /2) - c5) / pow(s_arrayh[0].hresy,2);
                    }    
                    /**
					* For debugging
					* printf("x %.2lf y %.2lf\n", fx, fy);
					*/	
                    fx2 = pow(fx, 2);
					fy2 = pow(fy, 2);
                    rise = sqrt((fx2 + fy2));
                    if(slope_mode < 7 || slope_mode > 7)
                    {
						slprad = atan (rise);                           /*!< Calculating s-gradient in radians */
					}	
					if(slope_mode == 7)
					{
						slprad = atan (fy / fx);                        /*!< Calculating s-gradient in radians */
					}
					slpdeg = slprad * (180 / PI);                       /*!< Calculating s-gradient in degrees  */ 
					rast_degree[i][j] = slpdeg;                         /*!< Saving s-gradient in degrees  */
					rast_porcen[i][j] = slpdeg * 100 / 45;              /*!< Saving s-gradient in percentages  */
					/**
					* Classifying s-gradient and saving results
					*/
					if(slpdeg < 2)                   rast_class[i][j] = 1;    /*!< Very Slightly leaned/sloped - Muy ligeramente inclinado */
					if(slpdeg >= 2  && slpdeg < 5)   rast_class[i][j] = 2;    /*!< Slightly leaned/sloped - ligeramente inclinado */
					if(slpdeg >= 5  && slpdeg < 10)  rast_class[i][j] = 3;    /*!< leaned/sloped - inclinado */
					if(slpdeg >= 10 && slpdeg < 15)  rast_class[i][j] = 4;    /*!< heavily leaned/sloped - Fuertemente inclinado */
					if(slpdeg >= 15 && slpdeg < 30)  rast_class[i][j] = 5;    /*!< Moderately steep - Moderadamente escarpado */
					if(slpdeg >= 30 && slpdeg < 60)  rast_class[i][j] = 6;    /*!< Steep - Escarpado */
					if(slpdeg >= 60)                 rast_class[i][j] = 7;    /*!< Very steep Muy escarpado */	
					/**
					* Classifying s-gradient histogram frequency
					*/
					if(slpdeg == 0)                  degslo[0]++;    
					if(slpdeg  > 0  && slpdeg < 5)   degslo[1]++;    
					if(slpdeg >= 5  && slpdeg < 10)  degslo[2]++;     
					if(slpdeg >= 10 && slpdeg < 15)  degslo[3]++;     
					if(slpdeg >= 15 && slpdeg < 20)  degslo[4]++;     
					if(slpdeg >= 20 && slpdeg < 25)  degslo[5]++;    
					if(slpdeg >= 25 && slpdeg < 30)  degslo[6]++;     
					if(slpdeg >= 30 && slpdeg < 35)  degslo[7]++;     
					if(slpdeg >= 35 && slpdeg < 40)  degslo[8]++;     
					if(slpdeg >= 40 && slpdeg < 45)  degslo[9]++;     
					if(slpdeg >= 45 && slpdeg < 50)  degslo[10]++;   
					if(slpdeg >= 50 && slpdeg < 55)  degslo[11]++;      
					if(slpdeg >= 55 && slpdeg < 60)  degslo[12]++;      
					if(slpdeg >= 60 && slpdeg < 65)  degslo[13]++;      
					if(slpdeg >= 65 && slpdeg < 70)  degslo[14]++;      
					if(slpdeg >= 70 && slpdeg < 75)  degslo[15]++;      
					if(slpdeg >= 75 && slpdeg < 80)  degslo[16]++;         
					if(slpdeg >= 80 && slpdeg < 85)  degslo[17]++; 
					if(slpdeg >= 85 )                degslo[18]++;	
				}	
				if(c5 == dem_nulval)                                    /*!< if z = null value */
				{
					rast_degree[i][j] = dem_nulval;                     /*!< save null value */
					rast_porcen[i][j] = dem_nulval;
					rast_class[i][j]  = dem_nulval;
					degslo[18]++;
				}
			}
		}
	}
	/**
	* Calling raster writting functions
	*/
	s_arrayh[1].hn_cx      = s_arrayh[0].hn_cx;
	s_arrayh[1].hn_fy      = s_arrayh[0].hn_fy;
	s_arrayh[1].hxlo       = s_arrayh[0].hxlo;
	s_arrayh[1].hxhi       = s_arrayh[0].hxhi;
	s_arrayh[1].hylo       = s_arrayh[0].hylo;
	s_arrayh[1].hyhi       = s_arrayh[0].hyhi;
	s_arrayh[1].hzlo       = 0;
	s_arrayh[1].hzhi       = 90;                                        /*!< Assign z max value for raster head data */
	s_arrayh[1].hresx      = s_arrayh[0].hresx ;
	s_arrayh[1].hresy      = s_arrayh[0].hresy;
	
	sprintf(subfix,"_slope-degrees_mode-%i.grd", slope_mode);
	file_outrast = get_pathnam(dir_out, namfile, subfix, 1);
	write_grdrasterF(file_outrast, rast_degree, s_arrayh, 1, 1);        /*!< Write S-slope degree */
	
	s_arrayh[1].hzhi       = 7;                                         /*!< Assign z max value for raster head data */	
	sprintf(subfix,"_slope-class_mode-%i.grd", slope_mode);
	file_outrast = get_pathnam(dir_out, namfile, subfix, 1);
	write_grdrasterF(file_outrast, rast_class, s_arrayh, 1, 1);         /*!< Write S-slope class */
	
	s_arrayh[1].hzhi       = 100;                                       /*!< Assign z max value for raster head data */	
	sprintf(subfix,"_slope-percen_mode-%i.grd", slope_mode);
	file_outrast = get_pathnam(dir_out, namfile, subfix, 1);
	write_grdrasterF(file_outrast, rast_porcen, s_arrayh, 1, 1);        /*!< Write S-slope percentaje */
	
	/**
	* Write histograms
	*/
	printf("\nWrite histograms - slope-gradient %i\n", slope_mode);
	i = 0;
	sprintf(subfix,"_slope-histogram_mode-%i.csv", slope_mode);
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
		for(i=0;i<19;i++)
		{
			/**
			* Calc normalized values
			*/
			if(degslo[i] == 0)degpor[i]=0;
			if(degslo[i]  > 0)degpor[i]= (double)degslo[i] / (double)s_arrayh[0].htotdemcel;
			superf = (degslo[i] * s_arrayh[0].hresx * s_arrayh[0].hresy)/1000000;					/*!< Calc total surface by s-gradient class in km2 */
			printf("Total cells %i, surface %.2lf km2 with s-grad class %i\n", degslo[i], superf, i);
			/**
			* Write final values
			*/
			fprintf(file,"%i %i %lf %lf\n",
				i,
				degslo[i],
				degpor[i],
				superf
				);
		}	
	} 
	fclose(file);
	free_RasterF("Rast slope degrees", rast_degree, s_arrayh[1].hn_fy);
	free_RasterF("Rast slope class", rast_class, s_arrayh[1].hn_fy);
	free_RasterF("Rast slope percen", rast_porcen, s_arrayh[1].hn_fy);
}



#endif /* _SLP_RAS */
