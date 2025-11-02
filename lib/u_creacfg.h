/*!
* Copyright (C) 2019, 2021  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Programa: Create cfg files for different main programs
* Modulo: Confing files
* 
* Autores: Jose M. Marrero.
* 
* Version: 1.1
* Fecha de creacion: 2023-06-30
* Fecha de ultima modificacion: 2023-06-30
* 
* Nombre aplicacion: u_creacfg.h
* 
* Descripcion:
* 
* Realiza diversos calculos geométricos o aritmétricos
* Evalua relaciones de tamaño entre datos
* 
* *********************************************************************/
#ifndef _CREA_H
#define _CREA_H
/*!
* INCLUDES
***********************************************************************/ 
/*! 
* Calculos geometricos, aritmetricos y comparaciones de datos
***********************************************************************/

/**
* Funciones externas
***********************************************************************/

/*! 
* Librerias estandar
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//**********************************************************************
/*! SEAE method config files */
//**********************************************************************


void seae_maincfg(const char *namcfg, int strate, int mode, float unifam, const char *codenam,
		const char *baserar, const char *basecv, const char *suprar, const char *supcv, 
		const char *dirout, const char *cfgbase, const char *cfgsup)
{
FILE *file;
	if ((file = fopen(namcfg,"wt"))== NULL)
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
	
		fprintf(file,"#****ESTRATEGIA_ANALISIS****\n");
		fprintf(file,"ESTRATEGIA:%i\n",strate);
		fprintf(file,"MODE:%i\n", mode);
		fprintf(file,"UNIFAM:%f\n", unifam);
		fprintf(file,"CODE:%s\n", codenam);
		fprintf(file,"DIR_INBASER:%s\n", baserar);
		fprintf(file,"DIR_INBASEC:%s\n", basecv);
		fprintf(file,"DIR_INSUPR:%s\n", suprar);
		fprintf(file,"DIR_INSUPC:%s\n", supcv);
		fprintf(file,"DIR_OUTGEN:%s\n", dirout);
		fprintf(file,"BAS_INFILES:%s\n", cfgbase);
		fprintf(file,"SUP_INFILES:%s\n", cfgsup);
		fprintf(file,"#****ENDFILE****");
	}
	fclose(file);
}

//**********************************************************************
/*! Sum method config files */
//**********************************************************************

/*! 
 * Create main cfg file for sum
 * mode: type of approach (1) sum (2) Highest (3) substract
 * eqrast: if raster are equal in dimension (1) or different (0)
 * measur: unit used in raster file for z values (0 nothing, 1 mm, 2 cm, 3 dm, 4 m, etc..)
 * 	used to calculate volume in meters when relevant. 
 * outsubix: define outfile name
 * infile: secondary cfg file
 * dirin: path showing where to find raster files described in secondary cfg file
 * dirout: output file path
 * smothmod: if smooth mode is activate
 * smothmin: Min percentaje allowed to detect changes
 * smothmax: Max percetaje allowed to detect changes
 * smothcell: searching 
 */
void sum_maincfg(const char *namcfg, int mode, int eqrast, int measur, const char *outsubix, 
		const char *infile, const char *dirin, const char *dirout, 
		int smothmod, float smothmin, float smothmax, int smothcell)
{
FILE *file;
	printf("Writing main sum cfg file %s\n", namcfg);
	if ((file = fopen(namcfg,"wt"))== NULL)
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
		fprintf(file,"#*****INI********\n");
		fprintf(file,"R_MODE:%i\n", mode);
		fprintf(file,"EQ_RAST:%i\n", eqrast);
		fprintf(file,"MEASUR:%i\n", measur);
		fprintf(file,"OUTNAME:%s\n", outsubix);
		fprintf(file,"#*****CONFIG********\n");
		fprintf(file,"INFILE:%s\n", infile);
		fprintf(file,"#****DIRECTORIES********\n");
		fprintf(file,"DIR_IN:%s\n", dirin);
		fprintf(file,"DIR_OUT:%s\n", dirout);
		fprintf(file,"#**SMOOTHING PARAMETERS**\n");
		fprintf(file,"SOFTMODE:%i\n", smothmod);
		fprintf(file,"SOFTMIN:%f\n", smothmin);
		fprintf(file,"SOFTMAX:%f\n", smothmax);
		fprintf(file,"SOFTCEL:%i\n", smothcell);
		fprintf(file,"#****ENDFILE****");
	}
	fclose(file);
}	


/*! 
 * Create secondary cfg file with the list of raster to be summed
 * writemode: w write file first time
 * writemode: a add line at the end of the file
 * First file has to be the most relevant and the one used as base
 * If all raster are equal, order is not relevant
*/
void sum_listcfg(const char *namcfg, int pros, int ftype, int vnull, double vmin, double vmax,
		int esfloat, int orden, double displa_x, double displa_y, const char *name, const char *code, 
		const char *writemode, int lastdat)
{
FILE *file;
	printf("Writing second sum cfg file %s\n", namcfg);
	if ((file = fopen(namcfg,writemode))== NULL)
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
		if (strcmp (writemode, "w") ==0 )fprintf(file,"pros ftype vnull vmin vmax esfloat orden displa_x displa_y name code\n");
			
			
		if (lastdat == 0) fprintf(file,"%i %i %i %lf %lf %i %i %lf %lf %s %s\n",  
			pros, ftype, vnull, vmin, vmax, esfloat, orden, displa_x, displa_y, name, code);     
		else fprintf(file,"%i %i %i %lf %lf %i %i %lf %lf %s %s",  
			pros, ftype, vnull, vmin, vmax, esfloat, orden, displa_x, displa_y, name, code);  
	}
	fclose(file);
}






#endif /* _CREA_H */
