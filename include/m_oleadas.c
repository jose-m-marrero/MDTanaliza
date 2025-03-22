/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Program name: MDTanaliza
* Modulo:  m_oleadas.h
* 
* Autores: Ramon Ortiz 
*          Jose M. Marrero
*          Based in Felpeto.....
* 
* Version: 0.1.0
* Fecha de creacion: 2007-04-06
* Fecha de ultima modificacion: 2022-09-30
* 
* Nombre aplicacion: m_oleadas.h
* 
* Descripcion:
* 
* Model Pyroclastic Flow and Block and Ash volcanic hazards
* Based in energy cone model
* 
* *********************************************************************/
#ifndef _OLE_RAS
#define _OLE_RAS

/*! 
* Variables globales y estructuras
***********************************************************************/
#include "../include/m_global.h"
/*! 
* Almacena nombres de directorios y archivos de las capas a analizar
***********************************************************************/
#include "../include/struc_mdtanaliza.h"

/*!
* INCLUDES - MODULOS USO GENERAL
***********************************************************************/ 
/*! 
* Lectura de raster y manejo de arrays*
***********************************************************************/
#include "u_arrays.h"
/*! 
* Almacena datos de cabecera, el numero de arrays se define en el main
***********************************************************************/
#include "struct_array-cabecera.h"
/*! 
* Calculos geometricos, aritmetricos y comparaciones de datos
***********************************************************************/
#include "u_calculos.h"
/*! 
* Utilidades varias *
***********************************************************************/
#include "u_strings.h"


#include <string.h>
/*! Almacena los datos de la cabecera de arrays-raster grd */
struct CabeceraR s_arrayh[NUMARR];

extern double **rast_volhaz;



//**********************************************************************
/*! Oleadas Functions */
//**********************************************************************


/*!
 * Article (kaneko1992) Kaneko Takayuki anb Suzuki-Kamata, K. 
 * Energy Line/Cone Simulations of the 1991 pyroclastic flows of Unzen Volcano 
 * Bulletin of the Volcanological Society of Japan, 1992, 37, 35-45
*/

void oleada(double **raster, struct CabeceraR s_arrayh[], double alfa, double ht, int idxyf, int idxcx)
{
int i, j;
double s, d;
	
	for (j=0;j<s_arrayh[0].hn_fy;j++)
	{
        for(i=0;i<s_arrayh[0].hn_cx;i++)
		{
			d = sqrt((s_arrayh[0].hresx * 2)*(i-idxcx)*(i-idxcx) + (s_arrayh[0].hresy * 2)*(j-idxyf)*(j-idxyf)); 
			s = ht - d*alfa - raster[j][i];
			if(s >= 0)rast_volhaz[j][i] = s;
			else      rast_volhaz[j][i] = 0;
		}
    }
}

int barre(double **raster, struct CabeceraR s_arrayh[])
{
int i, j, tot;
double s;
    tot = 0;

	for (j=1;j<s_arrayh[0].hn_fy-1;j++)
	{
        for(i=1;i<s_arrayh[0].hn_cx-1;i++)
		{
		    if(rast_volhaz[j][i]<0)
		    {
		        s = -rast_volhaz[j][i]; //change in positive
		        if(rast_volhaz[j+1][i]>0   && raster[j][i]<= raster[j+1][i]  ) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j-1][i]>0   && raster[j][i]<= raster[j-1][i]  ) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j][i+1]>0   && raster[j][i]<= raster[j][i+1]  ) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j][i-1]>0   && raster[j][i]<= raster[j][i-1]  ) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j+1][i+1]>0 && raster[j][i]<= raster[j+1][i+1]) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j+1][i-1]>0 && raster[j][i]<= raster[j+1][i-1]) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j-1][i+1]>0 && raster[j][i]<= raster[j-1][i+1]) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j-1][i-1]>0 && raster[j][i]<= raster[j-1][i-1]) rast_volhaz[j][i] = s;
		        if(rast_volhaz[j][i]>0)tot = 1;
		    }
        }
    }
    return tot;
}

#endif /* _OLE_RAS */
