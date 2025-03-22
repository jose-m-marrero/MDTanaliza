/*!
* Copyright (C) 2019, 2025  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Manin module: u_arrays.h
* Module:  struct_head-array.sh
* 
* Author: Jose M. Marrero 
* 
* Version: 0.1.1
* Creation Date: 2019-11-01
* Last update: 2022-03-26
* 
* Description: 
* 
* Store header from *.grd raster format
*/

#ifndef _HEAD_A
#define _HEAD_A


typedef struct HeadR
{
	int hn_cx;         //numero de columnas
	int hn_fy;         //numereo de filas
	double hxlo;       //coordenada x minima
	double hxhi;       //coordenada x maxima
	double hylo;       //coordenada y minima
	double hyhi;       //coordenada y maxima
	double hzlo;       //coordenada z minima
	double hzhi;       //coordenada z maxima
	double hresx;      //resolucion en x
	double hresy;      //resolucion en y
	double hinvresx;
	double hinvresy;
	double hresxx;
	double hresyy;
	int htotdemcel;    //total de celdas
	double hperim;     //perimetro
	double harea;      //area
	
}HeadR;


#endif /* _HEAD_A */
