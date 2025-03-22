/*!
* Copyright (C) 2019, 2021  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Programa: Analisis raster
* Modulo:  m_topolevel.h
* 
* Autores: Jose M. Marrero 
* 
* Version: 0.1.0
* Fecha de creacion: 2022-11-23
* Fecha de ultima modificacion: 2022-09-30
* 
* Nombre aplicacion: m_topolevel.h
* 
* Descripcion:
* 
* Calculate topographic profiles
* Get heigh value and heigh differences along a line between two points
* 
* *********************************************************************/
#ifndef _TOPL_RAS
#define _TOPL_RAS

Calcula la recta entre dos puntos
Anda por la recta cada cierta distancia
Captura la altura
Captura la diferencia de cota entre la altura y la linea que une los dos puntos
Genera dos columnas con ambos datos
Representa los resultados del perfil topografico y la linea entre ambos
Paraña recta hace falta la altura inicial y la final 
