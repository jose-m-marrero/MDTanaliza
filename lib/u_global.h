/*! 
* **********************************************************************
* Programa: Creacion de nombres de archivos
* Modulo:  General
* 
* Autores: Jose Manuel (a partir  de ejemplos tomados de internet)
* 
* Version: 1.0-2021-12-03
* Fecha de creacion: 2019-11-01
* Fecha de ultima modificacion: 2021-12-03
* 
* Nombre aplicacion: crea_nombre.h
* 
* Descripcion:
* 
* Crea y unifica cadenas de texto para los nombres de archivos
* 
* *********************************************************************/
#ifndef _UTILS_H
#define _UTILS_H
/*!
* INCLUDES
***********************************************************************/  
#include <stdio.h>
#include <string.h>
#include <time.h>


/*! TIME */
int getime(void)
{
    time_t mytime = time(NULL);
    char * time_str = ctime(&mytime);
    time_str[strlen(time_str)-1] = '\0';
    printf("Current Time : %s\n", time_str);
    return 0;
}



#endif /* _UTILS_H */
