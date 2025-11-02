/*!
* Copyright (C) 2019, 2021  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Programa: Lectura de arrays
* Modulo:  General
* 
* Autores: Jose M. Marrero. Ejemplos tomados de internet y Ramón Ortiz
* 
* Version: 1.2.3
* Fecha de creacion: 2019-11-01
* Fecha de ultima modificacion: 2022-03-26
* 
* Nombre aplicacion: u_geograph.h
* 
* Descripcion:
* 
* Lectura de archivos raster y manejo dinámico de memoria
* Calcula de indices y coordenadas en raster
* Importante: El numero de raster comienza en 0
* 
* *********************************************************************/
#ifndef _GEO_H
#define _GE0_H
/*!
* INCLUDES
***********************************************************************/  
/*! 
* Calculos geometricos, aritmetricos y comparaciones de datos
***********************************************************************/
#include "u_calculos.h"

/*! 
* Librerias estandar
***********************************************************************/
#include <stdio.h>
#include <math.h>

#define PI 3.14159265



/*! CONVERSION GEOGRAFICAS A UTM */
int convert_geograf(double tlat, double tlong)
{
int thuso, mercentgrad;
double a, b, a2, b2, e, eb, c, alfa;	
double mercentrad, distangrad, distanrad;
double A, xi, eta, v, dseta, A1, A2, J2, J4, J6;
double alfab2, beta, gammag, B;
double tlatrad;
double tlonrad;

	//info en http://www.gabrielortiz.com/index.asp?Info=058a
	
	//Elipsoide GRS80
	//Semieje mayor
	a = 6378137.0000000000;  //6378388.0  elipsoide de Hayford
	//Semieje menor
	b = 6356752.3142500000;  //6356911.946130  elipsoide de Hayford
	a2 = pow(a, 2);
	b2 = pow(b, 2);
	//----------------------
	//Excentricidad
	//e = sqrt((a2 - b2)/a) = 0.08199189; 
	e =1-pow((b / a),2);
	//Segunda Excentricidad
	//eb = sqrt((a2 - b2)/b) = 0.08226889;
	eb = pow((a / b),2) - 1;
	//radio polar de Curvatura
	c =	a2 / b; //6399936.608
	//aplanamiento o achatamiento
	alfa = (a - b) / b; //0.003367003
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	thuso   = 28;
	tlatrad = to_degrees(tlat);
	tlonrad = to_degrees(tlong);
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//Calcular el Huso
	if(thuso == -9999) thuso = (tlong / 6) + 31;
	//Meridiano central
	mercentgrad = (thuso * 6) -183;
	mercentrad  = (mercentgrad * PI) / 180;
	//Distancia angular
	distangrad = tlong - mercentgrad;
	distanrad = (distangrad * PI) / 180;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//Parametros ecuaciones Coticchia-Surace
	A        = cos(tlatrad) * sin(distanrad);
	xi       = 0.5 * log((1 + A) / (1 - A));
	eta      = atan(tan(tlatrad)/cos(distanrad)) - tlatrad;
	v        = c * 0.9996 / sqrt(1 + eb * pow(cos(tlatrad),2));
	dseta    = (eb / 2) * pow(xi, 2) * pow(cos(tlatrad),2);
	A1       = sin(2 * tlatrad);
	A2       = A1 * pow(cos(tlatrad),2);
	J2       = tlatrad + (A1 / 2);
	J4       = ((3 * J2) + A2) / 4;
	J6       = ((5 * J4) + (A2 * pow(cos(tlatrad),2))) / 3;
	alfab2   = 3.0 / 4.0 * eb;
	beta     = (5.0 / 3.0) * pow(alfab2, 2);
	gammag   = (35.0 / 27.0) * pow(alfab2, 3); //como una y
	B        = 0.9996 * c * (tlatrad-(alfab2 * J2) + (beta * J4) - (gammag * J6));
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//Solucion coordenadas UTM
	xcoor = xi * v * (1 + (dseta / 3)) + 500000;
	//y depende del hemisferio
	if(tlat > 0)  ycoor = ((eta * v) * (1 + dseta)) + B;
	if(tlat < 0)  ycoor = (((eta * v) * (1 + dseta)) + B) + 10000000;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	printf("a = %.8lf\n", a);
	printf("b = %.8lf\n", b);
	printf("a2 = %.8lf\n", a2);
	printf("b2 = %.8lf\n", b2);
	printf("e = %.8lf\n", e);
	printf("eb = %.8lf\n", eb);
	printf("c = %.8lf\n", c);
	printf("alfa = %.8lf\n", alfa);
	printf("thuso = %i\n", thuso);
	printf("mercentgrad = %i\n", mercentgrad);
	printf("mercentrad = %.8lf\n", mercentrad);
	printf("distangrad = %.8lf\n", distangrad);
	printf("distanrad = %.8lf\n", distanrad);
	printf("A = %.8lf\n", A);
	printf("xi = %.8lf\n", xi);
	printf("eta = %.8lf\n", eta);
	printf("v = %.8lf\n", v);
	printf("dseta = %.8lf\n", dseta);
	printf("A1 = %.8lf\n", A1);
	printf("A2 = %.8lf\n", A2);
	printf("J2 = %.8lf\n", J2);
	printf("J4 = %.8lf\n", J4);
	printf("J6 = %.8lf\n", J6);
	printf("alfab2 = %.8lf\n", alfab2);
	printf("beta = %.8lf\n", beta);
	printf("gammag = %.8lf\n", gammag);
	printf("B = %.8lf\n", B);
	printf("X = %.8lf\n", xcoor);
	printf("Y = %.8lf\n", ycoor);
	printf("rad en log %lf\n", tlonrad);
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// Imprimimos los datos por pantalla para verificar que la lectura es correcta
	printf("Coor %.6lf; %.6lf Huso %i- Coorgeo %.8lf; %.8lf\n", xcoor, ycoor, thuso, tlat, tlong);
	printf("--------\n");
    return 1;
    
}



#endif /* _GE0_H */
