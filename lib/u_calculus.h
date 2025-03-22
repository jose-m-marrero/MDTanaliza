/*!
* Copyright (C) 2019, 2021  Jose M. Marrero <josemarllin@gmail.com> 
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
* 
* **********************************************************************
* Programa: Calculos geométricos y aritmétricos
* Modulo:  General
* 
* Autores: Jose M. Marrero. Ejemplos tomados de internet y Ramón Ortiz
* 
* Version: 1.1
* Fecha de creacion: 2019-11-01
* Fecha de ultima modificacion: 2022-01-13
* 
* Nombre aplicacion: calculos.h
* 
* Descripcion:
* 
* Realiza diversos calculos geométricos o aritmétricos
* Evalua relaciones de tamaño entre datos
* 
* *********************************************************************/
#ifndef _CACL_H
#define _CACL_H
/*!
* INCLUDES
***********************************************************************/ 
/*! 
* Calculos geometricos, aritmetricos y comparaciones de datos
***********************************************************************/

/**
* Funciones externas
***********************************************************************/
/*! Define array una dimension en integer para operaciones intermedias, modulo u_arrays.h */
int* Crea_1DIarray(int arraySize); 
/*! 
* Librerias estandar
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.14159265

//**********************************************************************
/*! Definicion Funciones */
//**********************************************************************

/*! RADIANES A GRADOS */
double to_degrees(double radians) 
{
    return radians * (180.0 / PI);
}

/*! GRADOS A RADIANES */
double to_radians(double degree) 
{
    return degree * (PI /180.0);
}

/*! Calcula distancia entre dos puntos */
double calc_dist(double x1, double x2, double y1, double y2)
{
double distance;
	distance = 0;
	distance = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	return distance;
}	

/*! Calcula hipotenusa para distancia en 3D */
double calc_hipo(double distan, double alt)
{
double hipotenusa;
	hipotenusa = sqrt(pow(distan, 2)+ pow(alt, 2));
	return hipotenusa;
}

/*! Compara dos valores */
int compara_val(const double a, const double b)
{
int result;
	result=0;
	if (a  < b) result = 1;
	if (a == b) result = 2;
	if (a  > b) result = 3;
	return result;
}

/*! Busca elementos repetidos por valor en un array y altera valor */
double* calc_repetidos(double *val_prin, int totnum)
{
int i, j, get_com, conrep, repe;
double valini, valcom;
	conrep=1;
	for (i=0;i<totnum;i++)
	{
		repe=0;
		valini = val_prin[i];
		for (j=0;j<totnum;j++)
		{
			valcom = val_prin[j];
			get_com = compara_val(valini, valcom);
			if (get_com == 2)repe++;
		}
		if (repe > 1) 
		{
			val_prin[i] = val_prin[i] + (float)(conrep * 0.000001); //sumamos una cantidad small para cambiar el valor repetido
			conrep++;
		}		
	}
	return val_prin;
}

/*! Ordena valores por tamaño, creciente o decreciente */
int* ordena_size(double *val_prin, int totnum, int reorden)	
{
int i, j, get_com;
int *val_sort;
double valini, valcom;
	
	val_sort =  Crea_1DIarray(totnum);
	
	for (i=0;i<totnum;i++)
	{
		//con esta aproximacion el mayor es el idx mas grande
		val_sort[i]=0;
		valini     = val_prin[i];
		for (j=0;j<totnum;j++)
		{
			valcom  = val_prin[j];
			get_com = compara_val(valini, valcom);
			if (get_com == 3)val_sort[i]+=1;
			if (val_sort[i] > totnum) printf("Edimen: Superado numero permitido de valores %i:%i\n", val_sort[i], totnum);
		}	
	}
	if (reorden == 1)
	{
		for (j=0;j<totnum;j++)
		{
			valini      = val_sort[j];
			val_sort[j] = totnum - valini; //hay que cambiarlo al mas small	
		}
	}
	return val_sort;
}

/*! Calcula cuadrante, angulo y angulo perpendacular */
double *get_angulo(double x1, double x2, double y1, double y2, int mode)
{
int cuadrante;
double diffy, diffx, dist, tang, angfin;
static double datos[4];
	//longitud
	diffx = x2 -x1;
	diffy = y2 -y1;
	//cuadrante
	if (diffy > 0 && diffx > 0)cuadrante = 1;                           /*! Determina cuadrante segun signo */
	if (diffy < 0 && diffx > 0)cuadrante = 2;
	if (diffy > 0 && diffx < 0)cuadrante = 4;
	if (diffy < 0 && diffx < 0)cuadrante = 3;
	//angulo que forma la pendiente de la recta
	if (mode == 1)                                                      /*! Calcula arco tangente */
	{
		dist  = calc_dist(x1, x2, y1, y2);                              /*! Calcula distancia entre vertices */
		tang = atan(diffy/dist)*(180/PI);  
	}
	if (mode == 2) tang = atan(diffy/diffx)*(180/PI); 
	//ajuste del angulo al norte y sentido horario        
	if (cuadrante == 1) angfin = 90.0 - tang;                           /*! Calcula angulo perpendacular */
	if (cuadrante == 2) angfin = 90.0 + -(tang);
	if (cuadrante == 3) angfin = 180.0 + tang;
	if (cuadrante == 4) angfin = 360.0 + tang;
	//Almacena resultados
	datos[0] = tang;       //con el 0 al este
	datos[1] = angfin;     //con el 0 al norte
	datos[2] = cuadrante;  //el cuadrante sigue siendo el mismo pero cambia de posicion
	//perpendicular al 0 norte
	if (cuadrante <= 3) datos[3] = angfin + 90;
	if (cuadrante == 4) datos[3] = 90 - (360 - angfin);
	
	//printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", x1, x2, y1, y2, datos[0], datos[1], datos[2], datos[3]);
	return (datos);

}

/*! Calcula si el angulo entre vertices de un poligono cambia o no */
int ang_cambio(double *xver, double *yver, int nver, float ang_chg)
{
int i, cambio, cuadrante, getcua;
double tx, ty, tx2, ty2, tang, getang;
double difang;
double *datos;
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	cambio=0;
	for(i=0;i<nver;i++)
	{
		tx  = xver[i];                                                  /*! Captura primer vertice */
		ty  = yver[i];
		if (i >= 0 && i < nver-1)                                       /*! Segun posicion captura el siguiente */
		{
			tx2 = xver[i+1];
			ty2 = yver[i+1];
		}
		if (i  == nver-1)                                               /*! Si llego al ultimo, captura el anterior */
		{
			tx2 = xver[0];
			ty2 = yver[0];
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		datos     = get_angulo(tx, tx2, ty, ty2, 2);					/*! Calcula angulo, ang. perpendicular y cuadrante */
		tang      = datos[0];
		cuadrante = (int)datos[2];
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (i == 0)
		{
			getang = tang;                                              /*! Inicia valores primer transecto */
			getcua = cuadrante;
		}
		if (i > 0)
		{
			difang = 0;
			if (getcua == cuadrante)                                    /*! Si los cuadrantes son iguales */
			{
				difang = fabs(getang - tang);                           /*! Calcula diferencia de angulos */
				if (difang > ang_chg) cambio++;                         /*! Si es mayor de lo previsto, hay cambio */
			}
			if (getcua != cuadrante)cambio++;                           /*! Si los cuadrantes son distintos, hay cambio */
			getcua = cuadrante;                                         /*! Actualiza valores para comparar con el siguiente vertice */
			getang = tang;
		}
	}
	return cambio;
}


/*! Calcula valor random float entre 0 y el valor maximo dado en float */
float get_random(float max_interval)
{
float ran_value;
	srand(time(NULL));
	ran_value = ((float)rand()/(float)(RAND_MAX)) * max_interval;
	return ran_value;
}

/*! Calcula valor random float entre valor minimo y el valor maximo dado */
float random_float(const float min, const float max)
{
    srand(time(NULL));
    if (max == min) return min;
    else if (min < max) return (max - min) * ((float)rand() / RAND_MAX) + min;

    // return 0 if min > max
    return 0;
}

int get_randomI(int min_val, int max_val)
{
int rand_value;
	srand(time(NULL));
	rand_value = min_val + rand() % (max_val - min_val);
	return rand_value;
}


/*! Calcula incrementos para obtener puntos perpendiculares en secciones transversales
 * Necesita el punto origen, asi como las coordenadas de punto siguiente o previo
 * en funcion de como este dibujado el la linea y la direccion esperada */
double *calc_vector(double txori, double tyori, double nextx, double pretx, double nexty, double prety, double slope, double brect, int modo)
{
double xii, iix;
double yii, iiy;
double a, b, c, d;
double e, f;
double H, I, g;
static double nwcoor[2];
//double mod;
	if(modo == 1)
	{
		xii   = nextx; //zpoint[nex].xcoor;
		yii   = (slope * xii) + brect;
		iix   = pretx; //zpoint[pre].xcoor;
		iiy   = (slope * iix) + brect;
		printf("conver %lf %lf %lf %lf %lf %lf\n", txori, tyori, xii, yii, iix, iiy);
	}	
	if(modo == 2)
	{
		xii     = nextx; //zpoint[nex].xcoor;
		yii     = nexty; //zpoint[nex].ycoor;
		iix     = pretx; //zpoint[pre].xcoor;
		iiy     = prety; //zpoint[pre].ycoor;
	}
	a   =  iix - txori;  //diff in x with next
	b   =  iiy - tyori;  //diff in y with next
	c   =  xii - txori;  //diff in x with pre
	d   =  yii - tyori;  //diff in y with pre
	e   =  sqrt(pow(a,2) + pow(b,2));  //horizontal distance with next
	f   =  sqrt(pow(d,2) + pow(c,2));  //horizontal distance with pre
	H   =  d / f - b / e;   //(dif en y with pre / horizontal dist with pre) - (dif en x with next / horizontal dist with next)
	I   =  a / e - c / f;   //(dif en x with next / horizontal dist with next) - (dif en x with pre / horizontal dist with pre)
	g   =  sqrt(pow(H,2) + pow(I,2));
	nwcoor[0] = H / g; //vnewx
	nwcoor[1] = I / g; //vnewy
	//printf("init en vec %lf %lf %lf %lf %lf %lf mode %i\n", txori, tyori, nextx, nexty, pretx, prety, modo);
	
	//printf("vec a %lf b %lf c %lf d %lf e %lf f %lf H %lf I %lf g %lf nwx %lf nwy %lf\n", a, b, c, d, e, f, H, I, g, nwcoor[0], nwcoor[1]);
		
	//mod   = sqrt(pow(vnewx,2) + pow(vnewy,2));
	return (nwcoor);
}

#endif /* _CACL_H */
