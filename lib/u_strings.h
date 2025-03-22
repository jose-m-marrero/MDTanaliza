/*! 
* **********************************************************************
* Programa: Manejo de string y paths nombres de archivos
* Modulo:  General
* 
* Autores: Jose Manuel (a partir  de ejemplos tomados de internet adaptados)
* 
* Version: 1.2
* Fecha de creacion: 2019-11-01
* Fecha de ultima modificacion: 2029-01-19
* 
* Nombre aplicacion: u_strings.h
* 
* Descripcion:
* 
* Crea y unifica cadenas de texto para los nombres de archivos
* 
* *********************************************************************/
#ifndef _WKSTRING_H
#define _WKSTRING_H

/*!
* DEFINES
***********************************************************************/ 
#define MAX_PATH    260
/*!
* INCLUDES
***********************************************************************/  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libgen.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>



/*! VERIFICA SI ES HDD EXTERNO MEDIANTE PATH */
int es_enmedia(const char *pwd)
{
	if (strstr(pwd, "media") != NULL)
	{
		return 1;
	}
	else return 0;
}	
	
/*! SUSTITUYE HOME POR MEDIA */
char* new_home(char *pwd)
{
int totsize;
char *pth1, *pth2, *pth3;
	
	char *splt_pth = strtok(pwd, "/");                                  /*!< captura primer elemento */
	pth1 = splt_pth;                                                    /*!< almacena primer elemento */
	splt_pth = strtok(NULL, "/");                                       /*!< captura segundo elemento */
	pth2 = splt_pth;                                                    /*!< almacena segundo elemento */
	splt_pth = strtok(NULL, "/");                                       /*!< captura tercer elemento */
	pth3 = splt_pth;                                                    /*!< almacena tercer elemento */
	
	totsize = strlen(pth1)+strlen(pth3)+strlen(pth3)+10;     /*!< Calcula size de toda la cadena completa */
	char *dirhome = (char *)malloc(totsize+1);                            /*!< define dirhome */
	if (dirhome== NULL)
	{
		printf("Atencion, no se pudo generar alocate para %s\n", pwd);
		free(dirhome);
		exit(0);
	}
	else
	{
		sprintf(dirhome, "/%s/%s/%s", pth1, pth2, pth3);                    /*!< Almacena cadena */
		return dirhome;
	}
}

/*! RESERVA ESPACIO DE MEMORIA PARA UN STRING ARRAY
 * ASINGA VALOR */
char *string_array(const char *string)
{
char *new_string;
	new_string = (char *)malloc(strlen(string)+1);
	if (new_string == NULL)
	{
		printf("Erro in malloc for name %s\n", string);
		free(new_string);
		exit(0);
	}
	else strcpy(new_string, string);	
	return new_string;
}

/*! CAPTURA PATH */
char* get_path(char* path)
{
	char *pthfin = (char *)malloc(strlen(path)+1); 
	if (pthfin== NULL)
	{
		printf("Atencion, no se pudo generar alocate para %s\n", path);
		free(pthfin);
		exit(0);
	}
	else
	{
		pthfin = dirname(path);
		strcat(pthfin,"/");
		return pthfin;
	}
}

/*! CAPTURA NOMBRE ARCHIVO */
const char* get_name(const char* path)
{
	if(path != NULL)
	{
		for(size_t i = strlen(path);  i > 0; --i)
		{
			//if (path[i-1] == separator)
			if(path[i-1]=='\\' || path[i-1]=='/' )
			{
				return &path[i];
			}
		}
	}
	return path;
	/*
	char *namfin = (char *)malloc(strlen(path)+1); 
	if (namfin== NULL)
	{
		printf("Atencion, no se pudo generar alocate para %s\n", path);
		free(namfin);
		exit(0);
	}
	else
	{
		namfin = basename(path);
		return namfin;
	}
	*/
}

/*! CAPTURA EXTENSION */
char* get_exten(char* glob_namfile)
{
	printf("nombre %s\n", glob_namfile);
	char* ext = strrchr(glob_namfile, '.');
	if (!ext) {
		/* no extension */
		printf("Atencion, no se capturo extension de %s\n", glob_namfile);
		exit(0);
	} 
	else {
			return ext;
	}
}

/*! ELIMINA EXTENSION DEL ARCHIVO */
char* remove_extension(const char* glob_namfile)
{
int i, size, nstr, ndot, endot, posdot[10];
char *newname, arr[2];

	size = strlen(glob_namfile);
	newname = (char *)malloc(size+1);       /*! define size memoria */
	if (newname== NULL)
	{
		printf("Atencion, no se pudo generar alocate para %s\n", glob_namfile);
		free(newname);
		exit(0);
	}
	else
	{
		strcpy(newname, glob_namfile);        /*! Copia valor de nombre */
		for (i = 0; i < 10; i++)posdot[i]=0;  /*! pon a cero valores */
		nstr=ndot=endot=0;
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		for (i = 0; i < size; i++)
		{
			arr[0] = *(newname+i);        /*! Captura caracter del nombre */
			arr[1] = '\0';			
			if (strcmp(arr, ".") == 0)    /*! Compara caracter para ver si es punto */
			{
				posdot[ndot] = nstr;      /*! Almacena la posicion del caracter */
				ndot++;                   /*! Cuenta numero de puntos */
			}
			nstr++;                       /*! Cuenta numero de caracteres */
		}
		//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		endot = posdot[ndot-1];           /*! Localiza la posicion del ultimo punto */
		newname[endot] = '\0';            /*! Modifica la posicion con terminacion de string */
		//printf("TEST: size %i dot %i endot %i fin %s\n", size, ndot, endot, newname);
		return newname;
	}
}

/*! CALCULA DIMENSION 2 STRING */
int calc_size2(const char* str1, const char* str2)
{
	int tot_size;
	
	tot_size = strlen(str1)+ strlen(str2)+1;
	return tot_size;	
	
}

/*! CONSTRUYE PATH: CALCULA DIMENSION 3 STRING */
int calc_size3(const char* str1, const char* str2, const char* str3)
{
	int tot_size;
	tot_size = strlen(str1)+ strlen(str2)+strlen(str3)+1;
	return tot_size;	
	
}

/*! CONSTRUYE PATH: CALCULA DIMENSION 4 STRING */
int calc_size4(const char* str1, const char* str2, const char* str3, const char* str4)
{
	int tot_size;
	tot_size = strlen(str1)+ strlen(str2)+strlen(str3)+strlen(str4)+1;
	return tot_size;	
	
}

/*! CREA NOMBRE ARCHIVO */
const char* calc_name(const char* pathom, const char* pathdir, const char* filname, int numstr, int totsize)
{
int txtsize;
	if (numstr == 2)txtsize = calc_size2(pathom, filname);              /*!< Calcula longitud string para dos palabras */
	if (numstr == 3)txtsize = calc_size3(pathom, pathdir, filname);     /*!< Calcula longitud string para tres palabras */
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	char filout[txtsize];                                               /*!< Crea cadena con la longitud string necesaria */
	
	if(txtsize > totsize) 
	{
		printf("Se ha excedido size de ruta de archivo en %s\n", filname);
		exit(0);
	}	
	if (numstr == 2)sprintf(filout,"%s%s", pathom, filname);            /*!< asigna valores en cadena */
	if (numstr == 3)sprintf(filout,"%s%s%s", pathom, pathdir, filname); /*!< asigna valores en cadena */
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	if (txtsize < MAX_PATH)
	{
		/*!< no se puede devolver el string porque es una variable local que se destruye, hay que pasarlo a un puntero */
		char *name = (char *)malloc(txtsize+1);                               /*!< crea puntuero con longitud de cadena */
		if (name== NULL)
		{
			printf("Atencion, no se pudo generar alocate para %s\n", pathom);
			free(name);
			exit(0);
		}
		else
		{
			strncpy(name, filout, txtsize);                                     /*!< asigna valor de cadena a puntero */
			return name;                                                        /*!< devuelve puntero */
		}
	}
	else
	{
		printf("Atencion la cadea excede el espacio permitido\n");
		printf("%i: %s %s\n", txtsize,  pathdir, filname);
		printf("Revise ruta de entrada\n");
		exit(0);
	}
}

/*! Construye nombre de archivo de entrada dinamicamente */
char *get_pathnam(const char* pathom, const char* pathdir, const char *namfile, int tipo)
{
int totsize;
char *nomefin;
	if (tipo == 1) 
	{
		totsize = calc_size3(pathom, pathdir, namfile);
		nomefin = (char *)malloc(totsize+1);
		if (nomefin== NULL)
		{
			printf("Atencion, no se pudo generar alocate para %s\n", pathom);
			free(nomefin);
			exit(0);
		}
		else
		{
			sprintf(nomefin,"%s%s%s", pathom, pathdir, namfile);
		}
	}
	if (tipo == 2) 
	{
		totsize = calc_size2(pathom, namfile);
		nomefin = (char *)malloc(totsize+1);
		if (nomefin== NULL)
		{
			printf("Atencion, no se pudo generar alocate para %s\n", pathom);
			free(nomefin);
			exit(0);
		}
		else
		{
			sprintf(nomefin,"%s%s", pathom, namfile);
		}
	}
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//printf("TEST: size nome %i %s\n", totsize, nomefin);
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	return nomefin;
}

/*! Determina si el archivo de entrada csv ha sido correctamente leido
 * Ligado a estructuras de almacenamiento. Contempla distintos tipos de estrategias */
void chk_readfiles(const char *tipofile, int totelem, int totespe, int totpros, int tipopros, int ncol, int especol, int tipocol, int totprev, int tipopre)
{
	printf("Verificacion de archivo %s\n", tipofile);
	printf("Numero total de elementos = %i\n", totelem);
	printf("Numero total de elementos permitidos = %i\n", totespe);
	/*!< si la lectura ha sido correcta */
	if (totelem == 0)
	{
		printf("Atencion error de lectura en archivo %s, se ha leido %i elementos\n", tipofile, totelem);
		exit(0);
	}
	/*!< Cuando el archivo distingue entre procesados y no procesados */
	if (tipopros==1)
	{
		printf("Numero total de elementos a procesar = %i con eval de proceso %i\n", totpros, tipopros);
		if (totpros == 0)  
		{
			printf("Atencion: No se han detectado capas a procesar en %s, revise configuracion\n", tipofile);
			printf("Ponga el valor de la primera columna en 1 para las que quiera procesar\n");
			exit(0);
		}
	}
	/*!< Cuando el archivo evalua numero de columnas posibles */
	if (tipocol==1)
	{
		printf("Numero total de columnas leidas %i con eval columnas %i\n", ncol, especol);
		if (ncol != especol && ncol > 0)
		{
			printf("Atencion el numero de columnas es distinto al esperado\n");
			printf ("Fila: %i lee: %i esperada: %i\n", totelem, ncol, especol); 
			exit(0);
		}
	}
	/*!< si se ha excedido el size permitido de almacenamiento */
	if (totelem > totespe)
	{
		printf("Atencion, se han excedido el numero de elementos permitidos %i en %i\n", totespe, totelem);
		printf("Incremente valor en #define correspondiente para %s\n", tipofile);
		exit(0);
	}
	/*!< Cuando existe otro archivo que debe tener la misma dimension que este, por ejeplo edi y edixyz */
	if (tipopre==1)
	{
		if (totelem != totprev)
		{
			printf("Atencion las filas actuales e iniciales no coinciden en %s\n", tipofile);
			exit(0); 
		}
	}
}

/*! Evalua si existe archivo */
int if_filexist(const char *fname)
{
	if( access( fname, F_OK ) == 0 ) {
		return 1;  //existe
	} else {
		return 0;
	}
}

/*! Evalua si existe directorio */
int if_direxist(const char *dir)
{
	printf("chk %s\n", dir);
	struct stat stats;
	if(stat(dir,&stats) == 0)
	{
		return 1;
	}
	return 0;
}


/*! Crea directorio */
void crea_mkdir(char *path, mode_t mode)
{
char *restpath, *sumpath; /*! Disminuye path y suma path, en uno elimina elemento en el otro los va añadiendo */
int cont;
	
	//intenta la primera vez si existe
	if (mkdir(path, mode) == 0)                                         /*! Intenta crear directorio */
	{
	   printf("Directorio creado\n");
	}
	else
	{
		/*! Si no puede, evalua si existe o no */
		if (errno == EEXIST) printf("El directorio Ya existe\n");   
		if (errno != EEXIST) 
		{
			/*! Si no existe, descompona y evalua toda la estructura */
			printf("Atencion: Directorio no creado, se analiza el arbol completo\n");
			restpath = (char *)malloc(strlen(path)+1);                  /*! Reserva memoria */
			sumpath = (char *)malloc(strlen(path)+1);
			strcpy(restpath, path);                                     /*! Copia el path completo en var alternativa */
			char *ptr = strtok(restpath, "/");                          /*! Caputra primer nivel path */
			strcpy(sumpath, "/");                                       /*! Crea / como primer elemento del path */
			cont=0;                                                     /*! Contador del numero de niveles */
			while (ptr != NULL)                                         /*! Mientras existan niveles */
			{
				if (cont > 0)strcat(sumpath, "/");                      /*! A partir del segundo nivel isertamos separador */
				strcat(sumpath, ptr);                                   /*! Concatenamos nuevo nivel */
				if (cont > 1)                                           /*! Si estamos mas alla del /home/usuario/ */
				{
					if (!if_direxist(sumpath))                          /*! Si el directorio NO existe */
					{
						if (mkdir(sumpath, mode) == -1)                 /*! Lo creamos, si devuelve error */
						{
							if (errno == EEXIST) printf("Porque el directorio Ya existe\n");   
							if (errno != EEXIST) 
							{
								printf("Hubo un error al crear la ruta '%s'\n%m\n", sumpath); 
								exit(0);
							}
						}
						else printf("Nuevo directorio creado en %i nivel\n", cont);  /*! Lo creamos, sin error */
					}
				}
				ptr = strtok(NULL, "/");                                /*! Capturamos siguiente nivel */
				cont++;                                                 /*! Adicionamos nivel */
			}
		}
	}
}

#endif /* _WKSTRING_H */
