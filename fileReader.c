#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"

int readFile(char * files, int flag){
	
	char * pathFile;
	char file_str[MAX_file_len]="\0";
	FILE   *nomi_file;
	
	int ii;


	if (flag == PRE) pathFile = "./fileList/listaLFSIDE.txt";
	else pathFile = "./fileList/BV_res.txt";

	
     	if ( (nomi_file = fopen(pathFile, "r")) == NULL)
           printf("cazzo 1\n");
 	ii = 0 ; 
  	while (fscanf(nomi_file, "%s\n", file_str)!=EOF) {
           sscanf(file_str,"%s\0",(files+ii*MAX_file_len));
          ii++;
        }
        printf ("%s\n",pathFile);
	int n_righe_file = ii;
        fclose(nomi_file);

	return n_righe_file;
	
}
