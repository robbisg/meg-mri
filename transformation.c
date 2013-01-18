#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"


/*Main coordinates transformation routine*/
void transform (double * coordIn, double * coordOut, double * conf){

    int i;
    
    double * angoli = calloc(3, sizeof (double));
    double * trasl = calloc(3, sizeof(double));
    double * zoom = calloc(3, sizeof(double));
    double * shear = calloc(3, sizeof(double));
    
    double vettZ[3], vettS[3], vettR[3]; 
    
    for (i = 0; i < 3; i++){
           
	trasl[i] = conf[i+0];
	angoli[i] = conf[i+3];
	zoom[i] = conf[i+6];
	shear[i] = conf[i+9];
    }
      
      zooming(coordIn, vettZ, zoom);
      shearing(vettZ, vettS, shear);
      rotationForward(vettS,vettR, angoli);
      //rotationRewind(vettS, vettR, angoli);
      
    for (i = 0; i < 3; i++){
	
	coordOut[i] = vettR[i] + trasl[i];
    
    }
    
    free(zoom); free(shear); free(angoli); free(trasl);
    
    
    return;
    
}

void transformInv (double * coordIn, double * coordOut, double * conf){
    
    int i;
    
    double * angoli = calloc(3, sizeof (double));
    double * trasl = calloc(3, sizeof(double));
    double * zoom = calloc(3, sizeof(double));
    double * shear = calloc(3, sizeof(double));
    
    double vettZ[3], vettS[3], vettR[3]; 
    
    for (i = 0; i < 3; i++){
           
	trasl[i] = conf[i+0];
	angoli[i] = conf[i+3];
	zoom[i] = conf[i+6];
	shear[i] = conf[i+9];
    }
      
      /*      |1  a  b|         |1  -a  ac-b|   */
      /*  S = |0  1  c|  S^-1 = |0   1   -c |   */
      /*      |0  0  1|         |0   0    1 |   */
     
      shear[0] = -1*conf[9];
      shear[1] = conf[9]*conf[11]-conf[10];
      shear[2] = -1*conf[11];
      
      
      /*      |a  0  0|         |1/a  0   0 |   */
      /*  Z = |0  b  0|  Z^-1 = | 0  1/b  0 |   */
      /*      |0  0  c|         | 0   0  1/c|   */
      
      
     
      zoom[0] = 1/zoom[0];  
      zoom[1] = 1/zoom[1];
      zoom[2] = 1/zoom[2];
      
            
      rotationRewind(coordIn, vettZ, angoli);
      shearing(vettZ, vettS, shear);
      zooming(vettS,vettR, zoom);
      
      
    for (i = 0; i < 3; i++){
	
	coordOut[i] = vettR[i] - trasl[i];
    
    }
    
    free(zoom); free(shear); free(angoli); free(trasl);
    
    
    return;
}



/*Transformation second approach*/
void transformLight (double * vettOrig, double * matrixRot, double * vettPost, double * trasl){
	
	int i;
	righe_per_colonne(vettOrig, matrixRot, vettPost);
	for (i = 0; i < 3; i++){
	  vettPost[i] = vettPost[i] + trasl[i];
	}
	
	
	
}



void invTransformLight (double * vettOrig, double * matrixRot, double * vettPost, double * trasl){
	int i;
	righe_per_colonne(vettOrig, matrixRot, vettPost);
	for (i = 0; i < 3; i++){
	  vettPost[i] =  vettPost[i] - trasl[i];
	}

}



void forwardMatrix (double * matrixOut, double * conf){

  
  int i;
  
    double * angoli = calloc(3, sizeof (double));
    double * zoom = calloc(3, sizeof(double));
    double * shear = calloc(3, sizeof(double));

    for (i = 0; i < 3; i++){
           
	angoli[i] = conf[i+3];
	zoom[i] = conf[i+6];
	shear[i] = conf[i+9];
    }
    
    double * matrixRot = calloc(9, sizeof(double));
    double * matrixS = calloc(9, sizeof(double));
    double * matrixZ = calloc(9, sizeof(double));
    double * matrixRS = calloc(9, sizeof(double));
    
    matrix_XYZ(angoli, matrixRot);
    matrix_shear(shear, matrixS);
    matrix_zoom(zoom, matrixZ);
    
    
    /*Given a vector x the transformation of zoom, shear and rotation in this
      order, is given by: 
    
      x' = Z*x;
      x'' = S*x' = S*Z*x;
      x''' = R*x'' = R*S*Z*x;
    
      So we have to calculate R*S*Z matrix to have an 
      output transformation matrix*/
    
    
    /*Product R*S*/
    matrixProduct(matrixRot, matrixS, matrixRS, 3);
     
    /*Product (R*S)*Z associative property*/
    matrixProduct(matrixRS, matrixZ, matrixOut, 3);
    
    
    free(matrixRS);
    free(matrixS);
    free(matrixZ);
    free(matrixRot);
        
    free(zoom); free(shear); free(angoli);
}

void inverseMatrix (double * matrixOut, double * conf){

  
  int i;
  
    double * angoli = calloc(3, sizeof (double));
    double * zoom = calloc(3, sizeof(double));
    double * shear = calloc(3, sizeof(double));
    for (i = 0; i < 3; i++){
           
	angoli[i] = conf[i+3];
	zoom[i] = conf[i+6];
	shear[i] = conf[i+9];
    }

           
	      
      /*      |1  a  b|         |1  -a  ac-b|   */
      /*  S = |0  1  c|  S^-1 = |0   1   -c |   */
      /*      |0  0  1|         |0   0    1 |   */
     
      shear[0] = -1*conf[9];
      shear[1] = conf[9]*conf[11]-conf[10];
      shear[2] = -1*conf[11];
      
      
      /*      |a  0  0|         |1/a  0   0 |   */
      /*  Z = |0  b  0|  Z^-1 = | 0  1/b  0 |   */
      /*      |0  0  c|         | 0   0  1/c|   */
      
      
     
      zoom[0] = 1/zoom[0];  
      zoom[1] = 1/zoom[1];
      zoom[2] = 1/zoom[2];
    
    
    double * matrixRot = calloc(9, sizeof(double)) ;
    double* matrixS = calloc(9, sizeof(double))  ;
    double* matrixZ = calloc(9, sizeof(double))  ;
    double *matrixRS  = calloc(9, sizeof(double)) ;
    
    
    matrix_zoom(zoom, matrixZ);
    matrix_shear(shear, matrixS);
    matrix_XYZ(angoli, matrixRot);
    
    
    /*Given a vector x the transformation of rotation, shear and zoom in this
      order, is given by: 
    
      x' = R*x;
      x'' = S*x' = S*R*x;
      x''' = Z*x'' = Z*S*R*x;
    
      So we have to calculate Z*S*R matrix to have an 
      output transformation matrix**/
    
    
    /*Product Z*S*/
    matrixProduct(matrixZ, matrixS, matrixRS, 3);
     
    /*Product (Z*S)*R associative property*/
    matrixProduct(matrixRS, matrixRot, matrixOut, 3);
    
    
    free(matrixRS);
    free(matrixS);
    free(matrixZ);
    free(matrixRot);
        
    free(zoom); free(shear); free(angoli);
}

void righe_per_colonne(double *vect_in, double *matrix, double *vect_out) {

    *vect_out = *(matrix) * (*vect_in)
            + *(matrix+1) * (*(vect_in+1))
            + *(matrix+2) * (*(vect_in+2));

    *(vect_out+1) = *(matrix+3) * (*vect_in)
                + *(matrix+4) * (*(vect_in+1))
                + *(matrix+5) * (*(vect_in+2));

    *(vect_out+2) = *(matrix+6) * (*vect_in)
                + *(matrix+7) * (*(vect_in+1))
                + *(matrix+8) * (*(vect_in+2));
   return;
}

void prod_vect(double *vect_1, double *vect_2, double *vect_out){

  *(vect_out)=*(vect_1+1) * *(vect_2+2) - *(vect_1+2) * *(vect_2+1) ;
  *(vect_out+1)=*(vect_1+2) * *(vect_2) - *(vect_1) * *(vect_2+2) ;
  *(vect_out+2)=*(vect_1) * *(vect_2+1) - *(vect_1+1) * *(vect_2) ;


  //printf("Ruoto attorno a Z di un angle =%f\n",theta/3.1416*180.);
    return;
}

void matrix_X (double phi, double *matrix){

    *(matrix)=1.;
    *(matrix+1)=0;
    *(matrix+2)=0.;
    *(matrix+3)=0;
    *(matrix+4)=cos(phi);
    *(matrix+5)=-sin(phi);
    *(matrix+6)=0.;
    *(matrix+7)=sin(phi);
    *(matrix+8)=cos(phi);

//    printf("Ruoto attorno a X di un angle =%f\n",phi/3.1416*180.);
    return;
}

void matrix_Y(double psi, double *matrix){

    *(matrix)=cos(psi  );
    *(matrix+1)=0.;
    *(matrix+2)=sin(psi  );
    *(matrix+3)=0.;
    *(matrix+4)=1.;
    *(matrix+5)=0.;
    *(matrix+6)=-sin(psi  );
    *(matrix+7)=0.;
    *(matrix+8)=cos(psi  );

//    printf("Ruoto attorno a Y di un angle =%f\n",psi/3.1416*180.);
    return;
}

void matrix_Z(double theta, double *matrix){

    *(matrix)=cos(theta);
    *(matrix+1)=-sin(theta);
    *(matrix+2)=0.;
    *(matrix+3)=sin(theta);
    *(matrix+4)=cos(theta);
    *(matrix+5)=0.;
    *(matrix+6)=0.;
    *(matrix+7)=0.;
    *(matrix+8)=1.;

//    printf("Ruoto attorno a Z di un angle =%f\n",theta/3.1416*180.);
    return;
}

void matrix_ZYX(double * tre_angoli, double * matrix){
    
/*Matrice di rotazione su ZYX: Effettuata come prodotto fra (Mx*My)*Mz
  dove Mx, My, Mz sono le matrici di rotazione sugli assi X, Y, Z.*/
	



	double alpha, theta, phi;
	alpha =  *(tre_angoli+2);	//Angolo di rotazione attorno a X
	phi =-1 * *(tre_angoli+1);	//attorno a Y
	theta =*(tre_angoli+0);		//attorno a Z

/*Le rotazioni sono espresse considerando il senso anti-orario*/	

    *(matrix+0) = cos(phi)*cos(theta);
    *(matrix+1) = sin(theta)*cos(phi);
    *(matrix+2) =-sin(phi);
    *(matrix+3) = sin(alpha)*sin(phi)*cos(theta) - cos(alpha)*sin(theta);
    *(matrix+4) = sin(alpha)*sin(phi)*sin(theta) + cos(alpha)*cos(theta);
    *(matrix+5) = sin(alpha)*cos(phi);
    *(matrix+6) = cos(alpha)*sin(phi)*cos(theta) + sin(alpha)*sin(theta);
    *(matrix+7) = cos(alpha)*sin(phi)*sin(theta) - sin(alpha)*cos(theta);
    *(matrix+8) = cos(alpha)*cos(phi);
    
    return;
}


void matrix_XYZ (double * tre_angoli, double * matrix){
  
/*Matrice di rotazione su XYZ: Effettuata come prodotto fra (Mz*My)*Mx
  dove Mx, My, Mz sono le matrici di rotazione sugli assi X, Y, Z.*/

 
 
	double alpha, theta, phi;
	alpha =-1 *  *(tre_angoli+2);//Angolo di rotazione attorno a X
	phi = *(tre_angoli+1);//attorno a Y
	theta =-1 * *(tre_angoli+0);//attorno a Z

/*Le rotazioni sono espresse considerando il senso anti-orario*/	

    *(matrix+0) = cos(theta)*cos(phi);
    *(matrix+3) =-sin(theta)*cos(phi);
    *(matrix+6) = sin(phi);
    *(matrix+1) = sin(alpha)*sin(phi)*cos(theta) + cos(alpha)*sin(theta);
    *(matrix+4) =-sin(alpha)*sin(phi)*sin(theta) + cos(alpha)*cos(theta);
    *(matrix+7) =-sin(alpha)*cos(phi);
    *(matrix+2) =-cos(alpha)*sin(phi)*cos(theta) + sin(alpha)*sin(theta);
    *(matrix+5) = cos(alpha)*sin(phi)*sin(theta) + sin(alpha)*cos(theta);
    *(matrix+8) = cos(alpha)*cos(phi);
    
    return;
}

void matrix_zoom (double * tre_zooming, double * matrix){
  
    double zoomX = *(tre_zooming+0);
    double zoomY = *(tre_zooming+1);
    double zoomZ = *(tre_zooming+2);
    
    *(matrix+0) = zoomX;
    *(matrix+4) = zoomY;
    *(matrix+8) = zoomZ;
    
    *(matrix+1) = 0.;
    *(matrix+2) = 0.;
    *(matrix+3) = 0.;
    *(matrix+5) = 0.;
    *(matrix+6) = 0.;
    *(matrix+7) = 0.;
    
    return;
}

void matrix_shear (double * tre_shear, double * matrix){
  
    double shearX = *(tre_shear+0);
    double shearY = *(tre_shear+1);
    double shearZ = *(tre_shear+2);
    
    *(matrix+1) = shearX;
    *(matrix+2) = shearY;
    *(matrix+5) = shearZ;
    
    *(matrix+0) = 1.;
    *(matrix+4) = 1.;
    *(matrix+8) = 1.;
    
    *(matrix+3) = 0.;
    *(matrix+6) = 0.;
    *(matrix+7) = 0.;
    
    return;
}
void rotationRewind (double * vett_orig, double * vett_roto, double * tre_angoli){

	/*Effettua una rotazione su XYZ */

	double * matrixRot = calloc(9, sizeof(double));

	matrix_XYZ(tre_angoli, matrixRot);
    	righe_per_colonne(vett_orig, matrixRot, vett_roto);

	free(matrixRot);

}

void rotationForward (double * vett_orig, double * vett_roto, double * tre_angoli){
	
	/*Effettua una rotazione su ZYX */

	double * matrixRot = calloc(9, sizeof(double));
	
	matrix_ZYX(tre_angoli, matrixRot);
    	righe_per_colonne(vett_orig, matrixRot, vett_roto);

	free(matrixRot);
}


void matrixProduct (double * matrix1, double * matrix2, double * result, int dim){
    /*Effettua il prodotto fra matrici quadrate*/
    
    /*dim è la dimensione della matrice quadrata*/
    
    int i, j, c, r;
    

    for (r = 0; r < dim; r++) {
      for (c = 0; c < dim; c++){
	for (i = 0; i < dim; i++){
	  //printf("%lf * %lf -- ",matrix1[r*dim+i],matrix2[i*dim+c]);
	  result[r*dim+c] += matrix1[r*dim+i]*matrix2[i*dim+c];
	}
	//printf("\n");
      }
    }
}


void shearing (double * vettOrig, double * vettPost, double * tre_shear){
     
  
      double * matrixS = calloc(9, sizeof(double));
     
      
      matrix_shear(tre_shear,matrixS);
      righe_per_colonne(vettOrig, matrixS, vettPost);
      
      free(matrixS);
      
      return;
      
}

void zooming (double * vettOrig, double * vettPost, double * tre_zoom){
	
  
	double * matrixZ = calloc(9, sizeof(double));


		
	matrix_zoom(tre_zoom, matrixZ);
    	righe_per_colonne(vettOrig, matrixZ,vettPost);
	
	free(matrixZ);
	
	return;
}

void testMatrixProduct (){
	  
	double  matrix1[] = {1.,2.,3.,4.,5.,6.,7.,8.,9.};
	double  matrix2[] = {9.,8.,7.,6.,5.,4.,3.,2.,1.};
	double * result = calloc (9, sizeof(double));
	
	matrixProduct(matrix1,matrix2,result,3);
	
	printa(result, 9);

}

void  ruotamento_on_X_Y_Z(double *v_orig_1, double *v_orig_2, double *v_orig_3, double *tre_angoli) {

/* Ruoto autovettori pre ut orientarli su assi X,Y,Z */

   double cos_alpha, sin_alpha, alpha, cos_phi, sin_phi, phi, cos_psi, sin_psi,psi;
   double *v_rotz_1, *v_rotz_2, *v_rotz_3;
   double *v_roty_1, *v_roty_2, *v_roty_3;
   double *v_rotx_1, *v_rotx_2, *v_rotx_3;

   double matrix_roZ[9]={0};
   double matrix_roX[9]={0};
   double matrix_roY[9]={0};

   v_rotz_1 = calloc (3, sizeof(double));
   v_rotz_2 = calloc (3, sizeof(double));
   v_rotz_3 = calloc (3, sizeof(double));
   v_roty_1 = calloc (3, sizeof(double));
   v_roty_2 = calloc (3, sizeof(double));
   v_roty_3 = calloc (3, sizeof(double));
   v_rotx_1 = calloc (3, sizeof(double));
   v_rotx_2 = calloc (3, sizeof(double));
   v_rotx_3 = calloc (3, sizeof(double));

    /* Angolo sul piano X, Y (X -> Y rotazione positiva) */
    cos_alpha=*v_orig_1 / sqrt( (*v_orig_1 * *v_orig_1 + *(v_orig_1+1) * *(v_orig_1+1)));
    sin_alpha=*(v_orig_1+1) / sqrt( (*v_orig_1 * *v_orig_1 + *(v_orig_1+1) * *(v_orig_1+1)));

    if (sin_alpha >= 0) alpha = acos(cos_alpha);
    else alpha = -1. * acos(cos_alpha);      //Angolo rispetto a X
    
	*(tre_angoli+0) = alpha;

    /* Ruoto attorno a Z di un angolo opposto per andare nel piano X-Z */
    matrix_Z(-1.*alpha  ,matrix_roZ);
    righe_per_colonne(v_orig_1, matrix_roZ, v_rotz_1);
  //  printf("v_rotz_1[0]=%f v_rotz_1[1]=%f v_rotz_1[2]=%f \n",*v_rotz_1,*(v_rotz_1+1),*(v_rotz_1+2));
    righe_per_colonne(v_orig_2, matrix_roZ, v_rotz_2);
  //  printf("v_rotz_2[0]=%f v_rotz_2[1]=%f v_rotz_2[2]=%f \n",*v_rotz_2,*(v_rotz_2+1),*(v_rotz_2+2));
    righe_per_colonne(v_orig_3, matrix_roZ, v_rotz_3);
  //  printf("v_rotz_3[0]=%f v_rotz_3[1]=%f v_rotz_3[2]=%f \n",*v_rotz_3,*(v_rotz_3+1),*(v_rotz_3+2));

    /* Angolo sul piano Z, X del vettore appena ruotato (Z -> X rotazione positiva)*/
    cos_phi=*(v_rotz_1) / sqrt( (*v_rotz_1 * *v_rotz_1 + *(v_rotz_1+2) * *(v_rotz_1+2)));
    sin_phi=*(v_rotz_1+2) / sqrt( (*v_rotz_1 * *v_rotz_1 + *(v_rotz_1+2) * *(v_rotz_1+2)));

    if (sin_phi >= 0) phi = acos(cos_phi);
    else phi = -1. * acos(cos_phi);      //Angolo rispetto a X
    *(tre_angoli+1) = phi;

    /* Ruoto attorno a Y di un angolo pari a phi per andare sull'asse X */
    matrix_Y(phi   ,matrix_roY);
    righe_per_colonne(v_rotz_1, matrix_roY, v_roty_1);
   // printf("v_roty_1[0]=%f v_roty_1[1]=%f v_roty_1[2]=%f \n",*v_roty_1,*(v_roty_1+1),*(v_roty_1+2));
    righe_per_colonne(v_rotz_2, matrix_roY, v_roty_2);
   // printf("v_roty_2[0]=%f v_roty_2[1]=%f v_roty_2[2]=%f \n",*v_roty_2,*(v_roty_2+1),*(v_roty_2+2));
    righe_per_colonne(v_rotz_3, matrix_roY, v_roty_3);
   // printf("v_roty_3[0]=%f v_roty_3[1]=%f v_roty_3[2]=%f \n",*v_roty_3,*(v_roty_3+1),*(v_roty_3+2));

    /* Angolo sul piano Y, Z del vettore che deve andaer su Y (Y -> Z rotazione positiva)*/
    cos_psi=*(v_roty_2+1) / sqrt( (*(v_roty_2+1) * *(v_roty_2+1) + *(v_roty_2+2) * *(v_roty_2+2)));
    sin_psi=*(v_roty_2+2) / sqrt( (*(v_roty_2+1) * *(v_roty_2+1) + *(v_roty_2+2) * *(v_roty_2+2)));

    if (sin_psi >= 0) psi = acos(cos_psi);
    else psi = -1. * acos(cos_psi);      //Angolo rispetto a Y
    *(tre_angoli+2) = psi;

    /* Ruoto attorno a X di un angolo opposto a psi per andare sull'asse Y */
    matrix_X(-1.*psi ,matrix_roX);
    righe_per_colonne(v_roty_1, matrix_roX, v_rotx_1);
   // printf("v_rotx_1[0]=%f v_rotx_1[1]=%f v_rotx_1[2]=%f \n",*v_rotx_1,*(v_rotx_1+1),*(v_rotx_1+2));
    righe_per_colonne(v_roty_2, matrix_roX, v_rotx_2);
   // printf("v_rotx_2[0]=%f v_rotx_2[1]=%f v_rotx_2[2]=%f \n",*v_rotx_2,*(v_rotx_2+1),*(v_rotx_2+2));
    righe_per_colonne(v_roty_3, matrix_roX, v_rotx_3);
 //   printf("v_rotx_3[0]=%f v_rotx_3[1]=%f v_rotx_3[2]=%f \n",*v_rotx_3,*(v_rotx_3+1),*(v_rotx_3+2));


	printf("\n%f %f %f\n",*(tre_angoli+0),*(tre_angoli+1),*(tre_angoli+2));


	free(v_rotz_1);
	free(v_rotz_2);
	free(v_rotz_3);
	free(v_roty_1);
	free(v_roty_2);
	free(v_roty_3);
	free(v_rotx_1);
	free(v_rotx_2);
	free(v_rotx_3);
    return ;
}

void  trans_bis_roto(double *v_orig, double *vect_out, double *tre_angoli_pre, double *tre_angoli_post) {

/* Ruoto autovettori pre ut orientarli su assi X,Y,Z */

   double cos_alpha, sin_alpha, alpha, cos_phi, sin_phi, phi, cos_psi, sin_psi,psi;
   double *v_rotz_a, *v_roty_a, *v_rotx_a;
   double *v_rotz_b, *v_roty_b, *v_rotx_b;
   double matrix_roZ[9]={0};
   double matrix_roX[9]={0};
   double matrix_roY[9]={0};

   v_rotz_a = calloc(3,sizeof(double));
   v_roty_a = calloc(3,sizeof(double));
   v_rotx_a = calloc(3,sizeof(double));
   v_rotz_b = calloc(3,sizeof(double));
   v_roty_b = calloc(3,sizeof(double));
   v_rotx_b = calloc(3,sizeof(double));

   *(v_orig+0) = *(v_orig+0);
   *(v_orig+1) = *(v_orig+1);
   *(v_orig+2) = *(v_orig+2);

    /* Angolo sul piano X, Y (X -> Y rotazione positiva) */
    alpha = *(tre_angoli_post+0);
    /* Ruoto attorno a Z di un angolo opposto per andare nel piano X-Z */
    matrix_Z(-1.*alpha  ,matrix_roZ);
    righe_per_colonne(v_orig, matrix_roZ, v_rotz_a);
 //   printf("v_rotz_1[0]=%f v_rotz_1[1]=%f v_rotz_1[2]=%f \n",*v_rotz_a,*(v_rotz_a+1),*(v_rotz_a+2));

    /* Angolo sul piano Z, X del vettore appena ruotato (Z -> X rotazione positiva)*/
    phi = *(tre_angoli_post+1);
    /* Ruoto attorno a Y di un angolo pari a phi per andare sull'asse X */
    matrix_Y(phi   ,matrix_roY);
    righe_per_colonne(v_rotz_a, matrix_roY, v_roty_a);
//    printf("v_roty_1[0]=%f v_roty_1[1]=%f v_roty_1[2]=%f \n",*v_roty_a,*(v_roty_a+1),*(v_roty_a+2));

    /* Angolo sul piano Y, Z del vettore che deve andaer su Y (Y -> Z rotazione positiva)*/
    psi = *(tre_angoli_post+2);
    /* Ruoto attorno a X di un angolo opposto a psi per andare sull'asse Y */
    matrix_X(-1.*psi ,matrix_roX);
    righe_per_colonne(v_roty_a, matrix_roX, v_rotx_a);
//    printf("v_rotx_1[0]=%f v_rotx_1[1]=%f v_rotx_1[2]=%f \n",*v_rotx_a,*(v_rotx_a+1),*(v_rotx_a+2));

/*******************************************************************/
/*  ORA GIRO gli angoli al contrario secondo la struttura orig *****/
/*******************************************************************/

    /* Angolo sul piano Y, Z del vettore che deve andaer su Y (Y -> Z rotazione positiva)*/
    psi =  *(tre_angoli_pre+2);
    /* Ruoto attorno a X di un angolo opposto a psi per andare sull'asse Y */
    matrix_X( 1.*psi ,matrix_roX);
    righe_per_colonne(v_rotx_a, matrix_roX, v_rotx_b);
//    printf("v_rotx_1[0]=%f v_rotx_1[1]=%f v_rotx_1[2]=%f \n",*v_rotx_b,*(v_rotx_b+1),*(v_rotx_b+2));

    /* Angolo sul piano Z, X del vettore appena ruotato (Z -> X rotazione positiva)*/
    phi = *(tre_angoli_pre+1);
    /* Ruoto attorno a Y di un angolo pari a phi per andare sull'asse X */
    matrix_Y(-1.*phi   ,matrix_roY);
    righe_per_colonne(v_rotx_b, matrix_roY, v_roty_b);
//    printf("v_roty_1[0]=%f v_roty_1[1]=%f v_roty_1[2]=%f \n",*v_roty_b,*(v_roty_b+1),*(v_roty_b+2));


    /* Angolo sul piano X, Y (X -> Y rotazione positiva) */
    alpha = *(tre_angoli_pre+0);
    /* Ruoto attorno a Z di un angolo opposto per andare nel piano X-Z */
    matrix_Z( 1.*alpha  ,matrix_roZ);
    righe_per_colonne(v_roty_b, matrix_roZ, v_rotz_b);
//    printf("v_rotz_1[0]=%f v_rotz_1[1]=%f v_rotz_1[2]=%f \n",*v_rotz_b,*(v_rotz_b+1),*(v_rotz_b+2));

	

    *(vect_out + 0) = *(v_rotz_b + 0);
    *(vect_out + 1) = *(v_rotz_b + 1);
    *(vect_out + 2) = *(v_rotz_b + 2);

   free(v_rotz_a); 
   free(v_roty_a); 
   free(v_rotx_a); 
   free(v_rotz_b); 
   free(v_roty_b); 
   free(v_rotx_b); 
    return ;
}

