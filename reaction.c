#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double T = 1.;
double tol = .5;
double dt = .005;



void updateF(double * z, double * f){
	f[0] = -z[0]*z[1];
	f[1] = z[0]*z[1]-z[1]*z[2];
	f[2] = z[1]*z[2]-z[2];
}

void updateFprime(double * z, double ** fprime){
	fprime[0][0] = -z[1];
	fprime[0][1] = -z[0];
	fprime[0][2] = 0.;
	fprime[1][0] = z[1];
	fprime[1][1] = z[0]-z[2];
	fprime[1][2] = -z[1];
	fprime[2][0] = 0.;
	fprime[2][1] = z[2];
	fprime[2][2] = z[1]-1;
}

void inverse33( double ** a , double ** inv ){
   int i;
   double det = 0.0; 
   for( i=0 ; i<3 ; ++i )
      det += a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]) ;
 
   int j;
   for(i=0;i<3;i++){
      for(j=0;j<3;j++)
           inv[j][i] = ((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/det ;
   }
}

void updateInv(double ** fprime, double ** inverse){
	double ** preInv = (double **) malloc(3*sizeof(double *));
	int k;
	for(k=0;k<3;k++){
		preInv[k] = (double *) malloc(3*sizeof(double));
	}
	int i,j;
	for(i = 0; i<3;i++){
		for(j = 0; j<3;j++){
			if(i != j){preInv[i][j] = (-1.)*dt*fprime[i][j];}
		}
		preInv[i][i] = 1 - dt*fprime[i][i];
	}

	inverse33(preInv, inverse);
	free(preInv);
}

void beuler(double * Z, double * F, double ** Fprime, double ** inverse, double time){
	
	int i;
	for(i=0;i<3;i++){
		Z[i] = Z[i] + dt*(inverse[i][1]*F[1]+inverse[i][2]*F[2]+inverse[i][0]*F[0]);
	}
	
}

void feuler(double * Z, double * F, double ** Fprime, double ** inverse, double time){
	
	int i;
	for(i=0;i<3;i++){
		Z[i] = Z[i] + dt*F[i];
	}
	
}

double get_stepFE(double * Z, double * F, double ** Fprime, double ** inverse, double time){
	double y1[3];
	double y2[3];
	int i;
	for(i=0;i<3;i++){
		y1[i]=Z[i];
		y2[i]=Z[i];
	}
	
	feuler(y1,F,Fprime,inverse,time);
	updateF(y1,F);
	feuler(y1,F,Fprime,inverse,time);
	updateF(Z,F);
	dt = dt*2.;
	feuler(y2,F,Fprime,inverse,time);
	dt = dt/2.;
	
	double dy = 0.0;
	for(i=0;i<3;i++) dy += fabs(y1[i]-y2[i]);
	
	double h_new = pow(tol/dy,.5)*dt;
	if(h_new>.01) h_new = .01;
	return(h_new);
}	

double get_stepBE(double * Z, double * F, double ** Fprime, double ** inverse, double time){
	double y1[3];
	double y2[3];
	int i;
	for(i=0;i<3;i++){
		y1[i]=Z[i];
		y2[i]=Z[i];
	}
	
	beuler(y1,F,Fprime,inverse,time);
	updateF(y1,F);
	updateFprime(y1, Fprime);
	updateInv(Fprime, inverse);
	beuler(y1,F,Fprime,inverse,time);
	updateF(Z,F);
	updateFprime(Z, Fprime);
	updateInv(Fprime, inverse);
	dt = dt*2.;
	beuler(y2,F,Fprime,inverse,time);
	dt = dt/2.;
	
	double dy = 0.0;
	for(i=0;i<3;i++) dy += fabs(y1[i]-y2[i]);
	
	double h_new = pow(tol/dy,.5)*dt;
	if(h_new>.01) h_new = .01;
	return(h_new);
}
	
	

void main(){

	double time = 0.;
	double * Z = (double *) malloc(3*sizeof(double));
	double * F = (double *) malloc(3*sizeof(double));
	double ** Fprime = (double **) malloc(3*sizeof(double *));
	double ** inverse = (double **) malloc(3*sizeof(double *));
	int j;
	for(j=0;j<3;j++){
		Fprime[j] = (double *) malloc(3*sizeof(double));
		inverse[j] = (double *) malloc(3*sizeof(double));
	}
	
	Z[0] = 150.;
	Z[1] = 10.;
	Z[2] = 10.;
	
	updateF(Z,F);
	updateFprime(Z, Fprime);
	updateInv(Fprime, inverse);
	
	//Forward Euler
	
	FILE * feData = fopen("feuler.txt","w");
	fprintf(feData,"%15.15f    %15.15f     %15.15f        %15.15f\n", Z[0],Z[1],Z[2],time);
	int countFE = 0;
	while(time<T){
		fprintf(feData,"%15.15f    %15.15f     %15.15f        %15.15f\n", Z[0],Z[1],Z[2],time);
		dt = get_stepFE(Z,F,Fprime,inverse,time);
		feuler(Z,F,Fprime,inverse,time);
		updateF(Z,F);
		time += dt;
		countFE++;
	}
	fclose(feData);
	
	//Reinitialize starting conditions to do BE
	Z[0] = 150.;
	Z[1] = 10.;
	Z[2] = 10.;
	time = 0;
	
	updateF(Z,F);
	updateFprime(Z, Fprime);
	updateInv(Fprime, inverse);
	
	//Backward Euler
	FILE * beData = fopen("beuler.txt","w");
	fprintf(beData,"%15.15f    %15.15f     %15.15f        %15.15f\n", Z[0],Z[1],Z[2],time);
	int countBE = 0;
	while(time<T){
		fprintf(beData,"%15.15f    %15.15f     %15.15f        %15.15f\n", Z[0],Z[1],Z[2],time);
		dt = get_stepBE(Z,F,Fprime,inverse,time);
		beuler(Z,F,Fprime,inverse,time);
		updateF(Z,F);
		updateFprime(Z, Fprime);
		updateInv(Fprime, inverse);
		time += dt;
		countBE++;
	}
	fclose(beData);
	
	printf("%d\n",countFE);
	printf("%d\n",countBE);
	
	//Freeing all the mallocs
	for(j=0;j<3;j++){
		free(Fprime[j]);
		free(inverse[j]);
	}
	free(Z);
	free(F);
	free(Fprime);
	free(inverse);
}
