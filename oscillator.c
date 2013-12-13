#include <stdio.h>
#include <math.h>

//double gamma = 1.;
//double omega = 30.;
double T = 10.;
int N = 100000;
double dt;
double tol = .0001;

double dxdt(double x,double v,double t){
	return v;
}

double dvdt(double x,double v,double t){
//I commented out the more general way because it was giving an inexplicable error.
	//double a = -gamma*v-pow(omega,2)*x;
	double a = -v-900*x;
	return a;
}

void feuler(double * x,double * v,double * t){
	double xnew = *x+dxdt(*x,*v,*t)*dt;
	double vnew = *v+dvdt(*x,*v,*t)*dt;
	*x = xnew;
	*v = vnew;
	*t = *t+dt;
}

void midpoint(double * x,double * v,double * t){
	double xmid = *x+dxdt(*x,*v,*t)*dt/2;
	double vmid = *v+dvdt(*x,*v,*t)*dt/2;
	double xnew = *x+dxdt(xmid,vmid,*t)*dt;
	double vnew = *v+dvdt(xmid,vmid,*t)*dt;
	*x = xnew;
	*v = vnew;
	*t = *t+dt;
}

void beuler(double * x,double * v,double * t){
	double xnew = *x+dxdt(*x,*v,*t)*dt;
	double vnew = *v+dvdt(*x,*v,*t)*dt;
	double tnew = *t+dt;
	
	double errorx = xnew-*x-dxdt(xnew,vnew,tnew)*dt;
	double errorv = vnew-*v-dvdt(xnew,vnew,tnew)*dt;
	while(fabs(errorx)>tol){
		xnew = xnew-errorx;
		errorx = xnew-*x-dxdt(xnew,vnew,tnew)*dt;
	}
	while(fabs(errorv)>tol){
		vnew = vnew-errorv;
		errorv = vnew-*v-dvdt(xnew,vnew,tnew)*dt;
	}
	
	
	*x = xnew;
	*v = vnew;
	*t = tnew;
}

void cn(double * x,double * v,double * t){
	dt = dt/2;
	
	feuler(x,v,t);
	
	beuler(x,v,t);
	
	dt = dt*2;
}


void main(){
	double x = 0;
	double v = 1;
	double t = 0;
	dt = T/N;
	
	FILE* tcn = fopen("cn.txt","w");
	int i;
	for(i=0; i<N; i++){
		fprintf(tcn,"%f    %f    %f\n", x,v,t);
		cn(&x,&v,&t);
	}
	fclose(tcn);
	
	FILE* tfeuler = fopen("feuler.txt","w");
	x = 0;
	v = 1;
	t = 0;
	for(i=0; i<N; i++){
		fprintf(tfeuler,"%f    %f    %f\n", x,v,t);
		feuler(&x,&v,&t);
	}
	fclose(tfeuler);
	
	FILE* tbeuler = fopen("beuler.txt","w");
	x = 0;
	v = 1;
	t = 0;
	for(i=0; i<N; i++){
		fprintf(tbeuler,"%f    %f    %f\n", x,v,t);
		beuler(&x,&v,&t);
	}
	fclose(tbeuler);
	
	FILE* tmidpoint = fopen("midpoint.txt","w");
	x = 0;
	v = 1;
	t = 0;
	for(i=0; i<N; i++){
		fprintf(tmidpoint,"%f    %f    %f\n", x,v,t);
		midpoint(&x,&v,&t);
	}
	
	fclose(tmidpoint);
}
