#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
const int E=100, UN=400, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel;

//void ecfcal(int,double[]);

int  N, trim, check, pf, nxk, nyk, nzk, fc, dumb, mum, tum,  gord, tric, drum, sj, jdt, dum, lc, l, m, n, xk, yk, zk, ip, nlies, mo, no, wbuq, titi, td, rdir, pk, i, j, k, chum, gum, jt, mt, it, nf, block1, block2, block3, block[3], blockd, jtt, ri, ittt, lsd, sitt, iltt, wgdum, ad1, ad2, ad3, wt, dr, cont, dummy, kt, tri, puc, ns, ie, trick, sini, ini, it, cow, jet, lo, loo, ix, iy, iz, lct, rul, sla, veei, veec, SHA, ptr, ptrr, bn, bnll, sst, da, sitts,  rpid1, rpid2, rpid, ntt;

double xi, yi, zi, xf, yf, zf, dx, dy, dz, mxi, myi, mzi, dxi, dyi, dzi, u1, u2, u3, u4, u5, u6, u7, u8, u9, v1, v2, v3, w1, w2, w3, w4, w5, w6, w7, w8, w9, re1, re2, re3, r1, r2, r3, yx1, yx2, yx3, xx1, xx2, xx3, mr, mrp, mr1, zt1, zt2, dzt1, dzt2, theta, theta1, theta2, phi, rdphi, q1, q2, q3, q4, q5, q6, q7, q11, q12, q13, q21, q22, q23, q31, q32, q33, pdx, pdy, pdz, x1per, x2per, y1per, y2per, z1per, z2per, ra1, ra2, cdtp1, cdtp2, cdtp, dtp1, dtp2, dtp, dom0, dom1, dom2, cp, cph, cpsh, lt[3][4], lts[3][3],  dump[4], rs,  sol[3], del[9], dli[9], rdom[3], rdli[9], f1, f2, f3, f[3], bl1, bl2, bl3, bl[3], ar, per; 

double step, omega, red, boss,  com, kd, qua, quao, quai, quaoi, rijs, xt, yt, lem, lemm, lemms, lems, lembda, cm, ep, ec, di, dive, rd, rdd, beta, duc,  tpar, sq1per, sq2per, tpart,  tper,  ar, ivf, fvf, dumm, rdnp, asp, piee, fe, tpiee, dis, rtm[3][3], abrtm[3][3], ra, rb, ree[3], del1[9], del2[9], a, b, c, dd, trig, s[2], cfun[3], EPSILON, eb, teb, tbl, pbc,  dif, difu, difus, difl, difusl, adefus, adefusl, ble, vees, rst, rsmt, rdth, istep, prtd, uti, bdm;  



int main(void)
{

	FILE* bu;
	bu=fopen("Base_data.txt", "r");
	
	printf("\n\n\n                                ******:::::::   Prameters given to the system by using Bash File   :::::::******\n\n\n");
	fscanf(bu,"%lf",&bdm);
	N=(int)round(bdm);
	//if (N==0) printf("raam tere desh men \n");
	// Definition of universal constants starts!
	
	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;
	
	// Definition of universal constants ends!
	
	
	
	// Definition of global constants starts!
	
	//step=0.01;
	//fscanf(bu,"%lf",&step);
	//fscanf(bu,"%ld",&tsim);
	fscanf(bu,"%lf",&bdm);
	step=bdm;
	fscanf(bu,"%lf",&bdm);
	tsim=(int)round(bdm);
	
	//tsim=100000;
	ivf=0.01;
	
	// Definition of global constants ends       !!!!!!!!!!!!!!!!
	
					
	// Correlation function        !!!!!!!!!!!!!!!!
	
	
	
	
	// Correlation function declaration end         !!!!!!!!!!!!!!!!			
	
	
	
	
  	// Definition of simulation constants starts        !!!!!!!!!!!!!!!!
	
	
	//fvf= 0.02 ;
	fscanf(bu,"%lf",&bdm);
	fvf=bdm;
	//fscanf(bu,"%lf",&fvf);
	
	int war;
	
	//fscanf(bu,"%d",&war);
	fscanf(bu,"%lf",&bdm);
	war=(int)round(bdm);
	
	//war=3;
	
	 
	double* ar=malloc(war*sizeof(double));
	
	double* per=malloc(war*sizeof(double));
	
	double* sp=malloc(war*sizeof(double));
	
	double* spo=malloc(war*sizeof(double));
	
	int* nt=malloc(war*sizeof(int));
	
	int* sco=malloc(war*sizeof(int));
	
	
	int** pop=malloc(war*sizeof(int*));
	for (i=0; i<war; i++)
	{
		pop[i]=malloc(war*sizeof(int));
	}
	
	// printf("raam tere desh men \n");
	
	
	
	
	for (i=0; i<war; i++)
	{
		sco[i]=0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			sco[i]=(int)round(bdm);
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		ar[i]=1.0;
	}
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			ar[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	//per[0]=1.0;
	for (i=0; i<war-1; i++)
	{
		per[i]=0.0;
	}
	
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war-1; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			per[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	red=0.0;
	for (i=0; i<war-1; i++)
	{
		red=red+per[i];
	}
	
	per[war-1]=1.0-red;
	
	
	for (i=0; i<war  ; i++)
	{
		nt[i]=(int)round((float)N*per[i]);
		//printf("%d  \n", nt[i]);
	}
	
	
	//Last fraction        !!!!!!!

	
	
	
	// Shape parameters 0 and < 1.0 ; Oblate and 1 and > 1.0 ; Prolate          !!!!!!!!!!!!!! 
	
	SHA = 1 ;
	
	
	
	
	// Definition of simulation constants ends        !!!!!!!!!!!
	
	
	
	// Definition of compression parameter starts          !!!!!!!!!!
	
	veec=200;
	vees=0.005;    // Compression steps!
	vel=(int)(10000.0*(fvf/0.5));
	//vel=10;
	
	
	// Definition of compression parameter ends            !!!!!!!!!!
	
	

	// Declaration of core bond variable starts          !!!!!!!!!!
	
	double epsi;
	
	// Declaration of core bond variable ends         !!!!!!!!!!
	
	
	
	// Evaluation of bond variable starts           !!!!!!!!!
	
	
	
	bn=40;
	bnll=40;
	sst=1;
	
	
	//double epsi1, epsi2, epsi3, omega2, beta1, beta2, beta3, betai;
	
	// Declaration of core bond variable ends!
	
	
	// Evaluation of bond variable starts!
	

	//pop =  0 ;              // 0: hard core;   1: HC+Jenus;    2: HC+Patchy;    3: HC+Isotropic;     4: HC+jenus+Isotropic;   5: HC+patchy+Isotropic; 

	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			pop[i][j]=0;
			pop[j][i]=0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				pop[i][j]=(int)round(bdm);
				pop[j][i]=(int)round(bdm);
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}



	
	
	// Potential properties define here!
	
	
	
	
	
	
	
	
	
	double* poten=malloc((((war*(war-1))/2)+war)*sizeof(double*));
	
	// Jenus potential properties;
	
	double** epsi1=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi1[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta1=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta1[i]=malloc(war*sizeof(double)); 
	
	}
	
	
	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi1[i][j]=0.0;
			epsi1[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi1[i][j]=bdm;
				epsi1[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta1[i][j]=0.0;
			beta1[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta1[i][j]=bdm;
				beta1[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	
	
	
	
	// Patchy potential properties;
	
	double** epsi2=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi2[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta2=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta2[i]=malloc(war*sizeof(double)); 
	
	}
	
	double* omega2=malloc(war*sizeof(double));
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi2[i][j]=0.0;
			epsi2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi2[i][j]=bdm;
				epsi2[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta2[i][j]=0.0;
			beta2[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta2[i][j]=bdm;
				beta2[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	for (i=0; i<war; i++)
	{
		omega2[i]=0.0;
	}
	
	fscanf(bu,"%lf",&bdm);
	for (i=0; i<war; i++)
	{
		fscanf(bu,"%lf",&bdm);
		if (bdm==50.0) 
		{
			break;
		}
		else
		{
			omega2[i]=bdm;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	// Isotropic potential properties;
	
	double** epsi3=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		epsi3[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** beta3=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta3[i]=malloc(war*sizeof(double)); 
	
	}
	

	
	
	jt=0;
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			epsi3[i][j]=0.0;
			epsi3[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				epsi3[i][j]=bdm;
				epsi3[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			beta3[i][j]=0.0;
			beta3[j][i]=0.0;
		}
	}
	
	
	fscanf(bu,"%lf",&bdm);
	for(i=0; i<war; i++)
	{
		for(j=i; j<war; j++)
		{
			fscanf(bu,"%lf",&bdm);
			if (bdm==50.0) 
			{
				break;
			}
			else
			{
				beta3[i][j]=bdm;
				beta3[j][i]=bdm;
			}
			
		}
		if (bdm==50.0) 
		{
			break;
		}
	}
	for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	

	fclose(bu);
	remove("Base_data.txt");



	
	
	// Declaired!
	
	double** beta=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		beta[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** betai=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		betai[i]=malloc(war*sizeof(double)); 
	
	}
	
	epsi=0.0;
	for(i=0; i<war; i++)
	{
		for(j=0; j<war; j++)
		{
			
	
			if (pop[i][j]==0)
			{
				epsi=0.0;
			}
			
			else if (pop[i][j]==1)
			{
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi1[i][j]>epsi) epsi=epsi1[i][j];
						beta[i][j]=beta1[i][j];
						beta[j][i]=beta1[i][j];
						
					}
				}
			}
			else if (pop[i][j]==2)
			{
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi2[i][j]>epsi) epsi=epsi2[i][j];
						beta[i][j]=beta2[i][j];
						beta[j][i]=beta2[i][j];
						
					}
				}
			}
			else if (pop[i][j]==3)
			{
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
				
				
			}
			else if (pop[i][j]==4)
			{
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi1[i][j]>epsi) epsi=epsi1[i][j];
						beta[i][j]=beta1[i][j];
						beta[j][i]=beta1[i][j];
						
					}
				}
				
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
				
				
			}	
			else if (pop[i][j]==5)
			{
				
				
				//epsi=0.0;
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi2[i][j]>epsi) epsi=epsi2[i][j];
						beta[i][j]=beta2[i][j];
						beta[j][i]=beta2[i][j];
						
					}
				}
				
				for(i=0; i<war; i++)
				{
					for(j=i; j<war; j++)
					{
						if (epsi3[i][j]>epsi) epsi=epsi3[i][j];
						betai[i][j]=beta3[i][j];
						betai[j][i]=beta3[i][j];
						
					}
				}
			}
			
		}
	}
	
	
	
	
	
	
	// Evaluation of bond variable ends!
	
	printf("Number of particles-    %d \n", N);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Step size-   %lf \n", step);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of simulation step-   %ld \n", tsim);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Starting volume fraction-   %lf\n", ivf);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Volume fraction-    %lf \n", fvf);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Type of particles-    %d \n", war);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Shape factor (aspect ratios) are given by \n");
	for (i=0; i<war; i++)
	{
		printf("Aspect ratio for particle type %d is given by-   %lf\n", i+1, ar[i]);
	}
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Particles fractions are given by \n");
	for (i=0; i<war; i++)
	{
		printf("Fraction for particle type %d is given by-   %lf\n", i+1, per[i]);
		
	}
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Number of particles for each type (aspect ratios) are given by \n");
	for (i=0; i<war; i++)
	{
		printf("Number of particles for particle type %d is given by-   %d\n", i+1, nt[i]);
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("Potential parameters are given by:- \n");
	printf("\n");
			 printf("\n");
			 printf("\n");
	
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			 printf("Kind of potential between family of particles %d and %d is defined as-   ", i,j);
			 if (pop[i][j]==0)
			 {
			 	printf("Only hard core.\n");
			 	printf("Epsi-   Not applicable.\n");
			 	printf("Beta-   Not applicable.\n");
			 	printf("Omega-   Not applicable.\n");
			 	
			 }
			 else if (pop[i][j]==1)
			 {
			 	printf("Hard core + Jenus.\n");
			 	printf("Epsi-   %lf \n", epsi1[i][j]);
			 	printf("Beta-   %lf \n", beta1[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	 
			 
			  else if (pop[i][j]==2)
			 {
			 	printf("Hard core + Patchy.\n");
			 	printf("Epsi-   %lf \n", epsi2[i][j]);
			 	printf("Beta-   %lf \n", beta2[i][j]);
			 	printf("Omega for particle %d, for the patchy potential defined between type of particle %d and %d-   %lf \n", i, i, j, omega2[i]);
			 }	
			 
			  else if (pop[i][j]==3)
			 {
			 	printf("Hard core + Jenus.\n");
			 	printf("Epsi-   %lf \n", epsi3[i][j]);
			 	printf("Beta-   %lf \n", beta3[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	
			  else if (pop[i][j]==4)
			 {
			 	printf("Hard core + Jenus + Isotropic.\n");
			 	printf("Epsi and Epsi-isotropic-   %lf  %lf \n", epsi1[i][j], epsi3[i][j]);
			 	printf("Beta and Beta-isotropic-   %lf  %lf \n", beta1[i][j], beta3[i][j]);
			 	printf("Omega-   Not applicable\n");
			 }	
			 
			 else if (pop[i][j]==5)
			 {
			 	printf("Hard core + Jenus + Isotropic.\n");
			 	printf("Epsi and Epsi-isotropic-   %lf  %lf \n", epsi2[i][j], epsi3[i][j]);
			 	printf("Beta-\n");
			 	printf("Beta and Beta-isotropic-   %lf  %lf \n", beta2[i][j], beta3[i][j]);
			 	printf("Omega for particle %d, for the patchy potential defined between type of particle %d and %d-   %lf \n", i, i, j, omega2[i]);
			 }	
			 printf("\n");
			 printf("\n");
			 printf("\n");
		 
		}
	
	}
	printf("Values corresponding to the particular parameter is given by\n");
	
	
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	/*
	printf("Applicable parameters is given as:   \n");
	
	printf("Epsi %lf \n", epsi);
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			if (pop[i][j]>0)
			{
				printf("Applied beta is given between particles %d and %d is given by- %lf   \n", i, j, beta[i][j]);
			
				if (pop[i][j]>2) printf("Applied betai is given between particles %d and %d is given by- %lf   \n", i, j, betai[i][j]);
			}
			else
			{
				printf("Not applicable!!!\n");
			}
		}
		
	}
	// Printing things starts now !!!!!!!
	
	*/
	
	
	
	
	
	
	
	
	
	
	// Printing things ends now !!!!!!!
	
	// Declaration of simulation variables starts!
	
	
	// Declaration of simulation variables ends!
	
	
	
	
	// Evaluation of simulation variables starts! 
	
	double* d=malloc(war*sizeof(double));
	
	double* r=malloc(war*sizeof(double));
	
	double* rsm=malloc(war*sizeof(double));
	
	double* rm=malloc(war*sizeof(double));
	
	double* rs=malloc(war*sizeof(double));
	
	double** e=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		e[i]=malloc(3*sizeof(double)); 
	
	}
	
	double** ds=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		ds[i]=malloc(war*sizeof(double)); 
	
	}
	
	double** dss=malloc(war*sizeof(double*));
	for(i=0; i<war; i++)
	{
		dss[i]=malloc(war*sizeof(double)); 
	
	}
	
	for (i=0; i<war; i++)
	{
		if (sco[i]==0)
		{
			d[i]=pow(ar[i],(2.0/3.0));
			r[i]=d[i]/2.0;
			rs[i]=r[i]*r[i];
			rsm[i]=rs[i]/(ar[i]*ar[i]);
			rm[i]=r[i]/ar[i];
			
			e[i][0]=rm[i]+(epsi/2.0);
			e[i][1]=rm[i]+(epsi/2.0);
			e[i][2]=r[i]+(epsi/2.0);
		}
		else
		{
			if (ar[i]>1.0)
			{
				d[i]=1.0;
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=rm[i]+(epsi/2.0);
				e[i][1]=rm[i]+(epsi/2.0);
				e[i][2]=r[i]+(epsi/2.0);
			}
			else
			{
				d[i]=ar[i];
				r[i]=d[i]/2.0;
				rs[i]=r[i]*r[i];
				rsm[i]=rs[i]/(ar[i]*ar[i]);
				rm[i]=r[i]/ar[i];
				
				e[i][0]=rm[i]+(epsi/2.0);
				e[i][1]=rm[i]+(epsi/2.0);
				e[i][2]=r[i]+(epsi/2.0);
			
			}
		}
		
	}	
	for (i=0; i<war; i++)
	{
		for (j=0; j<war; j++)
		{
			if (rm[i]<r[i]) 
			{ 
				ds[i][j]=rm[i];
				dss[i][j]=r[i];
			}
			else
			{
				ds[i][j]=r[i];
				dss[i][j]=rm[i];
			}
			
			if (rm[j]<r[j]) 
			{ 
				ds[i][j]=ds[i][j]+rm[j];
				dss[i][j]=dss[i][j]+r[j];
			}
			else
			{
				ds[i][j]=ds[i][j]+r[j];
				dss[i][j]=dss[i][j]+rm[j];
			}
			
			dss[i][j]=dss[i][j]+epsi;
			
			dss[i][j]=dss[i][j]*dss[i][j];
			ds[i][j]=ds[i][j]*ds[i][j];
			
			
		}
		
		
	}
	//printf("%lf\n", ds[1][0]);
	fe=0.0;
	for (i=0; i<war; i++)
	{
		fe=fe+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*fvf));
	}
	fe=pow(fe,(1.0/3.0));
	
	f1=0.0;
	
	for (i=0; i<war; i++)
	{
		f1=f1+((4.0*piee*r[i]*r[i]*r[i]*(float)nt[i])/(3.0*ar[i]*ar[i]*ivf));
	}
	
	f1=pow(f1,(1.0/3.0));
	f2=f1;
	f3=f1;
	
	ble=1.0;
	for (i=0; i<war; i++)
	{
		if (rm[i]<r[i])
		{
			if (d[i]>ble) ble=d[i];
		}
		else
		{
			if ((d[i]/ar[i])>ble) ble=d[i]/ar[i];
		}
	}
	
	ble=ble+epsi;
	
	
	
	
	
	
	bl1=ble;
	bl2=ble;
	bl3=ble;
	
	block1=(int)floor(f1/bl1);
	bl1=f1/(float)block1;
	f1=bl1*(float)block1;
	
	block2=(int)floor(f2/bl2);
	bl2=f2/(float)block2;
	f2=bl2*(float)block2;
	
	block3=(int)floor(f3/bl3);
	bl3=f3/(float)block3;
	f3=bl3*(float)block3;
	
	f[0]=f1;
	f[1]=f2;
	f[2]=f3;
	bl[0]=bl1;
	bl[1]=bl2;
	bl[2]=bl3;
	block[0]=block1;
	block[1]=block2;
	block[2]=block3;
	blockd=block1+20;
	
	
	printf("                                ******:::::::   Calculated parameters are given by   :::::::******\n");
	printf("\n");
	printf("Starting edge is given by-  %lf\n", f1);
	printf("\n");
	printf("\n");
	printf("Aimed edge is given by-  %lf\n", fe);
	printf("\n");
	printf("\n");
	printf("Starting number of shell at each edge is given by-  %d\n", block[0]);
	printf("\n");
	printf("\n");
	printf("Shell length is given by-   %lf \n\n\n", ble);
	
	pbc=ble+0.5;
	
	for (i=0; i<war; i++)
	{
		printf("For family of particle %d-\n", i);
		printf("Diameter-  %lf\n", d[i]);
		printf("Length of semi-symmetry-axis-  %lf\n", r[i]);
		printf("Length of semi-axis perpendicular to the symmetry-axis-  %lf \n \n \n \n", rm[i]);
	
	}

	// Evaluation of simulation variables ends!
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Declaration of core variables starts!
	
	int pit[400], wpi;
	
	
	double** li=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		li[i]=malloc(9*sizeof(double));
	}
	
	double** fli=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		fli[i]=malloc(9*sizeof(double));
	}
	
	double** el=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		el[i]=malloc(9*sizeof(double));
	}
	
	double** ami=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		ami[i]=malloc(3*sizeof(double));
	}
	int* pid=malloc(N*sizeof(int));
	
	double** fau=malloc(N*sizeof(double*));
	
	for (i=0;i<N;++i)
	{
		fau[i]=malloc(3*sizeof(double));
	}
	
	
	
	
	double ldcxyz[N][3], ldcxy[N][2], ldcz[N], lom[N][3], mdcxyz[N][3], mdcxy[N][2], mdcz[N], mom[N][3], sdcxyz[N][3], sdcxy[N][2], sdcz[N], som[N][3];
	
	
	
	int** nrdi=malloc(N*sizeof(int*));
	int* wnrd=malloc(N*sizeof(int));
	
	for (i=0; i<N; i++)
	{
		nrdi[i]=malloc(UN*sizeof(int));
	}
	
	
	
	int**** dri=malloc(blockd*sizeof(int***));
	int*** wdr=malloc(blockd*sizeof(int**));
		
	for (j=0; j<blockd; j++)
	{
		dri[j]=malloc(blockd*sizeof(int**));
		wdr[j]=malloc(blockd*sizeof(int*));
		for (i=0; i<blockd; i++)
		{
			dri[j][i]=malloc(blockd*sizeof(int*));
			wdr[j][i]=malloc(blockd*sizeof(int));
	
			for (k=0; k<blockd; k++)
			{
				dri[j][i][k]=malloc(80*sizeof(int));
			}
		}	
	}
	
	// Declaration of compresser variables starts!
	
	int vee;
	
	// Declaration of compresser variables ends!
	
	// Declaration of core variables ends!
	
	
	
	
	
	
	// Declaration of bonding variable starts!
	
	
	
	
		int** bobb=malloc(N*sizeof(int*));
		int* wbobb=malloc(N*sizeof(int));
		for (k=0; k<N; k++)
		{
			bobb[k]=malloc(E*sizeof(int));
		
		}	
			
		
		
		int** cbobb=malloc(N*sizeof(int*));
		int* wcbobb=malloc(N*sizeof(int));
		for (k=0; k<N; k++)
		{
			cbobb[k]=malloc(E*sizeof(int));
		
		}	
			
		
		
		int wpbobb;
		int* psbobb=malloc(E*sizeof(int));

		int hold[N][E], whold[N], bumbque[E];
		
	
		
		
		
		
		
		
		
		int** bobbi=malloc(N*sizeof(int*));
		int* wbobbi=malloc(N*sizeof(int));
		for (k=0; k<N; k++)
		{
			bobbi[k]=malloc(E*sizeof(int));
		
		}	
			
		
		
		int** cbobbi=malloc(N*sizeof(int*));
		int* wcbobbi=malloc(N*sizeof(int));
		for (k=0; k<N; k++)
		{
			cbobbi[k]=malloc(E*sizeof(int));
		
		}	
			
		
		
		int wpbobbi;
		int* psbobbi=malloc(E*sizeof(int));

		int holdi[N][E], wholdi[N], bumbquei[E], wbuqi;
		
		
		
	

	
	// Declaration of bonding variable ends! 
	
	
	
	
	
	
	
	// Initialization of core variables start!
	
	
	adefus=0.0;
	adefusl=0.0;
	for (i=0; i<N; i++)
	{
		
		sdcxyz[i][0]=0.0;
		sdcxyz[i][1]=0.0;
		sdcxyz[i][2]=0.0;
		sdcxy[i][0]=0.0;
		sdcxy[i][1]=0.0;
		sdcz[i]=0.0;
		som[i][0]=0.0;
		som[i][1]=0.0;
		som[i][2]=0.0;
		
		mdcxyz[i][0]=0.0;
		mdcxyz[i][1]=0.0;
		mdcxyz[i][2]=0.0;
		mdcxy[i][0]=0.0;
		mdcxy[i][1]=0.0;
		mdcz[i]=0.0;
		mom[i][0]=0.0;
		mom[i][1]=0.0;
		mom[i][2]=0.0;
		
		ldcxyz[i][0]=0.0;
		ldcxyz[i][1]=0.0;
		ldcxyz[i][2]=0.0;
		ldcxy[i][0]=0.0;
		ldcxy[i][1]=0.0;
		ldcz[i]=0.0;
		lom[i][0]=0.0;
		lom[i][1]=0.0;
		lom[i][2]=0.0;
		
		
		wnrd[i]=0;
		for (j=0; j<UN; j++)
		{
			nrdi[i][j]=6000;
		
		}			
		
	}
	
	for (i=0; i<blockd; i++)
	{
		
		for (j=0; j<blockd; j++)
		{
			for (k=0; k<blockd; k++)
			{
				wdr[i][j][k]=0;
				for (l=0; l<80; l++)
				{
					dri[i][j][k][l]=6000;
				}
			}	
		}
	}
	
	

	
		
	
	
	
	for (i=0; i<N; i++)
	{
		jt=0;
		jtt=0;
		for (j=0; j<war; j++)
		{
			jtt=jtt+(int)(per[j]*(float)N);
			if (i>jtt)    jt=jt+1;		
		}
		pid[i]=jt;
	}
	
	
	
	
	double* rdtheta=malloc(war*sizeof(double));
	double* grxy=malloc(war*sizeof(double));
	double* grz=malloc(war*sizeof(double));
	double* grth=malloc(war*sizeof(double));
	
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Dynamical parameters are given as- \n");
	
	for (jt=0; jt<war; jt++)
	{
		if(rm[jt]<r[jt])
		{
			asp=ar[jt];
			if (asp==1.0)
			{
				
				grz[jt]=1.0;
				grxy[jt]=1.0;
				grth[jt]=1.0;
			
			}
			else
			{
				q1=asp+pow((pow(asp,2.0)-1.0),0.5);
				q2=asp-pow((pow(asp,2.0)-1.0),0.5);
				q3=pow(asp,2.0)-1.0;
				q4=pow(q3,(3.0/2.0));
				
				q5=1.0-(asp*asp);
				q6=2.0*(asp*asp);
				q7=q6-1.0;
				grz[jt]=(8.0/3.0)*(1.0/(((2.0*asp)/q5)+((q7/q4)*log(q1/q2))));
				grxy[jt]=(8.0/3.0)*(1.0/((asp/q3)+(((q6-3.0)/q4)*log(q1))));
				grth[jt]=(2.0/3.0)*((pow(asp,4.0)-1.0)/(asp*(((q7/pow(q3,0.5))*log(q1))-asp)));
				
				q11=rm[jt]/0.5;
				grz[jt]=pow((1.0/(grz[jt]*q11)),0.5);
				grxy[jt]=pow((1.0/(grxy[jt]*q11)),0.5);
				grth[jt]=pow((1.0/grth[jt]),0.5);
			}
			
		}
		else
		{
			asp=ar[jt];
			if (asp==1.0)
			{
				grz[jt]=1.0;
				grxy[jt]=1.0;
				grth[jt]=1.0;
			
			}
			else
			{
				q1=1.0-(2.0*asp*asp);
				q2=1.0-(asp*asp);
				q3=pow(q2,0.5);
				q4=3.0-(2.0*asp*asp);
				q5=-q1;
				q6=(1.0/q3)*acos(asp);
				grz[jt]=(6.0/(8.0*q2))*(((q1*acos(asp))/q3)+asp);
				grxy[jt]=(3.0/(8.0*q2))*(((q4*acos(asp))/q3)-asp);
				grth[jt]=(3.0/2.0)*((((q5*q6)-asp)*asp)/(pow(asp,4.0)-1.0));
				
				q11=rm[jt]/0.5;
				grz[jt]=pow((grz[jt]/q11),0.5);
				grxy[jt]=pow((grxy[jt]/q11),0.5);
				grth[jt]=pow(grth[jt],0.5);
			}
			
		}
		printf("For family of particles %d:     grz= %lf,  grxy=%lf, and  grth=%lf\n", jt, grz[jt], grxy[jt], grth[jt]); 
	}
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	for (jt=0; jt<war; jt++)
	{
		
		rdtheta[jt]=(pow(2.0,0.5))*step*grth[jt];
		printf("For family of particles %d:     S_xy = %lf,  S_z = %lf and  S_r = %lf\n", jt, grxy[jt]*step, grz[jt]*step, rdtheta[jt]);
		
	}
	
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	// Initialization of core variables ends!
	
	
	
	
	
	
	
	
	
	// Initialization of bond variables starts!
	
	
	for (i=0; i<N; i++)
	{
		wbobb[i]=0;
		wcbobb[i]=0;
		for (j=0; j<E; j++)
		{
			bobb[i][j]=6000;
		}
		for (j=0; j<E; j++)
		{
			cbobb[i][j]=6000;
		}				
	}
		
	
	
	
	
	
	wpbobb=0;
	for(i=0; i<E; i++)
	{
		psbobb[i]=6000;
	}
		
	

	for (i=0; i<N; i++)
	{
		wbobbi[i]=0;
		wcbobbi[i]=0;
		for (j=0; j<E; j++)
		{
			bobbi[i][j]=6000;
		}
		for (j=0; j<E; j++)
		{
			cbobbi[i][j]=6000;
		}				
	}
		
	
	
	
	
	
	wpbobbi=0;
	for(i=0; i<E; i++)
	{
		psbobbi[i]=6000;
	}
	
	
	
	
	
	
	// Initialization of bond variables ends!
	
	
	
	
	
	
	
	// Setting initial coordinates and direction starts!
	
	ix=block1;
	iy=block2;
	iz=block3;
	
	FILE* cfu;
	
	FILE* vu;
	vu=fopen("Relaxed_system.txt","r");
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	
	if (vu==NULL)
	{
		printf("Initial configuration within the simulation box have been started generating: \n \n");
		for(i=0; i<N; i++)	
	    	{	
			//di[i]=rand() % 100002;
			/*if (SHA==1)
			{
				li[i][6]=1.0;
				li[i][7]=0.0;
				li[i][8]=0.0;
				
			}
			else
			{
				li[i][6]=0.0;
				li[i][7]=0.0;
				li[i][8]=1.0;
			}*/
			
			
			ra1=drand48();
	        	theta=acos(1.0-(2.0*ra1));
			ra2=drand48();
	       	phi=2.0*piee*ra2;
			
			
			pdx=sin(theta)*cos(phi);
			pdy=sin(theta)*sin(phi);
			pdz=cos(theta);
			
			li[i][6]=pdx;
			li[i][7]=pdy;
			li[i][8]=pdz;
			
			cow=0;
			while (cow==0)	
			{	
				
				xi=((float)(rand() % ix))*bl[0]+(bl[0]/2.0);
		    		yi=((float)(rand() % iy))*bl[1]+(bl[1]/2.0);
		    		zi=((float)(rand() % iz))*bl[2]+(bl[2]/2.0);
		    		xk=(int)floor(xi/bl[0]);
		    		yk=(int)floor(yi/bl[1]);
		    		zk=(int)floor(zi/bl[2]);
		    		//printf("i, iz %d %lf %lf %lf %d \n",i, xi, yi, zi, wdr[xk][yk][zk]);
				if (wdr[xk][yk][zk]==0)
				{
					dum=wdr[xk][yk][zk];
					dri[xk][yk][zk][dum]=i;
					wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
					ami[i][0]=xi;
					ami[i][1]=yi;
					ami[i][2]=zi;
					fau[i][0]=xi;
					fau[i][1]=yi;
					fau[i][2]=zi;
					
					
					cow=1;
					//fprintf(gu,"%lf  %lf  %lf \n", fau[i][0], fau[i][1], fau[i][2]);
				}
			}
			//printf("number filled %d\n", i);
		}
		//printf("Go yourself\n");
		
		cow=0;
		sitts=0;
		quao=0.0;
		quaoi=0.0;
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		printf("Initial configuration generated successfully!!!!!!! \n \n");
	}
	else
	{
		printf("A pre-prepared system found and intial configuration will be set accordingly. \n");
		printf("Initial configuration is being read from the available file!!! \n \n  ");
		cow=1;
		
		fscanf(vu, "%d", &sitts);
		
		fscanf(vu, "%lf", &quao);
		
		fscanf(vu, "%lf", &quaoi);
		
		char filename[25] = {0};

		sprintf(filename, "%d.txt", (int)floor((float)(sitts+1)/100.0));

		
		cfu=fopen(filename,"w");
		
		
		
		
		for (i=0; i<N; i++)
		{
			
			for(l=0;l<3;l++)
			{
				fscanf(vu,"%lf", &ami[i][l]);
			}
			
			fscanf(vu,"%d", &wnrd[i]);
			for (j=0; j<wnrd[i]; j++) 
			{
				fscanf(vu,"%d", &nrdi[i][j]);
			}
			
			fscanf(vu,"%d", &wbobb[i]);
			for (j=0; j<wbobb[i]; j++) 
			{
				fscanf(vu,"%d", &bobb[i][j]);
			}
			
			fscanf(vu,"%d", &wbobbi[i]);
			for (j=0; j<wbobbi[i]; j++) 
			{
				fscanf(vu,"%d", &bobbi[i][j]);
			}
			
			for(l=6;l<9;l++)
			{
				fscanf(vu,"%lf", &li[i][l]);
				
			}
			for(l=0;l<3;l++)
			{
				fscanf(vu,"%lf", &fau[i][l]);
			}
		}
		
		
		for ( i=0; i<3; i++)
		{
			f[i]=fe;
			block[i]=(int)(floor(f[i]/ble));
			bl[i]=f[i]/(float)block[i];
			f[i]=bl[i]*(float)block[i];
			//printf("ptr %d \n", ptr);
		}
		
			
		f1=f[0];
		f2=f[1];
		f3=f[2];
		block1=block[0];
		block2=block[1];
		block3=block[2];
		bl1=bl[0];
		bl2=bl[1];
		bl3=bl[2];
		for (i=0; i<blockd; i++)
		{
			
			for (j=0; j<blockd; j++)
			{
				for (k=0; k<blockd; k++)
				{
					wdr[i][j][k]=0;
					for (l=0; l<80; l++)
					{
						dri[i][j][k][l]=6000;
					}
				}	
			}
		}
		
		for(i=0; i<N; i++)	
	    	{	
			
			xi=ami[i][0];
			yi=ami[i][1];
			zi=ami[i][2];
	    		xk=(int)floor(xi/bl[0]);
	    		yk=(int)floor(yi/bl[1]);
	    		zk=(int)floor(zi/bl[2]);
	    		fau[i][0]=xi;
			fau[i][1]=yi;
			fau[i][2]=zi;
	    		
			dum=wdr[xk][yk][zk];
			dri[xk][yk][zk][dum]=i;
			wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
			if (xk>block1 || yk>block2 || zk>block3) printf("YOU ARE \n");
		
			
		}
		cow=1;
		
		printf("Initial configuration have been copied from the file and implemented into the system successfully!!!!!!! \n \n  ");
		
			//printf("Go and yourself\n");
		
		fclose(vu);
	}
	
	//printf("raam\n");
	
	
	for (k=0; k<N; k++)
	{
		theta=acos(li[k][8]);
		//printf("%lf\n", theta);
		if (fabs(li[k][8]) > 0.99999)
		{
		dx=1.0;
		dy=0.0;
		dz=0.0;
		}
		else
		{
		dx=-li[k][7]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
		dy=li[k][6]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
		dz=0.0;
		}
		q1=cos(theta/2.0);
		q2=dx*sin(theta/2.0);
		q3=dy*sin(theta/2.0);
		q4=dz*sin(theta/2.0);
		q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
		q21=2.0*((q2*q3)+(q4*q1));
		q31=2.0*((q2*q4)-(q3*q1));
		q12=2.0*((q2*q3)-(q4*q1));
		q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
		q32=2.0*((q3*q4)+(q2*q1));
		q13=2.0*((q2*q4)+(q3*q1));
		q23=2.0*((q3*q4)-(q2*q1)); 
		q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));

		li[k][0]=1.0*q11;
		li[k][1]=1.0*q21;
		li[k][2]=1.0*q31;
		li[k][3]=1.0*q12;
		li[k][4]=1.0*q22;
		li[k][5]=1.0*q32;
		li[k][6]=1.0*q13;
		li[k][7]=1.0*q23;
		li[k][8]=1.0*q33;
		
		fli[k][0]=1.0*q11;
		fli[k][1]=1.0*q21;
		fli[k][2]=1.0*q31;
		fli[k][3]=1.0*q12;
		fli[k][4]=1.0*q22;
		fli[k][5]=1.0*q32;
		fli[k][6]=1.0*q13;
		fli[k][7]=1.0*q23;
		fli[k][8]=1.0*q33;
		
		
		
		//printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  \n", q11, q21, q31, q12, q22, q32, q13, q23, q33);
	}
	
	
	for (i=0; i<N; i++)
	{
		rpid=pid[i];
		if (rpid>=war) printf("raaam naam\n");
		el[i][0]=(rsm[rpid]*li[i][0]*li[i][0])+(rsm[rpid]*li[i][3]*li[i][3])+(rs[rpid]*li[i][6]*li[i][6]);
		el[i][1]=(rsm[rpid]*li[i][0]*li[i][1])+(rsm[rpid]*li[i][3]*li[i][4])+(rs[rpid]*li[i][6]*li[i][7]);
		el[i][2]=(rsm[rpid]*li[i][0]*li[i][2])+(rsm[rpid]*li[i][3]*li[i][5])+(rs[rpid]*li[i][6]*li[i][8]);
		el[i][3]=(rsm[rpid]*li[i][0]*li[i][1])+(rsm[rpid]*li[i][3]*li[i][4])+(rs[rpid]*li[i][6]*li[i][7]);
		el[i][4]=(rsm[rpid]*li[i][1]*li[i][1])+(rsm[rpid]*li[i][4]*li[i][4])+(rs[rpid]*li[i][7]*li[i][7]);
		el[i][5]=(rsm[rpid]*li[i][1]*li[i][2])+(rsm[rpid]*li[i][4]*li[i][5])+(rs[rpid]*li[i][7]*li[i][8]);
		el[i][6]=(rsm[rpid]*li[i][0]*li[i][2])+(rsm[rpid]*li[i][3]*li[i][5])+(rs[rpid]*li[i][6]*li[i][8]);
		el[i][7]=(rsm[rpid]*li[i][1]*li[i][2])+(rsm[rpid]*li[i][4]*li[i][5])+(rs[rpid]*li[i][7]*li[i][8]);
		el[i][8]=(rsm[rpid]*li[i][2]*li[i][2])+(rsm[rpid]*li[i][5]*li[i][5])+(rs[rpid]*li[i][8]*li[i][8]);	
	}
	
	
	// Setting initial coordinates and direction ends!
	


	
	
     	// Randomization starts here!
     	
     	srand48(time(NULL));
	srand(time(NULL));
	//srand48(6);
	//srand(6);

	static unsigned long long  x=123456789,y=987654321,z=43219876,cc=6543217; 
	unsigned int JKISS()
	{ 
		unsigned long long t;
		x=314527869*x+1234567; 
		y^=y <<5;y^=y>>7;y^=y<<22;
		t = 4294584393ULL*z+cc; cc = t>>32; z= t;
		return x+y+z; 
	}

	// Randomization ends here!






	//Preparing for simulation!
     
     	FILE* swu;
     	swu=fopen("Diffusion_coefficient_el_s.txt","w");
     	
     	FILE* mwu;
     	mwu=fopen("Diffusion_coefficient_el_m.txt","w");
     	
     	FILE* lwu;
     	lwu=fopen("Diffusion_coefficient_el_l.txt","w");
     	
     	//FILE* dwu;
     	//dwu=fopen("Data_cord_dire.txt","w");
     	
     
     
     
     
     
     
     
     	int* trac=malloc(war*sizeof(int));
     	
     	for (i=0; i<war; i++)
     	{
     		trac[i]=0;
     	}
     	
     	int* tl=malloc(war*sizeof(int));
     	
     	for (i=0; i<war; i++)
     	{
     		tl[i]=(int)((float)nt[i]*2.0);
     	}
     
     	
     	
     	
     	
     	
     	
     	
	lsd=(int)floor((1.0/(step*step)));
	iltt=(int)floor((float)tsim/(float)lsd);
	rdd=(float)lsd/10.0;
	kd=0.0;
	istep=1.0/(step*step);
	
	//cow=0;
	titi=0;
	ri=0;
	nlies=2*N;
	teb=f1-0.01;
	eb=teb;
	rul=1;
	sla=0;
	ve=0;
	ptr=1;
	veei=1;
	vee=0;
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("t_sim required to achive unit physical time for the given step length-    ");
	printf("%d \n", lsd);
	printf("\n");
	printf("\n");

	printf("Maximum possible physical time with given step length and tsim is given by-   ");
	printf("%d\n",iltt);
	printf("\n");
	printf("\n");

	printf("\n");
	printf("\n");
	printf("\n");
	
	printf("                                       *** Shree Ganeshaay Namah ***  :::\n");
	//printf("aimed edge, f1 %lf %lf %d %d %d\n", f1, fe, block1, blockd, iz);
	printf("\n");
	printf("\n");
	printf("\n");
	printf("      *  ::  Be aware of yourself. Rest I will take care. By the way, I am running for you  ::  *\n");
	printf("\n");
	printf("\n");
	printf("\n");
	
	//Preparations done!
	
	//printf("%d\n", pid[299]);
	
	
	
	
	
	
	
	
	for (i=0; i<war; i++)
	{
		spo[i]=0.0;
	}
	
	
	
	
	//Simulation starts here!
	
	for(sitt=sitts; sitt<iltt; sitt++)
	{	
		
		
		
		if (sitt==20)
		{
			
			for (i=0; i<N; i++)
			{
				ldcxyz[i][0]=0.0;
				ldcxyz[i][1]=0.0;
				ldcxyz[i][2]=0.0;
				ldcxy[i][0]=0.0;
				ldcxy[i][1]=0.0;
				ldcz[i]=0.0;
				lom[i][0]=0.0;
				lom[i][1]=0.0;
				lom[i][2]=0.0;
				
				fli[i][0]=li[i][0];
				fli[i][1]=li[i][1];
				fli[i][2]=li[i][2];
				fli[i][3]=li[i][3];
				fli[i][4]=li[i][4];
				fli[i][5]=li[i][5];
				fli[i][6]=li[i][6];
				fli[i][7]=li[i][7];
				fli[i][8]=li[i][8];
				cow=1;
				kd=0.0;
				titi=0;
			}
		
		}
		
		for(ittt=0; ittt<lsd; ittt++)
		{	
		
		
			jet=0;
			//printf("raam %d\n",ittt);
			/*if(cow==1)
			{
				printf("raam%d\n", ittt);
			}*/
			
			/*for (i=0; i<N; i++)
			{
				
				wcbobb[i]=0;
				for (j=0; j<E; j++)
				{
					cbobb[i][j]=6000;
				}				
			}*/
			
			
			
			
			
			if (cow==0)
			{
				ittt=0;
				cow=1;
				
				//printf("raam\n");
				if (teb >= fe)
				{
					//printf("raam\n");
					ptr=4;
					for (jt=0; jt<N; jt++)
					{
						if (ami[jt][0]>=eb-0.01)
						{
							ptr=1;
							//titi=titi+1;
							cow=0;
						}
						else if (ami[jt][1]>=eb-0.01)
						{
							if (ptr>=2) ptr=2;
							cow=0;
						}
						else if (ami[jt][2]>=eb-0.01)
						{
							if (ptr>=3)ptr=3;
							cow=0;
						}
						
					}
					//printf("raam ptr %d\n", ptr);
					for ( i=0; i<(ptr-1); i++)
					{
						f[i]=teb;
						block[i]=(int)(floor(f[i]/ble));
						bl[i]=f[i]/(float)block[i];
						f[i]=bl[i]*(float)block[i];
						//printf("ptr %d \n", ptr);
					}
					
					if (f[0]!=f1 || f[1]!=f2 || f[2]!=f3 )
					{
						
						f1=f[0];
						f2=f[1];
						f3=f[2];
						block1=block[0];
						block2=block[1];
						block3=block[2];
						bl1=bl[0];
						bl2=bl[1];
						bl3=bl[2];
						for (i=0; i<blockd; i++)
						{
							
							for (j=0; j<blockd; j++)
							{
								for (k=0; k<blockd; k++)
								{
									wdr[i][j][k]=0;
									for (l=0; l<80; l++)
									{
										dri[i][j][k][l]=6000;
									}
								}	
							}
						}
						
						for(i=0; i<N; i++)	
					    	{	
							
							xi=ami[i][0];
							yi=ami[i][1];
							zi=ami[i][2];
					    		xk=(int)floor(xi/bl[0]);
					    		yk=(int)floor(yi/bl[1]);
					    		zk=(int)floor(zi/bl[2]);
					    		fau[i][0]=xi;
							fau[i][1]=yi;
							fau[i][2]=zi;
					    		
							dum=wdr[xk][yk][zk];
							dri[xk][yk][zk][dum]=i;
							wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
							if (xk>block1 || yk>block2 || zk>block3) printf("YOU ARE doomed\n");
						
							
						}
					
					}
					if (rul == 0) 
					{	
						cow=0;
						ptr=1;
					}
					if (veei == 0) 
					{
						ve=ve+1;
						cow=0;
						ptr=1;
					}
					
					if (cow==1)
					{
						
						//vee=0;
						rul=0;
						teb=teb-0.001;
						eb=teb;
						
						cow=0;
						ptr=1;
						
						if (teb<((fe+((float)veec*vees)-vees)-((float)vee*vees))) 
						{
							veei=0;
						}
						
						//printf("teb %lf \n", teb);
						
					}
					
					if (ve > vel) 
					{
						ve=0; veei=1; vee=vee+1; //printf("raaam\n");
					}
					
				}
				
				else
				{
					
					
					if (ri==0)
					{
						veei = 0; ve = 0; ri = 1; vel=20000;
					}	
					
					if (veei == 0) 
					{
						ve = ve+1;
						cow = 0;
				
					}					
					if (ve > vel) 
					{
						 veei = 1; //printf("raaam\n");
					}			
				
				}
				
				
				//if(rul==0) printf("rule %d\n",rul);
				
				//printf("%d \n",titi);
				if(cow==1)
				{
					
					for (i=0; i<N; i++)
					{
						sdcxyz[i][0]=0.0;
						sdcxyz[i][1]=0.0;
						sdcxyz[i][2]=0.0;
						sdcxy[i][0]=0.0;
						sdcxy[i][1]=0.0;
						sdcz[i]=0.0;
						som[i][0]=0.0;
						som[i][1]=0.0;
						som[i][2]=0.0;
						
						mdcxyz[i][0]=0.0;
						mdcxyz[i][1]=0.0;
						mdcxyz[i][2]=0.0;
						mdcxy[i][0]=0.0;
						mdcxy[i][1]=0.0;
						mdcz[i]=0.0;
						mom[i][0]=0.0;
						mom[i][1]=0.0;
						mom[i][2]=0.0;
						
						ldcxyz[i][0]=0.0;
						ldcxyz[i][1]=0.0;
						ldcxyz[i][2]=0.0;
						ldcxy[i][0]=0.0;
						ldcxy[i][1]=0.0;
						ldcz[i]=0.0;
						lom[i][0]=0.0;
						lom[i][1]=0.0;
						lom[i][2]=0.0;
						
						xi=ami[i][0];
						yi=ami[i][1];
						zi=ami[i][2];
				    		fau[i][0]=xi;
						fau[i][1]=yi;
						fau[i][2]=zi;
						
					}
					
					
					FILE* ku;
     					ku=fopen("Compressed_system.txt","w");
     					
					for (i=0; i<N; i++)
					{
						fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
						
						fprintf(ku, "%d\n", wnrd[i]);
						for (j=0; j<wnrd[i]; j++) 
						{
							fprintf(ku,"   %d", nrdi[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku, "%d\n", wbobb[i]);
						for (j=0; j<wbobb[i]; j++) 
						{
							fprintf(ku,"   %d", bobb[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku, "%d\n", wbobbi[i]);
						for (j=0; j<wbobbi[i]; j++) 
						{
							fprintf(ku,"   %d", bobbi[i][j]);
						}
						fprintf(ku, "\n");
						
						fprintf(ku,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
						
						fprintf(ku,"%lf  %lf  %lf\n", fau[i][0], fau[i][1], fau[i][2]);
					}
					
					
					
					fclose(ku);
					printf("Compression happend!!! Now the actual simulation will start!!! %lf %lf %lf\n", f1, f2, f3);
				}
						 
			}
			//printf("raammmmmmmmmdsfdfdsfdsfmm\n");
			
			sla=sla+1;
			
			//printf("raam\n");
			
			for (i=0; i<war; i++)
		     	{
		     		trac[i]=0;
		     	}
			
			for(mo=0; mo<nlies; mo++)
			{
				
				
				
            			//k=rand() % N;
            			boss=drand48();
            			
            			ntt=N;
            			for (i=0; i<war; i++)
            			{
            				if(trac[i]>tl[i])	ntt=ntt-nt[i];
            			}
            			k=rand()  % ntt;
            			
            			
            			loo=0;
            			for (i=0; i<war; i++)
            			{
            				if(trac[i]>tl[i])
            				{
            					if (k>loo) k=k+nt[i];
            				}
            				loo=loo+nt[i];
            			}
            			
            			
            			rpid1=pid[k];
            			trac[rpid1]=trac[rpid1]+1;
            			
            		
            			
            			pf=0;
            			loo=0;
            			
            	               /*
            			
            			if (cow==0)
				{
					
					
				}
				
				*/
				 			
            			mxi=ami[k][0];
	    			myi=ami[k][1];
	    			mzi=ami[k][2];
				
				
				
	    			l=(int)floor(mxi/bl1);
	    			m=(int)floor(myi/bl2);
	    			n=(int)floor(mzi/bl3);
				//printf("raamtlllllllllllllllll%d\n",k);
				
				
				
				wpbobb=0;
				for(i=0; i<E; i++)
				{
					psbobb[i]=6000;
				}
		    		wpbobbi=0;
				for(i=0; i<E; i++)
				{
					psbobbi[i]=6000;
				}
				
				
				
				
            			if (boss >= 0.5)
				{
					
					//printf("raamtt\n");
					
					rdphi=tpiee*drand48();
					rdth=rdtheta[rpid1];
					
					xf=sin(rdth)*cos(rdphi);
					yf=sin(rdth)*sin(rdphi);
					zf=cos(rdth);
					dx=(-yf)*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
					dy=xf*(1.0/(pow(((xf*xf)+(yf*yf)),0.5)));
					dz=0.0;
					rdom[0]=dx;
					rdom[1]=dy;
					rdom[2]=dz;
					q1=cos(rdth/2.0);
					q2=dx*sin(rdth/2.0);
					q3=dy*sin(rdth/2.0);
					q4=dz*sin(rdth/2.0);
					q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
					q12=2.0*((q2*q3)+(q4*q1));
					q13=2.0*((q2*q4)-(q3*q1));
					q21=2.0*((q2*q3)-(q4*q1));
					q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
					q23=2.0*((q3*q4)+(q2*q1));
					q31=2.0*((q2*q4)+(q3*q1));
					q32=2.0*((q3*q4)-(q2*q1)); 
					q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
					rdli[0]=q11;
					rdli[1]=q12;
					rdli[2]=q13;
					rdli[3]=q21;
					rdli[4]=q22;
					rdli[5]=q23;
					rdli[6]=q31;
					rdli[7]=q32;
					rdli[8]=q33;
					
					
					theta=acos(li[k][8]);
					if (fabs(li[k][8]) > 0.99999)
					{
						dx=1.0;
						dy=0.0;
						dz=0.0;
					}
					else
					{
						dx=-li[k][7]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
						dy=li[k][6]*(1.0/pow(((li[k][6]*li[k][6])+(li[k][7]*li[k][7])),0.5));
						dz=0.0;
					}
					q1=cos(theta/2.0);
					q2=dx*sin(theta/2.0);
					q3=dy*sin(theta/2.0);
					q4=dz*sin(theta/2.0);
					q11=1.0-2.0*(pow(q3,2.0)+pow(q4,2.0));
					q21=2.0*((q2*q3)+(q4*q1));
					q31=2.0*((q2*q4)-(q3*q1));
					q12=2.0*((q2*q3)-(q4*q1));
					q22=1.0-2.0*(pow(q2,2.0)+pow(q4,2.0));  
					q32=2.0*((q3*q4)+(q2*q1));
					q13=2.0*((q2*q4)+(q3*q1));
					q23=2.0*((q3*q4)-(q2*q1)); 
					q33=1.0-2.0*(pow(q2,2.0)+pow(q3,2.0));
					
					dli[0]=rdli[0]*q11+rdli[1]*q12+rdli[2]*q13;
					dli[1]=rdli[0]*q21+rdli[1]*q22+rdli[2]*q23;
					dli[2]=rdli[0]*q31+rdli[1]*q32+rdli[2]*q33;
					dli[3]=rdli[3]*q11+rdli[4]*q12+rdli[5]*q13;
					dli[4]=rdli[3]*q21+rdli[4]*q22+rdli[5]*q23;
					dli[5]=rdli[3]*q31+rdli[4]*q32+rdli[5]*q33;
					dli[6]=rdli[6]*q11+rdli[7]*q12+rdli[8]*q13;
					dli[7]=rdli[6]*q21+rdli[7]*q22+rdli[8]*q23;
					dli[8]=rdli[6]*q31+rdli[7]*q32+rdli[8]*q33;
					
					dom0=rdom[0]*q11+rdom[1]*q12+rdom[2]*q13;
					dom1=rdom[0]*q21+rdom[1]*q22+rdom[2]*q23;
					dom2=rdom[0]*q31+rdom[1]*q32+rdom[2]*q33;
					
					rst=rs[rpid1];
					rsmt=rsm[rpid1];
					del[0]=(rsmt*dli[0]*dli[0])+(rsmt*dli[3]*dli[3])+(rst*dli[6]*dli[6]);
					del[1]=(rsmt*dli[0]*dli[1])+(rsmt*dli[3]*dli[4])+(rst*dli[6]*dli[7]);
					del[2]=(rsmt*dli[0]*dli[2])+(rsmt*dli[3]*dli[5])+(rst*dli[6]*dli[8]);
					del[3]=(rsmt*dli[0]*dli[1])+(rsmt*dli[3]*dli[4])+(rst*dli[6]*dli[7]);
					del[4]=(rsmt*dli[1]*dli[1])+(rsmt*dli[4]*dli[4])+(rst*dli[7]*dli[7]);
					del[5]=(rsmt*dli[1]*dli[2])+(rsmt*dli[4]*dli[5])+(rst*dli[7]*dli[8]);
					del[6]=(rsmt*dli[0]*dli[2])+(rsmt*dli[3]*dli[5])+(rst*dli[6]*dli[8]);
					del[7]=(rsmt*dli[1]*dli[2])+(rsmt*dli[4]*dli[5])+(rst*dli[7]*dli[8]);
					del[8]=(rsmt*dli[2]*dli[2])+(rsmt*dli[5]*dli[5])+(rst*dli[8]*dli[8]);
					
					if (cow==0)
					{
						if (rul==0 && sla==5)
						{
							for (jt=0; jt<9; jt++)
							{
								dli[jt]=li[k][jt];
								del[jt]=el[k][jt];
								
								
							}
						}
					}			
						//printf("raamtt\n");
					check=0;
					
						
 					for (ip=0; ip<wnrd[k]; ip++)
					{
						//printf("raamtttttttttttttt\n");
						
						w1=dli[0];
						w2=dli[1];
						w3=dli[2];
						w4=dli[3];
						w5=dli[4];
						w6=dli[5];
						w7=dli[6];
						w8=dli[7];
						w9=dli[8];

						
						puc=nrdi[k][ip];
						rpid2=pid[puc];
						
						u1=li[puc][0];
						u2=li[puc][1];
						u3=li[puc][2];
						u4=li[puc][3];
						u5=li[puc][4];
						u6=li[puc][5];
						u7=li[puc][6];
						u8=li[puc][7];
						u9=li[puc][8];
						re1=ami[puc][0]-mxi;
						re2=ami[puc][1]-myi;
						re3=ami[puc][2]-mzi;
						if (re1 < -pbc) re1=re1+f1;
						else if (re1 > pbc) re1=re1-f1;
						if (re2 < -pbc) re2=re2+f2;
						else if (re2 > pbc) re2=re2-f2;
						if (re3 < -pbc) re3=re3+f3;
						else if (re3 > pbc) re3=re3-f3;
						
						
						mr=(re1*re1+re2*re2+re3*re3);
						mrp=pow(mr,0.5);
						
						ree[0]=(re1*w1+re2*w2+re3*w3);
						ree[1]=(re1*w4+re2*w5+re3*w6);
						ree[2]=(re1*w7+re2*w8+re3*w9);
						
						rtm[0][0]=(u1*w1+u2*w2+u3*w3);
						rtm[1][0]=(u1*w4+u2*w5+u3*w6);
						rtm[2][0]=(u1*w7+u2*w8+u3*w9);
						rtm[0][1]=(u4*w1+u5*w2+u6*w3);
						rtm[1][1]=(u4*w4+u5*w5+u6*w6);
						rtm[2][1]=(u4*w7+u5*w8+u6*w9);
						rtm[0][2]=(u7*w1+u8*w2+u9*w3);
						rtm[1][2]=(u7*w4+u8*w5+u9*w6);
						rtm[2][2]=(u7*w7+u8*w8+u9*w9);
						
						for ( i = 0; i < 3; i++)
						{
							for ( j = 0; j < 3; j++)
								abrtm[i][j] = fabs(rtm[i][j]) + EPSILON;
						}						
						
						tri=0;
						
						for (i = 0; i < 3; i++) 
						{
							ra = e[rpid1][i];
							rb = e[rpid2][0] * abrtm[i][0] + e[rpid2][1] * abrtm[i][1] + e[rpid2][2] * abrtm[i][2];
							if (fabs(ree[i]) > ra + rb) tri=1;
						}
						
						
						for (i = 0; i < 3; i++) 
						{
							ra = e[rpid1][0] * abrtm[0][i] + e[rpid1][1] * abrtm[1][i] + e[rpid1][2] * abrtm[2][i];
							rb = e[rpid2][i];
							if (fabs(ree[0] * rtm[0][i] + ree[1] * rtm[1][i] + ree[2] * rtm[2][i]) > ra + rb) tri=1;
						}
						
						
						ra= e[rpid1][1] * abrtm[2][0] + e[rpid1][2] * abrtm[1][0];
						rb= e[rpid2][1] * abrtm[0][2] + e[rpid2][2] * abrtm[0][1];
						if (fabs(ree[2] * rtm[1][0] - ree[1] * rtm[2][0]) > ra + rb) tri=1;
						
						ra= e[rpid1][1] * abrtm[2][1] + e[rpid1][2] * abrtm[1][1];
						rb= e[rpid2][0] * abrtm[0][2] + e[rpid2][2] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][1] - ree[1] * rtm[2][1]) > ra + rb) tri=1;
						
						ra= e[rpid1][1] * abrtm[2][2] + e[rpid1][2] * abrtm[1][2];
						rb= e[rpid2][0] * abrtm[0][1] + e[rpid2][1] * abrtm[0][0];
						if (fabs(ree[2] * rtm[1][2] - ree[1] * rtm[2][2]) > ra + rb) tri=1;
						
						ra = e[rpid1][0] * abrtm[2][0] + e[rpid1][2] * abrtm[0][0];
						rb = e[rpid2][1] * abrtm[1][2] + e[rpid2][2] * abrtm[1][1];
						if (fabs(ree[0] * rtm[2][0] - ree[2] * rtm[0][0]) > ra + rb) tri=1;
					
						ra= e[rpid1][0] * abrtm[2][1] + e[rpid1][2] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[1][2] + e[rpid2][2] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][1] - ree[2] * rtm[0][1]) > ra + rb) tri=1;
						
						ra= e[rpid1][0] * abrtm[2][2] + e[rpid1][2] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[1][1] + e[rpid2][1] * abrtm[1][0];
						if (fabs(ree[0] * rtm[2][2] - ree[2] * rtm[0][2]) > ra + rb) tri=1;
						
						ra= e[rpid1][0] * abrtm[1][0] + e[rpid1][1] * abrtm[0][0];
						rb= e[rpid2][1] * abrtm[2][2] + e[rpid2][2] * abrtm[2][1];
						if (fabs(ree[1] * rtm[0][0] - ree[0] * rtm[1][0]) > ra + rb) tri=1;
							
						ra= e[rpid1][0] * abrtm[1][1] + e[rpid1][1] * abrtm[0][1];
						rb= e[rpid2][0] * abrtm[2][2] + e[rpid2][2] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][1] - ree[0] * rtm[1][1]) > ra + rb) tri=1;
						
						ra= e[rpid1][0] * abrtm[1][2] + e[rpid1][1] * abrtm[0][2];
						rb= e[rpid2][0] * abrtm[2][1] + e[rpid2][1] * abrtm[2][0];
						if (fabs(ree[1] * rtm[0][2] - ree[0] * rtm[1][2]) > ra + rb) tri=1;
						
						if (tri==0)
						{
								
							w1=dli[6];
							w2=dli[7];
							w3=dli[8];
							
							
							
							u1=li[puc][6];
							u2=li[puc][7];
							u3=li[puc][8];
							del2[0]=el[puc][0];
							del2[1]=el[puc][1];
							del2[2]=el[puc][2];
							del2[3]=el[puc][3];
							del2[4]=el[puc][4];
							del2[5]=el[puc][5];
							del2[6]=el[puc][6];
							del2[7]=el[puc][7];
							del2[8]=el[puc][8];
							
							a=0.1;
							b=1.0;
							c=(b+a)/2.0;
							
							trig=1.0;
							
							cfun[0]=0.0;
							cfun[1]=0.0;
							tri=0;
							while (trig>0.0001)
							{
								tri=tri+1;
								
								
								lemm=1.0-c;
								lem=c;
								lemms=lemm*lemm;
								lems=lem*lem;
								wt=0;
								
								for (i=0; i<3; i++)
								{
									for (j=0; j<3; j++)
									{
										lt[i][j]=del[wt]*lemm+del2[wt]*lem;
										wt=wt+1;
									}
								}
								
								wt=0;
								
								for (i=0; i<3; i++)
								{
									for (j=0; j<3; j++)
									{
										lts[i][j]=del[wt]*lemms-del2[wt]*lems;
										wt=wt+1;
									}
								}
								lt[0][3]=re1;
								lt[1][3]=re2;
								lt[2][3]=re3;
								cont=1;
								for (it=0; it<2; it++) 
								{
								    	dummy=0;
								    	if  (lt[it][it]==0.0)
									{
										for (jt=0; jt<3; jt++) 
										{
									    		if (dummy == 0)	
											{
												if (lt[jt][it] != 0.0) 
												{
										    			mt=jt;
										    			dummy=dummy+1;
										    			dump[0]=lt[it][0];
													dump[1]=lt[it][1];
													dump[2]=lt[it][2];
													dump[3]=lt[it][3];
										    			lt[it][0]=lt[mt][0];
													lt[it][1]=lt[mt][1];
													lt[it][2]=lt[mt][2];
													lt[it][3]=lt[mt][3];
										    			lt[mt][0]=dump[0];
													lt[mt][1]=dump[1];
													lt[mt][2]=dump[2];
													lt[mt][3]=dump[3];
												}
											}
										}
										if (dummy==0) continue; 
									}
								    	for (jt=cont; jt<3; jt++)
									{
										dumm=lt[jt][it];
										for (mt=0; mt<4; mt++)
										{
									    		lt[jt][mt]=lt[jt][mt]-((lt[it][mt]/lt[it][it])*dumm);
										}
									}
								    	cont=cont+1;
								}
								//NOW CHECK EXISTENCE OF SOLUTION OF EQUATION
								sol[0]=0.0; sol[1]=0.0; sol[2]=0.0; kt=3;
								for (it=2; it>=0; it--)
								{    	sol[it]=lt[it][kt];
								    	for (jt=0; jt<3; jt++)
									{
										if (it==jt) continue;
										sol[it]=sol[it]-(lt[it][jt]*sol[jt]);
									}
								   	sol[it]=sol[it]/lt[it][it];
								}
								ec=0.0;
								for (i=0; i<3; i++)
								{
									ec=ec+(lts[i][0]*sol[0]+lts[i][1]*sol[1]+lts[i][2]*sol[2])*sol[i];
								}
													
								if (ec>0.0)
								{
									a=c;
									cfun[1]=ec;
								}
			    					else
			    					{
									b=c;
									cfun[2]=ec;
								}
								c=(a+b)/2.0;	
			    					trig=fabs(b-a);
			    					
			    					ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
				    				//if (ec>1.0) trig=0.000001;	
								
							}
							
							ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
							if (ec<1.0)
							{
								check=1;
								jet=jet+1;
							}
							else
							{	
							
								
								
								if (pop[rpid1][rpid2]==1)
								
								{
									ec=ec/mr;
									ec=ec*(mrp-epsi1[rpid1][rpid2])*(mrp-epsi1[rpid1][rpid2]);
									if (ec<1.0)
									{
									
										
										
										zt1=(re1*w1+re2*w2+re3*w3)/mrp;
										zt2=(re1*u1+re2*u2+re3*u3)/mrp;
										if ( zt1 > 0.0 && zt2 < 0.0 )
										{
											psbobb[wpbobb]=puc;
											wpbobb=wpbobb+1;	
										}
										
									}
								
								}
								
								
								else if (pop[rpid1][rpid2]==2)
								
								{
									ec=ec/mr;
									ec=ec*(mrp-epsi2[rpid1][rpid2])*(mrp-epsi2[rpid1][rpid2]);
									if (ec<1.0)
									{
									
										//For patch
									
										zt1=(fabs(re1*w1+re2*w2+re3*w3))/mrp;
										zt2=(fabs(re1*u1+re2*u2+re3*u3))/mrp;
										if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
										{
											psbobb[wpbobb]=puc;
											wpbobb=wpbobb+1;	
										}
										
									}
								
								}
								
								
								else if (pop[rpid1][rpid2]==3)
								{
									ec=ec/mr;
									ec=ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2]);
									if (ec<1.0)
									{
									
										psbobbi[wpbobbi]=puc;
										wpbobbi=wpbobbi+1;	
										
									}
								}
								
								else if (pop[rpid1][rpid2]==4)
								
								{
									ec=ec/mr;
									
									
									
									if ((ec*(mrp-epsi1[rpid1][rpid2])*(mrp-epsi1[rpid1][rpid2])) < 1.0) 
									{
									
										
										zt1=(re1*w1+re2*w2+re3*w3)/mrp;
										zt2=(re1*u1+re2*u2+re3*u3)/mrp;
										
										if ( zt1 > 0.0 && zt2 < 0.0 )
										{
											psbobb[wpbobb]=puc;
											wpbobb=wpbobb+1;
											
										
										}
									}
									
									
									if ((ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2])) < 1.0)
									
									{
										psbobbi[wpbobbi]=puc;
										wpbobbi=wpbobbi+1;
											
									}
									
								
								}
								
								
								else if (pop[rpid1][rpid2]==5)
								
								{
									ec=ec/mr;
										
									if ((ec*(mrp-epsi2[rpid1][rpid2])*(mrp-epsi2[rpid1][rpid2])) < 1.0) 
									{
									
										
										zt1=(fabs(re1*w1+re2*w2+re3*w3))/mrp;
										zt2=(fabs(re1*u1+re2*u2+re3*u3))/mrp;
										if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
										{
											psbobb[wpbobb]=puc;
											wpbobb=wpbobb+1;
											
										
										}
									}
									
									if ((ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2])) < 1.0)
									
									{
										psbobbi[wpbobbi]=puc;
										wpbobbi=wpbobbi+1;
											
									}
								
								}
											
									
							}
							
								
						}
						
						if (check==1)
						{
							break;
						}
						
						
					}
						
					
						
					if (check==0) 
					{	
						
						for (it=0; it<wbobb[k]; it++)
						{
							check=1;
							for (i=0; i < wpbobb; i++)
							{
								if (bobb[k][it]==psbobb[i]) check=0;
								
							}
							if (check==1)  break; 
						}
						
						
						
						
						if (check==0)
						{
							for (it=0; it<wbobbi[k]; it++)
							{
								check=1;
								for (i=0; i < wpbobbi; i++)
								{
									if (bobbi[k][it]==psbobbi[i]) check=0;
									
								}
								if (check==1)  break; 					
							}
						}
							
						
						
						
						
						
						
						
					}	
							
						
						
					
					//printf("out\n");
	
				
			
					if (check==0) 
					{
						
						if(cow==1)
						{
							mum=wcbobb[k];
								
							for (jt=0; jt<mum; ++jt)
		    					{
								chum=cbobb[k][jt];
								drum=wcbobb[chum];
								for (sj=0; sj<drum; ++sj)
								{
									if (cbobb[chum][sj]==k)
									{
										for (tric=sj; tric<drum; tric++)
										{
											cbobb[chum][tric]=cbobb[chum][(tric+1)];
										}
										break;
									}
								
								}
								wcbobb[chum]=drum-1;
								cbobb[k][jt]=6000;
							}
								
							wcbobb[k]=0;	
							for (jdt=0; jdt < wpbobb; jdt++)
							{	
								//printf("raam is ");
								td=psbobb[jdt];
								cbobb[k][jdt]=td;
									tum=wcbobb[td];
									cbobb[td][tum]=k;
								wcbobb[td]=tum+1;
								wcbobb[k]=wcbobb[k]+1;
							}
							
							
							
								
								mum=wcbobbi[k];
								
								for (jt=0; jt<mum; ++jt)
			    					{
									chum=cbobbi[k][jt];
									drum=wcbobbi[chum];
									for (sj=0; sj<drum; ++sj)
									{
										if (cbobbi[chum][sj]==k)
										{
											for (tric=sj; tric<drum; tric++)
											{
												cbobbi[chum][tric]=cbobbi[chum][(tric+1)];
											}
											break;
										}
									
									}
									wcbobbi[chum]=drum-1;
									cbobbi[k][jt]=6000;
								}
									
								wcbobbi[k]=0;	
								for (jdt=0; jdt < wpbobbi; jdt++)
								{	
									//printf("raam is ");
									td=psbobbi[jdt];
									cbobbi[k][jdt]=td;
										tum=wcbobbi[td];
										cbobbi[td][tum]=k;
									wcbobbi[td]=tum+1;
									wcbobbi[k]=wcbobbi[k]+1;
								}
							
							
							
							
							
							
							
							
						}
						
						
						for (jt=0; jt<9; jt++)
						{
							li[k][jt]=dli[jt];
							el[k][jt]=del[jt];
							
							
						}
						som[k][0]=som[k][0]+dom0;
						som[k][1]=som[k][1]+dom1;
						som[k][2]=som[k][2]+dom2;
						
						mom[k][0]=mom[k][0]+dom0;
						mom[k][1]=mom[k][1]+dom1;
						mom[k][2]=mom[k][2]+dom2;
						
						lom[k][0]=lom[k][0]+dom0;
						lom[k][1]=lom[k][1]+dom1;
						lom[k][2]=lom[k][2]+dom2;
						
					}
					//printf("raamtttt\n");  	
				}
                           	
                		else
				{
					//printf("raamtttt\n");            			
					trim=2;
					ra1=drand48();
	        			theta=acos(1.0-(2.0*ra1));
					//ra2=drand48();
					ra2 = JKISS() / 4294967296.0;
	       			phi=2.0*piee*ra2;
					//printf("raaaam\n");
       				        /*dx=step*sin(theta)*cos(phi);
	        			dy=step*sin(theta)*sin(phi);
	        			dz=step*cos(theta);*/
					
					
					  						
					//x2per=step*pdx-dx;
					//y2per=step*pdy-dy;
					//z2per=step*pdz-dz;
					//printf("%lf  %lf  %lf\n", x2per, y2per, z2per);
					/*dx=step*sin(theta)*cos(phi);
					dy=step*sin(theta)*sin(phi);
					dz=step*cos(theta);*/
					
					if(cow==0)
					{
						if(veei == 0)
						{
							dx=0.01*sin(theta)*cos(phi);
							dy=0.01*sin(theta)*sin(phi);
							dz=0.01*cos(theta);
							xi=dx+ami[k][0];
							dxi=dx+fau[k][0];
							yi=dy+ami[k][1];
							dyi=dy+fau[k][1];
							zi=dz+ami[k][2];
							dzi=dz+fau[k][2];	
							
						}
						else
						{
							if (ptr==1)
							{
								dx=0.01*sin(theta)*cos(phi);
								dy=0.01*sin(theta)*sin(phi);
								dz=0.01*cos(theta);
								dx=dx-0.001;
								if (rul==0)
								{
									if (sla==5)
									{
										dx=0.0;
										dy=0.0;
										dz=0.0;
									}
								}
								xi=dx+ami[k][0];
								dxi=dx+fau[k][0];
								yi=dy+ami[k][1];
								dyi=dy+fau[k][1];
								zi=dz+ami[k][2];
								dzi=dz+fau[k][2];
								if (xi<0.0)
								{
									trim=1;
								}
								if (ami[k][0]<eb-0.01)
								{
									if (xi>eb-0.01) 
									{
										trim=1;
									}
								}
							}
							
							else if (ptr==2)
							{
								dx=0.01*sin(theta)*cos(phi);
								dy=0.01*sin(theta)*sin(phi);
								dz=0.01*cos(theta);
								dy=dy-0.001;
								if (rul==0)
								{
									if (sla==5)
									{
										dx=0.0;
										dy=0.0;
										dz=0.0;
									}
								}
								xi=dx+ami[k][0];
								dxi=dx+fau[k][0];
								yi=dy+ami[k][1];
								dyi=dy+fau[k][1];
								zi=dz+ami[k][2];
								dzi=dz+fau[k][2];
								if (yi<0.0)
								{
									trim=1;
								}
								if (ami[k][1]<eb-0.01)
								{
									if (yi>eb-0.01) 
									{
										trim=1;
									}
								}
							}
							else if (ptr==3)
							{
								dx=0.01*sin(theta)*cos(phi);
								dy=0.01*sin(theta)*sin(phi);
								dz=0.01*cos(theta);
								dz=dz-0.001;
								if (rul==0)
								{
									if (sla==5)
									{
										dx=0.0;
										dy=0.0;
										dz=0.0;
									}
								}
								xi=dx+ami[k][0];
								dxi=dx+fau[k][0];
								yi=dy+ami[k][1];
								dyi=dy+fau[k][1];
								zi=dz+ami[k][2];
								dzi=dz+fau[k][2];
								if (zi<0.0)
								{
									trim=1;
								}
								if (ami[k][2]<eb-0.01)
								{
									if (zi>eb-0.01) 
									{
										trim=1;
									}
								}
							}		
						}
					}
				
					else
					{
						//printf("zppppppppppp\n");
						pdx=sin(theta)*cos(phi);
						pdy=sin(theta)*sin(phi);
						pdz=cos(theta);
						
						tpart=(li[k][6]*pdx+li[k][7]*pdy+li[k][8]*pdz)*step*grz[rpid1];
						x1per=(li[k][0]*pdx+li[k][1]*pdy+li[k][2]*pdz)*step*grxy[rpid1];
						y1per=(li[k][3]*pdx+li[k][4]*pdy+li[k][5]*pdz)*step*grxy[rpid1];
						
						//tpart=(li[k][6]*pdx+li[k][7]*pdy+li[k][8]*pdz)*step;
						//x1per=(li[k][0]*pdx+li[k][1]*pdy+li[k][2]*pdz)*step;
						//y1per=(li[k][3]*pdx+li[k][4]*pdy+li[k][5]*pdz)*step;
						
						dx=li[k][6]*tpart+li[k][0]*x1per+li[k][3]*y1per;
						dy=li[k][7]*tpart+li[k][1]*x1per+li[k][4]*y1per;
						dz=li[k][8]*tpart+li[k][2]*x1per+li[k][5]*y1per; 
						
						xi=dx+ami[k][0];
						dxi=dx+fau[k][0];
						yi=dy+ami[k][1];
						dyi=dy+fau[k][1];
						zi=dz+ami[k][2];
						dzi=dz+fau[k][2];
					}
					
				
					
	        			if (xi >= f1) xi=xi-f1;
	        			else if (xi < 0.0) xi=xi+f1;
					if (yi >= f2) yi=yi-f2;
	        			else if (yi < 0.0) yi=yi+f2;
					if (zi >= f3) zi=zi-f3;
	        			else if (zi < 0.0) zi=zi+f3;
	        			
	        			
	        			if (xi==f1)
					{
						xi=xi-0.000001;
					}
					if (yi==f2)
					{
						yi=yi-0.000001;
					}
					if (zi==f3)
					{
						zi=zi-0.000001;
					}
	        			
	        			
	        			
	        			xk=(int)floor(xi/bl1);
	        			yk=(int)floor(yi/bl2);
	        			zk=(int)floor(zi/bl3);
					//printf("raaaam2 %d %d %d %lf %lf %lf\n",xk,yk,zk,dxi,dyi,bl);
						
						
					
		    				
	    				for (j=0; j<400; j++)
					{
						pit[j]=6000;
					}
					wpi=0;	
					dum=wdr[xk][yk][zk];
					//printf("raaaam3\n");
					for (jdt=0; jdt<dum; jdt++)
					{
						da=dri[xk][yk][zk][jdt]; rpid2=pid[da];
	                			xf=ami[da][0];
	                			yf=ami[da][1];
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
								
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[xk][yk][zk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
					//printf("raaaam\n");
	        			nxk=xk+1;
	        			nyk=yk;
	        			nzk=zk;
	        			if (nxk>=block1)
					{
	            				nxk=0;
	            				cp=f1;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
						da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];
	                			xf=ami[da][0]+cp;
	                			yf=ami[da][1];
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;
	                        				jet=jet+1;
	                        				break;
					
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
	        			nxk=xk-1;
	        			nyk=yk;
	        			nzk=zk;
	        			if (nxk==-1)
					{
	            				nxk=nxk+block1;
	            				cp=f1;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cp;
	                			yf=ami[da][1];
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk;
	        			nyk=yk+1;
	        			nzk=zk;
	        			if (nyk>=block2)
					{
	            				nyk=0;
	            				cp=f2;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]+cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  

	        			nxk=xk;
	        			nyk=yk-1;
	        			nzk=zk;
	        			if (nyk==-1)
					{
	            				nyk=nyk+block2;
	           				cp=f2;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]-cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk;
	        			nyk=yk;
	        			nzk=zk+1;
	        			if (nzk>=block3)
	            			{
						nzk=0;
					        cp=f3;
					}
	        			else cp=0.0;
					
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1];
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	       
	        			nxk=xk;
	        			nyk=yk;
	        			nzk=zk-1;
	        			if (nzk==-1)
					{
	            				nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1];
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  

	        			nxk=xk+1;
	        			nyk=yk-1;
	        			nzk=zk;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nyk == -1)
	            			{
						nyk=nyk+block2;
	            				cp=f2;
	        			}
					else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cph;
	                			yf=ami[da][1]-cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
					nxk=xk+1;
	        			nyk=yk+1;
	        			nzk=zk;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nyk >= block2)
	            			{
						nyk=0;
	            				cp=f2;
	        			}
					else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cph;
	                			yf=ami[da][1]+cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk-1;
	        			nyk=yk+1;
	        			nzk=zk;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nyk >= block2)
	            			{
						nyk=0;
	            				cp=f2;
	        			}
					else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cph;
	                			yf=ami[da][1]+cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			nxk=xk-1;
	        			nyk=yk-1;
	        			nzk=zk;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nyk == -1)
	            			{
						nyk=nyk+block2;
	            				cp=f2;
	        			}
					else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cph;
	                			yf=ami[da][1]-cp;
	                			zf=ami[da][2];
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk+1;
	        			nyk=yk;
	        			nzk=zk-1;
	        			if (nxk >= block1)
	        			{ 
		   				nxk=0;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
					{
	            				nzk=nzk+block3;
	            				cp=f3;
					}
				        else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cph;
	                			yf=ami[da][1];
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  

					nxk=xk+1;
	        			nyk=yk;
	        			nzk=zk+1;
	        			if (nxk >= block1)
	        			{ 
		   				nxk=0;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
					{
	            				nzk=0;
	            				cp=f3;
					}
				        else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cph;
	                			yf=ami[da][1];
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	    
	       			nxk=xk-1;
	        			nyk=yk;
	        			nzk=zk+1;
	        			if (nxk == -1)
	        			{ 
		   				nxk=nxk+block1;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
					{
	            				nzk=0;
	            				cp=f3;
					}
				        else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cph;
	                			yf=ami[da][1];
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	       				nxk=xk-1;
	        			nyk=yk;
	        			nzk=zk-1;
	        			if (nxk == -1)
	        			{ 
		   				nxk=nxk+block1;
	            				cph=f1;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
					{
	            				nzk=nzk+block3;
	            				cp=f3;
					}
				        else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cph;
	                			yf=ami[da][1];
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk;
	        			nyk=yk+1;
	        			nzk=zk-1;
	        			if (nyk >= block2)
	        			{
			 			nyk=0;
	                 			cph=f2;
					}
	       	 			else cph=0.0;
	        			if (nzk == -1)
	        			{
		    				nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
					
					nxk=xk;
	        			nyk=yk+1;
	        			nzk=zk+1;
	        			if (nyk >= block2)
	        			{
			 			nyk=0;
	                 			cph=f2;
					}
	       	 			else cph=0.0;
	        			if (nzk >= block3)
	        			{
		    				nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	
	        			nxk=xk;
	        			nyk=yk-1;
	        			nzk=zk+1;
	        			if (nyk == -1)
	        			{
			 			nyk=nyk+block2;
	                 			cph=f2;
					}
	       	 			else cph=0.0;
	        			if (nzk >= block3)
	        			{
		    				nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			nxk=xk;
	        			nyk=yk-1;
	        			nzk=zk-1;
	        			if (nyk == -1)
	        			{
			 			nyk=nyk+block2;
	                 			cph=f2;
					}
	       	 			else cph=0.0;
	        			if (nzk == -1)
	        			{
		    				nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0];
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        
	        			nxk=xk+1;
	        			nyk=yk+1;
	        			nzk=zk+1;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk >= block2)
	       				{    
						nyk=0;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
	        			{    
						nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cpsh;
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk+1;
	        			nyk=yk+1;
	        			nzk=zk-1;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk >= block2)
	       				{    
						nyk=0;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
	        			{    
						nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cpsh;
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk+1;
	        			nyk=yk-1;
	        			nzk=zk+1;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk == -1)
	       				{    
						nyk=nyk+block2;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
	        			{    
						nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cpsh;
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk-1;
	        			nyk=yk+1;
	        			nzk=zk+1;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk >= block2)
	       				{    
						nyk=0;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
	        			{    
						nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cpsh;
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk-1;
	        			nyk=yk-1;
	        			nzk=zk+1;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk == -1)
	       				{    
						nyk=nyk+block2;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk >= block3)
	        			{    
						nzk=0;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cpsh;
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]+cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk-1;
	        			nyk=yk+1;
	        			nzk=zk-1;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk >= block2)
	       				{    
						nyk=0;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
	        			{    
						nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cpsh;
	                			yf=ami[da][1]+cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk-1;
	        			nyk=yk-1;
	        			nzk=zk-1;
	        			if (nxk == -1)
					{
	            				nxk=nxk+block1;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk == -1)
	       				{    
						nyk=nyk+block2;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
	        			{    
						nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]-cpsh;
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					}  
	        			
					nxk=xk+1;
	        			nyk=yk-1;
	        			nzk=zk-1;
	        			if (nxk >= block1)
					{
	            				nxk=0;
	           	 			cpsh=f1;
					}
	        			else cpsh=0.0;
	        			if (nyk == -1)
	       				{    
						nyk=nyk+block2;
	           			 	cph=f2;
					}
	        			else cph=0.0;
	        			if (nzk == -1)
	        			{    
						nzk=nzk+block3;
	            				cp=f3;
					}
	        			else cp=0.0;  
					dum=wdr[nxk][nyk][nzk];
					for (jdt=0; jdt<dum; jdt++)
					{
	                			da=dri[nxk][nyk][nzk][jdt];  rpid2=pid[da];         xf=ami[da][0]+cpsh;
	                			yf=ami[da][1]-cph;
	                			zf=ami[da][2]-cp;
	                			if (da != k)
						{
	                    				rijs=((xf-xi)*(xf-xi))+((yf-yi)*(yf-yi))+((zf-zi)*(zf-zi));  
	                                                if (rijs <= ds[rpid1][rpid2])
							{
	                        				trim=1;jet=jet+1;
	                        				break;
							}   
	                   				else if (rijs <= dss[rpid1][rpid2])
							{
	                       					pf=1;
	                        				pit[wpi]=(dri[nxk][nyk][nzk][jdt]);
								wpi=wpi+1;
							}
						}
					} 
					 
					
					 
					if (wpi>19)
					{
						//printf("wpi %d\n", wpi);
					}
					
					
					
					
					if (trim == 2)
					{
						
						
						del1[0]=el[k][0];
						del1[1]=el[k][1];
						del1[2]=el[k][2];
						del1[3]=el[k][3];
						del1[4]=el[k][4];
						del1[5]=el[k][5];
						del1[6]=el[k][6];
						del1[7]=el[k][7];
						del1[8]=el[k][8];
						
						
						
	 					for (ip=0; ip<wpi; ip++)
						{
							//printf("wpi %d\n", wpi);
						
							w1=li[k][0];
							w2=li[k][1];
							w3=li[k][2];
							w4=li[k][3];
							w5=li[k][4];
							w6=li[k][5];
							w7=li[k][6];
							w8=li[k][7];
							w9=li[k][8];
						
							puc=pit[ip];
							rpid2=pid[puc];	
										
							u1=li[puc][0];
							u2=li[puc][1];
							u3=li[puc][2];
							u4=li[puc][3];
							u5=li[puc][4];
							u6=li[puc][5];
							u7=li[puc][6];
							u8=li[puc][7];
							u9=li[puc][8];
							re1=ami[puc][0]-xi;
							re2=ami[puc][1]-yi;
							re3=ami[puc][2]-zi;
							if (re1 < -pbc) re1=re1+f1;
							else if (re1 > pbc) re1=re1-f1;
							if (re2 < -pbc) re2=re2+f2;
							else if (re2 > pbc) re2=re2-f2;
							if (re3 < -pbc) re3=re3+f3;
							else if (re3 > pbc) re3=re3-f3;
							
							
							mr=(re1*re1+re2*re2+re3*re3);
							mrp=pow(mr,0.5);
							
							
							ree[0]=(re1*w1+re2*w2+re3*w3);
							ree[1]=(re1*w4+re2*w5+re3*w6);
							ree[2]=(re1*w7+re2*w8+re3*w9);
							
							rtm[0][0]=(u1*w1+u2*w2+u3*w3);
							rtm[1][0]=(u1*w4+u2*w5+u3*w6);
							rtm[2][0]=(u1*w7+u2*w8+u3*w9);
							rtm[0][1]=(u4*w1+u5*w2+u6*w3);
							rtm[1][1]=(u4*w4+u5*w5+u6*w6);
							rtm[2][1]=(u4*w7+u5*w8+u6*w9);
							rtm[0][2]=(u7*w1+u8*w2+u9*w3);
							rtm[1][2]=(u7*w4+u8*w5+u9*w6);
							rtm[2][2]=(u7*w7+u8*w8+u9*w9);
							
							for (i = 0; i < 3; i++)
								for (j = 0; j < 3; j++)
									abrtm[i][j] = fabs(rtm[i][j]) + EPSILON;
							
							
							tri=0;
							
							for (i = 0; i < 3; i++) 
							{
								ra = e[rpid1][i];
								rb = e[rpid2][0] * abrtm[i][0] + e[rpid2][1] * abrtm[i][1] + e[rpid2][2] * abrtm[i][2];
								if (fabs(ree[i]) > ra + rb) tri=1;
							}
							
							
							for (i = 0; i < 3; i++) 
							{
								ra = e[rpid1][0] * abrtm[0][i] + e[rpid1][1] * abrtm[1][i] + e[rpid1][2] * abrtm[2][i];
								rb = e[rpid2][i];
								if (fabs(ree[0] * rtm[0][i] + ree[1] * rtm[1][i] + ree[2] * rtm[2][i]) > ra + rb) tri=1;
							}
							
							
							
							
							ra= e[rpid1][1] * abrtm[2][0] + e[rpid1][2] * abrtm[1][0];
							rb= e[rpid2][1] * abrtm[0][2] + e[rpid2][2] * abrtm[0][1];
							if (fabs(ree[2] * rtm[1][0] - ree[1] * rtm[2][0]) > ra + rb) tri=1;
							
							ra= e[rpid1][1] * abrtm[2][1] + e[rpid1][2] * abrtm[1][1];
							rb= e[rpid2][0] * abrtm[0][2] + e[rpid2][2] * abrtm[0][0];
							if (fabs(ree[2] * rtm[1][1] - ree[1] * rtm[2][1]) > ra + rb) tri=1;
							
							ra= e[rpid1][1] * abrtm[2][2] + e[rpid1][2] * abrtm[1][2];
							rb= e[rpid2][0] * abrtm[0][1] + e[rpid2][1] * abrtm[0][0];
							if (fabs(ree[2] * rtm[1][2] - ree[1] * rtm[2][2]) > ra + rb) tri=1;
							
							ra = e[rpid1][0] * abrtm[2][0] + e[rpid1][2] * abrtm[0][0];
							rb = e[rpid2][1] * abrtm[1][2] + e[rpid2][2] * abrtm[1][1];
							if (fabs(ree[0] * rtm[2][0] - ree[2] * rtm[0][0]) > ra + rb) tri=1;
						
							ra= e[rpid1][0] * abrtm[2][1] + e[rpid1][2] * abrtm[0][1];
							rb= e[rpid2][0] * abrtm[1][2] + e[rpid2][2] * abrtm[1][0];
							if (fabs(ree[0] * rtm[2][1] - ree[2] * rtm[0][1]) > ra + rb) tri=1;
							
							ra= e[rpid1][0] * abrtm[2][2] + e[rpid1][2] * abrtm[0][2];
							rb= e[rpid2][0] * abrtm[1][1] + e[rpid2][1] * abrtm[1][0];
							if (fabs(ree[0] * rtm[2][2] - ree[2] * rtm[0][2]) > ra + rb) tri=1;
							
							ra= e[rpid1][0] * abrtm[1][0] + e[rpid1][1] * abrtm[0][0];
							rb= e[rpid2][1] * abrtm[2][2] + e[rpid2][2] * abrtm[2][1];
							if (fabs(ree[1] * rtm[0][0] - ree[0] * rtm[1][0]) > ra + rb) tri=1;
								
							ra= e[rpid1][0] * abrtm[1][1] + e[rpid1][1] * abrtm[0][1];
							rb= e[rpid2][0] * abrtm[2][2] + e[rpid2][2] * abrtm[2][0];
							if (fabs(ree[1] * rtm[0][1] - ree[0] * rtm[1][1]) > ra + rb) tri=1;
							
							ra= e[rpid1][0] * abrtm[1][2] + e[rpid1][1] * abrtm[0][2];
							rb= e[rpid2][0] * abrtm[2][1] + e[rpid2][1] * abrtm[2][0];
							if (fabs(ree[1] * rtm[0][2] - ree[0] * rtm[1][2]) > ra + rb) tri=1;
							
							
							
							//printf("tri %d\n", tri);
							if (tri==0)
							{
								//printf("tri %d\n", tri);
								w1=li[k][6];
								w2=li[k][7];
								w3=li[k][8];
								u1=li[puc][6];
								u2=li[puc][7];
								u3=li[puc][8];
								del2[0]=el[puc][0];
								del2[1]=el[puc][1];
								del2[2]=el[puc][2];
								del2[3]=el[puc][3];
								del2[4]=el[puc][4];
								del2[5]=el[puc][5];
								del2[6]=el[puc][6];
								del2[7]=el[puc][7];
								del2[8]=el[puc][8];
								
								a=0.1;
								b=1.0;
								c=(b+a)/2.0;
								
								trig=1.0;
								
								cfun[0]=0.0;
								cfun[1]=0.0;
								tri=0;
								lc=0;
								while (trig>0.0001)
								{
									tri=tri+1;
									lc=lc+1;
									
									lemm=1.0-c;
									lem=c;
									lemms=lemm*lemm;
									lems=lem*lem;
									wt=0;
									
									for (i=0; i<3; i++)
									{
										for (j=0; j<3; j++)
										{
											lt[i][j]=del1[wt]*lemm+del2[wt]*lem;
											wt=wt+1;
										}
									}
									
									wt=0;
									
									for (i=0; i<3; i++)
									{
										for (j=0; j<3; j++)
										{
											lts[i][j]=del1[wt]*lemms-del2[wt]*lems;
											wt=wt+1;
										}
									}
									lt[0][3]=re1;
									lt[1][3]=re2;
									lt[2][3]=re3;
									cont=1;
									for (it=0; it<2; it++) 
									{
									    	dummy=0;
									    	if  (lt[it][it]==0.0)
										{
											for (jt=0; jt<3; jt++) 
											{
										    		if (dummy == 0)	
												{
													if (lt[jt][it] != 0.0) 
													{
											    			mt=jt;
											    			dummy=dummy+1;
											    			dump[0]=lt[it][0];
														dump[1]=lt[it][1];
														dump[2]=lt[it][2];
														dump[3]=lt[it][3];
											    			lt[it][0]=lt[mt][0];
														lt[it][1]=lt[mt][1];
														lt[it][2]=lt[mt][2];
														lt[it][3]=lt[mt][3];
											    			lt[mt][0]=dump[0];
														lt[mt][1]=dump[1];
														lt[mt][2]=dump[2];
														lt[mt][3]=dump[3];
													}
												}
											}
											if (dummy==0) continue; 
										}
									    	for (jt=cont; jt<3; jt++)
										{
											dumm=lt[jt][it];
											for (mt=0; mt<4; mt++)
											{
										    		lt[jt][mt]=lt[jt][mt]-((lt[it][mt]/lt[it][it])*dumm);
											}
										}
									    	cont=cont+1;
									}
									//NOW CHECK EXISTENCE OF SOLUTION OF EQUATION
									sol[0]=0.0; sol[1]=0.0; sol[2]=0.0; kt=3;
									for (it=2; it>=0; it--)
									{    	sol[it]=lt[it][kt];
									    	for (jt=0; jt<3; jt++)
										{
											if (it==jt) continue;
											sol[it]=sol[it]-(lt[it][jt]*sol[jt]);
										}
									   	sol[it]=sol[it]/lt[it][it];
									}
									ec=0.0;
									for (i=0; i<3; i++)
									{
										ec=ec+(lts[i][0]*sol[0]+lts[i][1]*sol[1]+lts[i][2]*sol[2])*sol[i];
									}
										
									/*if (tri<5)
									{	
										if (ec>0.0)
										{
											a=c;
											cfun[1]=ec;
										}
					    					else
					    					{
											b=c;
											cfun[2]=ec;
										}
										c=(a+b)/2.0;
										printf("lccccccccc    %lf %lf\n",c,b); 
									}
									else
									{
										cfun[0]=cfun[1];
										cfun[1]=cfun[2];
										cfun[2]=ec;
										mrp=((cfun[1]*cfun[2]*a)/((a-b)*(a-c)));
										mr=((cfun[0]*cfun[1]*c)/((c-a)*(c-b))) ;
										red=mrp+mr+  ((cfun[0]*cfun[2]*b)/((b-a)*(b-c)));
		
										a=b;
										b=c;
										c=red;
									}*/
									
									
									
									
									if (ec>0.0)
									{
										a=c;
										cfun[1]=ec;
									}
				    					else
				    					{
										b=c;
										cfun[2]=ec;
									}
									c=(a+b)/2.0;
										
										
				    					trig=fabs(b-a);
				    					ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
				    					//if (ec>1.0) trig=0.000001;
				    						
									
								}
								
								ec=(sol[0]*re1+sol[1]*re2+sol[2]*re3)*lemm*lem;
								//printf("lccccccccc%lf\n",cfun[2]);
								
								//printf("tri %d\n", tri);
								
								if (ec<1.0)
								{
									trim=1;
									jet=jet+1;
								}
								else
								{	
								
									
									
									if (pop[rpid1][rpid2]==1)
									
									{
										ec=ec/mr;
										ec=ec*(mrp-epsi1[rpid1][rpid2])*(mrp-epsi1[rpid1][rpid2]);
										if (ec<1.0)
										{
										
											
											
											zt1=(re1*w1+re2*w2+re3*w3)/mrp;
											zt2=(re1*u1+re2*u2+re3*u3)/mrp;
											if ( zt1 > 0.0 && zt2 < 0.0 )
											{
												psbobb[wpbobb]=puc;
												wpbobb=wpbobb+1;	
											}
											
										}
									
									}
									
									
									else if (pop[rpid1][rpid2]==2)
									
									{
										ec=ec/mr;
										ec=ec*(mrp-epsi2[rpid1][rpid2])*(mrp-epsi2[rpid1][rpid2]);
										if (ec<1.0)
										{
										
											//For patch
										
											zt1=(fabs(re1*w1+re2*w2+re3*w3))/mrp;
											zt2=(fabs(re1*u1+re2*u2+re3*u3))/mrp;
											if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
											{
												psbobb[wpbobb]=puc;
												wpbobb=wpbobb+1;	
											}
											
										}
									
									}
									
									
									else if (pop[rpid1][rpid2]==3)
									{
										ec=ec/mr;
										ec=ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2]);
										if (ec<1.0)
										{
										
											psbobbi[wpbobbi]=puc;
											wpbobbi=wpbobbi+1;	
											
										}
									}
									
									else if (pop[rpid1][rpid2]==4)
									
									{
										ec=ec/mr;
										
										
										
										if ((ec*(mrp-epsi1[rpid1][rpid2])*(mrp-epsi1[rpid1][rpid2])) < 1.0) 
										{
										
											
											zt1=(re1*w1+re2*w2+re3*w3)/mrp;
											zt2=(re1*u1+re2*u2+re3*u3)/mrp;
											
											if ( zt1 > 0.0 && zt2 < 0.0 )
											{
												psbobb[wpbobb]=puc;
												wpbobb=wpbobb+1;
												
											
											}
										}
										
										
										if ((ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2])) < 1.0)
										
										{
											psbobbi[wpbobbi]=puc;
											wpbobbi=wpbobbi+1;
												
										}
										
									
									}
									
									
									else if (pop[rpid1][rpid2]==5)
									
									{
										ec=ec/mr;
											
										if ((ec*(mrp-epsi2[rpid1][rpid2])*(mrp-epsi2[rpid1][rpid2])) < 1.0) 
										{
										
											
											zt1=(fabs(re1*w1+re2*w2+re3*w3))/mrp;
											zt2=(fabs(re1*u1+re2*u2+re3*u3))/mrp;
											if ( zt1 > omega2[rpid1] && zt2 > omega2[rpid2] )
											{
												psbobb[wpbobb]=puc;
												wpbobb=wpbobb+1;
												
											
											}
										}
										
										if ((ec*(mrp-epsi3[rpid1][rpid2])*(mrp-epsi3[rpid1][rpid2])) < 1.0)
										
										{
											psbobbi[wpbobbi]=puc;
											wpbobbi=wpbobbi+1;
												
										}
									
									}
									
									
									
									
								}
								
								
							}
							if (trim==1)
							{
								break;
							}
							
							
							
							
							
							
							//printf("lc       %d\n",lc);
						}
						
						
						
						
					}
					
					if (trim==2)
					{
					
							for (it=0; it<wbobb[k]; it++)
							{
								trim=1;
								for (i=0; i < wpbobb; i++)
								{
									if (bobb[k][it]==psbobb[i]) trim=2;
									
								}
								if (trim==1)  break; // if the loop will keep on without breaking then may be in the next loop we will get the particle and the trim will become 2 so that we needed to break the loop at the moment						
							}
							
							
							
							if (trim==2)
							{
								for (it=0; it<wbobbi[k]; it++)
								{
									trim=1;
									for (i=0; i < wpbobbi; i++)
									{
										if (bobbi[k][it]==psbobbi[i]) trim=2;
										
									}
									if (trim==1)  break; 					
								}
							}
								
							
							
							
		
						
						
					}
					
					if (trim==2)
					{
						dum=wdr[l][m][n];
						for (j=0;j<dum;j++)
						{	
							
							if (dri[l][m][n][j]==k)
							{
								for (tric=j; tric<dum; tric++)
								{
									dri[l][m][n][tric]=dri[l][m][n][tric+1];
									
								}
								break;
							}
									
						}
						wdr[l][m][n]=wdr[l][m][n]-1;
						dum=wdr[xk][yk][zk];
						dri[xk][yk][zk][dum]=k;
						wdr[xk][yk][zk]=wdr[xk][yk][zk]+1;
						ami[k][0]=xi;
		    				fau[k][0]=dxi;
		    				ami[k][1]=yi;
		    				fau[k][1]=dyi;
		    				ami[k][2]=zi;
		    				fau[k][2]=dzi;
						
						sdcxy[k][0]=sdcxy[k][0]+x1per;
						sdcxy[k][1]=sdcxy[k][1]+y1per;
						sdcz[k]=sdcz[k]+tpart;
						sdcxyz[k][0]=sdcxyz[k][0]+dx;
						sdcxyz[k][1]=sdcxyz[k][1]+dy;
						sdcxyz[k][2]=sdcxyz[k][2]+dz;
						
						mdcxy[k][0]=mdcxy[k][0]+x1per;
						mdcxy[k][1]=mdcxy[k][1]+y1per;
						mdcz[k]=mdcz[k]+tpart;
						mdcxyz[k][0]=mdcxyz[k][0]+dx;
						mdcxyz[k][1]=mdcxyz[k][1]+dy;
						mdcxyz[k][2]=mdcxyz[k][2]+dz;
						
						ldcxy[k][0]=ldcxy[k][0]+x1per;
						ldcxy[k][1]=ldcxy[k][1]+y1per;
						ldcz[k]=ldcz[k]+tpart;
						ldcxyz[k][0]=ldcxyz[k][0]+dx;
						ldcxyz[k][1]=ldcxyz[k][1]+dy;
						ldcxyz[k][2]=ldcxyz[k][2]+dz;
						
							
						mum=wnrd[k];
						
						for (jt=0; jt<mum; ++jt)
	    					{
							chum=nrdi[k][jt];
							drum=wnrd[chum];
							for (sj=0; sj<drum; ++sj)
							{
								if (nrdi[chum][sj]==k)
								{
									for (tric=sj; tric<drum; tric++)
									{
										nrdi[chum][tric]=nrdi[chum][(tric+1)];
									}
									break;
								}
							
							}
							wnrd[chum]=drum-1;
							nrdi[k][jt]=6000;
						}
							
						wnrd[k]=0;	
						for (jdt=0; jdt < wpi; jdt++)
						{	
							//printf("raam is ");
							td=pit[jdt];
							nrdi[k][jdt]=td;
								tum=wnrd[td];
								nrdi[td][tum]=k;
							wnrd[td]=tum+1;
							wnrd[k]=wnrd[k]+1;
						}	
						
						
						
							if (cow==1)
							{
								mum=wcbobb[k];
									
								for (jt=0; jt<mum; ++jt)
			    					{
									chum=cbobb[k][jt];
									drum=wcbobb[chum];
									for (sj=0; sj<drum; ++sj)
									{
										if (cbobb[chum][sj]==k)
										{
											for (tric=sj; tric<drum; tric++)
											{
												cbobb[chum][tric]=cbobb[chum][(tric+1)];
											}
											break;
										}
									
									}
									wcbobb[chum]=drum-1;
									cbobb[k][jt]=6000;
								}
									
								wcbobb[k]=0;	
								for (jdt=0; jdt < wpbobb; jdt++)
								{	
									//printf("raam is ");
									td=psbobb[jdt];
									cbobb[k][jdt]=td;
										tum=wcbobb[td];
										cbobb[td][tum]=k;
									wcbobb[td]=tum+1;
									wcbobb[k]=wcbobb[k]+1;
								}
								
								
								
									
								mum=wcbobbi[k];
								
								for (jt=0; jt<mum; ++jt)
			    					{
									chum=cbobbi[k][jt];
									drum=wcbobbi[chum];
									for (sj=0; sj<drum; ++sj)
									{
										if (cbobbi[chum][sj]==k)
										{
											for (tric=sj; tric<drum; tric++)
											{
												cbobbi[chum][tric]=cbobbi[chum][(tric+1)];
											}
											break;
										}
									
									}
									wcbobbi[chum]=drum-1;
									cbobbi[k][jt]=6000;
								}
									
								wcbobbi[k]=0;	
								for (jdt=0; jdt < wpbobbi; jdt++)
								{	
									//printf("raam is ");
									td=psbobbi[jdt];
									cbobbi[k][jdt]=td;
										tum=wcbobbi[td];
										cbobbi[td][tum]=k;
									wcbobbi[td]=tum+1;
									wcbobbi[k]=wcbobbi[k]+1;
								}
								
								
								
								
								
								
								
							}
							
						
						
					}
	
				}
				
				
			//printf("snn750 %d\n",mo);
			}
			
			
			
			for (j=0; j<N; j++)
			{
				for (k=0; k<E;  k++)
				{
					hold[j][k]=6000;
					whold[j]=0;
				}
				
				
				for (k=0; k<wcbobb[j];  k++)
				{
					tri=0;
					for (i=0; i<wbobb[j]; i++)
					{
						if (bobb[j][i] == cbobb[j][k])
						{
							tri=1;
						}
				
					}
					if (tri==0)
					{
						hold[j][whold[j]]=cbobb[j][k];
						whold[j]=whold[j]+1;
					}
					//printf("just one time\n");
				}
				
				
				
				
				
				
					for (k=0; k<E;  k++)
					{
						holdi[j][k]=6000;
						wholdi[j]=0;
					}
					
					
					for (k=0; k<wcbobbi[j];  k++)
					{
						tri=0;
						for (i=0; i<wbobbi[j]; i++)
						{
							if (bobbi[j][i] == cbobbi[j][k])
							{
								tri=1;
							}
					
						}
						if (tri==0)
						{
							holdi[j][wholdi[j]]=cbobbi[j][k];
							wholdi[j]=wholdi[j]+1;
						}
						//printf("just one time\n");
					}
				
				
				
			}	
				
			for (j=0; j<N; j++)
			{
				wbuq=0;
				for (jt=0; jt<E; jt++)
				{
					bumbque[jt]=6000;
				}
				rpid1=pid[j];
            			for (ini=0; ini<wbobb[j]; ini++)
				{
                                   	rd=drand48();
                			bumbque[wbuq]=bobb[j][ini];
					wbuq=wbuq+1;
					rpid2=pid[bobb[j][ini]];
                			if (bobb[j][ini]>j)
					{
						if (rd <= beta[rpid1][rpid2])
						{
                        				bumbque[wbuq-1]=6000;
							wbuq=wbuq-1;
						}
					}
					chum=bobb[j][ini];
					drum=wbobb[chum];
					for (sj=0; sj<drum; ++sj)
					{
						if (bobb[chum][sj]==j)
						{
							for (tric=sj; tric<drum; tric++)
							{
								bobb[chum][tric]=bobb[chum][(tric+1)];
							}
							break;
						}
					
					}
					wbobb[chum]=drum-1;
					bobb[j][ini]=6000;
                        	}
				wbobb[j]=0;	
				for (jt=0; jt < wbuq; jt++)
				{	
					//printf("raam is ");
					td=bumbque[jt];
					bobb[j][jt]=td;
					tum=wbobb[td];
					bobb[td][tum]=j;
					wbobb[td]=tum+1;
					wbobb[j]=wbobb[j]+1;
				}	
				
				
				
				
				
				
				
				
				
					wbuqi=0;
					for (jt=0; jt<E; jt++)
					{
						bumbquei[jt]=6000;
					}
		    			for (ini=0; ini<wbobbi[j]; ini++)
					{
		                           	rd=drand48();
		        			bumbquei[wbuqi]=bobbi[j][ini];
		        			rpid2=pid[bobbi[j][ini]];
						wbuqi=wbuqi+1;
		        			if (bobbi[j][ini]>j)
						{
							if (rd <= betai[rpid1][rpid2])
							{
		                				bumbquei[wbuqi-1]=6000;
								wbuqi=wbuqi-1;
							}
						}
						chum=bobbi[j][ini];
						drum=wbobbi[chum];
						for (sj=0; sj<drum; ++sj)
						{
							if (bobbi[chum][sj]==j)
							{
								for (tric=sj; tric<drum; tric++)
								{
									bobbi[chum][tric]=bobbi[chum][(tric+1)];
								}
								break;
							}
						
						}
						wbobbi[chum]=drum-1;
						bobbi[j][ini]=6000;
		                	}
					wbobbi[j]=0;	
					for (jt=0; jt < wbuqi; jt++)
					{	
						//printf("raam is ");
						td=bumbquei[jt];
						bobbi[j][jt]=td;
						tum=wbobbi[td];
						bobbi[td][tum]=j;
						wbobbi[td]=tum+1;
						wbobbi[j]=wbobbi[j]+1;
					}	
				
				
				
				
				
				
				
				
				
				
				
			
			}
			
			for (j=0; j<N; j++)
			{
				for (k=0; k<whold[j];  k++)
				{
					if (wbobb[j] < bnll)
					{
					
						bobb[j][wbobb[j]]=hold[j][k];
						wbobb[j]=wbobb[j]+1;
					}
					//printf("just one time\n");
				}
				
				
				
				
				for (k=0; k<wholdi[j];  k++)
				{
					if (wbobbi[j] < bnll)
					{
					
						bobbi[j][wbobbi[j]]=holdi[j][k];
						wbobbi[j]=wbobbi[j]+1;
					}
					//printf("just one time\n");
				}
				
				
			}
			
			
			
			
			
			
			//if (pr==1) prtd=prtd+istep;
			
			if (jet==0)
			{
				rul=1;
			}
			
			if (sla==5)
			{
				sla=0;
			}			
			
			
			
		
		
			
						
		}
		printf("Simulation is going on si %d \n", sitt);
		
		
		
		
		
		jt=0;
		
		for (j=0; j<war; j++)
		{
			dif=0.0;
			difu=0.0;
			difus=0.0;
			difusl=0.0;
			jtt=jt+nt[j];
			for (i=jt; i<jtt; i++)
			{	
				dif=dif+(sdcxy[i][0]*sdcxy[i][0]+sdcxy[i][1]*sdcxy[i][1]);
				difu=sdcz[i]*sdcz[i]+difu;
				difus=difus+(sdcxyz[i][0]*sdcxyz[i][0]+sdcxyz[i][1]*sdcxyz[i][1]+sdcxyz[i][2]*sdcxyz[i][2]);
				difusl=difusl+((som[i][0]*som[i][0])+(som[i][1]*som[i][1])+(som[i][2]*som[i][2]));
				
				
				
				
			}
			
			dif=dif/(float)nt[j];
			
			difu=difu/(float)nt[j];

			difus=difus/(float)nt[j];
			
			difusl=(difusl*rdtheta[j]*rdtheta[j])/(float)nt[j];
			
			fprintf(swu,"%d  %lf  %lf  %lf  %lf\n", (sitt+1), dif, difu, difus, difusl);
			
			
			printf("small %d  %lf  %lf  %lf %lf %d\n", (sitt+1), dif, difu, difus, difusl, jet);
			
			for (i=jt; i<jtt; i++)
			{
				sdcxyz[i][0]=0.0;
				sdcxyz[i][1]=0.0;
				sdcxyz[i][2]=0.0;
				sdcxy[i][0]=0.0;
				sdcxy[i][1]=0.0;
				sdcz[i]=0.0;
				som[i][0]=0.0;
				som[i][1]=0.0;
				som[i][2]=0.0;
			}
			jt=jtt;
			
		}
		
		
		for (i=0; i<N; i++)
		{
			fprintf(cfu,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
			
			fprintf(cfu, "%d\n", wbobb[i]);
			for (j=0; j<wbobb[i]; j++) 
			{
				fprintf(cfu,"   %d", bobb[i][j]);
			}
			fprintf(cfu, "\n");
			
			
			fprintf(cfu, "%d\n", wbobbi[i]);
			for (j=0; j<wbobbi[i]; j++) 
			{
				fprintf(cfu,"   %d", bobbi[i][j]);
			}
			fprintf(cfu, "\n");
			
			
			fprintf(cfu,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
			fprintf(cfu,"%lf  %lf  %lf\n", fau[i][0], fau[i][1], fau[i][2]);
		}
		
		
		
		
		
		
		
		
			
		
		
		//fprintf(dwu,"\n");
		
		if (((float)(sitt+1)/10.0)==((float)floor((float)(sitt+1)/10.0)))
		{
		
		
			
			jt=0;
		
			for (j=0; j<war; j++)
			{
				dif=0.0;
				difu=0.0;
				difus=0.0;
				difusl=0.0;
				jtt=jt+nt[j];
				for (i=jt; i<jtt; i++)
				{	
					dif=dif+(mdcxy[i][0]*mdcxy[i][0]+mdcxy[i][1]*mdcxy[i][1]);
					difu=mdcz[i]*mdcz[i]+difu;
					difus=difus+(mdcxyz[i][0]*mdcxyz[i][0]+mdcxyz[i][1]*mdcxyz[i][1]+mdcxyz[i][2]*mdcxyz[i][2]);
					difusl=difusl+((mom[i][0]*mom[i][0])+(mom[i][1]*mom[i][1])+(mom[i][2]*mom[i][2]));
					
					
					
					
				}
				
				dif=dif/(float)nt[j];
				
				difu=difu/(float)nt[j];

				difus=difus/(float)nt[j];
				
				difusl=(difusl*rdtheta[j]*rdtheta[j])/(float)nt[j];
				
				dif=dif/10.0;
				
				difu=difu/10.0;

				difus=difus/10.0;
				
				difusl=difusl/10.0;
				
				fprintf(mwu,"%d  %lf  %lf  %lf  %lf\n", (sitt+1), dif, difu, difus, difusl);
				
				
				printf("small %d  %lf  %lf  %lf %lf %d\n", (sitt+1), dif, difu, difus, difusl, jet);
				
				for (i=jt; i<jtt; i++)
				{
					mdcxyz[i][0]=0.0;
					mdcxyz[i][1]=0.0;
					mdcxyz[i][2]=0.0;
					mdcxy[i][0]=0.0;
					mdcxy[i][1]=0.0;
					mdcz[i]=0.0;
					mom[i][0]=0.0;
					mom[i][1]=0.0;
					mom[i][2]=0.0;
				}
				jt=jtt;
				
			}
			
						
			
			
			
			
			
			
		}
		
		
		
		
		if (((float)(sitt+1)/100.0)==((float)floor((float)(sitt+1)/100.0)))
		{
		
		
		
			fclose(cfu);
			char filename[25] = {0};

    			sprintf(filename, "%d.txt", (int)floor((float)(sitt+1)/100.0));

			FILE* cfu;
			cfu=fopen(filename,"w");
			
			
			
			
			FILE* ku;
			ku=fopen("Relaxed_system.txt","w");
			
			fprintf(ku,"%d\n", sitt+1);
			
			fprintf(ku,"%lf\n", quao);
			
			fprintf(ku,"%lf\n", quaoi);
			
			for (i=0; i<N; i++)
			{
				fprintf(ku,"%lf  %lf  %lf\n", ami[i][0], ami[i][1], ami[i][2]);
				
				fprintf(ku, "%d\n", wnrd[i]);
				for (j=0; j<wnrd[i]; j++) 
				{
					fprintf(ku,"   %d", nrdi[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku, "%d\n", wbobb[i]);
				for (j=0; j<wbobb[i]; j++) 
				{
					fprintf(ku,"   %d", bobb[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku, "%d\n", wbobbi[i]);
				for (j=0; j<wbobbi[i]; j++) 
				{
					fprintf(ku,"   %d", bobbi[i][j]);
				}
				fprintf(ku, "\n");
				
				fprintf(ku,"%lf  %lf  %lf\n", li[i][6], li[i][7], li[i][8]);
				
				fprintf(ku,"%lf  %lf  %lf\n", fau[i][0], fau[i][1], fau[i][2]);
			}
			
			fclose(ku);
			
			
			
			
			
			
		}
		
		if (((float)(sitt+1)/1000.0)==((float)floor((float)(sitt+1)/1000.0)))
		{
		
		
			qua=0.0;
	
			for (i=0; i<N; i++)
			{
				qua=qua+(float)wbobb[i];
			
			}
			qua=qua/(float)N;
			
			quai=0.0;
	
			for (i=0; i<N; i++)
			{
				quai=quai+(float)wbobbi[i];
			
			}
			quai=quai/(float)N;

			
				if (fabs(qua-quao) < 0.05 && fabs(quai-quaoi) < 0.05) { printf("please don't run leave me alone\n"); break;}
				quao=qua;
				quaoi=quai;
		
		
		
		}
		
		
		
		

	}
	
	
	jt=0;
		
	for (j=0; j < war; j++)
	{
		dif=0.0;
		difu=0.0;
		difus=0.0;
		difusl=0.0;
		jtt=jt+nt[j];
		for (i=jt; i<jtt; i++)
		{	
			dif=dif+(ldcxy[i][0]*ldcxy[i][0]+ldcxy[i][1]*ldcxy[i][1]);
			difu=ldcz[i]*ldcz[i]+difu;
			difus=difus+(ldcxyz[i][0]*ldcxyz[i][0]+ldcxyz[i][1]*ldcxyz[i][1]+ldcxyz[i][2]*ldcxyz[i][2]);
			difusl=difusl+((lom[i][0]*lom[i][0])+(lom[i][1]*lom[i][1])+(lom[i][2]*lom[i][2]));
			
			
			
			
		}
		
		dif=dif/(float)nt[j];
		
		difu=difu/(float)nt[j];

		difus=difus/(float)nt[j];
		
		difusl=(difusl*rdtheta[j]*rdtheta[j])/(float)nt[j];
		
		dif=dif/10.0;
		
		difu=difu/10.0;

		difus=difus/10.0;
		
		difusl=difusl/10.0;
		
		dif=dif/((float)sitt-20.0-sitts);
	
		difu=difu/((float)sitt-20.0-sitts);

		difus=difus/((float)sitt-20.0-sitts);
		
		difusl=difusl/((float)sitt-20.0-sitts);
		
		
		fprintf(lwu,"%d  %lf  %lf  %lf  %lf\n", (sitt+1), dif, difu, difus, difusl);
		
		
		printf("small %d  %lf  %lf  %lf %lf %d\n", (sitt+1), dif, difu, difus, difusl, jet);
		
		for (i=jt; i<jtt; i++)
		{
			ldcxyz[i][0]=0.0;
			ldcxyz[i][1]=0.0;
			ldcxyz[i][2]=0.0;
			ldcxy[i][0]=0.0;
			ldcxy[i][1]=0.0;
			ldcz[i]=0.0;
			lom[i][0]=0.0;
			lom[i][1]=0.0;
			lom[i][2]=0.0;
		}
		jt=jtt;
		
	}
	
	
	
	
	
	
	
	fclose (swu);
	fclose (mwu);
	fclose (lwu);
	fclose (cfu);

	//SImulaton ends!





	//FINAL CORD Zone starts!

	
	
	//FINAL CORD Zone ends!







	//VTK Zone starts!
	
	FILE* fu;
     	fu=fopen("Cord.vtk","w");
     	
	fprintf(fu,"# vtk DataFile Version 3.0\n");
	fprintf(fu,"Random data to test tensors\n");
	fprintf(fu,"ASCII\n");
	fprintf(fu,"DATASET POLYDATA\n");
	fprintf(fu,"POINTS %d float\n", N);
	
	for (jt=0; jt<N; jt++)
	{
		for(i=0;i<3;i++)
		{
	  		fprintf(fu,"%lf   ",ami[jt][i]);
		}
		fprintf(fu,"\n");
	}

	fprintf(fu,"\n");
	fprintf(fu,"POINT_DATA %d\n", N);
	fprintf(fu,"\n");
	fprintf(fu,"\n");
	fprintf(fu,"TENSORS non-spherical_ellipsoid float\n");
	fprintf(fu,"\n");
	
	for (i=0; i<N; i++)
	{
		rst=r[pid[i]];
		rsmt=rm[pid[i]];
		del[0]=(rsmt*li[i][0]*li[i][0])+(rsmt*li[i][3]*li[i][3])+(rst*li[i][6]*li[i][6]);
		del[1]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
		del[2]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
		del[3]=(rsmt*li[i][0]*li[i][1])+(rsmt*li[i][3]*li[i][4])+(rst*li[i][6]*li[i][7]);
		del[4]=(rsmt*li[i][1]*li[i][1])+(rsmt*li[i][4]*li[i][4])+(rst*li[i][7]*li[i][7]);
		del[5]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
		del[6]=(rsmt*li[i][0]*li[i][2])+(rsmt*li[i][3]*li[i][5])+(rst*li[i][6]*li[i][8]);
		del[7]=(rsmt*li[i][1]*li[i][2])+(rsmt*li[i][4]*li[i][5])+(rst*li[i][7]*li[i][8]);
		del[8]=(rsmt*li[i][2]*li[i][2])+(rsmt*li[i][5]*li[i][5])+(rst*li[i][8]*li[i][8]);
		
		fprintf(fu,"%lf    %lf    %lf\n",del[0], del[1], del[2]);
		fprintf(fu,"%lf    %lf    %lf\n",del[3], del[4], del[5]);
		fprintf(fu,"%lf    %lf    %lf\n",del[6], del[7], del[8]);
		fprintf(fu,"\n");
	}
	
	
	
	fclose(fu);

	//VTK Zone ends!




	
	//The war of Freedom started!
	

	for (i=0;i<N;++i)
	{
		free(li[i]);
	}
		
	free(li);
	
	
	for (i=0;i<N;++i)
	{
		free(el[i]);
	}

	free(el);
	
	for (i=0;i<N;++i)
	{
		free(ami[i]);
	}
		
	
	free(ami);
	
	
	for (i=0;i<N;++i)
	{
		free(fau[i]);
	}
	free(fau);
	
	
	
	
	
	
	/*for (i=0;i<N;++i)
	{
		free(som[i]);
	}
		
	free(som);
	

	for (i=0;i<N;++i)
	{
		free(sdcxyz[i]);
	}
	
	free(sdcxyz);
	
	for (i=0;i<N;++i)
	{
		free(sdcxy[i]);
	}
		
	free(sdcxy);
	
	
	free(sdcz);
	
	
	
	for (i=0;i<N;++i)
	{
		free(mom[i]);
	}
		
	free(mom);
	

	for (i=0;i<N;++i)
	{
		free(mdcxyz[i]);
	}
	
	free(mdcxyz);
	
	for (i=0;i<N;++i)
	{
		free(ldcxy[i]);
	}
		
	free(ldcxy);
	
	
	free(ldcz);
	
	
	
	for (i=0;i<N;++i)
	{
		free(lom[i]);
	}
		
	free(lom);
	

	for (i=0;i<N;++i)
	{
		free(ldcxyz[i]);
	}
	
	free(ldcxyz);
	
	for (i=0;i<N;++i)
	{
		free(ldcxy[i]);
	}
		
	free(mdcxy);
	
	
	free(mdcz);
	*/
	
	for (k=0; k<N; k++)
	{
		free(nrdi[k]);
	
	}
		
	free(nrdi);
	free(wnrd);
	
	
		
	for (i=0; i<blockd; i++)
	{
		
		for (j=0; j<blockd; j++)
		{
		
			for (k=0; k<blockd; k++)
			{
				free(dri[i][j][k]);
			}
			free(dri[i][j]);
			free(wdr[i][j]);
		}
		free(dri[i]);
		free(wdr[i]);
		
	}
	free(dri);
	free(wdr);
	

	//Now freedom is achieved!	


	
}





