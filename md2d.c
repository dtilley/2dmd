#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// partical typedef
typedef struct {
  double x_now,x_new,y_now,y_new;
  double vx_now,vx_new,vy_now,vy_new;
  double fx_now,fx_new,fy_now,fy_new;
  int id;
} atom; 


// things that should be in a header
int n; 
double m=1;
double epsilon=1;
double sigma=1;
double alpha=0.98;
double kb=1;
double lx,ly;
int state=0;
double eta1,eta2,ksi1,ksi2;
double drandmax= (double) RAND_MAX;
double pi=3.141592653589793;

// psudorandom double generator (0,1]
double drand(){
  return (double) (1+rand())/(1+drandmax);
}
// random number from guassian distribution
double gaussian(){
  if (state==0){
    ksi1=drand();
    ksi2=drand();
    eta1=sqrt(-2*log(ksi1))*cos(2*pi*ksi2);
    eta2=sqrt(-2*log(ksi1))*sin(2*pi*ksi2);
    state=1;
    return eta1;
  }else{
    state=0;
    return eta2;
  }
}

// in retrospect this is aweful 
void setpostion(double *x, double *y){
  int i = 0;
  int isodd=0;
  double gridx, gridy;
  if (n % 2 != 0){
    x[i]=lx/2.0;
    y[i]=ly/2.0;
    i++;
    isodd=1;
    if (n==1){
      return;
    }
  }
  if (isodd==1 && (n-1)%4==0){
    // test case atoms=5,9
    gridx=lx/(n-1);
    gridy=ly/(n-1);
  }else if (n%4==0){
    // test case atoms=4,8
    gridx=lx/n;
    gridy=ly/n;
  }else if (n%4!=0 && n%2==0){
    // populate the midline test case atoms=10,6,2
    //printf("right test case");
    if((n-2)!=0){
      gridx=lx/(n-2);
      gridy=ly/(n-2);
    }else{
      x[i]=lx/3;
      y[i]=ly/2;
      i++;
      x[i]=lx-lx/3;
      y[i]=ly/2;
      return;
    }
    x[i]=lx/3;
    y[i]=ly/2;
    i++;
    x[i]=lx-lx/3;
    y[i]=ly/2;
    i++;
  }else if (isodd==1 && (n-1)%4!=0 && (n-1)%2==0){
    // populate the midline in the test case atoms=9,7,3
    if((n-3)!=0){
      gridx=lx/(n-3);
      gridy=ly/(n-3);
    }else{
      gridx=lx/n;
      gridy=ly/n;
    }
    x[i]=gridx;
    y[i]=ly/2;
    i++;
    x[i]=lx-gridx;
    y[i]=ly/2;
    i++;
  }
  int j=n/4;
  int k=1;
  while (k<=j){
    x[i]=gridx*k;
    y[i]=gridy*k;
    i++;
    x[i]=lx-gridx*k;
    y[i]=gridy*k;
    i++;
    x[i]=gridx*k;
    y[i]=ly-gridy*k;
    i++;
    x[i]=lx-gridx*k;
    y[i]=ly-gridy*k;
    i++;
    k++;
  }
  return;
}

// Initiate velocities
void setvelocity(double v,double *vxy){
  // Velocity faction in x-direction
  double vsquared=v*v;
  *vxy= sqrt((vsquared)*drand());
  *(vxy+1)= sqrt(vsquared-(*vxy)*(*vxy));
  // +/- direction
  if(drand()<0.5){
    *(vxy)= *(vxy) *-1.0;
  }
  if(drand()<0.5){
    *(vxy+1)= *(vxy+1) * -1.0;
  }
}

void jump(double T, double *fxy, double dt){
  double bsquared,beta;
  beta=sqrt((2*alpha*kb*T)/dt)*gaussian();
  bsquared=beta*beta;
  *fxy=sqrt((bsquared)*drand());
  *(fxy+1)= sqrt(bsquared-(*fxy)*(*fxy));
  // +/- direction
  if(drand()<0.5){
    *(fxy)= *(fxy) *-1.0;
  }
  if(drand()<0.5){
    *(fxy+1)= *(fxy+1) * -1.0;
  }
}


// Calculate wall interaction
double wall_repulsion(double *xy, double *fxy){
  double U;
  U=epsilon*pow(sigma,12)*(pow(*(xy),-12)+pow((*(xy)-lx),-12));
  U=U+epsilon*pow(sigma,12)*(pow(*(xy+1),-12)+pow((*(xy+1)-ly),-12));
  *(fxy)=epsilon*pow(sigma,12)*12.0*(pow(*(xy),-13)+pow((*(xy)-lx),-13));
  *(fxy+1)=epsilon*pow(sigma,12)*12.0*(pow(*(xy+1),-13)+pow((*(xy+1)-ly),-13));
  return U;
}

// Calculate pair interactions
double lj_potential(double *atm1, double *atm2, double *fxy, int periodic){
  double r_mag,f_mag,x,y,absx,absy,U;
  x=*atm2-*atm1;
  y=*(atm2+1)-*(atm1+1);
  absx=sqrt(x*x);
  absy=sqrt(y*y);
  if (periodic==1 && (lx/2)<absx){
    if (x>0){
      x=x-lx;
    }else{
      x=lx+x;
    }
  }
  if (periodic==1 && (ly/2)<absy){
    if (y>0){
      y=y-ly;
    }else{
      y=ly+y;
    }
  }
  r_mag=sqrt(x*x + y*y);
  U=pow(sigma,6)*(pow(r_mag,-6));
  U= epsilon*(U*U - 2*U)+epsilon; 
  f_mag=12*epsilon*(pow(sigma,12)*(pow(r_mag,-13)-pow(sigma,6)*pow(r_mag,-7)));
  *fxy=f_mag*(x/r_mag);
  *(fxy+1)=f_mag*(y/r_mag);
  return U;
}

int main(int argc, char *argv[] ){
  int h,periodic,c; 
  double e,KE,dt;
  double T=0;
  double U=0;
  double t=0;
  double vx,vy;
  double xy_fpair[2];
  double fxy[2];
  double vxy[2];
  double xy1[2];
  double xy2[2];

  if (argc != 14){
    printf("useage: %s -nsteps (int steps) -natoms (int n) -dt (double dt) -e (double KE) -box (double lx ly ) -periodic (int 1 or 0) \n", argv[0]);
    return 0;
  }
  else {
    for (int i=1;i<argc;i++){
      if (strcmp(argv[i],"-nsteps")==0){
	h=atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-natoms")==0){
	n=atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-dt")==0){
	dt=atof(argv[i+1]);
      } else if(strcmp(argv[i],"-e")==0){
	e=atof(argv[i+1]);
      } else if(strcmp(argv[i],"-box")==0){
	lx=atof(argv[i+1]);
	ly=atof(argv[i+2]);
      } else if(strcmp(argv[i],"-periodic")==0){
	periodic=atoi(argv[i+1]);
      }
    }
  }

  // Create particles and set ID,v,r
  e=e/n;
  atom prtcls[n];
  double xs[n];
  double ys[n];
  double v_mag=0;
  setpostion(xs,ys);
  for(int i=0;i<n;i++){
    prtcls[i].id=i;
    setvelocity(e,vxy);
    prtcls[i].vx_now=vxy[0];
    prtcls[i].vy_now=vxy[1];
    prtcls[i].x_now=xs[i];
    prtcls[i].y_now=ys[i];
    v_mag=v_mag+sqrt(prtcls[i].vy_now*prtcls[i].vy_now + prtcls[i].vx_now*prtcls[i].vx_now);
  }
  T=(m/kb)*(1.0/3.0)*((v_mag*v_mag)/n);
  KE=0.5*m*(v_mag*v_mag);

  // Calculate initial forces

    if (periodic==0){ 
      for (int j=0;j<n;j++){
	// get coords for particle
	xy1[0]=prtcls[j].x_now;
	xy1[1]=prtcls[j].y_now;
	U=U+wall_repulsion(xy1,fxy);
	// store forces associated with the particle
	prtcls[j].fx_now=fxy[0];
	prtcls[j].fy_now=fxy[1];
      }
    }
    c=0;
    while (c<n){
      // get coords for particle 1
      xy1[0]=prtcls[c].x_now;
      xy1[1]=prtcls[c].y_now;
      for (int j=0;j<n;j++){
	if (j!=c){
	  // get coords for particle 2
	  xy2[0]=prtcls[j].x_now;
	  xy2[1]=prtcls[j].y_now;
	  // store forces associated with the particle 
	  U=U+0.5*lj_potential(xy1,xy2,xy_fpair,periodic);
	  if (periodic==0){
	    prtcls[c].fx_now=prtcls[c].fx_now+xy_fpair[0];
	    prtcls[c].fy_now=prtcls[c].fy_now+xy_fpair[1];
	  }else {
	    prtcls[c].fx_now=xy_fpair[0];
	    prtcls[c].fy_now=xy_fpair[1];
	    // periodic and more than 1 particle, so bath interations and dissipations
	    jump(T,fxy,dt);
	    prtcls[c].fx_now=prtcls[c].fx_now+fxy[0];
	    prtcls[c].fy_now=prtcls[c].fy_now+fxy[1];
	  }
	}
      }
      c++;
    }
    // Print Energies
    printf("Time KE U\n");
    printf("%.5f %.5f %.5f\n",t,KE,U);


  // run md2d
  for (int i=0;i<h;i++){

    // Time tracker
    t=t+dt;
    U=0;

    // Calculate new postion
    for (int j=0;j<n;j++){
      /*print postion time
      if (j< n-1){
	printf("%.5f %.5f ",prtcls[j].vx_now,prtcls[j].vy_now);
      }else{
	printf("%.5f %.5f\n",prtcls[j].vx_now,prtcls[j].vy_now);
      }
      */

      prtcls[j].x_new=prtcls[j].x_now+dt*prtcls[j].vx_now+((dt*dt)/(2*m))*prtcls[j].fx_now;
      prtcls[j].y_new=prtcls[j].y_now+dt*prtcls[j].vy_now+((dt*dt)/(2*m))*prtcls[j].fy_now;

      if (periodic==1){
	if (prtcls[j].x_new > lx){
	  prtcls[j].x_new=prtcls[j].x_new-lx;
	}else if (prtcls[j].x_new < 0){
	  prtcls[j].x_new=prtcls[j].x_new+lx;
	}
	if (prtcls[j].y_new > ly){
	  prtcls[j].y_new=prtcls[j].y_new-ly;
	}else if (prtcls[j].y_new < 0){
	  prtcls[j].y_new=prtcls[j].y_new+ly;
	}
      }
    }

    // Calculate new forces
    if (periodic==0){ 
      for (int j=0;j<n;j++){
	// get coords for particle
	xy1[0]=prtcls[j].x_new;
	xy1[1]=prtcls[j].y_new;
	U=U+wall_repulsion(xy1,fxy);
	// store forces associated with the particle
	prtcls[j].fx_new=fxy[0];
	prtcls[j].fy_new=fxy[1];
      }
    }
    c=0;
    while (c<n){
      // get coords for particle 1
      xy1[0]=prtcls[c].x_new;
      xy1[1]=prtcls[c].y_new;
      for (int j=0;j<n;j++){
	if (j!=c){
	  // get coords for particle 2
	  xy2[0]=prtcls[j].x_new;
	  xy2[1]=prtcls[j].y_new;
	  // store forces associated with the particle 
	  U=U+0.5*lj_potential(xy1,xy2,xy_fpair,periodic);
	  if (periodic==0){
	    prtcls[c].fx_new=prtcls[c].fx_new+xy_fpair[0];
	    prtcls[c].fy_new=prtcls[c].fy_new+xy_fpair[1];
	  }else {
	    prtcls[c].fx_new=xy_fpair[0];
	    prtcls[c].fy_new=xy_fpair[1];
	  }
	}
      }
      c++;
    }


    // Calculate new momentum, set now<-new, calculate energies
    v_mag=0;
    vx=0;
    vy=0;
    for (int j=0;j<n;j++){
      prtcls[j].vx_new=prtcls[j].vx_now+(dt/(2*m))*(prtcls[j].fx_now+prtcls[j].fx_new);
      prtcls[j].vy_new=prtcls[j].vy_now+(dt/(2*m))*(prtcls[j].fy_now+prtcls[j].fy_new);
      prtcls[j].vx_now=prtcls[j].vx_new;
      prtcls[j].vy_now=prtcls[j].vy_new;
      prtcls[j].x_now=prtcls[j].x_new;
      prtcls[j].y_now=prtcls[j].y_new;
      prtcls[j].fx_now=prtcls[j].fx_new;
      prtcls[j].fy_now=prtcls[j].fy_new;
      vx=prtcls[j].vx_now;
      vy=prtcls[j].vy_now;
      v_mag=v_mag+sqrt(vx*vx + vy*vy);
      vx=0;
      vy=0;
    }
    KE=0.5*m*(v_mag*v_mag);
    // Print energies and time
    if (i%100==0){
      printf("%.5f %.5f %.5f\n",t,KE,U);
    }
  }
}
