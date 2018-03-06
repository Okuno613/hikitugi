/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

float
  =======================================*/
#include "LIF.h"
#include <omp.h>



double dv(double v ,double inp){
  return (-v+inp)/TAU;
}

double runge(double *v, double *inp, int i){
  double kv1 = DT*dv(v[i], inp[i]);
  double kv2 = DT*dv(v[i]+kv1*0.5, inp[i]);
  double kv3 = DT*dv(v[i]+kv2*0.5, inp[i]);
  double kv4 = DT*dv(v[i]+kv3,inp[i]);
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;
}


void init(double *v,int *spike, double *inp,int *spikecnt)
{
#pragma omp parallel for
  for(int i=0;i<NUM;i++){
    v[i] = V0;
    inp[i] = 0;
    spike[i] = 0;
    spikecnt[i]=0;
  }
}

void calv(double *v,int *spike, double *inp, double t,int *spikecnt,double TH,FILE *fp1,FILE *fp2)
{
  int i0=(t/TEND)*NUM/2;

#pragma omp parallel for
  for(int i=0;i<NUM;i++){
    inp[i] = 20;// * exp( -( ((i-i0)*(i-i0)) / (2*pow(sigma,2))));//*sin(2.0*M_PI*t/100.0);

    runge(v,inp,i);
    
  if(v[i]>=20){
    v[i]=V0;
  }
  
  if(v[i] > TH){
    v[i] =30 ;
    spike[i] = spike[i]+1;
    spikecnt[i]=spike[i];
    fprintf(fp1,"%d\t %d\t %d\n",i,int(t),spikecnt[i] );
  }
  
    if(int(t)%10==0){
      spike[i]=0;
      }
  }
  // fprintf(fp2,"%lf\t \n",v[0] );
  fprintf(fp2,"%lf\t %lf \n",t,v[0] );
}



void Simulation::sim()
{
    int count = 0;


    double *v =new double[NUM];
    double *u =new double[NUM];
    double *inp =new double[NUM];
    //    int sq = sqrt(NUM);




    double t = 0.0;
    int *spike= new int[NUM];
    int *spikecnt = new int[NUM];

    FILE *fp1,*fp2;
    fp1=fopen("double_LIF_cpp_moved.txt","w");
    fp2=fopen("double_LIF_cpp_V.txt","w");
    
    init(v, spike, inp,spikecnt);



    for(;;){

      calv(v, spike, inp, t, spikecnt,TH,fp1,fp2);
      count++;
      t = count * DT;
      if( t > TEND){
	break;
      }

   
}
    fclose(fp1);
    fclose(fp2);
}

 
 
int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
