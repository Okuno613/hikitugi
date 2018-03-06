/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

float
  =======================================*/
#include "izhikevich.h"
#include <omp.h>



double dv(double v,double u, double inp){
  return 0.04*v*v +5*v +140 -u +inp;
}

double du(double v,double u){
  return a*(b*v-u);
}


double runge(double *v, double *u, double *inp, int i){
  double kv1 = DT*dv(v[i],u[i], inp[i]);
  double ku1 = DT*du(v[i],u[i]);
    
  double kv2 = DT*dv(v[i]+kv1*0.5,u[i]+ku1*0.5, inp[i]);
  double ku2 = DT*du(v[i]+kv1*0.5,u[i]+ku1*0.5);
    

  double kv3 = DT*dv(v[i]+kv2*0.5,u[i]+ku2*0.5, inp[i]);
  double ku3 = DT*du(v[i]+kv2*0.5,u[i]+ku2*0.5);

  double kv4 = DT*dv(v[i]+kv3, u[i]+ku3, inp[i]);
  double ku4 = DT*du(v[i]+kv3, u[i]+ku3);
    
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;
  u[i] += (ku1 + 2.0*ku2 + 2.0*ku3 + ku4)/6.0;
}


void init(double *v, double *u,int *spike, double *inp,int *spikecnt)
{
#pragma omp parallel for
  for(int i=0;i<NUM;i++){
    v[i] = V0;
    u[i] = U0;
    inp[i] = 0;
    spike[i] = 0;
    spikecnt[i]=0;
  }
}

void calv(double *v, double *u, int *spike, double *inp, double t,int *spikecnt,double TH,FILE *fp1,FILE *fp2)
{
  int i0=(t/TEND)*NUM/2;

#pragma omp parallel for
  for(int i=0;i<NUM;i++){
    inp[i] = I0*exp(0.01);// * exp( -( ((i-i0)*(i-i0)) / (2*pow(sigma,2))))*sin(2.0*M_PI*t/100.0);

    runge(v,u,inp,i);

    if(v[i] > TH){
      v[i] =c ;
      u[i] =u[i]+d;
      spike[i] = spike[i]+1;
      spikecnt[i]=spike[i];
      
    }
    if(int(t)%10==0 && i==0){
      fprintf(fp1,"%lf\t %d\n",t,spikecnt[0] );
      spike[0]=0;
      }
  }
  fprintf(fp2,"%lf\t \n",v[0] );
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
    fp1=fopen("double_izh_cpp_psth.txt","w");
    fp2=fopen("double_izh_cpp_v.txt","w");
    
    init(v,u, spike, inp,spikecnt);



    for(;;){

      calv(v,u, spike, inp, t, spikecnt,TH,fp1,fp2);
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
