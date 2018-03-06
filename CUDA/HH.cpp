/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

float
  =======================================*/
#include "HH.h"
#include <omp.h>

double dv(double v,double m,double h,double n,double inp){
  return( (inp+ g_na*(pow(m,3))*h*(E_na-v)+ g_k*pow(n,4)*(E_k-v) +g_l*(E_l-v)  )/Cm );
}

double dm(double m,double v){
  double am= -0.1*(v+40)/(exp((-v-40)/10)-1);
  double bm= 4*exp((-v-65)/18);
  double taum = 1/(am+bm);
  return( ( am*(1-m) - bm*m ) );
}

double dh(double h,double v){
  double ah = 0.07*exp((-v-65)/20);
  double bh = 1/(1+exp((-v-35)/10));
  double tauh = 1/(ah+bh);
  return( (ah*(1-h) - bh *h) );
}

double dn(double n,double v){
  double an = -0.01*(v+55) / (exp((-v-55)/10)-1) ;
  double bn = 0.125*exp((-v-65)/80);
  double taun = 1/(an+bn);
  return( (an*(1-n) - bn*n) );
}


double runge(double *v, double *inp, double *m,double *h, double *n,int i){

  double kv1 = DT*dv(v[i], m[i], h[i], n[i] ,inp[i]);
  double km1 = DT*dm(m[i], v[i]);
  double kh1 = DT*dh(h[i], v[i]);
  double kn1 = DT*dn(n[i], v[i]);
    
  double kv2 = DT*dv(v[i]+kv1*0.5, m[i]+km1*0.5, h[i]+kh1*0.5, n[i]+kn1*0.5 ,inp[i]);
  double km2 = DT*dm(m[i]+km1*0.5, v[i]+kv1*0.5);
  double kh2 = DT*dh(h[i]+kh1*0.5, v[i]+kv1*0.5);
  double kn2 = DT*dn(n[i]+kn1*0.5, v[i]+kv1*0.5);
  
  double kv3 = DT*dv(v[i]+kv2*0.5, m[i]+km2*0.5, h[i]+kh2*0.5, n[i]+kn2*0.5 ,inp[i]);
  double km3 = DT*dm(m[i]+km2*0.5, v[i]+kv2*0.5);
  double kh3 = DT*dh(h[i]+kh2*0.5, v[i]+kv2*0.5);
  double kn3 = DT*dn(n[i]+kn2*0.5, v[i]+kv2*0.5);
  
  double kv4 = DT*dv(v[i]+kv3, m[i]+km3, h[i]+kh3, n[i]+kn3 ,inp[i]);
  double km4 = DT*dm(m[i]+km3, v[i]+kv3);
  double kh4 = DT*dh(h[i]+kh3, v[i]+kv3);
  double kn4 = DT*dn(n[i]+kn3, v[i]+kv3);

    
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;
  m[i] += (km1 + 2.0*km2 + 2.0*km3 + km4)/6.0;
  h[i] += (kh1 + 2.0*kh2 + 2.0*kh3 + kh4)/6.0;
  n[i] += (kn1 + 2.0*kn2 + 2.0*kn3 + kn4)/6.0;
  
}

  

void init(double *v, double *m,double *h, double *n, int *spike, double *inp,double *THl,int *spikecnt,int *cnt)
{
  #pragma omp parallel for
  for(int i=0;i<NUM;i++){
    v[i] = V0;
    m[i] = m0;
    h[i] = h0;
    n[i] = n0;
    inp[i] = I0lif;
    spike[i] = 0;
    spikecnt[i]=0;
    cnt[i]=0;
    THl[i]=0;
  }
}

void calv(double *v, double *m,double *h, double *n,int *spike, double *inp, double t,int *cnt,double *THl,int *spikecnt,FILE *fp1,FILE *fp2)
{
  int i0=(t/TEND)*NUM/2;
  //  int j0=sq/2;
  fprintf(fp2,"%lf \t %lf\n",t,v[0] );
   //fprintf(fp1,"%lf\t %lf\n",t,n[0]);
  #pragma omp parallel for
  for(int i=0;i<NUM;i++){
    inp[i] = I0lif;//*exp(0.01);//*sin((2*M_PI*t*50/1000)) ;

    runge(v,inp, m,h,n,i);

    
    if(v[i] > 0 and cnt[i] == 0){
      //spike[i] = spike[i]+1;
      cnt[i]=1;
      fprintf(fp1,"%d\t %d\t %03d\n",i,int(t),(spike[i]-spikecnt[i]) );
      spikecnt[i]=int(t);
    }
    spike[i]=int(t);

      /*
    if(int(t)%10==0 && i==0){
      fprintf(fp1,"%lf\t %d\n",t,spikecnt[0] );
      spike[0]=0;
      spikecnt[0]=0;
      }
      */
    
    
    if(cnt[i]==1 and v[i]<=-55){
      cnt[i]=0;
    }
  }
  //fprintf(fp2,"%lf\t %lf\n",t,v[0] );
}


void Simulation::sim()
{
    int count = 0;


    double *v =new double[NUM];
    double *m =new double[NUM];
    double *h =new double[NUM];
    double *n =new double[NUM];
    
    double *inp =new double[NUM];
    
    int *cnt=new int[NUM];

    double t = 0.0;
    int *spike= new int[NUM];
    int *spikecnt = new int[NUM];

    FILE *fp1,*fp2;
    fp1=fopen("double_HH_cpp_h.txt","w");
    fp2=fopen("double_HH_cpp_v.txt","w");
    
    double *THl =new double[NUM];

    init(v,m,h, n, spike, inp, THl,spikecnt,cnt);
    for(;;){

      calv(v, m, h, n,spike,inp, t, cnt,THl,spikecnt,fp1,fp2);
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
