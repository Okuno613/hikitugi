/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

double
  =======================================*/
#include "HH.h"

__device__ double dV(double v,double m,double h,double n,double inp){
  return( ((g_na*m*m*m*h*(E_na-v)+g_k*n*n*n*n*(E_k-v)+g_l*(E_l-v) +inp)) /Cm );
}

__device__ double dm(double m,double v){
  double am= (0.1*(v+40))/(1-exp((-v-40)/10));
  double bm= 4*exp((-v-65)/18);
  double taum = 1/(am+bm);
  return( ( am*(1-m) - bm*m ) );
}

__device__ double dh(double h,double v){
  double ah =( (0.07*exp((-v-65)/20)));
  double bh = (1/(1+exp((-v-35)/10)));
  double tauh = 1/(ah+bh);
  return( (ah*(1-h) - bh *h) );
}

__device__ double dn(double n,double v){
  double an =(0.01*(v+55) / (1-exp((-v-55)/10) ) );
  double bn = ( 0.125*exp((-v-65)/80) );
  double taun = 1/(an+bn);
  return( (an*(1-n) - bn*n) );
}

__global__ void init(double *V, double *m,double *h, double *n, int *spike, double *inp,int *spikecnt,int *cnt)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
    V[i] = V0;
    m[i] = 0;
    h[i] = 0;
    n[i] = 0;
    inp[i] = I0lif;
    spike[i] = 0;
    spikecnt[i]=0;
    cnt[i]=0;
    //TH[i]=TH;

}

__global__ void calv(double *V, double *m,double *h, double *n,int *spike, double *inp, double *t,int *cnt,int NUM,double sigma)
{
  //int i0=(t[0]/TEND)*NUM/2;
  //  int j0=sq/2;
  int i = threadIdx.x + blockIdx.x * blockDim.x;

  inp[i] = I0lif*__expf(0.01);// * /* (sqrt(2*M_PI*pow(sigma,2))) * */ __expf( -( ((i-i0)*(i-i0)) / (2*sigma*sigma)))*sin(2.0*M_PI*t[0]/100.0);

  double kV1 = DT*dV(V[i], m[i], h[i], n[i] ,inp[i]);
  double km1 = DT*dm(m[i], V[i]);
  double kh1 = DT*dh(h[i], V[i]);
  double kn1 = DT*dn(n[i], V[i]);

  double kV2 = DT*dV(V[i]+kV1*0.5, m[i]+km1*0.5, h[i]+kh1*0.5, n[i]+kn1*0.5 ,inp[i]);
  double km2 = DT*dm(m[i]+km1*0.5, V[i]+kV1*0.5);
  double kh2 = DT*dh(h[i]+kh1*0.5, V[i]+kV1*0.5);
  double kn2 = DT*dn(n[i]+kn1*0.5, V[i]+kV1*0.5);

  double kV3 = DT*dV(V[i]+kV2*0.5, m[i]+km2*0.5, h[i]+kh2*0.5, n[i]+kn2*0.5 ,inp[i]);
  double km3 = DT*dm(m[i]+km2*0.5, V[i]+kV2*0.5);
  double kh3 = DT*dh(h[i]+kh2*0.5, V[i]+kV2*0.5);
  double kn3 = DT*dn(n[i]+kn2*0.5, V[i]+kV2*0.5);

  double kV4 = DT*dV(V[i]+kV3, m[i]+km3, h[i]+kh3, n[i]+kn3 ,inp[i]);
  double km4 = DT*dm(m[i]+km3, V[i]+kV3);
  double kh4 = DT*dh(h[i]+kh3, V[i]+kV3);
  double kn4 = DT*dn(n[i]+kn3, V[i]+kV3);


  V[i] += (kV1 + 2.0*kV2 + 2.0*kV3 + kV4)/6.0;
  m[i] += (km1 + 2.0*km2 + 2.0*km3 + km4)/6.0;
  h[i] += (kh1 + 2.0*kh2 + 2.0*kh3 + kh4)/6.0;
  n[i] += (kn1 + 2.0*kn2 + 2.0*kn3 + kn4)/6.0;


}


void Simulation::sim()
{
    int count = 0;
    int size_d = sizeof(double)*NUM;

    double *V,*d_V;
    double *m,*d_m;
    double *h,*d_h;
    double *n,*d_n;
    double *inp,*d_inp;
    V = (double *)malloc(size_d);
    m = (double *)malloc(size_d);
    h = (double *)malloc(size_d);
    n = (double *)malloc(size_d);
    inp = (double *)malloc(size_d);
    cudaMalloc((void **)&d_V, size_d);
    cudaMalloc((void **)&d_m, size_d);
    cudaMalloc((void **)&d_h, size_d);
    cudaMalloc((void **)&d_n, size_d);
    cudaMalloc((void **)&d_inp, size_d);


    double t = 0.0;
    double *d_t;
    cudaMalloc((void **)&d_t, sizeof(double));

    int *spike,*d_spike;
    int *spikecnt,*d_spikecnt;
    int *cnt,*d_cnt;

    int size_i = sizeof(int)*NUM;
    spike  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spike, size_i);
    spikecnt  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spikecnt, size_i);
    cnt  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_cnt, size_i);

    FILE *fp1;
    fp1=fopen("cuda_double_HH_volt.txt","w");


    init<<<NUM/Threads,Threads>>>(d_V,d_m,d_h, d_n, d_spike, d_inp,d_spikecnt,d_cnt);


    for(;;){

      calv<<<NUM/Threads,Threads>>>(d_V, d_m, d_h, d_n,d_spike,d_inp, d_t, d_cnt,NUM,sigma);

      cudaMemcpy(V, d_V, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(m, d_m, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(h, d_h, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(n, d_n, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(inp, d_inp, size_d, cudaMemcpyDeviceToHost);

      cudaMemcpy(spike, d_spike, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt, d_spikecnt, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(cnt, d_cnt, size_i, cudaMemcpyDeviceToHost);

      fprintf(fp1,"%lf \n",V[0]);


      count++;
      t = count * DT;
      if( t > TEND){
	break;
      }


    }
    free(V);
    free(n);
    free(h);
    free(m);
    cudaFree(d_V);
    cudaFree(d_n);
    cudaFree(d_h);
    cudaFree(d_m);
    free(cnt);
    cudaFree(d_cnt);
    free(spike);
    free(spikecnt);
    cudaFree(d_spike);
    cudaFree(d_spikecnt);

}



int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
