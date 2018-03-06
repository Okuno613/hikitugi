/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

double
  =======================================*/
#include "LIF.h"



__device__ double dv(double v ,double inp){
  return (-v+inp)/TAU;
}


__global__ void init(double *v,int *spike, double *inp,int *spikecnt)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  v[i] = V0;
  inp[i] = 0;
  spike[i] = 0;
  spikecnt[i]=0;

}

__global__ void calv(double *v,int *spike, double *inp, double *t,int *spikecnt,double TH,int NUM,int sigma)
{

  int i0=(t[0]/TEND)*NUM/2;

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  inp[i] = 20*__expf(0.01);// * __expf( -( ((i-i0)*(i-i0)) / (2*sigma*sigma)))*sin(2.0*M_PI*t[0]/100.0);

  double kv1 = DT*dv(v[i], inp[i]);
  double kv2 = DT*dv(v[i]+kv1*0.5, inp[i]);
  double kv3 = DT*dv(v[i]+kv2*0.5, inp[i]);
  double kv4 = DT*dv(v[i]+kv3,inp[i]);
  v[i] += (kv1 + 2.0*kv2 + 2.0*kv3 + kv4)/6.0;

  if(v[i]>=20){
    v[i]=V0;
  }

  if(v[i] > TH){
    v[i] =30 ;
    spike[i] = spike[i]+1;
    spikecnt[i]=spike[i];

  }
  if(int(t[0])%10==0){
    spike[i]=0;
  }
}




void Simulation::sim()
{
    int count = 0;

    int size_d = sizeof(double)*NUM;

    double *v,*d_v;
    double *inp,*d_inp;
    v = (double *)malloc(size_d);
    inp = (double *)malloc(size_d);
    cudaMalloc((void **)&d_v, size_d);
    cudaMalloc((void **)&d_inp, size_d);




    double t = 0.0;
    double *d_t;
    cudaMalloc((void **)&d_t, sizeof(double));

    int *spike,*d_spike;
    int *spikecnt,*d_spikecnt;

    FILE *fp1;
    fp1=fopen("cuda_double_lif_volt.txt","w");

    int size_i = sizeof(int)*NUM;
    spike  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spike, size_i);
    spikecnt  = (int *)malloc(size_i);
    cudaMalloc((void **)&d_spikecnt, size_i);



    init<<<NUM/Threads,Threads>>>(d_v, d_spike, d_inp,d_spikecnt);

    for(;;){

      calv<<<NUM/Threads,Threads>>>(d_v, d_spike, d_inp, d_t, d_spikecnt,TH,NUM,sigma);

      cudaMemcpy(v, d_v, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(inp, d_inp, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(spike, d_spike, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt, d_spikecnt, size_i, cudaMemcpyDeviceToHost);

      fprintf(fp1,"%lf \n",v[0]);

      count++;
      t = count * DT;
      if( t > TEND){
	break;
      }

    }
    free(v);
    free(spike);
    free(spikecnt);
    cudaFree(d_v);
    cudaFree(d_spike);
    cudaFree(spikecnt);


}



int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
