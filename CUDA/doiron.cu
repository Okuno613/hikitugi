/*=======================================
  Since : May/19/2008
  Update: <2016/02/26>

CUDA
  =======================================*/
#include "dendrite.h"

__device__ double dV_d(double V_d,double V_s,double h_d,double n_d,double p_d){
  double m_inf_d=1/(1+exp(-(V_d-V12_d)/k_m_inf_d));
  return( g_na_d* (pow(m_inf_d,2)) *h_d * (V_na - V_d) + g_dr_d * (pow(n_d,2)) *p_d *(V_k -V_d)+ (g_c/(1-kappa))*(V_s-V_d) + g_leak *(V_l - V_d));
}
      
__device__ double dh_d(double h_d, double V_d){
  double h_inf_d=1/(1+exp(-(V_d-V12_h_d)/k_h_d));
  return( (h_inf_d - h_d) /tau_h_d );
}

__device__ double dn_d(double n_d, double V_d){
  double n_inf_d=1/(1+exp(-(V_d-V12_n_d)/k_n_d));
  return( (n_inf_d - n_d) /tau_n_d );
}

__device__ double dp_d(double p_d, double V_d){
  double p_inf_d=1/(1+exp(-(V_d-V12_p_d)/k_p_d));
  return( (p_inf_d - p_d) /tau_p_d );
}
	    
__device__ double dV_s(double V_s, double inp, double V_d,double n_s){
  double m_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return (inp + g_na_s * (pow(m_inf_s,2) ) * (1-n_s )* (V_na - V_s) +g_dr_s * pow(n_s,2) *(V_k -V_s)+ (g_c/kappa)*(V_d-V_s) +g_leak *(V_l-V_s));
}

__device__ double dn_s(double n_s, double V_s){
  double n_inf_s=1/(1+exp(-(V_s-V12_s)/k_m_inf_s));
  return( (n_inf_s - n_s) /tau_n_s );
}

  

__global__ void init(double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d, int *spike_s,int *spike_d, double *inp,double *THl,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d){
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  V_s[i] = V0;
  n_s[i] = 0.5;
  V_d[i] = V0;
  h_d[i] = 0.1;
  n_d[i] = 0.1;
  p_d[i] = 0.1;
  inp[i] = 0;
  spike_s[i] = 0;
  spike_d[i] =0;
  spikecnt_s[i]=0;
  spikecnt_d[i]=0;
  count_s[i]=0;
  count_d[i]=0;
  THl[i]=TH;
}

__global__ void calv(double *V_s, double *n_s,double *V_d,double *h_d, double *n_d,double *p_d,  int *spike_s,int *spike_d, double *inp, double *t,int *spikecnt_s,int *spikecnt_d,int *count_s,int *count_d,double *THl,int sq,double sigma)
{

  
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  //  int i0=(t[0]/TEND)*sq/2;
  //int j0=sq/2;
  
  //inp[i] = I0lif * (sqrt(2*M_PI*(__powf(sigma,2)) *__expf( -( (__powf(i-i0,2) +__powf(j0,2))) / (2*(__powf(sigma,2))))));//*sin(2.0*M_PI*t[0]/100.0); //gauss+sin imp
  inp[i]=I0lif;//*sin(2.0*M_PI*t[0]/100.0); //gauss+sin imp

  double kV_s1 = DT*dV_s(V_s[i], inp[i], V_d[i],n_s[i]);
  double kn_s1 = DT*dn_s(n_s[i], V_s[i]);
  double kV_d1 = DT*dV_d(V_d[i], V_s[i] ,h_d[i], n_d[i], p_d[i]);
  double kh_d1 = DT*dh_d(h_d[i], V_d[i]);
  double kn_d1 = DT*dn_d(n_d[i], V_d[i]);
  double kp_d1 = DT*dp_d(p_d[i], V_d[i]);
  
  double kV_s2 = DT*dV_s(V_s[i]+kV_s1*0.5, inp[i], V_d[i]+kV_d1*0.5 ,n_s[i]+kn_s1*0.5);
  double kn_s2 = DT*dn_s(n_s[i]+kn_s1*0.5, V_s[i]+kV_s1*0.5);
  double kV_d2 = DT*dV_d(V_d[i]+kV_d1*0.5, V_s[i]+kV_s1*0.5 ,h_d[i]+kh_d1*0.5, n_d[i]+kn_d1*0.5, p_d[i]+kp_d1*0.5);
  double kh_d2 = DT*dh_d(h_d[i]+kh_d1*0.5, V_d[i]+kV_d1*0.5);
  double kn_d2 = DT*dn_d(n_d[i]+kn_d1*0.5, V_d[i]+kV_d1*0.5);
  double kp_d2 = DT*dp_d(p_d[i]+kp_d1*0.5, V_d[i]+kV_d1*0.5);

  double kV_s3 = DT*dV_s(V_s[i]+kV_s2*0.5, inp[i], V_d[i]+kV_d2*0.5 ,n_s[i]+kn_s2*0.5);
  double kn_s3 = DT*dn_s(n_s[i]+kn_s2*0.5, V_s[i]+kV_s2*0.5);
  double kV_d3 = DT*dV_d(V_d[i]+kV_d2*0.5, V_s[i]+kV_s2*0.5 ,h_d[i]+kh_d2*0.5, n_d[i]+kn_d2*0.5, p_d[i]+kp_d2*0.5);
  double kh_d3 = DT*dh_d(h_d[i]+kh_d2*0.5, V_d[i]+kV_d2*0.5);
  double kn_d3 = DT*dn_d(n_d[i]+kn_d2*0.5, V_d[i]+kV_d2*0.5);
  double kp_d3 = DT*dp_d(p_d[i]+kp_d2*0.5, V_d[i]+kV_d2*0.5);

  double kV_s4 = DT*dV_s(V_s[i]+kV_s3, inp[i], V_d[i]+kV_d2 ,n_s[i]+kn_s2);
  double kn_s4 = DT*dn_s(n_s[i]+kn_s3, V_s[i]+kV_s3);
  double kV_d4 = DT*dV_d(V_d[i]+kV_d3, V_s[i]+kV_s3 ,h_d[i]+kh_d3, n_d[i]+kn_d3, p_d[i]+kp_d3);
  double kh_d4 = DT*dh_d(h_d[i]+kh_d3, V_d[i]+kV_d3);
  double kn_d4 = DT*dn_d(n_d[i]+kn_d3, V_d[i]+kV_d3);
  double kp_d4 = DT*dp_d(p_d[i]+kp_d3, V_d[i]+kV_d3);

  V_s[i] += (kV_s1 + 2.0*kV_s2 + 2.0*kV_s3 + kV_s4)/6.0;
  n_s[i] += (kn_s1 + 2.0*kn_s2 + 2.0*kn_s3 + kn_s4)/6.0;
  V_d[i] += (kV_d1 + 2.0*kV_d2 + 2.0*kV_d3 + kV_d4)/6.0;
  h_d[i] += (kh_d1 + 2.0*kh_d2 + 2.0*kh_d3 + kh_d4)/6.0;
  n_d[i] += (kn_d1 + 2.0*kn_d2 + 2.0*kn_d3 + kn_d4)/6.0;
  p_d[i] += (kp_d1 + 2.0*kp_d2 + 2.0*kp_d3 + kp_d4)/6.0;

  
  if(V_s[i] > 20 and count_s[i] == 0){
    spike_s[i] = spike_s[i]+1;
    spikecnt_s[i]=spike_s[i];
    //THl[i]=THl[i]+THup;
    count_s[i]=1;
    //fprintf(fp1,"%d\t %d\t %d\n \n",i,int(t),spikecnt_s[i] );
  
  }
  
  if(V_d[i] > 20 and count_d[i] == 0){
    spike_d[i] = spike_d[i]+1;
    spikecnt_d[i]=spike_d[i];
    count_d[i]=1;
    //THl[i]=THl[i]+THup;
    //fprintf(fp2,"%d\t %d\t %d\n \n",i,int(t),spikecnt_d[i] );
  
  }

  if(int(t[0])%10==0){
    spike_s[i]=0;
    spike_d[i]=0;
  }
    
  if(count_d[i]==1 and V_d[i]<=-55){
    count_d[i]=0;
  }
  if(count_s[i]==1 and V_s[i]<=-55){
    count_s[i]=0;
  }


  //fprintf(fp3,"%lf \t %lf \n",t,V_s[0]);  
  //fprintf(fp4,"%lf \t %lf \n",t,V_d[0]);


}



void Simulation::sim()
{
    int count = 0;
    int size_d = sizeof(double)*NUM;
    int size_i = sizeof(int)*NUM;
    double *V_s, *n_s,*V_d,*h_d,*n_d,*p_d;
    double *d_V_s, *d_n_s,*d_V_d,*d_h_d,*d_n_d,*d_p_d;
    double *inp,*d_inp;

    V_s = (double *)malloc(size_d);
    n_s = (double *)malloc(size_d);
    V_d = (double *)malloc(size_d);
    h_d = (double *)malloc(size_d);
    n_d = (double *)malloc(size_d);
    p_d = (double *)malloc(size_d);    
    inp = (double *)malloc(size_d);

    cudaMalloc((void **)&d_V_s, size_d);
    cudaMalloc((void **)&d_n_s, size_d);
    cudaMalloc((void **)&d_V_d, size_d);
    cudaMalloc((void **)&d_h_d, size_d);
    cudaMalloc((void **)&d_n_d, size_d);
    cudaMalloc((void **)&d_p_d, size_d);
    cudaMalloc((void **)&d_inp, size_d);

    int *count_s,*d_count_s;
    int *count_d,*d_count_d;
    count_s = (int *)malloc(size_i);
    count_d = (int *)malloc(size_i);
    cudaMalloc((void **)&d_count_s, size_i);
    cudaMalloc((void **)&d_count_d, size_i);
    
    double t = 0.0;
    double *d_t;
    cudaMalloc((void **)&d_t, sizeof(double));
    
    int *spike_s,*d_spike_s;
    int *spikecnt_s,*d_spikecnt_s;
    int *spike_d,*d_spike_d;
    int *spikecnt_d,*d_spikecnt_d;
    double *THl,*d_THl;

    spike_s = (int *)malloc(size_i);
    spikecnt_s = (int *)malloc(size_i);
    spike_d = (int *)malloc(size_i);
    spikecnt_d = (int *)malloc(size_i);
    THl = (double *)malloc(size_d);

    cudaMalloc((void **)&d_spike_s, size_i);
    cudaMalloc((void **)&d_spikecnt_s, size_i);
    cudaMalloc((void **)&d_spike_d, size_i);
    cudaMalloc((void **)&d_spikecnt_d, size_i);
    cudaMalloc((void **)&d_THl, size_d);


    
    FILE *fp1,*fp2,*fp3,*fp4;    
    fp1=fopen("Vs_moved.txt","w");
    fp2=fopen("Vd_moved.txt","w");
    fp3=fopen("cuda_double_Vs0_volt.txt","w");
    fp4=fopen("cuda_double_Vd0_volt.txt","w");

    init<<<NUM/ Threads, Threads>>>(d_V_s,d_n_s, d_V_d, d_h_d, d_n_d, d_p_d, d_spike_s,d_spike_d, d_inp, d_THl,d_spikecnt_s,d_spikecnt_d,d_count_s,d_count_d);


    //fprintf(fp2,"%lf \t %lf \n",t,V_s[0]);     

    for(;;){
      //cudaMemcpy(d_t, &t, sizeof(double), cudaMemcpyHostToDevice);
      
      
 
      calv<<<NUM/ Threads,Threads>>>(d_V_s, d_n_s, d_V_d, d_h_d, d_n_d, d_p_d, d_spike_s,d_spike_d,d_inp, d_t, d_spikecnt_s,d_spikecnt_d,d_count_s,d_count_d,d_THl ,sq,sigma);


      cudaMemcpy(V_s,d_V_s, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(n_s,d_n_s, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(V_d,d_V_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(h_d,d_h_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(n_d,d_n_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(p_d,d_p_d, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(inp,d_inp, size_d, cudaMemcpyDeviceToHost);
      cudaMemcpy(count_s,d_count_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spike_s,d_spike_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt_s,d_spikecnt_s, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(count_d,d_count_d, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spike_d,d_spike_d, size_i, cudaMemcpyDeviceToHost);
      cudaMemcpy(spikecnt_d,d_spikecnt_d, size_i, cudaMemcpyDeviceToHost);

      for(int i=0;i<NUM;i++){
	if(V_s[0] > 20 and spikecnt_s[i] > 0){
	  // fprintf(fp1,"%d\t %d\t %d\n \n",i,int(t),spikecnt_s[i] );
	}
      }
      
      //  fprintf(fp3,"%lf \t %lf \n",t,V_s[0]); 
      //fprintf(fp4,"%lf \t %lf \n",t,V_d[0]);

      count++;
      t = count * DT;
      if( t > TEND){
	break;
      }
      /*
      cudaFree(d_V_s);
      cudaFree(d_n_s);
      cudaFree(d_V_d);
      cudaFree(d_h_d);
      cudaFree(d_n_d);
      cudaFree(d_p_d);

      cudaFree(d_count_s);
      cudaFree(d_count_d);
    
      cudaFree(spike_s);
      cudaFree(spike_d);
      cudaFree(spikecnt_s);
      cudaFree(spikecnt_d);
      */
    }
    
  
    cudaMemcpy(V_s,d_V_s, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(n_s,d_n_s, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(V_d,d_V_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_d,d_h_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(n_d,d_n_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(p_d,d_p_d, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(inp,d_inp, size_d, cudaMemcpyDeviceToHost);
    cudaMemcpy(count_s,d_count_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spike_s,d_spike_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spikecnt_s,d_spikecnt_s, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(count_d,d_count_d, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spike_d,d_spike_d, size_i, cudaMemcpyDeviceToHost);
    cudaMemcpy(spikecnt_d,d_spikecnt_d, size_i, cudaMemcpyDeviceToHost);
    
    
    free(V_s);
    free(n_s);
    free(V_d);
    free(h_d);
    free(n_d);
    free(p_d);
    cudaFree(d_V_s);
    cudaFree(d_n_s);
    cudaFree(d_V_d);
    cudaFree(d_h_d);
    cudaFree(d_n_d);
    cudaFree(d_p_d);
    free(count_s);
    free(count_d);
    cudaFree(d_count_s);
    cudaFree(d_count_d);
    free(spike_s);
    free(spike_d);
    free(spikecnt_s);
    free(spikecnt_d);
    cudaFree(THl);
    cudaFree(spike_s);
    cudaFree(spike_d);
    cudaFree(spikecnt_s);
    cudaFree(spikecnt_d);
    cudaFree(THl);

    
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);    
    

}

 
/* 
int main(int argc, char* argv[]){
 Simulation sim;
 sim.sim();
 return(0);
 }
*/

int main(int argc, char* argv[]){

 Simulation sim;
 sim.sim();
     
 return(0);
 }
