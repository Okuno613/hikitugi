
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
 
#define pi 3.1415926535         // 円周率

 // 離散フーリエ変換（読み込むファイル名, 書き込むファイル名）
int dft(char filename1[], char filename2[], char filename3[])
{  
        int k, n, N;
        int max = 100000;  // 読み込むデータ数の上限
        double f1[max+1];
	double f2[max+1];
	double ans;
	
        FILE *fp1, *fp2, *fp3;
    // ファイルオープン(フーリエ変換したいデータファイル, フーリエ変換後のデータ保存用ファイル)
        if((fp1=fopen(filename1,"r"))==NULL){
	  printf("FILE1 not open\n");
        return -1;
        }
        if((fp2=fopen(filename2,"r"))==NULL){
	  printf("FILE2 not open\n");
        return -1;
        }
	if((fp3=fopen(filename3,"w"))==NULL){
	  printf("FILE3 not open\n");
        return -1;
        }
        //データの読み込み
        for(N=0; N<max; N++) {
	  fscanf(fp1,"%lf", &f1[N]);
	  fscanf(fp2,"%lf", &f2[N]);
        }
        //実数部分と虚数部分に分けてフーリエ変換
        for(k=0; k<N; k++){
  ans=f1[k]-f2[k];
  if(ans>=0){
  fprintf(fp3,"%lf\t %lf \n",k*0.01, ans  );
  }
  if(ans<0){
  fprintf(fp3,"%lf\t %lf \n",k*0.01, -ans  );
  }
}
  
        fclose(fp1);
        fclose(fp2);
	fclose(fp3);
	
        return 0;
}

int main()
{

  //dft("Vs_volt.txt","float_Vs_volt.txt","hikaku.txt");
  dft("double_LIF_cpp_V.txt","float_LIF_cpp_V.txt","hikaku_LIF.txt");
  //dft("double_HH_cpp_v.txt","float_HH_cpp_v.txt","hikaku_HH.txt");
}

