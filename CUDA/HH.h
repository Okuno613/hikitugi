#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>

const double DT = 0.01;
const int TEND = 1000;
const double t= 0.0;

const int NUM = pow(2,0);
const int Threads = 128;

const double V0 =-55.3;
const double m0 = (-0.1*(V0+40)/(exp((-V0-40)/10)-1)) /( -0.1*(V0+40)/(exp((-V0-40)/10)-1) +4*exp((-V0-65)/18) );
const double h0 =(0.07*exp((-V0-65)/20))/( 0.07*exp((-V0-65)/20)+ (1/(1+exp((-V0-35)/10)) )  );
const double n0 =(-0.01*(V0+55) / (exp((-V0-55)/10)-1)) / ( (-0.01*(V0+55) / (exp((-V0-55)/10)-1) + ( 0.125*exp((-V0-65)/80) ) )) ;

const double TH=0;
const double g_na= 120;
const double g_k= 36;
const double g_l= 0.3;

const double E_na= 55.0;
const double E_k=-88.5;
const double E_l=-54.387;
const double Cm=1;


const double I0lif=20;

const int FT=500;
const double sq = sqrt(NUM);
const double sigma = NUM/4;

const int size_d = sizeof(double)*NUM;

class Simulation{
public:
    void sim();
};
#endif
