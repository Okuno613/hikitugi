#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>

const double DT = 0.01;
const int TEND = 1000;

const int NUM = pow(2,8);
const int Threads = 128;

const double TAU = 15.0;
// time constant of neuron (ms)

const double TH0 = 1.0;
const double a = 0.02;
const double b = 0.2;
const double c = -65.0;
const double d = 8.0;


const double INP = 5.0;
const double V0 = 0.0;// (mV) equal c
const double TH = 15.0;// (mV)
const double STD_N = 2.0;//standard deviation of threshold (mV)

const double I0=30;

const int FT=500;
const int sq = sqrt(NUM);
const double sigma = NUM/10;



class Simulation{
public:
    void sim();
};
#endif
