#ifndef CHAR_P_HPP
#define CHAR_P_HPP

const double Gamma = 0.0;

// const double rkappa = 1.0;
// const double rGamma = .05;
// const double sig1 = .2;			//s - vol
// const double sig2 = .2;			//r - vol

const double rkappa = .0349;
const double rGamma = .05;
const double sig1 = .2;			//s - vol
const double sig2 = .02;			//r - vol

const double rho13 = 0.3;		//relationship sr
// const double rho13 = 0.3;		//relationship sr


const double atemp = .5 * (pow(1 + rho13, .5) + pow(1 - rho13, .5));
const double btemp = .5 * (pow(1 + rho13, .5) - pow(1 - rho13, .5));

double fairfee;					//fair fee charged for insurance company

const double Rto = .75;			//interest rate domain
//changed from -1.0 to -0.1
const double Rfrom = - Rto;		//interest rate domain
const double Tfrom = 0.0;		//time domain
const double Tto = 10.0;			//time domain

const double withdrawalPen = .1;//withdrawal penalty
const double P0 = 1.0;		//initial premium
const double WF = 1.0;			//withdrawal frequence
const double G = P0 / (Tto * WF);	//guaranteed minimum withdrawal

const unsigned long Snum = 50;	//number of stock price, 50, 100
const unsigned long Rnum = 30;	//number of interest rate, 30, 60
const unsigned long Tnum = 1 * WF * Tto;	//number of time, 1
const unsigned long Bnum = WF * Tto;	//number of guaranteed account
const unsigned long Bnum_dy = WF * Tto;
//const unsigned long Bnum_dy = 100;

const double Sfrom = P0 / 100.0;		//stock price domain
//changed from 100 * Tto * P0 to 14 * P0
const double Sto = 100.0 * P0;	//stock price domain
const double Xfrom = log(Sfrom / P0);
const double Xto = log(Sto / P0);
const double Bfrom = 0.0;		//guaranteed account domain
const double Bto = P0;			//guaranteed account domain

// double a = 0;
// double alpha = 0;
// double b = 1.0;
// double Beta = 0.0;
// int kind = 6;
// const unsigned long orderS = 5;	//5, 9
// const unsigned long orderR = 3;	//3, 5
// double *wS;
// double *xS;
// double *wR;
// double *xR;

alglib::ae_int_t orderS = 5;	//5, 9
alglib::ae_int_t orderR = 3;	//3, 5
alglib::ae_int_t info = 1;
alglib::real_1d_array wS, xS, wR, xR;

#endif