#ifndef COEF_P_HPP
#define COEF_P_HPP

double Ic(const double s, const double B) {
	double minLeft = ((G > B)? B: G);
	double cf = ((1 - withdrawalPen) * B + withdrawalPen * minLeft);
	return ((cf > s)? cf: s);
	}

double Bfun(const double t, const double T) {return (1.0 / rkappa * (1.0 - exp(-rkappa * (T - t))));}
double Afun(const double t, const double T) { return ((rGamma - sig2 * sig2 / (2.0 * rkappa * rkappa)) * (Bfun(t, T) + t - T) - sig2 * sig2 * Bfun(t, T) * Bfun(t, T) / (4.0 * rkappa));}
double pfun(const double t, const double T, const double r) {return (exp(Afun(t, T) - r * Bfun(t, T)));}
double mur(const double r, const double deln) {return (r * exp(-rkappa * deln) + (rGamma - sig2 * sig2 / (rkappa * rkappa)) * (1.0 - exp(-rkappa * deln)) + sig2 * sig2 / (2.0 * rkappa * rkappa) * (1.0 - exp(-2.0 * rkappa * deln)));}
double sigr(const double deln) {return (sig2 * sig2 * (1.0 - exp(-2.0 * rkappa * deln)) / (2 * rkappa));}
double mux(const double x, const double r, const double deln) {return (x + (1.0 - exp(-rkappa * deln)) * (r + (1.0 - exp(-rkappa * deln)) * sig2 * sig2 / (2.0 * rkappa * rkappa)) / rkappa + (rGamma - sig2 * sig2 / (rkappa * rkappa)) * (deln - (1.0 - exp(-rkappa * deln)) / rkappa) - rho13 * sig1 * sig2 * (rkappa * deln - (1.0 - exp(-rkappa * deln))) / (rkappa * rkappa) - (fairfee + .5 * sig1 * sig1) * deln);}
double sigx(const double deln) {return (sig1 * sig1 * deln + sig2 * sig2 * (2.0 * rkappa * deln - 4.0 * (1.0 - exp(-rkappa * deln)) + (1.0 - exp(-2.0 * rkappa * deln))) / (2.0 * rkappa * rkappa * rkappa) + 2.0 * rho13 * sig1 * sig2 * (rkappa * deln - (1.0 - exp(-rkappa * deln))) / (rkappa * rkappa));}
double rhoxr(const double deln) {return ((rho13 * sig1 * sig2 * (1.0 - exp(-rkappa * deln)) / rkappa + sig2 * sig2 * (2.0 * (1.0 - exp(-rkappa * deln)) - (1.0 - exp(-2.0 * rkappa * deln))) / (2.0 * rkappa * rkappa)) / pow(sigx(deln) * sigr(deln), 0.5));}
//note that x is xm + ln(W(0))
double psitow(const double psi1, const double psi2, const double x, const double r, const double deln) {return (exp(pow(2.0, 0.5) * pow(sigx(deln), 0.5) * (atemp * psi1 + btemp * psi2) + mux(x, r, deln)));}
double psitor(const double psi1, const double psi2, const double x, const double r, const double deln) {return (pow(2.0, 0.5) * pow(sigr(deln), 0.5) * (btemp * psi1 + atemp * psi2) + mur(r, deln));}

#endif