#include "stdc++.h"
// #include "spline.hpp"		//from http://kluge.in-chemnitz.de/opensource/spline/
#include "alglib/interpolation.cpp"
// #include "hermite.cpp"
#include "characteristics_pavel.hpp"
#include "coefunc_pavel.hpp"
#include "range.cpp"
#include "ghqc.cpp"
#include "findfair_pavel.cpp"

int main() {
	// wS = new double[orderS];
	// xS = new double[orderS];
	// wR = new double[orderR];
	// xR = new double[orderR];
	// cgqf(orderS, kind, alpha, Beta, a, b, xS, wS);
	// cgqf(orderR, kind, alpha, Beta, a, b, xR, wR);

	alglib::gqgenerategausshermite(orderS, info, xS, wS);
	alglib::gqgenerategausshermite(orderR, info, xR, wR);

	Trange xrange(Xfrom, Xto);
	Trange rrange(Rfrom, Rto);
	Trange trange(Tfrom, Tto);
	Trange brange(Bfrom, Bto);

	double (*PsitoW) (const double, const double, const double, const double, const double) = &psitow;
	double (*PsitoR) (const double, const double, const double, const double, const double) = &psitor;
	double (*Pfun) (const double, const double, const double) = &pfun;

	double (*ic) (const double, const double) = &Ic;

	GHQC ghqc(PsitoW, PsitoR, Pfun, ic, xrange, rrange, trange, brange, Snum, Rnum, Tnum, Bnum, Bnum_dy, xS, wS, xR, wR, orderS, orderR, Dynamic);
	
	long maxiter = 10;
	double tol = 1e-2;

	double starget = P0;
	double rtarget = rGamma;

	FindFairP fairVal(ghqc, starget, rtarget, Bnum, maxiter, tol);
	fairVal.start();

	std::cout << fairfee;

	// ghqc.start();

	// std::vector<double> sarr = ghqc.getSARR();
	// std::vector<double> rarr = ghqc.getRARR();
	// std::vector<double> barr = ghqc.getBARR();
	// std::vector<double> final = ghqc.getres();

	// std::ofstream output_file("./test.txt");

	// output_file << std::setw(20) << "b: " << std::setw(20) << "r: " << std::setw(20) << "s: " << std::setw(20) << "value: " << "\n";

	// for (long bb = 0; bb< (Bnum + 1); bb++) {
	// 	for (long i = 0; i < (Rnum + 1); i++) {
	// 		for (long k = 0; k < (Snum + 1); k++) {
	// 			output_file << std::setw(20) << barr[bb] << std::setw(20) << rarr[i] << std::setw(20) << sarr[k] << std::setw(20) << final[bb * (Rnum + 1) * (Snum + 1) + i * (Snum + 1) + k] << "\n";
	// 		}
	// 	}
	// }

	return 0;
}