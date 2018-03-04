#ifndef GHQC_CPP
#define GHQC_CPP

#include "ghqc.hpp"

GHQC::GHQC() {}

GHQC::GHQC(double (*psiW) (const double, const double, const double, const double, const double), 
	double (*psiR) (const double, const double, const double, const double, const double),
	double (*pFun) (const double, const double, const double),
	double (*IC) (const double, const double),
	Trange xxrange, Trange rrrange, Trange ttrange, Trange bbrange,
	unsigned long ssnum, unsigned long rrnum, unsigned long ttnum, unsigned long bbnum, unsigned long bbnum_dy,
	double xS[], double wS[], double xR[], double wR[], unsigned long samples, unsigned long sampler, polType ptype) {
	
	psiTow = psiW;
	psiTor = psiR;
	pfun = pFun;
	ic = IC;

	xrange = xxrange;
	rrange = rrrange;
	trange = ttrange;
	brange = bbrange;

	current = trange.low();

	snum = ssnum;
	rnum = rrnum;
	tnum = ttnum;
	bnum = bbnum;
	sampleS = samples;
	sampleR = sampler;

	// psis = std::vector<double> (xS, xS + sizeof (* xS));
	// psiweights = std::vector<double> (wS, wS + sizeof (* wS));

	psiweights = std::vector<double> (sampleS);
	psis = std::vector<double> (sampleS);

	for (long i = 0; i < sampleS; i++) {
		psis[i] = xS[i];
		psiweights[i] = wS[i];
	}

	// psir = std::vector<double> (xR, xR + sizeof (* xR));
	// psiweightr = std::vector<double> (wR, wR + sizeof (* wR));

	psiweightr = std::vector<double> (sampleR);
	psir = std::vector<double> (sampleR);

	for (long i = 0; i < sampleR; i++) {
		psir[i] = xR[i];
		psiweightr[i] = wR[i];
	}

	ptyp = ptype;

	if (ptyp == Dynamic) bnum = bbnum_dy;

	// for (long p = 0; p < (sampleR); p++) std::cout << psir[p] << " ";

	init();

}

void GHQC::init() {
	h = (trange.high() - trange.low()) /  tnum;

	XARR = std::vector<double> (snum + 1);
	SARR = std::vector<double> (snum + 1);
	RARR = std::vector<double> (rnum + 1);
	TARR = std::vector<double> (tnum + 1);
	BARR = std::vector<double> (bnum + 1);
	tmp = std::vector<double> ((snum + 1) * (rnum + 1));
	res = std::vector<double> ((bnum + 1) * (snum + 1) * (rnum + 1));

	samplePointS = std::vector<double> ((snum + 1) * sampleS * (rnum + 1) * sampleR);
	samplePointR = std::vector<double> (sampleS * (rnum + 1) * sampleR);

	XARR = xrange.mesh(snum);
	RARR = rrange.mesh(rnum);
	TARR = trange.mesh(tnum);
	BARR = brange.mesh(bnum);

	for (long i = 0; i < (snum + 1); i++) {
		SARR[i] = P0 * exp(XARR[i]);
	}

	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (snum + 1); j++) {
			for (long k = 0; k < sampleR; k++) {
				for (long m = 0; m < sampleS; m++) {
					samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m] = psiTow(psis[m], psir[k], XARR[j] + log(P0), RARR[i], h);
					// std::cout << samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m] << " ";
					// std::cout << h << " ";
				}
			}
		}
	}

	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < sampleR; k++) {
			for (long m = 0; m < sampleS; m++) {
				samplePointR[i * sampleR * sampleS + k * sampleS + m] = psiTor(psis[m], psir[k], 0.0, RARR[i], h);
			}
		}
	}
}

void GHQC::start() {
	for (long bb = 0; bb < (bnum + 1); bb++) {
		for (long i = 0; i < (rnum + 1); i++) {
			for (long k = 0; k < (snum + 1); k++) {
				res[bb * (rnum + 1) * (snum + 1) + i * (snum + 1) + k] = ic(SARR[k], BARR[bb]);
			}
		}
	}

	while (!(finished())) {
		for (long bb = 0; bb < (bnum + 1); bb++) {
			for (long i = 0; i < (rnum + 1); i++) {
				for (long k = 0; k < (snum + 1); k++) {
					tmp[i * (snum + 1) + k] = res[bb * (rnum + 1) * (snum + 1) + i * (snum + 1) + k];
				}
			}

			advance();

			for (long i = 0; i < (rnum + 1); i++) {
				for (long k = 0; k < (snum + 1); k++) {
					res[bb * (rnum + 1) * (snum + 1) + i * (snum + 1) + k] = tmp[i * (snum + 1) + k];
				}
			}
		}

		current += h;

		if (((trunc(current + 1e-5) + 1e-5) >= current) && ((trunc(current + 1e-5) - 1e-5) <= current) && (!finished())) {
			if (ptyp == Static)  static_opti();
			if (ptyp == Dynamic)  dy_opti();
		}
	}
}

void GHQC::reset() {
	tmp = std::vector<double> ((snum + 1) * (rnum + 1));
	res = std::vector<double> ((bnum + 1) * (snum + 1) * (rnum + 1));

	samplePointS = std::vector<double> ((snum + 1) * sampleS * (rnum + 1) * sampleR);
	samplePointR = std::vector<double> (sampleS * (rnum + 1) * sampleR);
	
	current = trange.low();

	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (snum + 1); j++) {
			for (long k = 0; k < sampleR; k++) {
				for (long m = 0; m < sampleS; m++) {
					samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m] = psiTow(psis[m], psir[k], XARR[j] + log(P0), RARR[i], h);
					// std::cout << samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m] << " ";
					// std::cout << h << " ";
				}
			}
		}
	}

	for (long i = 0; i < (rnum + 1); i++) {
		for (long k = 0; k < sampleR; k++) {
			for (long m = 0; m < sampleS; m++) {
				samplePointR[i * sampleR * sampleS + k * sampleS + m] = psiTor(psis[m], psir[k], 0.0, RARR[i], h);
			}
		}
	}
}

void GHQC::advance() {
	std::vector<double> temp = tmp;
	std::vector<double> temptemp = tmp;

	alglib::real_1d_array x, y, f;

	x.setlength(snum + 1);
	y.setlength(rnum + 1);
	f.setlength((snum + 1) * (rnum + 1));

	for (long j = 0; j < (snum + 1); j++) {
		x[j] = SARR[j];	
	}
	for (long i = 0; i < (rnum + 1); i++) {
		y[i] = RARR[i];
	}

	for (long k = 0; k < ((snum + 1) * (rnum + 1)); k++) {
		f[k] = temp[k];
	}

    // double vx = 0.25;
    // double vy = 0.50;
	double v;

	alglib::spline2dinterpolant s;

    // build spline
	spline2dbuildbicubicv(x, snum + 1, y, rnum + 1, f, 1, s);

    // calculate S(0.25,0.50)
    // v = spline2dcalc(s, vx, vy);
    // printf("%.4f\n", double(v)); // EXPECTED: 1.0625

	// std:: cout << tmp[0] << " ";

	// std::cout << current << " ";

	for (long i = 0; i < (rnum + 1); i++) {
		for (long j = 0; j < (snum + 1); j++) {
			tmp[i * (snum + 1) + j] = 0;
			
			for (long k = 0; k < (sampleR); k++) {
				for (long m = 0; m < (sampleS); m++) {
					v = spline2dcalc(s, samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m], samplePointR[i * sampleR * sampleS + k * sampleS + m]);
					tmp[i * (snum + 1) + j] += (psiweights[m] * psiweightr[k] * double(v));
					// if ((i == 0) && (j == 10)) std::cout << tmp[i * (snum + 1) + j] << " " << v << " " << psiweights[m] << "\t";
				}
			}

					// tmp[i * (snum + 1) + j] *= (pfun(Tto - (current + h), Tto - current, RARR[i]) / 3.14159265358979323846);
			// if ((i == 0) && (j == 20)) std::cout << tmp[i * (snum + 1) + j] << " ";

			tmp[i * (snum + 1) + j] = tmp[i * (snum + 1) + j] * pfun(Tto - (current + h), Tto - current, RARR[i]) / 3.14159265358979323846;

			// if ((i == 0) && (j == 20)) std::cout << pfun(Tto - (current + h), Tto - current, RARR[i]) << " " << tmp[i * (snum + 1) + j] << " ";
		}
	}

	// std::cout << "\n";

	// for (long i = 0; i < (rnum + 1); i++) {
	// 	for (long j = 0; j < (snum + 1); j++) {
	// 		tmp[i * (snum + 1) + j] = 0;

	// 		for (long k = 0; k < (rnum + 1); k++) {
	// 			for (long m = 0; m < (snum + 1); m++) {
	// 				sp_s[m] = temp[k * (snum + 1) + m];
	// 			}
	// 			sp.set_points(SARR, sp_s);

	// 			for (long m = 0; m < (sampleS); m++) {
	// 				temptemp[k * (sampleS) + m] = sp(samplePointS[i * (snum + 1) * sampleR * sampleS + j * sampleR * sampleS + k * sampleS + m]);
	// 			}

	// 			// if (k == 0) std::cout << temptemp[0] << " " << i << " " << j << " " << current << "\n";
	// 		}

	// 		for (long k = 0; k < (sampleS); k++) {
	// 			for (long m = 0; m < (rnum + 1); m++) {
	// 				sp_r[m] = temptemp[m * (sampleS) + k];
	// 			}
	// 			sp.set_points(RARR, sp_r);

	// 			for (long m = 0; m < (sampleR); m++) {
	// 				samplePoint[m * (sampleS) + k] = sp(samplePointR[i * sampleR * sampleS + m * sampleS + k]);
	// 			}
	// 			// if ((k == 0) && (i == 0) && (j == 0)) std::cout << samplePoint[0] << " " << " " << current << " ";
	// 		}

	// 		for (long k = 0; k < (sampleR); k++) {
	// 			for (long m = 0; m < (sampleS); m++) {
	// 				// if ((i == 0) && (j == 0)) {
	// 				// 	std::cout << tmp[i * (snum + 1) + j] << " ";
	// 				// }
	// 				// if ((k == 0) && (i == 0) && (j == 0)) std::cout << samplePoint[k * (sampleS) + m] << " ";
	// 				tmp[i * (snum + 1) + j] += (psiweights[m] * psiweightr[k] * samplePoint[k * (sampleS) + m]);
	// 			}
	// 		}

	// 		// if ((i == 0) && (j == 0)) std::cout << "\n";

	// 		tmp[i * (snum + 1) + j] *= (pfun(Tto - (current + h), Tto - current, RARR[i]) / 3.14159265358979323846);

	// 		// if ((i == 0) && (j == 0)) std::cout << tmp[i * (snum + 1) + j] << "\n";
	// 	}
	// }

	// std::cout << temptemp[0] << " " << current << "\n";
}

void GHQC::static_opti() {
	std::vector<double> opttemp = res;
	for (long bb = 1; bb < (bnum + 1); bb++) {
		for (long i = 0; i < (rnum + 1); i++) {
			std::vector<double> bestvalue((snum + 1), 0);
			long bbsub = bb - 1;
			double withdrawn = G;
			tk::spline sp;
			std::vector<double> sp_x = SARR;
			std::vector<double> sp_y(snum + 1);
			for (long k = 0; k < (snum + 1); k++) {
				sp_x[k] = (sp_x[k] > withdrawn)? (sp_x[k] - withdrawn): SARR[0]; 
				sp_y[k] = opttemp[bbsub * (rnum + 1) * (snum + 1) + i * (snum + 1) + k];
			}
			sp.set_points(SARR, sp_y);
			for (long k = 0; k < (snum + 1); k++) {
				double valuenow = sp(sp_x[k]) + withdrawn - withdrawalPen * ((withdrawn > G)? (withdrawn - G): 0);
				bestvalue[k] = (bestvalue[k] > valuenow)? bestvalue[k]: valuenow;
			}

			for (long k = 0; k < (snum + 1); k++) {
				res[bb * (rnum + 1) * (snum + 1) + i * (snum + 1) + k] = bestvalue[k];
			}
		}
	}
}

void GHQC::dy_opti() {
	std::vector<double> opttemp = res;
	for (long bb = 1; bb < (bnum + 1); bb++) {
		for (long i = 0; i < (rnum + 1); i++) {
			std::vector<double> bestvalue((snum + 1), 0);
			for (long bbsub = 0; bbsub < (bb + 1); bbsub++) {
				double withdrawn = BARR[bb] - BARR[bbsub];
				tk::spline sp;
				std::vector<double> sp_x = SARR;
				std::vector<double> sp_y(snum + 1);
				for (long k = 0; k < (snum + 1); k++) {
					sp_x[k] = (sp_x[k] > withdrawn)? (sp_x[k] - withdrawn): SARR[0]; 
					sp_y[k] = opttemp[bbsub * (rnum + 1) * (snum + 1) + i * (snum + 1) + k];
				}
				sp.set_points(SARR, sp_y);
				for (long k = 0; k < (snum + 1); k++) {
					double valuenow = sp(sp_x[k]) + withdrawn - withdrawalPen * ((withdrawn > G)? (withdrawn - G): 0);
					bestvalue[k] = (bestvalue[k] > valuenow)? bestvalue[k]: valuenow;
				}
			}

			for (long k = 0; k < (snum + 1); k++) {
				res[bb * (rnum + 1) * (snum + 1) + i * (snum + 1) + k] = bestvalue[k];
			}
		}
	}
}

bool GHQC::finished() {
	if (current >= (trange.high() - 1e-5)) {return true;}
	else {return false;}
}

std::vector<double> GHQC::getSARR() {
	return SARR;
}

std::vector<double> GHQC::getRARR() {
	return RARR;
}

std::vector<double> GHQC::getBARR() {
	return BARR;
}

std::vector<double> GHQC::getres() {
	return res;
}

#endif