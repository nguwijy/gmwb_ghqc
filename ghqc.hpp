#ifndef GHQC_HPP
#define GHQC_HPP

enum polType {Static, Dynamic};

class GHQC {
private:
	double (*psiTow) (const double, const double, const double, const double, const double);
	double (*psiTor) (const double, const double, const double, const double, const double);
	double (*pfun) (const double, const double, const double);
	double (*ic) (const double, const double);

	unsigned long snum;
	unsigned long rnum;
	unsigned long tnum;
	unsigned long bnum;
	unsigned long sampleS;
	unsigned long sampleR;

	polType ptyp;

		double h;					//del T

		double current;				//time to maturity, start from maturity, current = 0

		std::vector<double> SARR;
		std::vector<double> XARR;
		std::vector<double> RARR;
		std::vector<double> TARR;
		std::vector<double> BARR;

		std::vector<double> psis;
		std::vector<double> psiweights;
		std::vector<double> psir;
		std::vector<double> psiweightr;

		std::vector<double> samplePointS;
		std::vector<double> samplePointR;

		Trange xrange;
		Trange rrange;
		Trange trange;
		Trange brange;

		std::vector<double> tmp;	//value in current time as specified by B
		std::vector<double> res;	//value in current time in overall

		void init();
		void advance();

	public:
		GHQC();
		GHQC(double (*psiW) (const double, const double, const double, const double, const double), 
			double (*psiR) (const double, const double, const double, const double, const double),
			double (*pFun) (const double, const double, const double),
			double (*IC) (const double, const double),
			Trange xxrange, Trange rrrange, Trange ttrange, Trange bbrange,
			unsigned long ssnum, unsigned long rrnum, unsigned long ttnum, unsigned long bbnum, unsigned long bbnum_dy,
			alglib::real_1d_array xS, alglib::real_1d_array wS, alglib::real_1d_array xR, alglib::real_1d_array wR, unsigned long samples, unsigned long sampler, polType ptype);
		void start();
		void reset();
		bool finished();

		void static_opti();
		void dy_opti();

		std::vector<double> getSARR();
		std::vector<double> getRARR();
		std::vector<double> getBARR();
		std::vector<double> getres();
	};

#endif
