#ifndef FFAIR_P_HPP
#define FFAIR_P_HPP

class FindFairP {
	private:
		GHQC scheme;

		std::vector<double> sarr;
		std::vector<double> rarr;
		std::vector<double> barr;

		std::vector<double> arr;

		double sTarget;
		double rTarget;

		double sWeight;
		double rWeight;

		long snum;
		long rnum;

		long sIndex;
		long rIndex;
		long bIndex;

		long maxIter;
		double Tol;

		double interestedV;
		double lastinterestedV;
		double lastlastinterestedV;

		double lastfair;
		double lastlastfair;

		long counter;

		void init();

	public:
		FindFairP(GHQC& input, double STARGET, double RTARGET, long BINDEX,
			long MAXITER, double TOL);
		void reset(GHQC& input);
		void start();

		//setter
		void setInput(GHQC& input);

		//getter
		double getsWeight();
		double getrWeight();

		long getsIndex();
		long getrIndex();

		double getInterestedV();
};

#endif