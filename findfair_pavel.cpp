#ifndef FFAIR_P_CPP
#define FFAIR_P_CPP

#include "findfair_pavel.hpp"

void FindFairP::init() {
	sarr = scheme.getSARR();
	rarr = scheme.getRARR();
	barr = scheme.getBARR();

	snum = sarr.size() - 1;
	rnum = rarr.size() - 1;

	bIndex = barr.size() - 1;

	fairfee = 0.0;
	lastlastfair = fairfee;

	scheme.reset();
	scheme.start();
	arr = scheme.getres();

	for (long ss = 0; ss < (snum + 1); ss++) {
		if (sarr[ss] <= sTarget) {
			sIndex = ss;
		} else {
			break;
		}
	}
	sWeight = (sTarget - sarr[sIndex]) / (sarr[sIndex + 1] - sarr[sIndex]);

	for (long rr = 0; rr < (rnum + 1); rr++) {
		if (rarr[rr] <= rTarget) {
			rIndex = rr;
		} else {
			break;
		}
	}
	rWeight = (rTarget - rarr[rIndex]) / (rarr[rIndex + 1] - rarr[rIndex]);

	lastlastinterestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + sIndex]
						+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + sIndex])
							+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + (sIndex + 1)]
							+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + (sIndex + 1)]);

	// lastlastinterestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
	// 					+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex])
	// 						+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
	// 						+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]);

	// std::cout << sIndex << " " << sWeight << " " << rIndex << " " << rWeight << " " << lastlastinterestedV << "\n";

	fairfee = 0.02;
	lastfair = fairfee;

	scheme.reset();
	scheme.start();
	arr = scheme.getres();

	for (long ss = 0; ss < (snum + 1); ss++) {
		if (sarr[ss] <= sTarget) {
			sIndex = ss;
		} else {
			break;
		}
	}
	sWeight = (sTarget - sarr[sIndex]) / (sarr[sIndex + 1] - sarr[sIndex]);

	for (long rr = 0; rr < (rnum + 1); rr++) {
		if (rarr[rr] <= rTarget) {
			rIndex = rr;
		} else {
			break;
		}
	}
	rWeight = (rTarget - rarr[rIndex]) / (rarr[rIndex + 1] - rarr[rIndex]);

	lastinterestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + sIndex]
						+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + sIndex])
							+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + (sIndex + 1)]
							+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + (sIndex + 1)]);

	// lastinterestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
	// 					+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex])
	// 						+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
	// 						+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]);
	
	interestedV = lastinterestedV;

	// std::cout << sIndex << " " << sWeight << " " << rIndex << " " << rWeight << " " << interestedV << "\n";

	counter = 0;
}

FindFairP::FindFairP(GHQC& input, double STARGET, double RTARGET, long BINDEX,
			long MAXITER, double TOL) {
	scheme = input;

	sTarget = STARGET;
	rTarget = RTARGET;
	// bIndex = BINDEX;

	maxIter = MAXITER;
	Tol = TOL;

	init();
}

void FindFairP::reset(GHQC& input) {
	scheme = input;
	
	init();
}

void FindFairP::start() {
	while ((counter < maxIter) && (((interestedV - sTarget) > Tol) || ((sTarget - interestedV) > Tol))) {
		fairfee = (lastlastfair * (lastinterestedV - sTarget) - lastfair * (lastlastinterestedV - sTarget)) / ((lastinterestedV - sTarget) - (lastlastinterestedV - sTarget));
		scheme.reset();
		scheme.start();
		arr = scheme.getres();
		interestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + sIndex]
						+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + sIndex])
							+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (snum + 1) + rIndex * (snum + 1) + (sIndex + 1)]
							+rWeight * arr[bIndex * (rnum + 1) * (snum + 1) + (rIndex + 1) * (snum + 1) + (sIndex + 1)]);

		// interestedV = (1 - sWeight) * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex]
		// 				+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + sIndex])
		// 					+ sWeight * ((1 - rWeight) * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + rIndex * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]
		// 					+ rWeight * arr[bIndex * (rnum + 1) * (vnum + 1) * (snum + 1) + (rIndex + 1) * (vnum + 1) * (snum + 1) + vIndex * (snum + 1) + (sIndex + 1)]);

		lastlastfair = lastfair;
		lastlastinterestedV = lastinterestedV;

		lastfair = fairfee;
		lastinterestedV = interestedV;
		
		counter ++;
	}
}

void FindFairP::setInput(GHQC& input) {
	scheme = input;
}

double FindFairP::getsWeight() {
	return sWeight;
}

double FindFairP::getrWeight() {
	return rWeight;
}

long FindFairP::getsIndex() {
	return sIndex;
}

long FindFairP::getrIndex() {
	return rIndex;
}

double FindFairP::getInterestedV() {
	return interestedV;
}

#endif