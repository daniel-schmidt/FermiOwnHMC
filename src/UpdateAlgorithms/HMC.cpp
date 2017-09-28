/*
 * HMC.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "HMC.h"

namespace FermiOwn {

HMC::HMC( double t, size_t nt, const BasicAction& naction, Field<Real>& nphi, std::ranlux48* rndGen ) :
	randomGenerator(rndGen),
	act(naction),
	phi(nphi),
	momentum(phi.getSize(), 1, randomGenerator, gaussianInit),
	integrator(phi, momentum, act, t, nt),
	dH(0.)
{
}

HMC::~HMC() {
}

bool HMC::update() {

	// integration
	momentum.setGaussian();
	double H_old = Hamiltonian();
	Field<Real> old = phi;
	integrator.integrate();
	double H_new = Hamiltonian();
	dH = H_old - H_new;
// 	std::cout << "difference of action: " << H_old-H_new << " exponential=" << dH << std::endl;
	// accept/reject step
	bool accepted = false;
	double r = uniformDistribution01(*randomGenerator);
	if( r < std::exp(dH) )
		accepted = true;
	else
		phi = old;
// 	std::cout << "H_old: " << H_old << " H_new: "<< H_new << " dH: " << dH << " r: " << r << " accepted: " << accepted << std::endl;

	return accepted;
}

// Private functions

double HMC::Hamiltonian() {
	double H = act.getAction(phi);
	H += 0.5 * momentum.cwiseMultAndSum(momentum);
	return H;
}

} // namespace FermiOwn
