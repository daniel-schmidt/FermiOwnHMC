/*
 * HMC_test.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "Action.h"
#include "HMC.h"

#include "Field.h"
#include "Lattice.h"

int main() {
	using namespace FermiOwn;
	std::cout << "Testing the HMC algorithm..." << std::endl;

	// initialize random number generator
	std::ranlux48 rndGen;

	// initialize lattice
	size_t Nt = 8, Ns = 8, dim = 3;
	Lattice lat(Nt, Ns, dim);

	// initialize action
	double kappa = 0.28, lambda = 1.145;
	Action act(lat, kappa, lambda);

	// initialize field
	Field<Real> phi(lat.getVol(), 1, &rndGen, randomInit);

	// initialize HMC

	for( size_t HMCnt = 10; HMCnt <= 1e6; HMCnt*= 10 )	{
		HMC updater( 1., HMCnt, act, phi, &rndGen );
		updater.update();
        std::cout << HMCnt<< "\t" << updater.getdH() << std::endl;
	}

	size_t updates = 1000;
	size_t thermal = 1000;
	size_t stepsPerUpdate=5;


        double HMC_len = 1.;
	size_t HMC_steps = 10;
        
        std::ofstream fMagAbs("magAbsHMC.dat");
        
        for( double kappa = 0.10; kappa < 0.30; kappa += 0.005 )
	{
            Action act( lat, kappa, lambda );
            HMC updater( HMC_len, HMC_steps, act, phi, &rndGen );
            
            size_t acceptance = 0.;
            double av_mag = 0.;
            
            for( size_t th = 0; th < thermal; th++ ) {
                    updater.update();
            }

            for( size_t u = 0; u<updates; u++ ) {
                    for( size_t spu = 0; spu < stepsPerUpdate; spu++ )
                    {
                            bool accepted = updater.update();
                            if(accepted) acceptance++;
                    }

                    double magnetization = 0.;
                    for( size_t x = 0; x < phi.getSize(); x++ )
                            magnetization += phi(x);
                    av_mag += std::fabs(magnetization)/phi.getSize();
            }
            double accRate = acceptance/double(updates*stepsPerUpdate);
            av_mag /= double(updates);

            std::cout << "Accepted "<< accRate << std::endl;
            std::cout << "Average magnetization: " << av_mag << std::endl;
            fMagAbs << kappa << "\t" << av_mag << "\t" << accRate << std::endl;
	}
	fMagAbs.close();
}
