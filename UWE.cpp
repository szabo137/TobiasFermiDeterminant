#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

#include <omp.h>
#include <stdio.h>
#include <sched.h>



#include <boost/program_options.hpp>


#include "matrix.hpp" // Class: Matrix


// GSL headers for determinant / inverse matrix
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>




#include "read_config.hpp"

const double pi = 3.1415926535897932384626433832;





double beta(double rs,double theta)
{
	return pow(( 9.0*pi*0.25 ), (-2.0/3.0)) *2.0*rs*rs/(theta);
}

double n(double rs)
{
	return 1.0 / ( 4.0/3.0 * pi * rs*rs*rs );
}


double L(rs,N)
{
	return pow( 4.0/3.0*pi*pow( rs, 3 )*N*2, 1.0/3.0 );
}

double Z_sum_index( double beta, double L, int max_index )
{
	double ans = 0.0;

 	double E = pi*pi*2.0/L/L;

	for(int x=1;x<1+max_index;x++)
	{
		double tmp = 2.0 * exp( -beta * (E * x*x ) );
		ans += tmp;
	}

 	ans += 1.0;

	return ans*ans*ans;
}







std::vector< std::vector < double > > UEG_N( int N, double beta, double L, int dim, int limit )
{
	std::vector<double> Fermi;
	std::vector<double> Bose;

	Matrix M_Fermi;
	Matrix M_Bose;

	M_Fermi.init( N, N );
	M_Bose.init( N, N );



	std::cout << "----N: " << N << "\n";


	for(int iX=0;iX<N;iX++)
	{
		double element_Fermi = 0.0;
		double element_Bose  = 0.0;

		for(int iY=0;iY<N;iY++)
		{

			if( iX == iY+1 )
			{
				element_Fermi = double( iX );
				element_Bose = -double( iX );
			}
			else if( iX <= iY )
			{
				int multi = iY+1-iX;

				element_Fermi = Z_sum_index( multi*beta, L, limit );
				element_Bose  = element_Fermi;

			}

			M_Fermi.set_value( iX, iY, element_Fermi );
			M_Bose.set_value(  iX, iY, element_Bose  );

		} // end loop iY

	} // end loop iX






	double Determinante_Fermi = M_Fermi.determinant();
	double Z_Fermi = Determinante_Fermi / tgamma( N+1 );

	double F_Fermi = - log(Z_Fermi) / beta;



	double Determinante_Bose = M_Bose.determinant();
	double Z_Bose = Determinante_Bose / tgamma( N+1 );

	double F_Bose = - log(Z_Bose) / beta;



	// Note: What we ultimately want is F_Bose, F_Fermi, and Z_Fermi/Z_Bose


// 	return sth here, whatever :)

}






int main(int ac, char* av[])
{


	read_config(ac, av);
	std::cout << "N: " << N << "\tbeta: " << inv_T << "\tdim: " << dim << "\n";


	// We assume that we have N, beta, dim from "read_config"










	return 1;



}
