#include <octave/oct.h>
#include <omp.h>

#define HELP ("Usage:  \n")

DEFUN_DLD(surface_force, args, ,HELP){
	octave_value_list retval;

	Matrix F(args(0).matrix_value());
	Matrix r(args(1).matrix_value());
	double z0 = args(2).scalar_value();
	double cf = args(3).scalar_value();
	double A = args(4).scalar_value();
	const unsigned npart = (unsigned)args(5).scalar_value(); 

	const double cf2 = cf*cf;
	
	for ( unsigned n=0; n<npart; n++ ){
		double dr = r(n, 2) - z0; 
		double dr2 = dr*dr;

		if ( dr2 < cf2 ){
			double rri = 1.0/dr2; 
			double rri3 = rri*rri*rri;
			double Fval = 48.0*rri3*(rri3 - 0.5*A)*rri;
			F(n,2) += dr*Fval;
		}
	}
	
	retval.append(F);	

	return retval;
}
