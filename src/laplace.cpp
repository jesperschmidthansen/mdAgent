#include <octave/oct.h>
#include <omp.h>

#define HELP ("Usage:  \n")

DEFUN_DLD(laplace, args, ,HELP){
	octave_value_list retval;

	NDArray u = args(0).array_value ();
	const double h = args(1).scalar_value();

	NDArray L( u.dims() );

	const int dims = u.dim1(); 
	const double hsq = h*h;

	for ( int i=1; i<dims-1; i++ ){ 
		for ( int j=1; j<dims-1; j++ ){
			for ( int k=1; k<dims-1; k++ ){
				double d2u_dx2 = (u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k)) / hsq;
				double d2u_dy2 = (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k)) / hsq;
				double d2u_dz2 = (u(i,j,k+1) - 2*u(i,j,k) + u(i,j,k-1)) / hsq;

				L(i,j,k) = d2u_dx2 + d2u_dy2 + d2u_dz2;
			}
		}
	}

	retval.append(L);

	return retval;
}	
