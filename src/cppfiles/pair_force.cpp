#include <octave/oct.h>
#include <omp.h>


#define HELP ("Usage:  \n")

void force(double *Fptr, const double *rptr, const unsigned nagents, const unsigned max_nagents, 
															const double cf, const double A){
	double dr[3], dist, Fval, fac;
	const unsigned npart = max_nagents;
	unsigned n, m, k, lvec = 3*npart;

#pragma omp parallel for if(npart>200) schedule(dynamic)	\
		private(dr, dist, Fval, fac, n, m, k)	\
		reduction(+:Fptr[:lvec]) 	
	for ( n=0; n<nagents-1; n++ ){
		for ( m=n+1; m<nagents; m++ ){
		
			for ( k=0; k<3; k++ ) 
				dr[k] = rptr[npart*k + n] - rptr[npart*k + m];

			dist = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
			if ( dist > 1e-3 && dist < cf ){
				Fval = -exp(-(dist-1.0)) + A*exp(-A*(dist-1.0));
				fac = Fval/dist;
				for ( k=0; k<3; k++ ){
					Fptr[npart*k + n] += dr[k]*fac; 
					Fptr[npart*k + m] -= dr[k]*fac;
				}
			}

		}
	}

}

DEFUN_DLD(pair_force, args, ,HELP){
	octave_value_list retval;

	Matrix F(args(0).matrix_value());
	Matrix r(args(1).matrix_value());
	double cf = args(2).scalar_value();
	double A = args(3).scalar_value();
	const unsigned nagents = (unsigned)args(4).scalar_value(); 
	const unsigned max_nagents = F.rows();

	const double *rptr = r.fortran_vec();
	double *Fptr = F.fortran_vec();

	force(Fptr, rptr, nagents, max_nagents, cf, A);

	retval.append(F);	

	return retval;


}

