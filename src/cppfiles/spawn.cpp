#include <octave/oct.h>
#include <omp.h>

#define _RAND() ( rand()/(RAND_MAX+1.0) )

#define HELP ("Usage:  \n")

DEFUN_DLD(spawn, args, ,HELP){
	octave_value_list retval;

	if ( args.length() != 7 ){
		octave_stdout << HELP ;
		return retval;
	}

	Matrix r(args(0).matrix_value());
	Matrix v(args(1).matrix_value());
	ColumnVector idx(args(2).vector_value());	
	RowVector random(args(3).vector_value());
	RowVector crit_prob(args(4).vector_value());
	unsigned nagents(args(5).int_value());
	const unsigned max_nagents(args(6).int_value());
	
	const unsigned ntype = idx.rows();

	unsigned nspawn = 0;
	for ( unsigned n=0; n<ntype; n++ ){
		unsigned midx = (unsigned)idx(n) - 1;
		if ( random(n)<crit_prob(midx) ){
			unsigned cidx = nagents;
			if ( cidx+1 == max_nagents ) break;

			for ( int k=0; k<3; k++ ){
				r(cidx,k) = r(midx,k);
				v(cidx,k) = 0.1*(_RAND() - 0.5);
			}
			nagents++;	nspawn++;
		}
	}	

	if ( nspawn>0 ){
		double sum_vel[3]={0.0f};
		for ( unsigned n=0; n<nagents; n++ )
			for ( int k=0; k<3; k++ ) sum_vel[k] += v(n, k);
		
		for ( int k=0; k<3; k++ ) sum_vel[k] = sum_vel[k]/nagents;

		for ( unsigned n=0; n<nagents; n++ )
			for ( int k=0; k<3; k++ ) v(n, k) = v(n, k) - sum_vel[k];
	}

	retval.append(r); 
	retval.append(v); 
	retval.append(nspawn);		

	return retval;
}
