#include <octave/oct.h>
#include <omp.h>

DEFUN_DLD(setnthreads, args, ,"Usage: setnthreads(<nthreads>)"){
	octave_value_list retval;

	int nthreads = args(0).scalar_value();
	omp_set_num_threads(nthreads);

	retval.append(1);
	return retval;
}
	

