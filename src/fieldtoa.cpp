#include <octave/oct.h>
#include <omp.h>

#define HELP ("Usage:  \n")

DEFUN_DLD(fieldtoa, args, ,HELP){

	octave_value_list retval;

	Matrix r(args(0).matrix_value());
	NDArray field(args(1).array_value());
	const double lbox = args(2).scalar_value();
	
	const unsigned ngrid = field.dim1();
	const unsigned nagents = r.rows();
	const double lgrid = lbox/ngrid;

	RowVector a(nagents);

	unsigned idx[3];
	for ( unsigned n=0; n<nagents; n++ ){
		for ( unsigned k=0; k<3; k++ ){
			unsigned idgrd = (unsigned)((r(n, k)+lbox*0.5)/lgrid);
			if ( idgrd >= ngrid ) idgrd = ngrid-1;
			idx[k] = idgrd;
		}
		a(n) = field(idx[0], idx[1], idx[2]);
	}
	
	retval.append(a);
	return retval;
}
