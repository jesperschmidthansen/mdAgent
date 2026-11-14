#include <octave/oct.h>
#include <omp.h>

#define HELP ("Usage:  \n")

DEFUN_DLD(atofield, args, ,HELP){
	octave_value_list retval;

	Matrix r(args(0).matrix_value());
	const unsigned ngrid = (unsigned)args(1).int_value();
	const double lbox = args(2).scalar_value();
	const unsigned npart = r.rows();
	const double lgrid = lbox/ngrid;

	dim_vector arraydim(ngrid, ngrid, ngrid);
	NDArray field(arraydim, 0.0f);

	unsigned idgrd[3];
	for ( unsigned n=0; n<npart; n++ ){
		for ( unsigned k=0; k<3; k++ ){
			idgrd[k] = (unsigned)((r(n, k)+lbox*0.5)/lgrid);
			if ( idgrd[k] >= ngrid )
				idgrd[k] = ngrid-1;
		}
		field(idgrd[0], idgrd[1], idgrd[2]) = field(idgrd[0], idgrd[1], idgrd[2]) + 1.0;
	}

	const double vol = lgrid*lgrid*lgrid;
	for ( unsigned n=0; n<ngrid; n++ )
		for ( unsigned m=0; m<ngrid; m++ )
			for ( unsigned k=0; k<ngrid; k++ )
				field(n,m,k) = field(n,m,k)/vol;


	retval.append(field);

	return retval;
}
