
/*

- Polynomial class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/

/*
 Variables description:

	values: y-coordinates.
	poly_ID: client's ID.
	val_size: number of y-coordinates.
*/
//**********************************************************************

#include "Hashtable.h"
//**********************************************************************

class Polynomial{

public:
	Polynomial(){};
	Polynomial(bigint* elem, string poly_ID, bigint * xpoints, int elem_size, int xpoints_size, bigint pubmoduli);
	bigint*  evaluate (bigint *, bigint * , int, int, bigint);
	bigint* evaluate_coeffs(bigint* coeff, bigint* x_points, int coeff_size, int xpoint_size, bigint pubmoduli);
	bigint* interpolate(int size, bigint* a, bigint* b, bigint N);
	bigint* get_values();
	void  blind_poly (bigint , bigint, int);
	void set_values(bigint* vals, int size);
	string get_poly_ID(){return poly_ID;}
	// variables
	bigint* values;
	string  poly_ID;
	int val_size;
};
//**********************************************************************
