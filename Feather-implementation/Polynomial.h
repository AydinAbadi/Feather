/*

- Polynomial class used by the clients and server in the Feather protocol.

*/

/*
 Variables description:

	values: y-coordinates.
	poly_ID: client's ID.
	val_size: number of y-coordinates.
*/
//**********************************************************************

#include "Hashtable.h"
//*********************************************************************

class Polynomial{

public:

	Polynomial(){};
	Polynomial(bigint* elem, string poly_ID, bigint * xpoints, int elem_size, int xpoints_size, bigint pubmoduli);
	Polynomial(bigint* elem, bigint * xpoints, int elem_size, int xpoints_size, bigint pubmoduli, int pubmoduli_size, unordered_map <string, int> map);
	bigint*  evaluate (bigint *, bigint * , int, int, bigint);
	bigint* evaluate_coeffs(bigint* coeff, bigint* x_points, int coeff_size, int xpoint_size, bigint pubmoduli);
	bigint* interpolate(int size, bigint* a, bigint* b, bigint N);
	bigint* get_values();
	void  blind_poly (bigint , bigint, int);
	void set_values(bigint* vals, int size);
	void blind_poly_(byte *key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod);
	void unblind_poly_(byte *key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod);
	bigint* unblind_poly_(bigint* blinded_vals, int num_of_vals, byte* key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod);
	bigint* blind_poly_(bigint* y_coordinates, int  val_size_, byte *key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod);
	//variables
	bigint* values;
	int val_size;
};
