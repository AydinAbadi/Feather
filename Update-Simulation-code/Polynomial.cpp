/*

- Polynomial class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/
//**********************************************************************

#include"Polynomial.h"
//**********************************************************************
// - Description: Constrcutor- given an array of elements and x-coordinates constructs a polynomial in the point-value form, i.e. calculates y-coordinates.

Polynomial::Polynomial(bigint* elem, string ID, bigint * xpoints, int elem_size, int xpoints_size, bigint pubmoduli){

	bigint ranx, ran, minus_one, *temp, *elem_temp;
	Random rd;
	Random rdx;
	gmp_randstate_t rand;
	gmp_randstate_t rand1;
	gmp_randinit_default(rand1);
	gmp_randstate_t randx;
	gmp_randstate_t rand1x;
	gmp_randinit_default(rand1x);
	rd.init_rand3(rand, ran, 8);
	gmp_randseed(rand1,ran);
	poly_ID = ID;
	mpz_init_set_str(minus_one, "-1", 10);
	temp = (mpz_t*)malloc(xpoints_size * sizeof(mpz_t));
	val_size = xpoints_size;
	rdx.init_rand3(randx, ranx, 8);
	gmp_randseed(rand1x,ranx);
	elem_temp = (mpz_t*)malloc(xpoints_size * sizeof(mpz_t));
	if(mpz_cmp(elem[0], minus_one) == 0){ // if eleme[0]==-1, then fills all elemements of the array: values, with random values.
		for(int i = 0; i < xpoints_size; i++){
			mpz_init(temp[i]);
			mpz_urandomb(temp[i], rand1, 50);
		}
		values = temp;
	}
	else{
		for(int j = 0; j < elem_size; j++){
			mpz_init_set(elem_temp[j], elem[j]);
			if(mpz_cmp(elem_temp[j],minus_one) == 0){
				mpz_urandomb(elem_temp[j], rand1x, 50);
			}
		}
		// computes y-coordibates by evaluating the polynomial of the form (x-elem[0])...(x-elem[elem_size-1]) at every x-coordinate.
		values = evaluate (elem_temp, xpoints, elem_size, xpoints_size, pubmoduli);
	}
}
//**********************************************************************
// - Function description: given an array of roots of a polynomial and array of x-coordinates, it evaluates
// the polynomial at each x-coordinate and returns an array of y-coordinates.

bigint* Polynomial::evaluate(bigint* elem, bigint* xp, int ele_size, int xp_size, bigint pubmod){

	bigint mult2, *val, temp, one;
	mpz_init(temp);
	mpz_init_set_str(one, "1", 10);
	val = (mpz_t*)malloc(xp_size * sizeof(mpz_t));
	for (int i = 0; i < xp_size; i++){
		mpz_init(val[i]);
		mpz_init_set(mult2, one);
		for (int j = 0; j < ele_size; j++){ //ele_size: # of elements in a bucket
			mpz_sub(temp, xp[i], elem[j]);
			mpz_mul(mult2, mult2, temp);
			mpz_mod(mult2, mult2, pubmod);
			mpz_set(val[i], mult2);
		}
	}
	return val;
}
//**********************************************************************
// - Function description: blinds the y-coordinates of a polynomial.

void Polynomial::blind_poly(bigint seed, bigint pubmod, int bit_size){

	bigint *pr_val, *ptr;
	ptr = get_values();
	pr_val = (mpz_t*)malloc(val_size * sizeof(mpz_t));
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	for(int i = 0; i < val_size; i++){
		mpz_init(pr_val[i]);
		mpz_urandomb(pr_val[i], rand, bit_size);
		mpz_add(pr_val[i], ptr[i], pr_val[i]); // blinds each y-coordinate
		mpz_mod(pr_val[i], pr_val[i], pubmod);
	}
	values = pr_val;
}

//**********************************************************************
// - Function description: returns x-coordinates.

bigint* Polynomial::get_values(){

	return values;
}
//**********************************************************************
// - Function description: given an array of coeeficients of a polynomial, and array of x-coordinates, it evaluates
// the polynomial at each x-coordinate and returns an array of y-coordinates.

bigint* Polynomial:: evaluate_coeffs(bigint* coeff, bigint* x_points, int coeff_size, int xpoint_size, bigint pubmoduli){

	bigint *res, temp;
	res = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	for(int i = 0; i < xpoint_size; i++){
		mpz_init(res[i]);
		for(int j = 0;j < coeff_size; j++){
			mpz_init(temp);
			mpz_powm_ui(temp, x_points[i],j, pubmoduli);
			mpz_mul(temp, temp, coeff[j]);
			mpz_add(res[i], temp, res[i]);
		}
		mpz_mod(res[i], res[i], pubmoduli);
	}
	mpz_clear(temp);
	return res;
}
//**********************************************************************
// - Function description: given an array of x-coordinates: a and array of y-coordinates: b, it interpolates
// a polynomial and returns an array containing the polynomial's coefficients.

bigint* Polynomial::interpolate(int size, bigint* a, bigint* b, bigint N){
	long m = size;
	bigint* prod;
	prod = (mpz_t*)malloc(size * sizeof(mpz_t));
	 for(int i = 0; i < size; i++){
		 mpz_init_set(prod[i], a[i]);
	 }
	 bigint t1, t2;
	mpz_init(t1);
	mpz_init(t2);
	int k, i;
	bigint* res;
	res = (mpz_t*)malloc(size * sizeof(mpz_t));
	bigint aa;
	for (k = 0; k < m; k++) {
		 mpz_init_set(aa , a[k]);
		 mpz_init_set_str(t1, "1", 10);
		 for (i = k-1; i >= 0; i--) {
			 mpz_mul(t1, t1, aa);
			 mpz_mod(t1, t1, N);
			 mpz_add(t1, t1, prod[i]);
		 }
		 mpz_init_set_str(t2, "0", 10);
		 for (i = k-1; i >= 0; i--) {
			 mpz_mul(t2, t2, aa);
			 mpz_mod(t2, t2,N);
			 mpz_add(t2, t2, res[i]);
		 }
		 mpz_invert(t1, t1, N);
		 mpz_sub(t2, b[k], t2);
		 mpz_mul(t1, t1, t2);
		 for (i = 0; i < k; i++) {
			 mpz_mul(t2, prod[i], t1);
			 mpz_mod(t2, t2,N);
			 mpz_add(res[i], res[i], t2);
			 mpz_mod(res[i], res[i], N);
		 }
		 mpz_init_set(res[k], t1);
		 mpz_mod(res[k], res[k], N);
		 if (k < m - 1) {
			 if (k == 0)
			 mpz_neg(prod[0], prod[0]);
			 else {
				 mpz_neg(t1, a[k]);
				 mpz_add(prod[k], t1, prod[k-1]);
				 for (i = k - 1; i >= 1; i--) {
					 mpz_mul(t2, prod[i], t1);
					 mpz_mod(t2, t2,N);
					 mpz_add(prod[i], t2, prod[i-1]);
				 }
				 mpz_mul(prod[0], prod[0], t1);
				 mpz_mod(prod[0], prod[0],N);
			 }
		 }
	 }
	 while (m > 0 && (res[m-1] == 0)) m--;
	 mpz_clear(t1);
	 mpz_clear(t2);
	 return res;
}
//**********************************************************************
