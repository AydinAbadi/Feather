/*

- Polynomial class used by the clients and server in the Feather protocol.

*/

#include"Polynomial.h"
//*********************************************************************
// - Description: Constrcutor- given an array of elements and x-coordinates constructs a polynomial in the point-value form, i.e. calculates y-coordinates.
Polynomial::Polynomial(bigint* elem, bigint * xpoints, int elem_size, int xpoints_size, bigint pubmoduli, int pubmoduli_size, unordered_map <string, int> map){

	bigint minus_one, *temp, *elem_temp;
	Random rd;
	mpz_init_set_str(minus_one, "-1", 10);
	val_size = xpoints_size;
	elem_temp = (mpz_t*)malloc(xpoints_size * sizeof(mpz_t));
	if(mpz_cmp(elem[0], minus_one) == 0){ // if eleme[0]==-1, then fills all elemements of the array: values, with random values.
		bigint* temp_= rd.get_nonconflict_randset(xpoints_size, pubmoduli, pubmoduli_size, map, xpoints, xpoints_size);
		values = temp_;
		return;
	}
	else if(mpz_cmp(elem[0], minus_one) != 0) {
		for(int j = 0; j < elem_size; j++){
			mpz_init_set(elem_temp[j], elem[j]);
			if(mpz_cmp(elem_temp[j], minus_one) == 0){
				bigint* temp_= rd.get_nonconflict_randset(1, pubmoduli, pubmoduli_size, map, xpoints, xpoints_size);
				mpz_init_set(elem_temp[j], temp_[0]);
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
	mpz_clear(mult2);
	mpz_clear(temp);
	mpz_clear(one);
	return val;
}
//**********************************************************************
// - Function description: given an array of y-coordinates, it blinds them using fresh pseudorandom values.
bigint* Polynomial::blind_poly_(bigint* y_coordinates, int val_size_, byte *key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod){

	//get y-coordinates
	bigint *ptr;
	ptr = y_coordinates;
	bigint *res_;
	res_ = (mpz_t*)malloc(val_size_ * sizeof(mpz_t));
	//-------- derive a key for indx
	string cipher;
	unsigned char prn_[byte_];
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(key, key_size, iv);
	StringSource s(to_string(indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
	memset(der_key_0, 0x00, key_size + 1);
	strcpy((char*)der_key_0, cipher.c_str());
	cipher.clear();
	//derive a key for counter_indx[indx].
	e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
	StringSource ss(to_string(counter_indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char derived_key_2[key_size]; // convert the ciphertext into a der_key
	memset(derived_key_2, 0x00, key_size + 1);
	strcpy((char*)derived_key_2, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key_2, key_size, iv);	// set key an iv.
	// use the derived key to generate a set of PRNs (for the bin at indx).
	for (int i = 0; i < val_size_; i++){
		cipher.clear();
		temp.clear();
		StringSource sss(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		//conver it to a bigint and store it somewhere.
		temp = cipher.substr (0, byte_);// truncate the ciphertext
		memset(prn_, 0x00, byte_ + 1);
		strcpy((char*)prn_, temp.c_str());
		mpz_init(res_[i]);
		mpz_import(res_[i], byte_, 1, 1, 0, 0, prn_);
		mpz_mod(res_[i], res_[i], pubmod);
		// use the bld factor to blind the y-coordinate
		mpz_add(res_[i], ptr[i], res_[i]); // blinds each y-coordinate
		mpz_mod(res_[i], res_[i], pubmod);
		mpz_clear(ptr[i]);
	}
	cipher.clear();
	temp.clear();
	free(ptr);
	return res_;
}
//**********************************************************************
// - Function description: blindes a polynomial, by fetching y-coordinates, and masking them with fresh PR values. It replaces the
// y-ccordinates with the blinded values.
void Polynomial::blind_poly_(byte *key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod){

	//get y-coordinates
	bigint *ptr;
	ptr = get_values();
	bigint *res_;
	res_ = (mpz_t*)malloc(val_size * sizeof(mpz_t));
	//-------- derive a key for indx
	string cipher;
	unsigned char prn_[byte_];
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(key, key_size, iv);
	StringSource s(to_string(indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
	memset(der_key_0, 0x00, key_size+1);
	strcpy((char*)der_key_0, cipher.c_str());
	cipher.clear();
	//derive a key for counter_indx[indx].
	e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
	StringSource ss(to_string(counter_indx), true,new StreamTransformationFilter(e,new StringSink(cipher)));
	unsigned char derived_key_2[key_size]; // convert the ciphertext into a der_key
	memset(derived_key_2, 0x00, key_size + 1);
	strcpy((char*)derived_key_2,cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key_2, key_size, iv);	// set key an iv.
	// use the derived key to generate a set of PRNs (for the bin at indx).
	for (int i = 0; i < val_size; i++){
		cipher.clear();
		temp.clear();
		StringSource sss(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		//conver it to a bigint and store it somewhere.
		temp = cipher.substr (0, byte_);// truncate the ciphertext
		memset(prn_, 0x00, byte_ + 1);
		strcpy((char*)prn_, temp.c_str());
		mpz_init(res_[i]);
		mpz_import(res_[i], byte_, 1, 1, 0, 0, prn_);
		mpz_mod(res_[i], res_[i], pubmod);
		// use the bld factor to blind the y-coordinate
		mpz_add(res_[i], ptr[i], res_[i]); // blinds each y-coordinate
		mpz_mod(res_[i], res_[i], pubmod);
	}
	cipher.clear();
	free(ptr);
	values = res_;
}
//**********************************************************************
// - Function description: given an array of blinded y-coordinates, it unblindes them and returns the y-coordinates.
bigint* Polynomial::unblind_poly_(bigint* blinded_vals, int num_of_vals, byte* key, byte* iv, int key_size, int indx, int counter_indx, int byte_, bigint pubmod){

	bigint *res_;
	res_ = (mpz_t*)malloc(num_of_vals * sizeof(mpz_t));
	//-------- derive a key for indx
	string cipher;
	unsigned char *prn_;
	prn_ = new  unsigned char[byte_];
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(key, key_size, iv);
	StringSource s(to_string(indx), true,new StreamTransformationFilter(e,new StringSink(cipher)));
	unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
	memset(der_key_0, 0x00, key_size + 1);
	strcpy((char*)der_key_0, cipher.c_str());
	cipher.clear();
	//derive a key for counter_indx[indx].
	e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
	StringSource ss(to_string(counter_indx), true,new StreamTransformationFilter(e,new StringSink(cipher)));
	unsigned char derived_key_2[key_size]; // convert the ciphertext into a der_key
	memset(derived_key_2, 0x00, key_size + 1);
	strcpy((char*)derived_key_2, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key_2, key_size, iv);	// set key an iv.
	// use the derived key to generate a set of PRNs (for the bin at indx).
	for (int i = 0;i < num_of_vals; i++){
		cipher.clear();
		temp.clear();
		StringSource sss(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		//conver it to a bigint and store it somewhere.
		temp = cipher.substr (0, byte_);// truncate the ciphertext
		memset(prn_, 0x00, byte_+1);
		strcpy((char*)prn_, temp.c_str());
		mpz_init(res_[i]);
		mpz_import(res_[i], byte_, 1, sizeof(prn_[0]), 0, 0, prn_);
		mpz_mod(res_[i], res_[i], pubmod);
		// use the bld factor to blind the y-coordinate
		mpz_sub(res_[i], pubmod, res_[i]);
		mpz_add(res_[i], res_[i], blinded_vals[i]); // unblinds the biginteger representing a bloom filters.
		mpz_mod(res_[i], res_[i], pubmod);
	}
	temp.clear();
	cipher.clear();
	delete[] prn_;
	return res_;
}
//**********************************************************************
// - Function description: blinds the y-coordinates of a polynomial.
bigint* Polynomial::get_values(){

	return values;
}
//*******************************************************************
// - Function description: Evaluates a polynomial at each element in an array: x_point, given a set of coeeficients and an array of points.
// This uses Horner's technique which is much faster than the naive approach that involves exponentiations.
bigint* Polynomial::evaluate_coeffs(bigint* coeff, bigint* x_points, int coeff_size, int xpoint_size, bigint pubmoduli){

	bigint* res;
	if(coeff_size == 1){
		res = (mpz_t*)malloc(1 * sizeof(mpz_t));
		mpz_init_set(res[0],coeff[0]);
		return res;
	}
	res = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	for(int i = 0; i < xpoint_size; i++){
		int n = coeff_size - 1;
		mpz_init(res[i]);
		mpz_mod(res[i], coeff[n], pubmoduli);
		while(n > 0){
			mpz_mul(res[i], res[i],  x_points[i]);
			n -= 1;
			mpz_add(res[i], res[i], coeff[n]);
			mpz_mod(res[i], res[i], pubmoduli);
		}
		mpz_mod(res[i], res[i], pubmoduli);
	}
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
		for (i = k - 1; i >= 0; i--) {
			mpz_mul(t1, t1, aa);
			mpz_mod(t1, t1, N);
			mpz_add(t1, t1, prod[i]);
		}
		mpz_init_set_str(t2, "0", 10);
		for (i = k - 1; i >= 0; i--) {
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
				mpz_add(prod[k], t1, prod[k - 1]);
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
	while (m > 0 && (res[m - 1] == 0)) m--;
	return res;
}
//**********************************************************************
