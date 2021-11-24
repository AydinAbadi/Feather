
/*

 - A simulation of update phase in the protocol proposed in:
   https://ieeexplore.ieee.org/abstract/document/7934388/

- More specifically, a client takes the following steps:
               (1) unblinds its entire outsourced data
		(2) updates a bin.
		(3) reencodes the elements of the bin.
		(4) re-blinds the entire dataset.

- To simulate the update, a set of random bigintegers are generated and
then the above operations are carried out sequentially on the set.
*/
//*********************************************************************

#include "Polynomial.h"

//*********************************************************************
// - Function description: given an array of polynomial's coefficients, it finds and returns the polynomials roots.

bigint* findroots(bigint *coeff, int coeff_size, int& number_of_roots, bigint pubmoduli){

	int counter_roots = 0;
	bigint *res;
	res = (mpz_t*)malloc(coeff_size * sizeof(mpz_t));
	char * tmp_mod = mpz_get_str(NULL, 10, pubmoduli);
	ZZ p = to_ZZ(tmp_mod);
	ZZ_p::init(p);
	ZZ_pX P;
	ZZ one(1);
	// finds the roots.
	for(int j = 0; j < coeff_size; j++){
		char * tmp = mpz_get_str(NULL, 10, coeff[j]);
		ZZ_p dd = to_ZZ_p(conv<ZZ> (tmp));
		SetCoeff(P, j, dd);
	}
		ZZ_p a = LeadCoeff(P);
		ZZ aa = rep(a);
		if(aa > one){MakeMonic(P);}
		Vec< Pair< ZZ_pX, long > > factors;
		CanZass(factors, P);
		vec_ZZ_p root;
		for(int j = 0;j < factors.length(); j++){
				if(factors[j].a.rep.length() == 2){
					root = FindRoots(factors[j].a);
					for(int k = 0; k < root.length(); k++){
						stringstream ss;
						ss << root[k];
						string tmpm = ss.str();
						char cv[tmpm.length()];
						strcpy(cv, tmpm.c_str());
						mpz_init_set_str(res[counter_roots], cv, 10);
						counter_roots++;
					}
				}
			}
			number_of_roots = counter_roots;
			return res;
	}
//**********************************************************************
// - Function description: given an array of x-coordinates: a and array of y-coordinates: b, it interpolates
//   a polynomial and returns an array containing the polynomial's coefficients.

bigint* interpolate(int size, bigint* a, bigint* b, bigint N){
	
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
	for (k = 0; k < m; k++){
		mpz_init_set(aa, a[k]);
		mpz_init_set_str(t1, "1", 10);
		for (i = k - 1; i >= 0; i--){
			mpz_mul(t1, t1, aa);
			mpz_mod(t1, t1, N);
			mpz_add(t1, t1, prod[i]);
		}
		mpz_init_set_str(t2, "0", 10);
		for (i = k - 1; i >= 0; i--){
			mpz_mul(t2, t2, aa);
			mpz_mod(t2, t2, N);
			mpz_add(t2, t2, res[i]);
		}
		mpz_invert(t1, t1, N);
		mpz_sub(t2, b[k], t2);
		mpz_mul(t1, t1, t2);
		for (i = 0; i < k; i++){
			mpz_mul(t2, prod[i], t1);
			mpz_mod(t2, t2, N);
			mpz_add(res[i], res[i], t2);
			mpz_mod(res[i], res[i], N);
		}
		mpz_init_set(res[k], t1);
		mpz_mod(res[k], res[k], N);
		if (k < m - 1){
			if (k == 0)
				mpz_neg(prod[0], prod[0]);
				else {
					mpz_neg(t1, a[k]);
					mpz_add(prod[k], t1, prod[k - 1]);
					for (i = k - 1; i >= 1; i--) {
						mpz_mul(t2, prod[i], t1);
						mpz_mod(t2, t2, N);
						mpz_add(prod[i], t2, prod[i - 1]);
					}
					mpz_mul(prod[0], prod[0], t1);
					mpz_mod(prod[0], prod[0], N);
				}
			}
		}
		while (m > 0 && (res[m - 1] == 0)) m--;
		mpz_clear(t1);
		mpz_clear(t2);
		return res;
	}
//**********************************************************************
// - Function description: unblids all blinded y-coordinates.
//   Also, it generates new blinding factors that will be used in step (4).

bigint** unblind(bigint** elem, bigint seed_, int table_size, int xpoint_size, bigint pubmoduli, bigint**& unbl_, int pub_moduli_bitsize){
	// generates fresh seeds.
	Random rd;
	gmp_randstate_t randx;
	bigint ranx, seedx, bin_key, bin_newkey;
	rd.init_rand3(randx, ranx, 8);
	mpz_init_set(seedx, ranx);
	gmp_randinit_default(randx);
	bigint**unbl, **new_blf;
	unbl = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	new_blf = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	gmp_randstate_t rand, rand3, rand4, rand_new;
	gmp_randinit_default(rand);
	gmp_randinit_default(rand_new);
	gmp_randseed(rand, seed_);
	gmp_randseed(rand_new, seedx);
	gmp_randinit_default(rand3);
	gmp_randinit_default(rand4);
	mpz_init(bin_key);
	mpz_init(bin_newkey);
	for(int i = 0; i < table_size; i++){
		unbl[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		new_blf[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_urandomb(bin_key, rand, pub_moduli_bitsize); // the key: k_j for each bin is re-generated.
		mpz_urandomb(bin_newkey, rand_new, pub_moduli_bitsize - 1); // a new key: k'_j for each bin is generated.
		gmp_randseed(rand4, bin_key);
		gmp_randseed(rand3, bin_newkey);
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(unbl[i][j]);
			mpz_init(new_blf[i][j]);
			mpz_urandomb(unbl[i][j], rand4, pub_moduli_bitsize);
			mpz_sub(unbl[i][j], pubmoduli, unbl[i][j]);	// generates additive inverse of the blinding factor.
			mpz_add(unbl[i][j], elem[i][j], unbl[i][j]); // unblinds the blinded y-coordinate.
			mpz_mod(unbl[i][j], unbl[i][j], pubmoduli);
			mpz_urandomb(new_blf[i][j], rand3, pub_moduli_bitsize - 1); // generates new blinding factor.
			mpz_mod(new_blf[i][j], new_blf[i][j], pubmoduli);
		}
	}
	unbl_ = unbl;
	return new_blf;
}
//**********************************************************************
// - Function description: given an element, it determines its bin's index in the hash table.

int gen_binIndx(bigint elem, int table_size){
	
	bigint b, zz;
	mpz_init(zz);
	mpz_set_ui(zz, table_size);
	string s_val;
	CryptoPP::SHA512 hash2;
	 s_val = mpz_get_str(NULL, 10, elem);
	unsigned int nDataLen = s_val.length();
	byte  digest[CryptoPP::SHA512::DIGESTSIZE];
	hash2.CalculateDigest(digest, (byte*)s_val.c_str(), nDataLen);
	s_val.clear();
	mpz_init(b);
	mpz_import(b, sizeof(digest), 1, sizeof(digest[0]), 0, 0, digest);
	mpz_mod(b, b, zz);
	int j = mpz_get_ui (b);
	return j;
}
//**********************************************************************
// - Function description: blinds y-coordniates of every bin in the hash table using the blinding
//   factors: blf that have already been generated.

bigint** blind_poly(bigint** elem, bigint** blf, int table_size, int xpoint_size, bigint pubmoduli){

	bigint** blinded_poly;
	blinded_poly = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	for(int i = 0; i < table_size; i++){
		blinded_poly[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(blinded_poly[i][j]);
			mpz_add(blinded_poly[i][j], elem[i][j], blf[i][j]);
			mpz_mod(blinded_poly[i][j], blinded_poly[i][j], pubmoduli);
		}
	}
	return blinded_poly;
}
//**********************************************************************

int main(){

	bigint seed;
	bigint elm, minus_one, **blf, *pubmoduli;
	int NoElem_in_bucket = 100;
	mpz_init_set_str(elm, "10", 10);
	mpz_init_set_str(minus_one, "-1", 10);
	int table_size = 30; // table_size can be set to any of the following values : 30, 61, 122, 245, 491, 983, 2621, 5242, 10485, 20971, 41943
	int xpoint_size = 201;
	int number_of_roots = 0;
	int pub_moduli_bitsize = 100;
	string outpoly_ID = "Client_ID";
	Random rd;
	gmp_randstate_t randx;
	bigint ranx, seedx;
	rd.init_rand3(randx, seed, 8); // generates a random seed to blind y-coordinates.
	cout<<"\n=============="<<endl;
	cout<<"\n\nTable_size:"<<table_size<<endl;
	cout<<"\n\nField size:"<<pub_moduli_bitsize<<endl;
	cout<<"\n=============="<<endl;
	int number_of_experiments = 100;
	double temp = 0;
	// generates a set of random bigintegers.
	// the total number of set elements is equals the total numner of elements is a hash table.
	for(int j = 0; j < number_of_experiments; j++){
		bigint **elem;
		elem = (mpz_t**)malloc(table_size * sizeof(mpz_t));
		for(int i = 0;i < table_size; i++){
			elem[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
			elem[i] = rd.gen_randSet(xpoint_size, 32);
		}
		blf = (mpz_t**)malloc(table_size * sizeof(mpz_t));
		pubmoduli = rd.gen_randSet(1, pub_moduli_bitsize); // generates a public modulus.
		bigint * xpoints = rd.gen_randSet(xpoint_size, pub_moduli_bitsize - 2);
		mpz_nextprime(pubmoduli[0], pubmoduli[0]);
		double start1 = clock();
		bigint** unblinded = unblind(elem, seed, table_size, xpoint_size, pubmoduli[0], blf, pub_moduli_bitsize); // unblinds the set elements.
		int inx = gen_binIndx(elm, table_size); // determines the index of the hashtable bin that will be updated.
		bigint* coeff = interpolate(xpoint_size, xpoints,unblinded[inx], pubmoduli[0]); // given x and y coordinates, interpolates a polynomial.
		bigint* roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli[0]); // finds roots of the polynomial interpolated.
		// inserts a new element and put back old set elements into the bin. Also, it padds the bin if needed.
		bigint* padded_bin;
		padded_bin = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
		mpz_init_set(padded_bin[0], elm);
		for(int i = 1; i < NoElem_in_bucket; i++){
			if(i < number_of_roots + 1){
				mpz_init_set(padded_bin[i], roots[i - 1]);
			}
			else{
				mpz_init_set(padded_bin[i], minus_one);
			}
		}
	Polynomial pol(padded_bin, outpoly_ID, xpoints, NoElem_in_bucket, xpoint_size, pubmoduli[0]); // constructs a poly representing the bin elements.
	unblinded[inx] = pol.values;
	bigint **updated_set = blind_poly(unblinded, blf, table_size, xpoint_size,pubmoduli[0]); // blindes all unblided polynomials in the hash table.
	double end1 = clock();
	temp += end1 - start1;
	}
	temp= temp/number_of_experiments;
	float dif1_f = temp/(double) CLOCKS_PER_SEC;
	cout<<"\n\nClient-side update- time:"<<dif1_f<<endl;
}
//**********************************************************************
