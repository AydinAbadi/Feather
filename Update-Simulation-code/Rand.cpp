/*

- Random class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/
//*********************************************************************

#include "Rand.h"

//*********************************************************************
// - Function description: given the source of randomness and the number of bytes required,
//  it writes random bytes to "buf".

void Random::get_rand_file(char* buf, int len, char* file){

	FILE* fp;
	char* p;
	fp = fopen(file, "r");
	p = buf;
	while(len){
		size_t s;
		s = fread(p, 1, len, fp);
		p += s;
		len -= s;
	}
	fclose(fp);
}
//**********************************************************************
// - Function description: provides the source of randomness: "/dev/urandom", to get_rand_file().

void Random::get_rand_devurandom(char* buf, int len){

	char* cg;
	string sg = "/dev/urandom";
	cg = new char[sg.length()];
	strcpy(cg,sg.c_str());
	get_rand_file(buf, len,cg);
}
//**********************************************************************
// - Function description: generates a truly random biginteger.

void Random::init_rand3(gmp_randstate_t& rand, bigint ran, int bytes){

	char* buf;
	mpz_t s;
	buf = new char[bytes];
	get_rand_devurandom (buf, bytes);
	gmp_randinit_default(rand);
	mpz_init(s);
	mpz_init(ran);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	mpz_init_set(ran, s);
	gmp_randseed(rand, s);
	mpz_clear(s);
	free(buf);
}
//**********************************************************************
// - Function description: generates an array of random bigintegers whose bit size is max_bitsize.

bigint* Random::gen_randSet (int size, int max_bitsize){ // the 2nd argument allows us 1-to set xpoints less than public moduli 2-to set public moduli smaller than N for each clients so they do not need to compare it with N.
	Random rd;
	mpz_t* pr_val;
	pr_val = (mpz_t*)malloc(size * sizeof(mpz_t));
	int max_bytesize = max_bitsize;
	gmp_randstate_t rand;
	bigint ran;
	rd.init_rand3(rand, ran, max_bytesize);
	for(int i = 0; i < size; i++){
		mpz_init(pr_val[i]);
		mpz_urandomb(pr_val[i], rand, max_bitsize);// The last argument is in bit
	}
	return pr_val;
}
//**********************************************************************
