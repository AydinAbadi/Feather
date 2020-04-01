/*

- Random class used in the Feather protocol.

*/
//**********************************************************************

#include "Rand.h"

//**********************************************************************
// - Function description: given the source of randomness and the number of bytes required,
// it writes random bytes to "buf".
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
	get_rand_file(buf, len, cg);
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
// - Function description: generates a set of distinct pseudorandom values
bigint * Random::get_nonconflict_randset(int size, bigint pub_mod, int pub_mod_size, unordered_map <string, int> map, bigint*x_points, int x_size){

	int rand_size = pub_mod_size/8;
	byte seed_[rand_size];
	AutoSeededRandomPool prng;
	bigint* res;
	res = (mpz_t*)malloc(size * sizeof(mpz_t));
	for(int i = 0; i < size; i++){
		mpz_init(res[i]);
		memset(seed_,0x00,rand_size+1);
		prng.GenerateBlock(seed_, rand_size);// seed_: master seed for PRF
		mpz_import(res[i], rand_size, 1, sizeof(seed_[0]), 0, 0, seed_);
	}
	return res;
}
//**********************************************************************
// - Function description: regenerates the most recent pseudorandom values of a bin.
bigint* Random::gen_bin_PRNs_(byte *key, byte* iv, int key_size, int indx, int counter_indx, int num_of_PRNs,int byte_){

	// regenertes the corresponding blinding factor.
	bigint *bld_factor;
	bld_factor = (mpz_t*)malloc(num_of_PRNs * sizeof(mpz_t));
	//-------- derive a BF_key for indx
	string cipher;
	unsigned char *prn_;
	prn_ = new  unsigned char[byte_];
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(key, key_size, iv);
	// encrypt indx to derive a key
	StringSource s(to_string(indx), true,new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
	memset(der_key_0, 0x00, key_size + 1);
	strcpy((char*)der_key_0, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
	StringSource ss(to_string(counter_indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char derived_key_2[key_size]; // convert the ciphertext into a der_key
	memset(derived_key_2, 0x00, key_size + 1);
	strcpy((char*)derived_key_2, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key_2, key_size, iv);	// set key an iv.
	for (int i = 0;i < num_of_PRNs; i++){
		cipher.clear();
		temp.clear();
		StringSource sss(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		//conver it to a bigint and store it somewhere.
		temp = cipher.substr (0, byte_);// truncate the ciphertext
		memset(prn_, 0x00, byte_ + 1);
		strcpy((char*)prn_, temp.c_str());
		mpz_init(bld_factor[i]);
		mpz_import(bld_factor[i], byte_, 1, sizeof(prn_[0]), 0, 0, prn_);
	}
	delete [] prn_;
	cipher.clear();
	temp.clear();
	return bld_factor;
}
//**********************************************************************
// - Function description: generates the most recent pseudorandom values for all bins of a hash table
bigint** Random::gen_HT_RRNs_(byte *key, byte* iv, int key_size, int HT_size, int* counter_indx, int num_of_PRNs, int byte_){

	bigint** res;
	res = (mpz_t**)malloc(HT_size * sizeof(mpz_t));
	for(int i = 0; i < HT_size; i++){
		res[i]= gen_bin_PRNs_(key, iv, key_size, i, counter_indx[i], num_of_PRNs, byte_);
	}
	return res;
}
//**********************************************************************
// - Function description: this is used by gen_PR_Labels to generate a label for indx-th bin.
bigint* Random::gen_PRN_(byte *key, byte* iv, int key_size, int indx, int byte_, int table_size){

	bigint *label,zero;
	mpz_init_set_str(zero, "0", 10);
	int counter = 1;
	label = (mpz_t*)malloc(1 * sizeof(mpz_t));
	string cipher;
	string temp;
	mpz_init_set_str(label[0], "0", 10);
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(key, key_size, iv);
	while(mpz_cmp(label[0], zero) == 0){
		cipher.clear();
		temp.clear();
		int indx_temp = indx * counter;
		StringSource s(to_string(indx_temp), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		mpz_init(label[0]);
		mpz_import(label[0], AES::BLOCKSIZE, 1, 1, 0, 0, cipher.c_str());
		counter=table_size * counter;
	}
	mpz_clear(zero);
	cipher.clear();
	temp.clear();
	return label;
}
//**********************************************************************
// - Function description: this is used to generate a set of pseudorandom labels.
bigint* Random::gen_PRNs_(byte *key, byte* iv, int key_size, int number_of_PRNs, int byte_, int table_size){

	bigint *labels, *temp;
	labels = (mpz_t*)malloc(number_of_PRNs * sizeof(mpz_t));
	temp = (mpz_t*)malloc(1 * sizeof(mpz_t));
	for(int i = 0;i < number_of_PRNs; i++){
		mpz_init(temp[0]);
		temp = gen_PRN_(key, iv, key_size, i, byte_, table_size);
		mpz_init_set(labels[i], temp[0]);
	}
	mpz_clear(temp[0]);
	free(temp);
	return labels;
}
//**********************************************************************
// - Function description: used to generated a set of pseudorandom masking values.
bigint** Random::gen_PRNs_forBins(int table_size, byte* key, byte* iv, int key_size, int bin_cap, int byte_, bigint pubmoduli_){

	// regenertes the corresponding blinding factor.
	bigint ** res_;
	res_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	//-------- derive a BF_key for indx
	string cipher;
	unsigned char *prn_;
	prn_ = new unsigned char[byte_];
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	// encrypt indx to derive a key
	for(int i = 0;i < table_size; i++){
		e.SetKeyWithIV(key, key_size, iv);
		cipher.clear();
		res_[i] = (mpz_t*)malloc(bin_cap * sizeof(mpz_t));
		StringSource s(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
		unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
		memset(der_key_0, 0x00, key_size + 1);
		strcpy((char*)der_key_0, cipher.c_str());
		e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
		for (int j = 0; j < bin_cap; j++){
			cipher.clear();
			temp.clear();
			StringSource sss(to_string(j), true, new StreamTransformationFilter(e, new StringSink(cipher)));
			temp = cipher.substr (0, byte_);// truncate the ciphertext
			memset(prn_, 0x00, byte_ + 1);
			strcpy((char*)prn_, temp.c_str());
			mpz_init(res_[i][j]);
			mpz_import(res_[i][j], byte_, 1, sizeof(prn_[0]), 0, 0, prn_);
		}
	}
	delete []prn_;
	cipher.clear();
	temp.clear();
	return res_;
}
//**********************************************************************
