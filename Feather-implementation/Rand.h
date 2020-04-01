/*

- Random class used in the Feather protocol.

*/
//**********************************************************************

#include <iostream>
#include "gmp.h"
#include <unordered_map>
#include <gmpxx.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <assert.h>
#include <time.h>
#include <bitset>
#include <vector>
#include <cryptopp/pwdbased.h>
#include <cryptopp/sha.h>
#include <cryptopp/filters.h>
#include <cryptopp/hex.h>
#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/RR.h>
#include <NTL/tools.h>
#include <NTL/vec_ZZVec.h>
#include <NTL/fileio.h>
#include <NTL/FacVec.h>
#include <NTL/new.h>
#include <boost/algorithm/string.hpp>
#include "cryptopp/secblock.h"
#include "cryptopp/cryptlib.h"
#include "bloom_filter/bloom_filter.hpp"
using namespace std;
using namespace NTL;
typedef mpz_t bigint;
typedef unsigned char byte;
using namespace CryptoPP;
using CryptoPP::HexEncoder;
using CryptoPP::HexDecoder;
using CryptoPP::StringSink;
using CryptoPP::StringSource;
using CryptoPP::StreamTransformationFilter;
#include "cryptopp/ccm.h"
using CryptoPP::CTR_Mode;
#include "assert.h"
using CryptoPP::CBC_Mode;
#include "cryptopp/osrng.h"
using CryptoPP::AutoSeededRandomPool;
#include "cryptopp/cryptlib.h"
using CryptoPP::Exception;
#include "cryptopp/hex.h"

//**********************************************************************

class Random{

public:
	Random(){};
	~Random(){};
	void init_rand3(gmp_randstate_t& rand, bigint ran, int bytes);
	bigint* gen_randSet(int size, int max_bitsize);
	bigint* gen_bin_PRNs_(byte *key, byte* iv, int key_size, int indx, int counter_indx, int num_of_PRNs, int byte_);
	bigint** gen_HT_RRNs_(byte *key, byte* iv, int key_size, int HT_size, int* counter_indx, int num_of_PRNs, int byte_);
	bigint* gen_PRN_(byte *key, byte* iv, int key_size, int indx, int byte_, int table_size);
	bigint* gen_PRNs_(byte *key, byte* iv, int key_size, int number_of_PRNs, int byte_, int table_size);
	bigint * get_nonconflict_randset(int size, bigint pub_mod, int pub_mod_size, unordered_map <string, int> map, bigint *x_points, int x_size);
	bigint** gen_PRNs_forBins(int table_size, byte* key, byte* iv, int key_size, int bin_cap, int byte_, bigint pubmoduli_);

private:
	void get_rand_devurandom(char* buf, int len);
	void get_rand_file(char* buf, int len, char* file);
};
//**********************************************************************
