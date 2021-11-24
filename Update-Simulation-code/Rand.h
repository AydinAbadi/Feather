/*

- Random class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/
//*********************************************************************

#include <iostream>
#include <gmp.h>
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

using namespace std;
using namespace NTL;
typedef mpz_t bigint;
typedef unsigned char byte;

//**********************************************************************

class Random{

public:
	Random(){};
	~Random(){};
	void init_rand3(gmp_randstate_t& rand, bigint ran, int bytes);
	bigint* gen_randSet(int size, int max_bitsize);

private:
	void get_rand_devurandom(char* buf, int len);
	void get_rand_file(char* buf, int len, char* file);
};
//**********************************************************************
