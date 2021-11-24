/*

- Hash table class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/

//*********************************************************************
#include"Rand.h"

//*********************************************************************

class Hashtable{

public:
	Hashtable ( int elem_in_bucket, bigint* elemen, int elem_size,int table_size);
	bigint* get_bucket(int index);
private:
	bigint **T; // hash table
};
//**********************************************************************
