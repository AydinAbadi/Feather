/*

- Hash table class used in the Feather protocol.

*/


#include"Rand.h"

//*********************************************************************

class Hashtable{

public:
	Hashtable ( int elem_in_bucket, bigint* elemen, int elem_size,int table_size);
	bigint* get_bucket(int index);
	void clear();
private:
	int NoElem_in_bucket_;
	int table_size_;
	bigint **T;// hash table
};
