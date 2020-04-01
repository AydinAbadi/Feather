/*

- Hash table class for a simulation of update phase in the protocol proposed in:
	https://ieeexplore.ieee.org/abstract/document/7934388/

*/
//**********************************************************************

#include"Hashtable.h"

//**********************************************************************
// - Description: Constructor - given the hash table parameters and an array of elements,
// it constructs a hash table and inserts the elements into it.

Hashtable::Hashtable(int NoElem_in_bucket, bigint* elem, int elem_size, int table_size){

	int indx[elem_size];
	// converts z to bigint zz, where zz will be used as a modulus.
	bigint zz, minus_one, *b;
	mpz_init(zz);
	mpz_set_ui(zz,table_size);
	T=(mpz_t**)malloc(table_size * sizeof(mpz_t));
	b=(mpz_t*)malloc(elem_size * sizeof(mpz_t));
	string s_val;
	CryptoPP::SHA512 hash2;
	byte  digest[CryptoPP::SHA512::DIGESTSIZE];
	for(int i = 0;i < table_size; i++){
		T[i] = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
		for(int k = 0; k < NoElem_in_bucket; k++){ // fills all bins with -1.
			mpz_init_set_str(T[i][k],"-1",10);
		}
	}

	int new_size;
	mpz_init_set_str(minus_one,"-1",10);
	for(int i = 0; i < elem_size; i++){
		s_val=mpz_get_str(NULL,10,elem[i]);
		unsigned int nDataLen = s_val.length();
		hash2.CalculateDigest(digest,( byte*)s_val.c_str(), nDataLen);
		s_val.clear();
		mpz_init(b[i]);
		mpz_import(b[i],sizeof(digest),1,sizeof(digest[0]),0,0,digest);
		mpz_mod(b[i],b[i],zz);
		indx[i] = mpz_get_ui (b[i]);
		for(int j = 0; j < NoElem_in_bucket;){
			if(mpz_cmp(T[indx[i]][j], minus_one) == 0){
					mpz_set(T[indx[i]][j], elem[i]);// stores the elements in the table' bins.
					break;
				}
				else{
					j++;
					if(j == NoElem_in_bucket){cout<<"\n===== OVERFLLOW; Value:"<<elem[i]<< ", was not inserted into HT"<<endl;
				}
			}
		}
	}
	mpz_clear(zz);
}
//**********************************************************************
// - Function description: given a bin's index, it returns the bin's content.

bigint* Hashtable::get_bucket(int index){

	return T[index];
}
//**********************************************************************
