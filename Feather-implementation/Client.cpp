/*

- Client-side computation of the Feather protocol

*/

#include"Client.h"
//**********************************************************************
// - Description:  Constructor- It generates the following private keys: seed,
// BF_key_, label_key, shuffle_key.
// It fetches the hashtable parameters, x-points, and field description from the server.
// It sets a moduli for blinding bloom filters and sets the bloom filter parameters.
Client::Client(Server*server, bigint *elements, int el_size){

  	key_size = AES::DEFAULT_KEYLENGTH;
	Random rd_;
	labels_bit_size = 64;
	int block_size = AES::BLOCKSIZE;
	elem = elements;
	serv = server;
	elem_size = el_size;
	//keys are generated.
	AutoSeededRandomPool prng;
	prng.GenerateBlock(seed_, key_size);// seed_: master seed for PRF
	prng.GenerateBlock(iv, block_size);// iv for PRF
	prng.GenerateBlock(BF_key_, key_size);// BF_key_: a mester seed for PRF to blind the BFs
	prng.GenerateBlock(BF_iv, block_size);
	prng.GenerateBlock(label_key_, key_size);// label_key_: a mester seed for PRF to generate a label for a bin
	prng.GenerateBlock(label_iv, block_size);
	Random rd;
	gmp_randstate_t rand4;
	bigint ran;
	mpz_init(ran);
	rd.init_rand3(rand4, ran, 8);
	mpz_init_set(shuffle_key, ran);
	int size;
	get_xpoints(size); //fetches x-coordonates from the server.
	xpoint_size = size;
	get_pubModuli(); //fetches the public moduli (description of the field) from the server.
	get_NoElem_in_bucket(); // gets the bins' capacity from the server.
	get_tablesize(); // gets the hash table length from the server.
	get_pubModuli_bitsize(); // gets the bit-size of the public moduli.
	// sets the bloom filters parameters.
	bf_parameters.projected_element_count = NoElem_in_bucket;
	bf_parameters.false_positive_probability = 0.0000000000009095;
	bf_parameters.random_seed = 0xA5A5A5A5;
	if (!bf_parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
	}
	bf_parameters.compute_optimal_parameters();
	pr_moduli_bitsize = 5780; // bit size of a bloom filter when it holds upto 100 elements, and the false positive rante is 2^{-40}.
	//generates a prime number: pr_moduli, of the above size.
	bigint* pr_mod;
	pr_mod = (mpz_t*)malloc(1 * sizeof(mpz_t));
	pr_mod = rd_.gen_randSet(1, pr_moduli_bitsize);
	mpz_init_set(pr_moduli, pr_mod[0]);
	//sets the counter.
	counter = new int [table_size];
	mpz_clear(pr_mod[0]);
  	free(pr_mod);
  	memset(counter, 0, table_size+1);
  	//inserts x_points to a hashtable.
  	string s_val;
	for(int i = 0;i < xpoint_size; i++){
		s_val.clear();
		s_val = mpz_get_str(NULL, 10, xpoints[i]);
		xPoint_map.insert(make_pair(s_val, 1));
	}
	gmp_randclear(rand4);
}
//**********************************************************************
void Client::free_client(){
	
  mpz_clear(pr_moduli);
  xPoint_map.clear();
}
//**********************************************************************
//Function description: blindes an array of shuffled bloom filters.
bigint** Client::blind_shuffled_bl_(bigint** s_bl, int table_size, byte* key, byte* iv, int key_size, int num_of_PRNs_, int byte_, bigint pubmoduli_){

	// regenertes the corresponding blinding factor.
	bigint ** res_;
  	res_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	//-------- derive a BF_key for indx
	string cipher;
  	unsigned char prn_[byte_];
	//unsigned char *prn_;
	string temp;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	// encrypt indx to derive a key
  	unsigned char der_key_0[key_size]; // convert the ciphertext into a der_key
	for(int i = 0;i < table_size; i++){
		e.SetKeyWithIV(key, key_size, iv);
    		cipher.clear();
    		res_[i] = (mpz_t*)malloc(num_of_PRNs_ * sizeof(mpz_t));
		StringSource s(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
    		memset(der_key_0, 0x00, key_size + 1);
		strcpy((char*)der_key_0,cipher.c_str());
		e.SetKeyWithIV(der_key_0, key_size, iv);	// set key an iv.
		for (int j = 0;j < num_of_PRNs_; j++){
			cipher.clear();
      			temp.clear();
      			StringSource sss(to_string(j), true,new StreamTransformationFilter(e, new StringSink(cipher)));
		  	temp = cipher.substr (0, byte_);// truncate the ciphertext
			memset(prn_, 0x00, byte_ + 1);
		 	strcpy((char*)prn_, temp.c_str());
		  	mpz_init(res_[i][j]);
		  	mpz_import(res_[i][j], byte_, 1, sizeof(prn_[0]), 0, 0, prn_);
			mpz_mod(res_[i][j], res_[i][j], pubmoduli_);
      			mpz_add(res_[i][j], res_[i][j], s_bl[i][j]);
      			mpz_mod(res_[i][j], res_[i][j], pubmoduli_);
		}
	}
	temp.clear();
  	cipher.clear();
	return res_;
}
//**********************************************************************
// - Function description: updates a client's outsourced data, i.e. deletes/inserts an element from/to the data.
string Client::update(bigint elem, string updt, bigint & label, string id){
	
	char * tmp_mod = mpz_get_str(NULL, 10, pubmoduli);
	ZZ p = to_ZZ(tmp_mod);
	ZZ_p::init(p);
	ZZ_pX P;
  	bool exits = false;
	string status;
  	bigint *un_bl;
	un_bl = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	bigint *temp_label, *bin, zero, minus_one, bf, *unblinded_bf, *temp_elem, *temp_res,zz, *blinded_bf, *temp_ar;
	temp_ar = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	mpz_init_set_str(zero, "0", 10);
	mpz_init_set_str(minus_one, "-1", 10);
	int number_of_roots = 0;
	int num_of_elements_found = 0;
	bin = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	temp_elem = (mpz_t*)malloc(1 * sizeof(mpz_t));
	temp_res = (mpz_t*)malloc(1 * sizeof(mpz_t));
	mpz_init_set(temp_elem[0], elem);
	int tmpx;
	int j = gen_binIndx(elem, table_size); // generates the bin index to which the element: elem, belongs to.
	Random rd_;
	temp_label = rd_.gen_PRN_(label_key_, label_iv, key_size, j, (pub_moduli_bitsize / 8), table_size);
	mpz_init_set(label, temp_label[0]);
	bin = serv->get_client_bin(label, id, bf, tmpx, exits, "update"); // retrieves bin: o_j, from the server and the corresponding bf.
	if(!exits){
		cout<<"\nIn update: the label does not exists!"<< endl;
		return "0";
	}
	unblinded_bf = unblind_BF_(bf, j, BF_key_, BF_iv, pr_moduli);
  	Polynomial poly;
  	un_bl = poly.unblind_poly_(bin, xpoint_size, seed_, iv, key_size, j, counter[j], pub_moduli_bitsize / 8, pubmoduli);//unblind bin
 	counter[j] += 1;
	int temp_counter = 0;
	temp_res = check_vals_in_BF(temp_elem, 1, unblinded_bf[0], bf_parameters, temp_counter); //checks if the elements (to be inserted/deleted) exists in the bin's bloom filters.
	bigint* new_Bigint_BF;
	new_Bigint_BF = unblinded_bf;
	if(updt == "insertion"){ // If the update is the insertion of the element: elem.
		number_of_roots = 0;
		num_of_elements_found = 0;
		if (mpz_cmp(temp_res[0], elem) != 0){
			Polynomial pol_;
			// exteracts the set elemenets from the bin by interpolating a polynomial and finding its roots.
			bigint* coeff = pol_.interpolate(xpoint_size, xpoints, un_bl, pubmoduli); // interpolate a polynomial, given x and y coordinates.
      			bigint* roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli, P); //finds the roots.
      			bigint* valid_roots = check_vals_in_BF(roots, number_of_roots, unblinded_bf[0], bf_parameters, num_of_elements_found); // extracts the valid roots, by checking which roots are in the bin's bloom filters.
      			bloom_filter filter(bf_parameters); // builds a new bloom filters for the bin.
			string s = mpz_get_str(NULL, 10, elem);
			filter.insert(s); // inserts the element, to be inserted, into the new bloom filters.
			// inserts the set elements of the bin into the new bloom flters.
			if(num_of_elements_found > 0){
				for(int i = 0; i < num_of_elements_found; i++){
					string s = mpz_get_str(NULL, 10, valid_roots[i]);
					filter.insert(s);
				}
			}
			new_Bigint_BF = convert_BF_to_bigint(filter); // converts the filter to a bigint.
      			bigint* padded_bin;
			padded_bin = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
			mpz_init_set(padded_bin[0],elem);
			for(int i = 1; i < NoElem_in_bucket; i++){
				if(i < num_of_elements_found + 1){
					mpz_init_set(padded_bin[i], valid_roots[i - 1]);}
					else{mpz_init_set(padded_bin[i], minus_one);}
			}
			Polynomial pol(padded_bin, xpoints, NoElem_in_bucket, xpoint_size, pubmoduli, pub_moduli_bitsize, xPoint_map);
      			un_bl = pol.values; // sets the y-coordinates.
			string ss = mpz_get_str(NULL, 10, elem);
			status="\nElement "+ss+" has been inserted";
		}
		else{cout<<"\nElement "<<elem<<" already exists in the set, no insertion was required"<<endl;}
	}
	if(updt == "deletion"){ // If the update is the deletion of the element: elem.
		number_of_roots = 0;
		num_of_elements_found = 0;
		if (mpz_cmp(temp_res[0], elem) == 0){
			Polynomial pol_;
			// exteracts the set elemenets from the bin by interpolating a polynomial and finding its roots.
			bigint* coeff = pol_.interpolate(xpoint_size, xpoints, un_bl, pubmoduli);
			bigint* roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli, P);//find the roots
			bigint* valid_roots = check_vals_in_BF(roots, number_of_roots, unblinded_bf[0], bf_parameters, num_of_elements_found);
			bloom_filter filter(bf_parameters); // builds a new bloom filters.
			if(num_of_elements_found > 0){
				for(int i = 0; i < num_of_elements_found; i++){
					// inserts all set eleements of the bin into the bloom filters, except the element to be deleted: elem.
					if(mpz_cmp(valid_roots[i], elem) != 0){
						string s = mpz_get_str(NULL, 10, valid_roots[i]);
						filter.insert(s);
					}
				}
			}
			new_Bigint_BF = convert_BF_to_bigint(filter);
			// inserts the valid roots (except elem) into the bin and pads it with -1's if needed.
			bigint* padded_bin;
			padded_bin = (mpz_t*)malloc(NoElem_in_bucket * sizeof(mpz_t));
			mpz_init_set(padded_bin[0], elem);
			for(int i = 0; i < NoElem_in_bucket; i++){
				if(i < num_of_elements_found){
					if(mpz_cmp(valid_roots[i], elem) != 0){
						mpz_init_set(padded_bin[i], valid_roots[i]);
					}
			else{ 
				mpz_init_set(padded_bin[i], minus_one);
			}
				}
				else{mpz_init_set(padded_bin[i], minus_one);}
			}
			Polynomial pol(padded_bin, xpoints, NoElem_in_bucket, xpoint_size, pubmoduli, pub_moduli_bitsize, xPoint_map);
			un_bl = pol.values; // sets y-coordinates
			string ss = mpz_get_str(NULL, 10, elem);
			status="Element "+ss+" has been deleted";
		}
		else{cout<<"Element "<<elem<<" does not exist in the set, no deletion was required"<<endl;}
	}
	double str_out_6_0 = clock();
	blinded_bf = blind_BF_(new_Bigint_BF[0], j, BF_key_, BF_iv,  pr_moduli);
	Polynomial poly_;
  	bigint *bl2;
  	bl2 = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
  	bl2 = poly_.blind_poly_(un_bl, xpoint_size, seed_, iv, key_size, j, counter[j], pub_moduli_bitsize/8, pubmoduli);
  	serv->update_client_bin(bl2, label, outpoly_ID, blinded_bf[0]); // asks the server to update the client's dataset.
	return status;
}
//**********************************************************************
// - Function description: generates a permutation map that includes a set of pairs (l_i,l_j)
// where each pair includes pseudorandom labels each blong to a different client.
bigint** Client::gen_map_(int size, byte* label_key1, byte* label_iv1, byte* label_key2, byte* label_iv2, int byte_){
	
	bigint **labels, *temp;
	labels= (mpz_t**)malloc(size * sizeof(mpz_t));
	Random rd;
	for(int i = 0; i < size; i++){
		labels[i] = (mpz_t*)malloc(2 * sizeof(mpz_t));
		mpz_init(labels[i][0]);
		mpz_init(labels[i][1]);
		temp = rd.gen_PRN_(label_key1, label_iv1, key_size, i, byte_, table_size);
		mpz_init_set(labels[i][0], temp[0]);
    		mpz_init(temp[0]);
		temp = rd.gen_PRN_(label_key2, label_iv2, key_size, i, byte_, table_size);
		mpz_init_set(labels[i][1], temp[0]);
	}
	mpz_clear(temp[0]);
	labels = R_shuffle(labels, size);
	free(temp);
	return labels;
}
//**********************************************************************
// - Function description: given two arrays permuted (under two different keys), it finds their matching indices.
// It is called by find_matched_bins().
 int* Client::find_matches(int* a, int* b, int size){

	 int* res;
	 res=new int[size];
	 for(int i = 0; i < size; i++){
		 for(int j = 0; j < size; j++){
			 if(a[i] == b[j]){
				 res[i] = j;
				 break;
			 }
		 }
	 }
	 return res;
 }
 //**********************************************************************
// - Function description: finds the match between two arrays permuted under two different keys.
// In particular, for an index "i" in vecotr a permuted using key k_1, it finds index "j" in anoher
// vector permuted under key k_2 such that "i" and "j" belong to the same index before they were
// permuted. It allocates values from 1 to "size" to each vector val_1 and va_2. Then,
// it permutes them and then uses find() to find the matched indices in the permuted vectors.
int* Client::find_matched_bins(bigint k_1, bigint k_2, int size){

	int val_1[size], val_2[size], *res;
	int *permuted_1, *permuted_2;
	for( int j = 0; j < size; j++){
		val_1[j] = j + 1;
		val_2[j] = j + 1;
	}
	permuted_1 = PR_shuffle(val_1, size, k_1);
	permuted_2 = PR_shuffle(val_2, size, k_2);
	res = find_matches(permuted_1, permuted_2, size);
	delete[]permuted_1;
  	delete[]permuted_2;
	return res;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of integers, using Fisher-Yates shuffle algorithm.
int* Client::PR_shuffle(int* elem, int size, bigint seed){

  	int *result;
	result = new int [size];
	int buffer;
	bigint r, big_j;
	mpz_init(r);
	mpz_init(big_j);
	for(int i = 0; i < size; i++){
		result[i] = elem[i];
	}
	if(size == 1) {
		return result;
	}
	int indx = 0;
	// use the seed to generate a random value between [1,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand,seed);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(r, rand, pub_moduli_bitsize);// Here, pub_moduli_bitsize is an arbitrary choice. It would fine as long as it's greater than the ceiling.
    		// gen random value in in range [0,j]
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buffer = result[j];
		result[j] = result[indx];
		result[indx] = buffer;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return result;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of big-integers, using Fisher-Yates shuffle algorithm.
bigint*  Client::PR_shuffle(bigint* elem, int size, bigint seed){
	
	bigint *res,buf, r, big_j;
	mpz_init(r);
	mpz_init(buf);
	mpz_init(big_j);
	res=(mpz_t*)malloc(size * sizeof(mpz_t));
	for(int i = 0; i < size; i++){
		mpz_init(res[i]);
		mpz_set(res[i], elem[i]);
	}
	if(size == 1) return res;
	int indx = 0;
	// use the seed to generate a random value between [1,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(r, rand, pub_moduli_bitsize);// Here, pub_moduli_bitsize is an arbitrary choice. It is fine as long as it's greater than the ceiling.
    		// gen random value in in range [0,j]
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		mpz_init_set(buf, res[j]);
		mpz_init(res[j]);
		mpz_set(res[j], res[indx]);
		mpz_set(res[indx], buf);
	}
	gmp_randclear(rand);
	mpz_clear(buf);
	mpz_clear(r);
	return res;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of bins, using Fisher-Yates shuffle algorithm.
bigint**  Client::PR_shuffle_bins(bigint** bins, int size, bigint seed_){
	
	bigint **s_bins, *buf,big_j, r;//blinding factors
	s_bins = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	buf = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	for(int i = 0; i < size; i++){
		s_bins[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(s_bins[i][j]);
			mpz_set(s_bins[i][j], bins[i][j]);
		}
	}
	if(size == 1) {return  s_bins;}
	mpz_init(r);
	mpz_init(big_j);
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed_);
	int indx = 0;
	// use the seed to generate a random value between [0,j]
	for (int  j = size - 1; j > 0; j--){
		// gen random value in in range [0,j]
		mpz_urandomb(r, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buf = s_bins[j];
		s_bins[j] = s_bins[indx];
		s_bins[indx] = buf;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return s_bins;
}
//**********************************************************************
// - Function description: pseudorandomly permutes an array of polynomials, using Fisher-Yates shuffle algorithm.
Polynomial*  Client::PR_shuffle_poly(Polynomial* pol, int size, bigint seed){

	Polynomial *ply_res, buf;
	ply_res = new Polynomial[size];
	for(int i = 0; i < size; i++){
		ply_res[i] = pol[i];
	}
	if(size == 1) return ply_res;
	int indx = 0;
	// use the seed to generate a random value between [0,j]
	gmp_randstate_t rand;
	gmp_randinit_default(rand);
	gmp_randseed(rand, seed);
	bigint r, big_j;
	mpz_init(r);
	mpz_init(big_j);
	for (int j = size - 1; j > 0; j--){
		// gen random value in in range [0,j]
		mpz_urandomb(r, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(r, r, big_j);
		indx = mpz_get_ui(r);
		//exchange
		buf = ply_res[j];
		ply_res[j] = ply_res[indx];
		ply_res[indx] = buf;
	}
	gmp_randclear(rand);
	mpz_clear(r);
	mpz_clear(big_j);
	return ply_res;
}
//**********************************************************************
// - Function description: "randomly" permutes an array of label pairs, using Fisher-Yates shuffle algorithm.
bigint**  Client::R_shuffle(bigint** elem, int size){

	bigint **res;
	res=(mpz_t**)malloc(size * sizeof(mpz_t));
	for (int i = 0; i < size; i++){
		res[i]=(mpz_t*)malloc(2 * sizeof(mpz_t));
	}
	for(int i = 0; i < size; i++){
		mpz_init(res[i][0]);
		mpz_init(res[i][1]);
		mpz_set(res[i][0], elem[i][0]);
		mpz_set(res[i][1], elem[i][1]);
	}
	if(size == 1){
		return res;
	}
	int indx = 0;
	bigint buf, buf1;
	Random rd;
	gmp_randstate_t rand;
	bigint ran;
	rd.init_rand3(rand, ran, 8);
	bigint temp, big_j;
	mpz_init(temp);
	mpz_init(big_j);
	for (int  j = size - 1; j > 0; j--){
		mpz_urandomb(temp, rand, pub_moduli_bitsize);
		mpz_set_ui(big_j, j + 1);
		mpz_mod(temp, temp, big_j);
		indx = mpz_get_ui(temp);
		mpz_init_set(buf, res[j][0]);
		mpz_init_set(buf1, res[j][1]);
		mpz_init(res[j][0]);
		mpz_init(res[j][1]);
		mpz_set(res[j][0], res[indx][0]);
		mpz_set(res[j][1], res[indx][1]);
		mpz_set(res[indx][0], buf);
		mpz_set(res[indx][1], buf1);
	}
	gmp_randclear(rand);
	mpz_clear(temp);
	mpz_clear(buf);
	mpz_clear(buf1);
	return res;
}
//**********************************************************************
// - Function description: given two arrays: v_a and v_b, permuted under two distinct keys: pk_1 and pk_2, it allows one to find a mateched element pairs: (e_{a,i},e_{b.j})
// such that both elements of each pair belong to the same index before the permutation. Then, it adds up the elements in each pair.
bigint** Client::combine_permuted_bins(bigint**& v_a, bigint**& v_b, bigint**& a, int v_size, int xpoint_size_, bigint pk_1, bigint pk_2, bigint pubmoduli){

	bigint** res;
	int* ar;
	res = (mpz_t**)malloc(v_size * sizeof(mpz_t));
	ar = new int[v_size];
	ar = find_matched_bins( pk_1, pk_2, v_size);
	for(int i = 0; i < v_size; i++){
		res[ar[i]] = (mpz_t*)malloc(xpoint_size_ * sizeof(mpz_t));
		// sums the elements of v_a at position i with its matched elements in v_b (at position ar[i]).
		// then sums the result with the elements of a at position i.
		// stores the result in res at position ar[i],  which is compatible with the permutation used in v_b.
		for(int j = 0; j < xpoint_size_; j++){
			mpz_init(res[ar[i]][j]);
			mpz_add(res[ar[i]][j], v_a[i][j], v_b[ar[i]][j]);
			mpz_add(res[ar[i]][j], res[ar[i]][j], a[i][j]);
			mpz_mod(res[ar[i]][j], res[ar[i]][j], pubmoduli);
			mpz_clear(v_a[i][j]);
			mpz_clear(v_b[ar[i]][j]);
			mpz_clear(a[i][j]);
		}
	}
	free(v_a);
  	free(v_b);
  	free(a);
	delete[] ar;
	return res;
}
//**********************************************************************
// - Function description: converts a bloom filter to a biginteger value.
bigint* Client::convert_BF_to_bigint(bloom_filter filter){
	
	int size = filter.bit_table_.size();
	bigint* res;
	res = (mpz_t*)malloc(1 * sizeof(mpz_t));
	unsigned char ar[size];
	for(int i = 0; i < size; i++){
		ar[i] = filter.bit_table_[i];
	}
	mpz_init(res[0]);
	mpz_import(res[0], sizeof(ar), 1, sizeof(ar[0]), 0, 0, ar);// converts the array of bytes to a biginteger.
	return res;
}
//**********************************************************************
// - Function description: converts a biginteger representing a bloom filter to a bloom filter.
// In the case where BF's very first bis are zero, when it is converted to a biginteger they would be lost.
// So, the function ensures they are put back (using a pad), otherwise BF would not reconstructed correctly.
 bloom_filter Client::convert_bigint_to_BF(bigint a, bloom_parameters parameters){
	 
	 // converts bigint to a bitstring and pad it if needed.
	 bloom_filter filter(bf_parameters);
	 filter.clear();
	 int size = filter.bit_table_.size();
	 int offset = 0;
	 int inx = 0;
	 string s_val, pad;
	 s_val=mpz_get_str(NULL, 2, a);
	 //padds the string if needed
	 if(s_val.length() < size * 8){
		 int dif = (size * 8) - s_val.length();
		 for (int j = 0;j < dif; j++){pad += "0";}
	 }
	 s_val = pad + s_val;
	 // stores each 8-bits of the string in the hex form in each index of a new bf vector.
	 while (offset < s_val.length() / 8 + 1){
	 	string tmp_binr = s_val.substr(offset * 8, 8);
		// converts tmp_binr to hex.
		const unsigned g_unMaxBits = 8;
		bitset <g_unMaxBits> bs(tmp_binr);
		unsigned n = bs.to_ulong();
		stringstream ss;
		ss << hex << n;
		string tmp_hex = "0x" + boost::to_upper_copy(ss.str());
		// stores hex in the filter.
	 	std::istringstream str(tmp_hex);
	 	int num;
	 	str >> std::hex >> num;
	 	filter.bit_table_[inx] = num;
	 	inx++;
	 	offset++;
		 str.clear();
	 }
	 s_val.clear();
	 pad.clear();
	 return filter;
}
//**********************************************************************
// - Function description: given an element, it determines its bin's index in the hash table.
int Client::gen_binIndx(bigint elem, int table_size){
	
	bigint b, zz;
	mpz_init(zz);
	mpz_set_ui(zz, table_size);
	string s_val;
  	s_val.clear();
	CryptoPP::SHA512 hash2;
 	s_val=mpz_get_str(NULL, 10, elem);
	unsigned int nDataLen = s_val.length();
	byte digest[CryptoPP::SHA512::DIGESTSIZE];
	hash2.CalculateDigest(digest, (byte*)s_val.c_str(), nDataLen); // hashes the element.
	s_val.clear();
  	mpz_init(b);
	mpz_import(b, sizeof(digest), 1, sizeof(digest[0]), 0, 0, digest);
	mpz_mod(b, b, zz);
	int j = mpz_get_ui(b); // converts the hash value to an integer.
  	mpz_clear(b);
  	mpz_clear(zz);
  	return j;
}
//**********************************************************************
// - Function description: given a hashtable containing set elements, it assigns a bloom filter to each bin of the table.
bigint* Client::assing_BFs2HT(Hashtable HT, int NoElem_in_bucket, int table_size, bloom_parameters parameters){
	
	bloom_filter filter(bf_parameters);
	bigint minus_one;
	mpz_init_set_str(minus_one, "-1", 10);
	bigint* temp_ar, *bigint_BF, *temp_bigint;
	bigint_BF = (mpz_t*)malloc(table_size * sizeof(mpz_t));
	temp_bigint = (mpz_t*)malloc(1 * sizeof(mpz_t));
	// retrives the set elements of each bin, and inserts them to bloom filters.
	for(int i = 0; i < table_size; i++){
		temp_ar = HT.get_bucket(i);
		for(int j = 0; j < NoElem_in_bucket; j++){
			if(mpz_cmp(temp_ar[j], minus_one) > 1){
				string s = mpz_get_str(NULL, 10, temp_ar[j]); // this is done to make query and insertion compatible.
				filter.insert(s);
				s.clear();
			}
		}
		temp_bigint = convert_BF_to_bigint(filter); // converts the bloom filter to a biginteger.
		mpz_init_set(bigint_BF[i], temp_bigint[0]); // stores the biginteger in an array.
		filter.clear();
	}
	mpz_clear(temp_bigint[0]);
  	free(temp_bigint);
  	mpz_clear(minus_one);
  	return bigint_BF;
}
//**********************************************************************
// - Function description: given an array of values and a biginteger representing a Bloom filter, it returns those
// elements of the array that exist in the bloom filter. In particular, given a polynomial's roots and Bloom filter belonging to the same bin,
// it returns the roots that have been encoded in the bloom filter.
bigint* Client::check_vals_in_BF(bigint* vals, int val_size, bigint bf, bloom_parameters parameters, int& counter){

	bigint *res, minus_one, zero;
	mpz_init_set_str(minus_one, "-1", 10);
	mpz_init_set_str(zero, "0", 10);
  	string s;
	if(mpz_cmp(bf,zero) == 0){
		res = (mpz_t*)malloc(1 * sizeof(mpz_t));
		mpz_init_set(res[0], minus_one);
		counter = 0;
		return res;
	}
	int indx = 0;
	res = (mpz_t*)malloc(val_size * sizeof(mpz_t));
	bloom_filter filter(bf_parameters);
	filter = convert_bigint_to_BF(bf, bf_parameters); //convers biginteger (representation of a bloom filters) to a bloom filters.
	for (int j = 0; j < val_size; j++){
		s.clear();
		s = mpz_get_str(NULL, 10, vals[j]);
		if(filter.contains(s)){
			mpz_init_set(res[indx], vals[j]);
			indx++;
			s.clear();
		}
	}
	filter.clear();
	mpz_clear(zero);
	mpz_clear(minus_one);
	counter = indx;
	return res;
}
//**********************************************************************
// - Function description: fetches bin's capacity: d, as apart of the hashtable parameter, from the server.
void Client::get_NoElem_in_bucket(){

	NoElem_in_bucket = serv->get_NoElem_in_bucket();
}
//**********************************************************************
// - Function description: fetches an array of x-coordinates from the server.
void Client::get_xpoints(int& size){

	xpoints = serv->get_xpoints(size);
	xpoint_size = size;
}
//**********************************************************************
// - Function description: retrives the public moduli bit-size from the server.
void Client::get_pubModuli_bitsize(){

	pub_moduli_bitsize = serv->get_pubModuli_bitsize();
}
//**********************************************************************
// - Function description: fetches the public moduli from the server.
void Client::get_pubModuli(){

	bigint *ptr = (mpz_t*)malloc(1 * sizeof(mpz_t));
	ptr = serv->get_pubModuli();
	mpz_init_set(pubmoduli, ptr[0]);
}
//**********************************************************************
// - Function description: fetches the hashtable length: h, as apart of the hashtable parameter, from the server.
void Client::get_tablesize(){

	table_size = serv->get_table_size();
}
//**********************************************************************
// - Function description: prepare the set elements and sends a blinded dataset to the server.
void Client::outsource_db(string& poly_ID){

	Client_Dataset db;
	bigint minus_one, *blinded_BF, *bigint_BF;
	mpz_init_set_str(minus_one, "-1", 10);
	Hashtable HT(NoElem_in_bucket, elem, elem_size, table_size); // contructs a hash table and inserts the element into it.
	if(poly_ID == "B_ID"){ // this is done only for test-- so it can be tested on devices with small memory-- it can be commented out.
    	bigint_BF = assing_BFs2HT(HT, NoElem_in_bucket, table_size, bf_parameters); // assigns a bloom filters to each bin. It returns an array of bigint representing bloom filters.
    	blinded_BF = blind_BFs_(bigint_BF, table_size , BF_key_, BF_iv, pr_moduli);
    	db.BF = PR_shuffle(blinded_BF, table_size, shuffle_key); // permutes the array of bigintegers and stores the result in Client_Dataset that will be sent to the server.
	}
	//sets parameters to represent each bin by a polynomial
	Polynomial *poly;
	poly = new Polynomial [table_size];
	outpoly_ID = poly_ID;
	// for every index in the hash table, it contructs a polynomial (in poly is decided whether dummy values shuold be used).
	bigint* labels;
	labels = (mpz_t*)malloc(table_size * sizeof(mpz_t));
	Random rd_;
	labels = rd_.gen_PRNs_(label_key_, label_iv, key_size, table_size, ((pub_moduli_bitsize + 6) / 8), table_size);
  	db.labels = PR_shuffle(labels, table_size, shuffle_key);
  	for(int i = 0; i < table_size; i++){
		if(poly_ID == "B_ID"){
			mpz_clear(bigint_BF[i]);
			mpz_clear(blinded_BF[i]);
		}
		mpz_clear(labels[i]);
		Polynomial pol(HT.get_bucket(i), xpoints, NoElem_in_bucket, xpoint_size, pubmoduli, pub_moduli_bitsize, xPoint_map);
		poly[i] = pol;
		// assigns a seed to every index of HT.	Each seed is used to blind corresponding poly.
		poly[i].blind_poly_(seed_, iv, AES::DEFAULT_KEYLENGTH, i, counter[i], pub_moduli_bitsize / 8, pubmoduli);
	}
	db.poly = PR_shuffle_poly(poly, table_size, shuffle_key);
	db.client_ID = poly_ID;
	serv->store_poly(db);
  	delete[] poly;
  	free(labels);
  	if(poly_ID == "B_ID"){
    	free(blinded_BF);
    	free(bigint_BF);
	}
	mpz_clear(minus_one);
	HT.clear();
}
//**********************************************************************
// - Function description: generates a request for PSI computation. This request is created by the result recipient client and
// sent to the other client for its permission.
CompPerm_Request * Client::gen_compPerm_req(byte (& tmp_key_)[AES::DEFAULT_KEYLENGTH], byte (& tmp_iv_)[AES::BLOCKSIZE]){
	
	AutoSeededRandomPool prng;
  	byte temp_key[AES::DEFAULT_KEYLENGTH];
  	byte temp_iv[AES::BLOCKSIZE];
  	prng.GenerateBlock(temp_key, sizeof(temp_key));
  	prng.GenerateBlock(temp_iv, sizeof(temp_iv));
	bigint ** bl, **s_bl;
	Random rd;
	CompPerm_Request* ptr;
	ptr = new CompPerm_Request;
	bl = rd.gen_HT_RRNs_(seed_, iv, key_size, table_size, counter, xpoint_size, pub_moduli_bitsize/8);// regenerates the blinding factors.
  	s_bl = PR_shuffle_bins(bl, table_size, shuffle_key); // pseudorandomly permutes the blinding factors.
  	ptr->r = blind_shuffled_bl_(s_bl, table_size, temp_key, temp_iv, key_size, xpoint_size, pub_moduli_bitsize/8, pubmoduli);
  	memcpy(tmp_key_, temp_key, AES::DEFAULT_KEYLENGTH);
  	memcpy(tmp_iv_, temp_iv, AES::BLOCKSIZE);
	// sets the values to be sent to the other client.
	mpz_init_set(ptr->shuffle_key_, shuffle_key);
	memcpy(ptr->label_key_, label_key_, AES::DEFAULT_KEYLENGTH);
	memcpy(ptr->label_iv, label_iv, AES::BLOCKSIZE);
  	free(s_bl);
  	free(bl);
	ptr->id = outpoly_ID;
	return ptr;
}
//**********************************************************************
// - Function description: given a request for PSI computation, it grants (or does not grant) the computation.
GrantComp_Info * Client::grant_comp(CompPerm_Request* com_req, bigint **&qq, bool accept){
	
	// if the client does not grant the computation, then it returns NULL.
  	Random rd_;
  	int byte_ = (pub_moduli_bitsize) / 8;
	GrantComp_Info * ptr;
	ptr = new GrantComp_Info;
	if(!accept){
		ptr = NULL;
		return ptr;
	}
	AutoSeededRandomPool prng;
  	byte temp_key[AES::DEFAULT_KEYLENGTH];
  	prng.GenerateBlock(temp_key, sizeof(temp_key));
  	byte temp_iv[AES::BLOCKSIZE];
  	prng.GenerateBlock(temp_iv, sizeof(temp_iv));
  	string cipher;
	CBC_Mode< AES >::Encryption e;
	bigint  **a, *Sw1, *Sw2, **q, **pm, **bl, **s_bl, **v_A, **v_B, *temp_1, *temp_w1, *temp_w2;
  	Sw1 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	Sw2 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
	v_A = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	v_B = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	a = (mpz_t**)malloc(table_size * sizeof(mpz_t));
  	qq = (mpz_t**)malloc(table_size * sizeof(mpz_t));
  	pm = gen_map_ (table_size, label_key_, label_iv, com_req->label_key_, com_req->label_iv, pub_moduli_bitsize/8);
  	bl = rd_.gen_HT_RRNs_(seed_, iv, key_size, table_size, counter, xpoint_size, pub_moduli_bitsize/8);// regenerates the blinding factors.
  	s_bl = PR_shuffle_bins(bl, table_size, shuffle_key); // permuted the blinding factors.
  	unsigned char der_key_0[key_size];
  	unsigned char der_key_1[key_size];
  	unsigned char der_key_2[key_size];
  	double start_grant_300 = clock();
  	for(int i = 0; i < table_size; i++){
		e.SetKeyWithIV(temp_key, key_size, temp_iv);
		v_A[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		v_B[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
    		a[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
    		qq[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
    		cipher.clear();
    		StringSource s(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));//encrypt i
    		//set the result as a key
    		unsigned char der_key[key_size]; // convert the ciphertext into a der_key
    		memset(der_key, 0x00, key_size + 1);
    		strcpy((char*)der_key, cipher.c_str());
    		e.SetKeyWithIV(der_key, key_size, temp_iv);	// set key an iv.
		for(int j = 0; j < 3; j++){
			cipher.clear();
      			StringSource s(to_string(j), true, new StreamTransformationFilter(e, new StringSink(cipher)));
			if(j == 0){
				memset(der_key_0, 0x00, key_size+1);
				strcpy((char*)der_key_0, cipher.c_str());
			}
			else if(j == 1){
				memset(der_key_1, 0x00, key_size+1);
				strcpy((char*)der_key_1, cipher.c_str());
			}
			else if(j == 2){
				memset(der_key_2, 0x00, key_size+1);
				strcpy((char*)der_key_2, cipher.c_str());
			}
		}
		cipher.clear();
		// generates two sets: Sw1 and Sw2, each of which contains d+1 random coefficients.
		Sw1 = rd_.gen_PRNs_(der_key_1, temp_iv, key_size, NoElem_in_bucket + 1, byte_, table_size);
    		Sw2 = rd_.gen_PRNs_(der_key_2, temp_iv, key_size, NoElem_in_bucket + 1, byte_, table_size);
    		Polynomial pol_1, pol_2;
    		temp_w1 = pol_1.evaluate_coeffs(Sw1, xpoints, NoElem_in_bucket + 1, xpoint_size, pubmoduli); // evaluates the random coefficients in Sw1 at x-coordinates.
    		temp_w2 = pol_2.evaluate_coeffs(Sw2, xpoints, NoElem_in_bucket + 1, xpoint_size, pubmoduli); // evaluates the random coefficients in Sw2 at x-coordinates.
		for(int j = 0; j < xpoint_size; j++){
			mpz_init(v_A[i][j]);
			mpz_init(v_B[i][j]);
      			mpz_mul(v_A[i][j], s_bl[i][j], temp_w1[j]);
      			mpz_mod(v_A[i][j], v_A[i][j], pubmoduli);
      			mpz_mul(v_B[i][j], com_req->r[i][j], temp_w2[j]);
      			mpz_mod(v_B[i][j], v_B[i][j], pubmoduli);
      			temp_1 = rd_.gen_PRN_(der_key_0, temp_iv, key_size, j, byte_, table_size);
      			mpz_mod(temp_1[0], temp_1[0], pubmoduli);
      			mpz_init_set(a[i][j], temp_1[0]);
      			mpz_clear(temp_w1[j]);
      			mpz_clear(temp_w2[j]);
		}
	}
	cipher.clear();
	qq = combine_permuted_bins(v_A, v_B, a, table_size, xpoint_size, shuffle_key, com_req->shuffle_key_, pubmoduli); // computes q in step d.8 in the protocol.
	// sets the values sent to the server.
  	memcpy(ptr->seed, temp_key, AES::DEFAULT_KEYLENGTH);
  	memcpy(ptr->iv, temp_iv, AES::DEFAULT_KEYLENGTH);
  	ptr->pm = pm;
	ptr->id = new string[2];
	ptr->id[0] = com_req->id;// Client B's ID
	ptr->id[1] = outpoly_ID;// Client A's ID
	free(Sw1);
  	free(Sw2);
  	free(temp_w1);
  	free(temp_w2);
  	return ptr;
}
//**********************************************************************
// - Function description: given an array of polynomial's coefficients, it finds and returns the polynomials roots.
bigint* Client::findroots(bigint *coeff, int coeff_size, int& number_of_roots, bigint pubmoduli, ZZ_pX P){

	int counter_roots = 0;
	bigint *res;
	res = (mpz_t*)malloc(coeff_size * sizeof(mpz_t));
	ZZ one(1);
	for(int j = 0; j < coeff_size; j++){
		char * tmp = mpz_get_str(NULL, 10, coeff[j]);
		ZZ_p dd = to_ZZ_p(conv<ZZ> (tmp));
		SetCoeff(P, j, dd);
	}
	ZZ_p a = LeadCoeff(P);
	ZZ aa = rep(a);
	if(aa > one){
		MakeMonic(P);
	}
	Vec< Pair < ZZ_pX, long > > factors;
	CanZass(factors, P);
	vec_ZZ_p root;
	for(int j = 0; j < factors.length(); j++){
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
				tmpm.clear();
			}
		}
	}
	number_of_roots = counter_roots;
	return res;
}
//**********************************************************************
// - Function description: given the server's response, it find the intersection. In particular,
// it  unblinds the result, finds the polynomials roots and retrives the valid ones.
vector <string> Client::find_intersection(Server_Result* res, int*& size, bigint*** Q, int number_of_clients){
	
	char * tmp_mod = mpz_get_str(NULL, 10, pubmoduli);
	ZZ p = to_ZZ(tmp_mod);
	ZZ_p::init(p);
	ZZ_pX P;
  	bigint zero, minus_one;
  	mpz_init_set_str(zero, "0", 10);
  	mpz_init_set_str(minus_one, "-1", 10);
	vector <string> all_valid_roots;
	bigint *un_bl, *unbl_BFs;
	bigint *roots, *coeff;
 	bigint* valid_roots;
	un_bl = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	coeff = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	roots = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	unbl_BFs = unblind_BFs_(res->BF, table_size, BF_key_, BF_iv, pr_moduli);
	int number_of_roots = 0;
  	int num_of_elements_found = 0;
  	bigint tmp_x, temp_q;
  	mpz_init(tmp_x);
	string tempstr;
	for(int i = 0; i < table_size; i++){
		if( mpz_cmp(unbl_BFs[i],zero)!= 0){
      		// removes the blinding factors
			for(int j = 0; j < xpoint_size; j++){
				mpz_init(un_bl[j]);
				mpz_init(temp_q);
				for(int t = 0; t < (number_of_clients - 1); t++){
					mpz_add(temp_q, temp_q, Q[t][i][j]);
					mpz_mod(temp_q, temp_q, pubmoduli);
				}
				mpz_sub(un_bl[j], pubmoduli, temp_q);
				mpz_add(un_bl[j], un_bl[j], res->result[i][j]);
				mpz_mod(un_bl[j], un_bl[j], pubmoduli);
			}
			number_of_roots = 0;
			num_of_elements_found = 0;
      			Polynomial pol;
      			coeff = pol.interpolate(xpoint_size, xpoints, un_bl, pubmoduli); // interpolates a polynomial given x and y coordinates
      			roots = findroots(coeff, xpoint_size, number_of_roots, pubmoduli, P); // finds the roots of the interpolated polynomial.
      			if(number_of_roots != 0){
				valid_roots = (mpz_t*)malloc(number_of_roots * sizeof(mpz_t));
        			mpz_init_set(valid_roots[0], minus_one);
        			valid_roots = check_vals_in_BF(roots, number_of_roots, unbl_BFs[i], bf_parameters, num_of_elements_found); // extracts the valid roots.
        			if(mpz_cmp(valid_roots[0], minus_one)!= 0 && num_of_elements_found > 0){
					for(int k = 0; k < num_of_elements_found; k++){ // stores the valid roots.
						tempstr.clear();
            					tempstr = mpz_get_str(NULL, 10, valid_roots[k]);
            					mpz_clear(valid_roots[k]);
            					all_valid_roots.push_back(tempstr);
            					tempstr.clear();
					}
				}
			}
			mpz_clear(unbl_BFs[i]);
		}
	}
	free(unbl_BFs);
	free(coeff);
	return all_valid_roots;
}
//**********************************************************************
// - Function description: given an array of bigintegers representing Bloom filters, it blinds and returns the blinded ones.
// Used in outsourcing phase.
bigint* Client::blind_BFs_(bigint* bf, int bf_size, byte* BF_key, byte* BF_iv, bigint pr_moduli){
	
	//genrates table size pseudorandom values using the BF_key.
	bigint* blinded_BF;
	blinded_BF = (mpz_t*)malloc(bf_size * sizeof(mpz_t));
	bigint *blinding_fac;
	for(int i = 0; i < bf_size; i++){
		mpz_init(blinded_BF[i]);
		blinding_fac = gen_BF_PRN_(i, counter[i], BF_key, BF_iv);
    		mpz_mod(blinding_fac[0], blinding_fac[0], pr_moduli);
		mpz_add(blinded_BF[i], blinding_fac[0], bf[i]); // blinds each biginteger representing a bloom filters.
		mpz_mod(blinded_BF[i], blinded_BF[i], pr_moduli);
	}
	mpz_clear(blinding_fac[0]);
  	free(blinding_fac);
	return blinded_BF;
}
//**********************************************************************
// - Function description: unblinds a biginteger representing a Bloom filters.
bigint* Client::gen_BF_PRN_(int indx, int counter_indx, byte* BF_key, byte* BF_iv){

	int key_size =  AES::DEFAULT_KEYLENGTH;
	// regenertes the corresponding blinding factor.
  	int Num_of_AES_Invocations = 45;
	bigint  *bld_factor;
	bld_factor = (mpz_t*)malloc(1 * sizeof(mpz_t));
	//-------- derive a BF_key for indx
	string cipher;
	cipher.clear();
	CBC_Mode< AES >::Encryption e;
	e.SetKeyWithIV(BF_key, key_size, BF_iv);
  	// encrypt indx to derive a key
	StringSource s(to_string(indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char derived_key[key_size]; // convert the ciphertext into a der_key
  	memset(derived_key, 0x00, key_size + 1);
	strcpy((char*)derived_key, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key, key_size, BF_iv);	// set key an iv.
	StringSource ss(to_string(counter_indx), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	unsigned char derived_key_2[key_size]; // convert the ciphertext into a der_key
  	memset(derived_key_2, 0x00, key_size+1);
	strcpy((char*)derived_key_2, cipher.c_str());
	cipher.clear();
	e.SetKeyWithIV(derived_key_2, key_size, BF_iv);	// set key an iv.
	for (int i = 0;i < (int)Num_of_AES_Invocations; i++){
		StringSource sss(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));
	}
	char *prn_;
  	int prn_size_ = 735;
	prn_ = new char[prn_size_];
	memset(prn_, 0x00, prn_size_ + 1);
	strcpy(prn_, cipher.c_str());
	mpz_init(bld_factor[0]);
	mpz_import(bld_factor[0], prn_size_, 1, sizeof(prn_[0]), 0, 0, prn_);
  	delete[]prn_ ;
  	cipher.clear();
	return bld_factor;
}
//**********************************************************************
// - Function description: unblinds an array of bigintegers representing Bloom filters.
bigint* Client::unblind_BFs_(bigint* BF, int bf_size, byte* BF_key, byte* BF_iv, bigint pr_moduli){
	
	bigint* unblinded_BFs, *temp, *bld_factor, *shuffled_bld;
	unblinded_BFs = (mpz_t*)malloc(bf_size * sizeof(mpz_t));
  	bld_factor = (mpz_t*)malloc(bf_size * sizeof(mpz_t));
  	for(int i = 0; i < bf_size; i++){
		temp = gen_BF_PRN_(i, counter[i], BF_key, BF_iv);
    		mpz_init(bld_factor[i]);
    		mpz_mod(bld_factor[i], temp[0], pr_moduli);
	}
	//suffle blinding factors
  	shuffled_bld = PR_shuffle(bld_factor, bf_size, shuffle_key);
  	for(int j = 0; j < bf_size; j++){
		mpz_init(unblinded_BFs[j]);
    		mpz_sub(unblinded_BFs[j], pr_moduli, shuffled_bld[j]);
    		mpz_add(unblinded_BFs[j], unblinded_BFs[j], BF[j]);
    		mpz_mod(unblinded_BFs[j], unblinded_BFs[j], pr_moduli);
    		mpz_clear(shuffled_bld[j]);
    		mpz_clear(bld_factor[j]);
	}
	free(shuffled_bld);
  	free(temp);
  	free(bld_factor);
  	return unblinded_BFs;
}
//**********************************************************************
// - Function description: unblinds a biginteger representing a Bloom filters.
bigint* Client::unblind_BF_(bigint BF, int  indx, byte* BF_key, byte* BF_iv, bigint pr_moduli){

	bigint* unblinded_BF, *bld_factor;
	unblinded_BF = (mpz_t*)malloc(1 * sizeof(mpz_t));
	mpz_init(unblinded_BF[0]);
	bld_factor = gen_BF_PRN_(indx, counter[indx], BF_key, BF_iv);
  	mpz_mod(bld_factor[0], bld_factor[0], pr_moduli);
	mpz_sub(unblinded_BF[0], pr_moduli, bld_factor[0]);
	mpz_add(unblinded_BF[0], unblinded_BF[0], BF); // unblinds the biginteger representing a bloom filters.
	mpz_mod(unblinded_BF[0], unblinded_BF[0], pr_moduli);
  	mpz_clear(bld_factor[0]);
  	free(bld_factor);
	return unblinded_BF;
}
//**********************************************************************
// - Function description: blinds a biginteger representing a Bloom filter.
bigint* Client::blind_BF_(bigint bf, int indx, byte* BF_key, byte* BF_iv, bigint pr_moduli){

	bigint* blinded_BF, *bld_factor;
	blinded_BF = (mpz_t*)malloc(1 * sizeof(mpz_t));
	mpz_init(blinded_BF[0]);
	bld_factor = gen_BF_PRN_(indx, counter[indx], BF_key, BF_iv);
  	mpz_mod(bld_factor[0], bld_factor[0], pr_moduli);
	mpz_add(blinded_BF[0], bf, bld_factor[0]);
	mpz_mod(blinded_BF[0], blinded_BF[0], pr_moduli);
  	mpz_clear(bld_factor[0]);
  	free(bld_factor);
	return blinded_BF;
}
//**********************************************************************
