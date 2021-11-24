/*

- Server-side computation of the Feather protocol

*/

#include"Server.h"

//*********************************************************************
// - Description: Constructor- It generates a set of x-coordinates and public moduli. Also, it sets a hash table parameters.
Server::Server(int num_xpoints, int dbs_size, int pub_mod_bitsize, int maxSetsize, int NoEl_bucket, int tb_size){

	max_setsize = maxSetsize;
	NoElem_in_bucket = NoEl_bucket;
	pub_moduli_bitsize = pub_mod_bitsize;
	xpoint_size = num_xpoints;
	Random rd;
	xpoints = rd.gen_randSet(num_xpoints, 28); // 20 is an arbitrary choioce. The main requirement is that x-coordinates must be non-zero.
	pu_moduli = rd.gen_randSet(1, pub_moduli_bitsize);
	mpz_nextprime(pu_moduli[0], pu_moduli[0]);
	db_size = dbs_size;
	count = 0;
	db = new Client_Dataset[dbs_size];
	table_size = tb_size;
	w_A_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	w_B_ = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	temp_w_A = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	temp_w_B = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	sum_a = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	a = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	ptr = new Server_Result;
	ptr->result = (mpz_t**)malloc(table_size * sizeof(mpz_t));
	temp_BFs = (mpz_t*)malloc(table_size * sizeof(mpz_t));
	for(int i = 0;i < table_size; i++){
		w_A_[i]= (mpz_t*)malloc((xpoint_size) * sizeof(mpz_t));
		w_B_[i]= (mpz_t*)malloc((xpoint_size) * sizeof(mpz_t));
		temp_w_A[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		temp_w_B[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		a[i]= (mpz_t*)malloc((xpoint_size) * sizeof(mpz_t));
		ptr->result[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		sum_a[i] = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
		mpz_init(temp_BFs[i]);
		for(int q = 0; q < xpoint_size; q++){
			mpz_init(ptr->result[i][q]);
			mpz_init(a[i][q]);
			mpz_init(sum_a[i][q]);
			mpz_init(w_A_[i][q]);
			mpz_init(w_B_[i][q]);
			mpz_init(temp_w_A[i][q]);
			mpz_init(temp_w_B[i][q]);
		}
	}
}
//**********************************************************************
// - Function description: to free some space at the server-side
void Server::free_server(){

	for(int i = 0;i < table_size; i++){
		for(int q = 0; q < xpoint_size; q++){
			mpz_clear(a[i][q]);
			mpz_clear(sum_a[i][q]);
			mpz_clear(temp_w_A[i][q]);
			mpz_clear(temp_w_B[i][q]);
		}
		free (sum_a[i]);
		free (a[i]);
		free (w_A_[i]);
		free (w_B_[i]);
	}
	free (sum_a);
	free (a);
	free (w_A_);
	free (w_B_);
}
//**********************************************************************
// - Function description: returns x_coordinates.
bigint* Server::get_xpoints(int& size){

	size = xpoint_size;
	bigint *ptr = xpoints;
	return ptr;
}
//**********************************************************************
// - Function description: returns public moduli.
bigint*  Server::get_pubModuli(){

	bigint *ptr;
	ptr = (mpz_t*)malloc(1 * sizeof(mpz_t));
	ptr = pu_moduli;
	return ptr;
}
//**********************************************************************
// - Function description: returns the upper bound on the set size.
int Server::get_maxSetsize(){

	return max_setsize;
}
//**********************************************************************
// - Function description: returns bin's capacity: d, as part of the hashtable parameters.
int Server::get_NoElem_in_bucket(){

	return NoElem_in_bucket;
}
//**********************************************************************
// - Function description: stores a client's dataset, given the avaialbe index at the server-side database.
// It's called by store_poly().
void Server::set_db(int index, Client_Dataset &p){
	//construct a map and insert lable and index to it.
	string s_val;
	for(int i = 0;i < table_size; i++){
		s_val.clear();
		s_val = mpz_get_str(NULL,10,p.labels[i]);
		if(i == 0){
			p.label_Index_map.insert(make_pair(s_val, 1000001));//xxxxxx why do we need make_pair???
		}//1000001 is used to distinguish between the case where the label does not
		// exists in the map (in this case map returns zero), and when the label exists and the corresponding value stored in the map is zero.
		else {
			p.label_Index_map.insert(make_pair(s_val, i));
		}
	}
	db[index] = p;
}
//**********************************************************************
// - Function description: stores a client's dataset at the server.
void Server::store_poly(Client_Dataset& p){
	if(count <= db_size){
		if(p.client_ID == "B_ID")
			set_db(0, p);
		else{
			count+=1;
			set_db(count, p);
		}
	}
	else {
		cout<<"\n Error: No space to store anymore Dataset."<<endl;
		return;
	}
}
//**********************************************************************
//Function description: computes the result, given all clients messages and their outsourced data.
Server_Result * Server::compute_result (GrantComp_Info** grantComp_info, byte* tmp_key_, byte *tmp_iv_){

	bigint  *o_A, *o_B, *temp_o_B, buf_1, buf_2;
	bigint  *Sw1, *Sw2, **q, **bl, **s_bl, *temp_1, *temp_2;
	string operation= "PSI_computation";
	Random rd_;
 	int byte_ = (pub_moduli_bitsize)/8;
 	bigint**tmp_bl_1, tmp_bf_A, tmp_bf_B;
 	tmp_bl_1= rd_.gen_PRNs_forBins(table_size, tmp_key_, tmp_iv_, 16, xpoint_size, byte_, pu_moduli[0]);// we might be include this in the loop below too, so we'll have only one loop
 	int key_size = AES::DEFAULT_KEYLENGTH;
 	string temp;
 	Sw1 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
 	Sw2 = (mpz_t*)malloc((NoElem_in_bucket + 1) * sizeof(mpz_t));
 	mpz_init(buf_1);
 	mpz_init(buf_2);
	bool exists;
	bigint temp_mul;
 	mpz_init(temp_mul);
 	unsigned char der_key_0[key_size];
 	unsigned char der_key_1[key_size];
 	unsigned char der_key_2[key_size];
 	CBC_Mode< AES >::Encryption e;
 	string cipher;
 	ptr->BF=(mpz_t*)malloc(table_size * sizeof(mpz_t));
 	int indx_A,indx_B;
	for (int t = 0; t < db_size - 1; t++){//db_size : number of clients
		for(int i = 0;i < table_size; i++){
			e.SetKeyWithIV(grantComp_info[t]->seed, key_size, grantComp_info[t]->iv);
			cipher.clear();
			StringSource s(to_string(i), true, new StreamTransformationFilter(e, new StringSink(cipher)));//encrypt i
			unsigned char der_key[key_size]; // convert the ciphertext into a der_key
			memset(der_key, 0x00, key_size+1);
			strcpy((char*)der_key, cipher.c_str());
			e.SetKeyWithIV(der_key, key_size, grantComp_info[t]->iv);	// set key an iv.
			for(int j = 0; j < 3; j++){
			cipher.clear();
			StringSource s(to_string(j), true, new StreamTransformationFilter(e, new StringSink(cipher)));
			if(j == 0){
				memset(der_key_0, 0x00, key_size+1);
				strcpy((char*)der_key_0, cipher.c_str());}
				else if(j == 1){
					memset(der_key_1, 0x00, key_size+1);
					strcpy((char*)der_key_1, cipher.c_str());}
					else if(j == 2){
						memset(der_key_2, 0x00, key_size+1);
						strcpy((char*)der_key_2, cipher.c_str());
					}
				}
				Sw1 = rd_.gen_PRNs_(der_key_1, grantComp_info[t]->iv, key_size, NoElem_in_bucket + 1, byte_, table_size);
				Sw2 = rd_.gen_PRNs_(der_key_2, grantComp_info[t]->iv, key_size, NoElem_in_bucket + 1, byte_, table_size);
				Polynomial pol_1, pol_2;
				w_A_[i] = pol_1.evaluate_coeffs(Sw1, xpoints, NoElem_in_bucket + 1, xpoint_size, pu_moduli[0]); // evaluates the random coefficients in Sw1 at x-coordinates.
				w_B_[i] = pol_2.evaluate_coeffs(Sw2, xpoints, NoElem_in_bucket + 1, xpoint_size, pu_moduli[0]); // evaluates the random coefficients in Sw2 at x-coordinates.
				for(int q = 0; q < xpoint_size; q++){
					temp_1 = rd_.gen_PRN_(der_key_0, grantComp_info[t]->iv, key_size, q, byte_, table_size);
					mpz_mod(temp_1[0], temp_1[0], pu_moduli[0] );
					mpz_init_set(a[i][q],temp_1[0]);
				}
			}
			for(int i = 0; i < table_size; i++){
				o_A = get_client_bin(grantComp_info[t]->pm[i][0], grantComp_info[t]->id[1], tmp_bf_A, indx_A, exists, "PSI_computation");
				if (t==db_size-2){
					o_B = get_client_bin(grantComp_info[t]->pm[i][1], grantComp_info[t]->id[0], tmp_bf_B, indx_B, exists, "PSI_computation");
					mpz_init_set(temp_BFs[indx_B], tmp_bf_B);
					mpz_clear(tmp_bf_B);
				}
				else{
					indx_B = get_client_bin_indx(grantComp_info[t]->pm[i][1], grantComp_info[t]->id[0], exists);
				}
				// given the map,  pseudorandom values, and clients' datasets, it computes the final result.
				for(int j = 0; j < xpoint_size; j++){
					if(t==0){
						mpz_init(temp_w_A[indx_B][j]);
						mpz_init(temp_w_B[indx_B][j]);
						mpz_init(sum_a[indx_B][j]);
					}
					mpz_mul(temp_mul, w_A_[indx_A][j], o_A[j]);
					mpz_add(temp_w_A[indx_B][j], temp_w_A[indx_B][j], temp_mul);
					mpz_add(temp_w_B[indx_B][j], temp_w_B[indx_B][j], w_B_[indx_B][j]);
					mpz_add(sum_a[indx_B][j], sum_a[indx_B][j], a[indx_A][j]);
					mpz_mod(sum_a[indx_B][j], sum_a[indx_B][j], pu_moduli[0]);
					if (t==db_size-2){
						mpz_add(buf_2, o_B[j], tmp_bl_1[indx_B][j]);
						mpz_mod(buf_2, buf_2, pu_moduli[0]);
						mpz_mul(buf_2, temp_w_B[indx_B][j], buf_2);
						mpz_mod(buf_2, buf_2, pu_moduli[0]);
						mpz_add(buf_2, temp_w_A[indx_B][j], buf_2);
						mpz_add(ptr->result [indx_B][j], sum_a[indx_B][j], buf_2);
						mpz_mod(ptr->result [indx_B][j], ptr->result [indx_B][j], pu_moduli[0]);
						mpz_clear(w_B_[indx_B][j]);
						mpz_clear(w_A_[indx_A][j]);
						mpz_clear(tmp_bl_1[indx_B][j]);
					}
				}
			}
		}
		ptr->BF = temp_BFs;
		free(tmp_bl_1);
		mpz_clear(buf_1);
		mpz_clear(buf_2);
		return ptr;
	}
//**********************************************************************
// - Function description: given client's ID, it finds the client's dataset index at the server-side.
 int Server::find_db_index(string id){

	 int i;
	 string s;
	 int temp = db_size;
	 bool found = false;
	 for(i = 0; i < temp; i++){ // in the server-side database containing  different clients datasets, find the clinet's index in there.
		 if(db[i].client_ID == id){
			 found=true;
			 return i;
			 break;
			 }
	 }
	 if(!found){
		 cout<<"There is exist no poly. in server with ID:"<<id<<endl;
		 return i = 1000000;// the value only is used to indicate NULL
	 }
}
//**********************************************************************
// - Function description: given client's update query, it applies the update to the client's dataset at the server-side
// It is called by update_client_bin();
void Server::update_db(bigint* vals,bigint label, string id, bigint bbf){

	bool found;
	int index = find_db_index(id); // finds the position in which the client's dataset is stored in the array of clients' datasets.
	int j=-1;
	string s_val = mpz_get_str(NULL, 10, label);
	int val_ = db[index].label_Index_map[s_val];
	if(val_!=0){
		// finds the client's bin tagged with the label.
		if (val_ == 1000001){
			j = 0;
		}
		else{
			j = val_;
		}
		mpz_set(db[index].BF[j], bbf);
		db[index].poly[j].values = vals;
		found = true;
	}
	if(!found){cout<<"\n Update didn't take place as the label does not exist!"<<endl;}
	s_val.clear();
}
//**********************************************************************
// - Function description: allows the client to send an update query to the server.
void Server::update_client_bin(bigint* vals, bigint label, string id, bigint bbf){

	update_db(vals, label, id, bbf);
}
//**********************************************************************
// - Function description: given a client's label, it returns the corresponding bin tagged with the label.
bigint* Server::get_client_bin(bigint label, string id, bigint& bf, int& indx, bool& exists, string operation){

	int j = -1;
	bigint *res;
	res = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	int index = find_db_index(id); // finds the position in which the client's dataset is stored in the array of clients' datasets.
	string s_val = mpz_get_str(NULL, 10, label);
	int val_ = db[index].label_Index_map[s_val];
	if(val_!= 0){
		 // finds the client's bin tagged with the label.
		 if (val_ == 1000001){
			 j = 0;
		 }
		 else{
			 j = val_;
		 }
		 exists = true;
	 }
	 if(j == -1) {
		 cout<<"\n InServer-get_client_bin-- The label for client ID: "<<id <<" does not exist!"<<endl;
		 return 0;
	 }
	 else{
		 res = db[index].poly[j].get_values();
		 if(id == "B_ID" || operation == "update"){
			 mpz_init_set(bf, db[index].BF[j]);
		 }
		 indx = j;
		 return res;
	 }
 }
//**********************************************************************
int Server::get_client_bin_indx(bigint label, string id, bool& exists){

	int j = -1;
	bigint *res;
	res = (mpz_t*)malloc(xpoint_size * sizeof(mpz_t));
	//improve this--it's finding time should be independent of the number of clients.
	int index = find_db_index(id); // finds the position in which the client's dataset is stored in the array of clients' datasets.
	string s_val = mpz_get_str(NULL, 10, label);
	int val_ = db[index].label_Index_map[s_val];
	if(val_ != 0){
		 // finds the client's bin tagged with the label.
		 if (val_ == 1000001){
			 j = 0;
		 }
		 else{
			 j = val_;
		 }
		 exists = true;
	 }
	 if(j == -1) {
		 cout<<"\n InServer-get_client_bin-- The label for client ID: "<<id <<" does not exist!"<<endl;
		 return 0;
	 }
		 return j;
}
//**********************************************************************
