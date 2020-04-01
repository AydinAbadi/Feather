
/*

- a test that runs both update and PSI computation in Feather protocol.

*/
//**********************************************************************

#include "Client.h"

//**********************************************************************
// - Function description: generates a set of random bigintegers,
// and ensures that the values are smaller than the public moduli and unequal to x-coordinates.
bigint* gen_randSet (int size, int max_bitsize, bigint* pubModuli, bigint* x_points, int xpoint_size){

	int counter = 0;
	Random rd;
	mpz_t *pr_val;
	pr_val=(mpz_t*)malloc(size * sizeof(mpz_t));
	unordered_map <string, int> map;
	string s_val;
	int max_bytesize = max_bitsize;
	gmp_randstate_t rand;
	bigint ran;
	rd.init_rand3(rand, ran, max_bytesize);
	bigint temp;
	mpz_init(temp);
	bool duplicated = false;
	for(int i = 0; i < size; i++){
		mpz_urandomb(temp, rand, max_bitsize);
		mpz_mod(temp, temp, pubModuli[0]);
		/*
		for(int k=0;k<counter; k++){
			if(mpz_cmp(pr_val[k],temp)==0)
			duplicated=true;
			}
			*/
		while (mpz_cmp(temp, pubModuli[0]) > 0 || duplicated == true){ // ensures the elements are smaller than the public moduli.
			mpz_init(temp);
			mpz_urandomb(temp, rand, max_bitsize);
			//extra checks-- ensures they are distinc.
			/*
			for(int k=0;k<counter; k++){
			if(mpz_cmp(pr_val[k],temp)==0){
			duplicated=true;break;
			}
			else{duplicated=false;
			}
			}
			*/
			for(int j = 0; j < xpoint_size; j++){ //checks the random element is not equal to any x_points.
				if(mpz_cmp(temp, x_points[j]) == 0){
					mpz_init(temp);
					mpz_urandomb(temp, rand, max_bitsize);
					for(int k=0;k<counter; k++){
						if(mpz_cmp(pr_val[k], temp) == 0){
							duplicated = true; break;
						}
						else{
							duplicated = false;
						}
					}
				}
			}
		}
		mpz_init_set(pr_val[i], temp);
		counter++;
		//s_val.clear();
		//s_val = mpz_get_str(NULL, 10, temp);
		//map.insert(make_pair(s_val, 1));
	}
	return pr_val;
}


int main(){

	int pub_mod_bitsize = 40;
	int max_setsize = 1024;//524288;//131072; //65536; // 1024; // 1048576;  // 4096;
	int table_length = 30;//5242 // 2621; // 30; // 41943; //122 ;
	int bucket_max_load = 100;
	int interSec_size = 1;//131072;
	if(interSec_size > max_setsize){cout<<"interSec_size > max_setsize"<<endl; return 0;}
	int xsize = 201;// Note that the number of x is determined by bucket_max_load
	if(xsize < (2 * bucket_max_load) + 1) {
		cout<<"\nxsize must be greater than 2*bucket_max_load)+1, reset it\n";
		return 0;
	}
	int number_of_experiments = 1;
	int number_of_clients = 2;
	double temp_req = 0;
	double temp_grant = 0;
	double temp_res = 0;
	double temp_intersect = 0;
	double temp_out = 0;
	double sum = 0;
	double start_out_1 = 0;
	double end_out_1 = 0;
	double sum_1 = 0;
	double diff_b = 0;
	for(int l = 0;l < number_of_experiments; l++){
		Server serv(xsize, number_of_clients, pub_mod_bitsize, max_setsize , bucket_max_load, table_length);
		Server * serv_ptr (& serv);
		int elem_bit_size = 100;
		bigint* pub_mod = serv.get_pubModuli();
		// Assigning random values to two sets a and b.
		cout<<"\n----------------------------------------------------------------\n";
		cout<<"\t Set_size:       "<<max_setsize<<endl;
		cout<<"\t Pub_mod_bitsize:  "<<pub_mod_bitsize<<endl;
		cout<<"\t Bucket_max_load:  "<<bucket_max_load<<endl;
		cout<<"\t Table_length:     "<<table_length<<endl;
		cout<<"\t Public modulous:   "<<pub_mod[0]<<endl;
		cout<<"\t Number of Clients :   "<<number_of_clients <<endl;
		cout<<"\n----------------------------------------------------------------\n";
		int t1 , t2;
		mpz_t *aa, *bb, **a;
		a = (mpz_t**)malloc((number_of_clients - 1) * sizeof(mpz_t));
		bb = (mpz_t*)malloc(max_setsize * sizeof(mpz_t));
		cout<<"\n------ Generating two random sets with distinct values and unequal to x_points"<<endl;
		bb = gen_randSet (max_setsize, elem_bit_size, serv.get_pubModuli(), serv.get_xpoints(t2), xsize);
		for(int i = 0;i < number_of_clients - 1; i++){
			a[i] = (mpz_t*)malloc(max_setsize * sizeof(mpz_t));
			a[i] = gen_randSet (max_setsize, elem_bit_size, serv.get_pubModuli(), serv.get_xpoints(t1), xsize);
		}
		bigint *x_p = serv.get_xpoints(t1);
		for(int j = 0; j < number_of_clients -1; j++){
			for(int i = 0; i < interSec_size; i++){
				mpz_set(a[j][i], bb[i]);
			}
		}
		//define the authorizers: A_i
		Client **A_;
		string A_IDs[number_of_clients - 1];
		A_ = new Client* [number_of_clients - 1];
		for(int j = 0; j < number_of_clients - 1; j++){
			A_[j] = new Client(serv_ptr, a[j], max_setsize);
			A_IDs[j]="A_"+to_string(j)+"ID";
		}
		Client B(serv_ptr, bb, max_setsize);
		string b_id = "B_ID";
		bigint label;
		cout<<"\n----------------- Client B outsourcing -----------------"<<endl;
		double start_out_b = clock();
		B.outsource_db(b_id);
		double end_out_b = clock();
		diff_b = end_out_b - start_out_b;
		cout<<"\n----------------- Clients are outsourcing -----------------"<<endl;
		free(bb);
		for(int j = 0; j < number_of_clients - 1; j++){
			cout<<"\n client "<<j<<" is outsourcing"<<endl;
			A_[j]->outsource_db(A_IDs[j]);
		}
		for(int j = 0; j < number_of_clients -1; j++){
			for(int i = 0; i < max_setsize; i++){
				mpz_clear(a[j][i]);
			}
			free(a[j]);
		}
		free (a);
		Random rd_;
		bigint* labels;
		int num_of_exper = 1;
		Random rd_1;
		bigint *temp_9;
		int size_;
		
		//-------Update----------
		bigint *temp;
		temp = (mpz_t*)malloc(1 * sizeof(mpz_t));
		temp= gen_randSet (1, elem_bit_size,serv.get_pubModuli(), serv.get_xpoints(t1), xsize);
		cout<<"\n inserting:"<<temp[0]<<endl;
		cout<<"\n++++++++++++++++++++++"<<endl;
		string ss = B.update(temp[0], "insertion", label, "B_ID");
 		string sss = B.update(temp[0], "deletion", label, "B_ID");
 		cout<<ss<<endl;
 		cout<<sss<<endl;
		//-----------Set Intersection------------
		bigint **q;
		int* sz;
		cout<<"\n---- Gennerate the Computation Request"<<endl;
		byte B_tk[AES::DEFAULT_KEYLENGTH];
		byte  B_tIV [AES::BLOCKSIZE];
		memset(B_tk, 0x00, (AES::DEFAULT_KEYLENGTH) + 1);
		memset(B_tIV, 0x00, (AES::BLOCKSIZE) + 1);
		double start_req = clock();
		CompPerm_Request* req = B.gen_compPerm_req(B_tk, B_tIV);
		double end_req = clock();
		temp_req += end_req - start_req;
		cout<<"\n---- Grant the Computation Done"<<endl;
		GrantComp_Info** ptr;
		ptr = new GrantComp_Info*[number_of_clients - 1];
		bigint ***Q;
		Q = (mpz_t***)malloc((number_of_clients - 1) * sizeof(mpz_t));
		double start_grant;
		double end_grant;
		for(int j = 0; j < number_of_clients - 1; j++){
			if(j == 0){start_grant = clock();
			}
			ptr[j] = A_[j]->grant_comp(req, Q[j], true);
			if(j == 0){end_grant = clock();
			}
		}
		for(int j = 0; j < number_of_clients - 1; j++){
			A_[j]->free_client();
			delete A_[j];
		}
		temp_grant += end_grant - start_grant;
		cout<<"\n***---- Server-side Result Computation."<<endl;
		double start_res=clock();
		Server_Result * res = serv.compute_result(ptr, B_tk, B_tIV);
		double end_res = clock();
		temp_res += end_res - start_res;
		//-----Just to free some memory----
		cout<<"\n cleanging the server"<<endl;
		serv.free_server();
		cout<<"\n---- Client-side Result Retirieval"<<endl;

		double start_intersect = clock();
		vector<string>  final_res = B.find_intersection(res, sz , Q, number_of_clients);
		double end_intersect = clock();
		temp_intersect += end_intersect - start_intersect;
		cout<<"\n\n\t======= Result ======="<<endl;
		for(int i = 0; i < final_res.size(); i++){
			cout<<"\n\nFinal_res "<<i + 1<<": "<<final_res[i]<<endl;
		}
	}
	cout<<"\n===================="<<endl;
	cout<<"\n\n\t============= Run time ==================="<<endl;
	double out = diff_b  /number_of_experiments;
	float out_time = out / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Outsourcing-- time:"<<out_time<<endl;
	double com_req = temp_req / number_of_experiments;
	float req_time = com_req / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Computation Request-- time:"<<req_time<<endl;
	double grant = temp_grant / number_of_experiments;
	float grant_time = grant / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Computation Grant-- time:"<<grant_time<<endl;
	double res_= temp_res / number_of_experiments;
	float res_time = res_ / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Server Computation-- time:"<<res_time<<endl;
	double inter = temp_intersect / number_of_experiments;
	float inter_time = inter / (double) CLOCKS_PER_SEC;
	cout<<"\n\n Find intersection-- time:"<<inter_time<<endl;
	cout<<"\n\n\t============================================"<<endl;

//-----------End of Set intersection------------

return 0;

}
//**********************************************************************
