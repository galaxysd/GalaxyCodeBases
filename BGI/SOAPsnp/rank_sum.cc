#include "soap_snp.h"

int Prob_matrix::rank_table_gen() {
	// When N <= 63, (so that n1<=31), use this table to test
	ubit64_t i, n1, N, T1;
	rate_t p_left, p_right;

	// Calculate the factorials
	double * fact = new double [64];
	fact[0]=(double)1.0;
	for(i=1;i!=64;i++) {
		fact[i] = fact[i-1]*i;
	}

	ubit64_t * rank_sum= new ubit64_t [64*64*2048]; // 6bit: N; 5bit: n1; 11bit; T1
	memset(rank_sum, 0, sizeof(ubit64_t)*64*64*2048);
	rank_sum[0]=1;
	for(N=1;N!=64;N++) {
		//cerr<<N<<endl;
		for(n1=0;n1<=N;n1++) {
			//cerr<<n1<<endl;
			for(T1=(1+n1)*n1/2;T1<=(N+N-n1+1)*n1/2;T1++) {
				//cerr<<T1<<endl;
				// Dynamic programming to generate the table
				rank_sum[N<<17|n1<<11|T1] = rank_sum[((N-1)<<17)|(n1<<11)|T1] + ((T1>=N && n1>0) ? rank_sum[((N-1)<<17)|((n1-1)<<11)|(T1-N)]:0);
				// Here, the p_rank is not cumulative
				p_rank[(N<<17)|(n1<<11)|T1] = rank_sum[N<<17|n1<<11|T1] / (fact[N]/(fact[n1]*fact[N-n1]));
				//cerr<<p_rank[(N<<17)|(n1<<11)|T1]<<endl;
			}
			p_left = 0.0, p_right =1.0;
			for(T1=(1+n1)*n1/2;T1<=(N+N-n1+1)*n1/2;T1++) {
				p_right = 1.0 - p_left;
				p_left += p_rank[(N<<17)|(n1<<11)|T1];
				p_rank[N<<17|n1<<11|T1] = (p_left<p_right?p_left:p_right);
				//std::cerr<<N<<"\t"<<n1<<"\t"<<T1<<"\t"<<p_rank[(N<<17)|(n1<<11)|T1]<<endl;
			}
		}
	}
	delete [] rank_sum;
	delete [] fact;
	return 1;
}

double Call_win::normal_test(int n1, int n2, double T1, double T2) {
	double u1, u2;
	u1 = (T1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
	u2 = (T2 - n2*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
	return normal_value(fabs(u1)>fabs(u2)?u1:u2);
}

double Call_win::table_test(double *p_rank, int n1, int n2, double T1, double T2) {
	if(n1<=n2) {
		return p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]+(T1-(int)T1)*(p_rank[(n1+n2)<<16|n1<<11|(int)(T1+1)]-p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]);
	}
	else {
		return p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]+(T2-(int)T2)*(p_rank[(n1+n2)<<16|n2<<11|(int)(T2+1)]-p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]);
	}
}

double Call_win::rank_test(Pos_info & info, char best_type, double * p_rank, Parameter * para) {
	if( (best_type&3) == ((best_type>>2)&3) ) {
		// HOM
		return 1.0;
	}
	if( info.count_uni[best_type&3]==0 || info.count_uni[(best_type>>2)&3]==0) {
		// HET with one allele...
		return 0.0;
	}
	//cerr<<"RankSum:"<<info.pos<<endl;
	//int * same_qual_count = new int [para->q_max-para->q_min+1];
	//memset(same_qual_count, 0, sizeof(int)*(para->q_max-para->q_min+1));
	//double * rank_array= new double [para->q_max-para->q_min+1];
	//memset(rank_array, 0, sizeof(double)*(para->q_max-para->q_min+1));
	int *same_qual_count = new int [64];
	double *rank_array = new double [64];
	memset(same_qual_count,0,sizeof(int)*64);
	memset(rank_array,0,sizeof(double)*64);

	int rank(0);
	double T[4]={0.0, 0.0, 0.0, 0.0};
	bool is_need[4] ={false,false,false,false};
	is_need[(best_type&3)]=true; is_need[((best_type>>2)&3)]=true;
	std::string::size_type o_base, strand;
	int  q_score, coord;
	for(o_base=0;o_base!=4;o_base++) {
		if(info.count_uni[o_base]==0 || !is_need[o_base]) continue;
		for(q_score=para->q_max-para->q_min;q_score>=0;q_score--) {
			for(coord=para->read_length-1;coord>=0;coord--) {
				for(strand=0;strand<2;strand++) {
					same_qual_count[q_score] += info.base_info[o_base<<15|strand<<14|q_score<<8|coord];
					//if(info.pos==1256 && info.base_info[o_base<<13|strand<<12|q_score<<6|coord]!=0) {
					//	cerr<<info.pos<<"\t"<<q_score<<"\t"<<same_qual_count[q_score]<<"\t"<<int(info.base_info[o_base<<13|strand<<12|q_score<<6|coord])<<endl;
					//}
				}
			}
		}
	}
	rank = 0;
	for(q_score=0;q_score<=(ubit64_t)(para->q_max-para->q_min+1);q_score++) {
		rank_array[q_score]= rank+(1+same_qual_count[q_score])/2.0;
		rank += same_qual_count[q_score];
	}
	for(o_base=0;o_base!=4;o_base++) {
		if(info.count_uni[o_base]==0 || !is_need[o_base]) continue;
		for(q_score=para->q_max-para->q_min;q_score>=0;q_score--) {
			for(coord=para->read_length-1;coord>=0;coord--) {
				for(strand=0;strand<2;strand++) {
					//cerr<<"ACTG"[o_base]<<endl;
					T[o_base] += (rank_array[q_score] * info.base_info[o_base<<15|strand<<14|q_score<<8|coord]);
					//cerr<<T[o_base]<<endl;
				}
			}
		}
	}
	delete [] same_qual_count;
	delete [] rank_array;
	//for(size_t a=0;a!=4;a++) {
	//	fprintf(stderr, "%lf\t", T[a]);
	//}
	//fprintf(stderr, "\n");
	//std::cerr<<"%d\t%d\t%lf\t%lf", info.count_uni[best_type&3],  info.count_uni[(best_type>>2)&3], T[best_type&3], T[(best_type>>2)&3]);
	if (info.count_uni[best_type&3]+info.count_uni[(best_type>>2)&3]<64) {
		return table_test(p_rank, info.count_uni[best_type&3], info.count_uni[(best_type>>2)&3], T[best_type&3], T[(best_type>>2)&3]);
	}
	else {
		return normal_test(info.count_uni[best_type&3], info.count_uni[(best_type>>2)&3],T[best_type&3], T[(best_type>>2)&3]);
	}
}
