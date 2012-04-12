#include "soap_snp.h"
Prob_matrix::Prob_matrix(){
	int i;
	p_matrix = new rate_t [256*256*4*4]; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	p_prior = new rate_t [8*4*4]; // 8(ref ACTGNNNN) * diploid(4x4)
	base_freq = new rate_t [4]; // 4 base
	type_likely = new rate_t [16+1]; //The 17th element rate_t[16] will be used in comparison
	type_prob = new rate_t [16+1];
	p_rank = new rate_t [64*64*2048]; // 6bit: N; 5bit: n1; 11bit; T1
	p_binom = new rate_t [256*256]; // Total * case
	for(i=0;i!=256*256*4*4;i++) {
		p_matrix[i] = 1.0;
	}
	for(i=0;i!=8*4*4;i++) {
		p_prior[i] = 1.0;
	}
	for(i=0;i!=4;i++) {
		base_freq[i] = 1.0;
	}
	for(i=0;i!=16+1;i++) {
		type_likely[i] = 0.0; // LOG10 Scale
		type_prob[i] = 0.0; // LOG10 Scale
	}
	for(i=0;i!=64*64*2048;i++) {
		p_rank[i] = 1.0;
	}
	for(i=0;i!=256*256;i++) {
		p_binom[i] = 1.0;
	}
}

Prob_matrix::~Prob_matrix(){
	delete [] p_matrix; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	delete [] p_prior; // 8(ref ACTGNNNN) * diploid(4x4)
	delete [] base_freq; // 4 base
	delete [] type_likely; //The 17th element rate_t[16] will be used in comparison
	delete [] type_prob;
	delete [] p_rank; // 6bit: N; 5bit: n1; 11bit; T1
	delete [] p_binom; // Total * case;
}

int Prob_matrix::matrix_gen(std::ifstream & alignment, Parameter * para, Genome * genome) {
	// Read Alignment files
	Soap_format soap;
	ubit64_t * count_matrix = new ubit64_t [256*256*4*4];
	map<Chr_name, Chr_info*>::iterator current_chr;
	current_chr = genome->chromosomes.end();
	ubit64_t ref(0);
	std::string::size_type coord;
	for(std::string line; getline(alignment, line);) {
		std::istringstream s(line);
		if( s>> soap ) {
			if(soap.get_pos() < 0) {
				continue;
			}
			//cerr<<soap<<endl;
			// In the overloaded "+" above, soap.position will be substracted by 1 so that coordiates start from 0
			if (current_chr == genome->chromosomes.end() || current_chr->first != soap.get_chr_name()) {
				current_chr = genome->chromosomes.find(soap.get_chr_name());
				if(current_chr == genome->chromosomes.end()) {
					for(map<Chr_name, Chr_info*>::iterator test = genome->chromosomes.begin();test != genome->chromosomes.end();test++) {
						cerr<<'!'<<(test->first)<<'!'<<endl;
					}
					cerr<<"Assertion Failed: Chromosome: !"<<soap.get_chr_name()<<"! NOT found"<<endl;
					exit(255);
				}
			}
			else {
				;
			}
			if (soap.get_pos()+soap.get_read_len()>=current_chr->second->length()) {
				continue;
			}
			if (soap.is_unique()) {
				for(coord = 0; coord != soap.get_read_len(); coord++) {
					if (soap.is_N(coord)) {
						;
					}
					else {
						if(! (soap.get_pos()+coord<current_chr->second->length())) {
							cerr<<soap<<endl;
							cerr<<"The program found the above read has exceed the reference length:\n";
							cerr<<"The read is aligned to postion: "<<soap.get_pos()<<" with read length: "<<soap.get_read_len()<<endl;
							cerr<<"Reference: "<<current_chr->first<<" FASTA Length: "<<current_chr->second->length()<<endl;
							exit(255);
						}
						ref = current_chr->second->get_bin_base(soap.get_pos()+coord);
						if ( (ref&12) !=0 ) {
							// This is an N on reference or a dbSNP which should be excluded from calibration
							;
						}
						else {
							if(soap.is_fwd()) {
								// forward strand
								count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | (coord<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
							}
							else {
								// reverse strand
								count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | ((soap.get_read_len()-1-coord)<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
							}
						}
					}
				}
			}
		}
	}
	ubit64_t o_base/*o_based base*/, t_base/*theorecical(supposed) base*/, type, sum[4], same_qual_count_by_type[16], same_qual_count_by_t_base[4], same_qual_count_total, same_qual_count_mismatch;
	char q_char/*fastq quality char*/;

	const ubit64_t sta_pow=10; // minimum number to say statistically powerful
	for(q_char=para->q_min; q_char<=para->q_max ;q_char++) {
		memset(same_qual_count_by_type, 0, sizeof(ubit64_t)*16);
		memset(same_qual_count_by_t_base, 0, sizeof(ubit64_t)*4);
		same_qual_count_total = 0;
		same_qual_count_mismatch = 0;
		for(coord=0; coord != para->read_length ; coord++) {
			for(type=0;type!=16;type++) {
				// If the sample is small, then we will not consider the effect of read cycle.
				same_qual_count_by_type[type] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				same_qual_count_by_t_base[(type>>2)&3] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				same_qual_count_total += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				//cerr<<(int)type<<'\t'<<same_qual_count_by_type[type]<<'\t'<<same_qual_count_by_t_base[0]<<endl;
				if(type % 5 != 0) {
					// Mismatches
					same_qual_count_mismatch += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				}
			}
		}
		for(coord=0; coord != para->read_length ; coord++) {
			//cerr<<(q_char)<<'\t'<<coord;
			memset(sum, (ubit64_t)0, sizeof(ubit64_t)*4);
			// Count of all ref base at certain coord and quality
			for(type=0;type!=16;type++) {
				sum[(type>>2)&3] += count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | type]; // (type>>2)&3: the ref base
				//cerr<<sum[type&3]<<endl;
			}
			for(t_base=0; t_base!=4; t_base++) {
				for(o_base=0; o_base!=4; o_base++) {
					//cerr<<'\t'<<count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base];
					if (count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base] > sta_pow) {
						// Statistically powerful
						p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base]) / sum[t_base];
					}
					else if (same_qual_count_by_type[t_base<<2|o_base] > sta_pow) {
						// Smaller sample, given up effect from read cycle
						p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] =  ((double)same_qual_count_by_type[t_base<<2|o_base]) / same_qual_count_by_t_base[t_base];
					}
					else if (same_qual_count_total > 0){
						// Too small sample, given up effect of mismatch types
						if (o_base == t_base) {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)(same_qual_count_total-same_qual_count_mismatch))/same_qual_count_total;
						}
						else {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)same_qual_count_mismatch)/same_qual_count_total;
						}
					}
					else {
						;
					}

					// For these cases like:
					// Ref: G o_base: G x10 Ax5. When calculate the probability of this allele to be A,
					// If there's no A in reference gives observation of G, then the probability will be zero,
					// And therefore exclude the possibility of this pos to have an A
					// These cases should be avoid when the dataset is large enough
					// If no base with certain quality is o_based, it also doesn't matter
					if( (p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]==0) || p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] ==1) {
						if (o_base == t_base) {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = (1-pow(10, -((q_char-para->q_min)/10.0)));
							if(p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]<0.25) {
								p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = 0.25;
							}
						}
						else {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = (pow(10, -((q_char-para->q_min)/10.0))/3);
							if(p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]>0.25) {
								p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = 0.25;
							}
						}
					}
				}
			}
		}
		//cerr<<'\t'<<same_qual_count_by_type[0]<<'\t'<<same_qual_count_by_t_base[0]<<endl;
	}
	delete [] count_matrix;

	//memset(&p_matrix[((ubit64_t)para->q_max-para->q_min+1)<<12], 0, 256*256*4*4 - (para->q_max-para->q_min+1)*256*4*4);
	// Note: from now on, the first 8 bit of p_matrix is its quality score, not the FASTQ char
	return 1;
}

int Prob_matrix::matrix_read(std::fstream &mat_in, Parameter * para) {
	int q_char, type;
	std::string::size_type coord;
	for(std::string line; getline(mat_in, line);) {
		std::istringstream s(line);
		s>>q_char>>coord;
		for(type=0;type!=16;type++) {
			s>>p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type];
			//cerr<<q_char<<"|"<<coord<<"|"<<p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type]<<endl;
			//exit(0);
		}
	}
	return 1;
}

int Prob_matrix::matrix_write(std::fstream &mat_out, Parameter * para) {
	for( char q_char = para->q_min; q_char <= para->q_max; q_char++ ) {
		for( std::string::size_type coord=0; coord != para->read_length; coord++) {
			mat_out<<((ubit64_t)q_char-para->q_min)<<'\t'<<coord;
			for(char type=0;type!=16;type++) {
				mat_out<<'\t'<<scientific<<showpoint<<setprecision(16)<<p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | type];
			}
			mat_out<<endl;
		}
	}
	return 1;
}
