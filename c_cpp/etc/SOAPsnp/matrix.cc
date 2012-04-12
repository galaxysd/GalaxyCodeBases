#include "soap_snp.h"

/**
 * DATE:  
 * FUNCTION:    initialize array will be used in the prob matrix
 * PARAMETER:   flag : whether need  the rank sum , if flag =1 means do the rank sum ,else don't 
 */
Prob_matrix::Prob_matrix(bool flag){
	int i;
	//update 11-11
	bool ran_sum_mode = flag;
	p_matrix = new rate_t [256*256*4*4]; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	p_prior = new rate_t [8*4*4]; // 8(ref ACTGNNNN) * diploid(4x4)
	base_freq = new rate_t [4]; // 4 base
	type_likely = new rate_t [16+1]; //The 17th element rate_t[16] will be used in comparison
	type_prob = new rate_t [16+1];
	//update 11-11
	if (ran_sum_mode == true)
	{
		p_rank = new rate_t [64*64*2048]; // 6bit: N; 5bit: n1; 11bit; T1
		p_binom = new rate_t [256*256]; // Total * case
	}
	else // if -u argument is not exit , don't need to allocate memory for the rank sum process. 
	{
		p_rank = NULL; 
		p_binom = NULL;
	}

	//Initialize  the array .
	for(i=0;i!=256*256*4*4;i++) 
	{
		p_matrix[i] = 1.0;
	}
	for(i=0;i!=8*4*4;i++) 
	{
		p_prior[i] = 1.0;
	}
	for(i=0;i!=4;i++) 
	{
		base_freq[i] = 1.0;
	}
	for(i=0;i!=16+1;i++)
	{
		type_likely[i] = 0.0; // LOG10 Scale
		type_prob[i] = 0.0; // LOG10 Scale
	}

	//update 11-11
	//if ran_sum_mode == false , p_rank and p_binom is NULL ; else initialize .
	if (ran_sum_mode == true)
	{
		for(i=0;i!=64*64*2048;i++) 
		{
			p_rank[i] = 1.0;
		}
		for(i=0;i!=256*256;i++) 
		{
			p_binom[i] = 1.0;
		}
	}
}

/**
 * DATE:  
 * FUNCTION:    delete the allocate memory   .
 * PARAMETER:   flag : whether need  the rank sum , if flag =1 means do the rank sum ,else don't 
 */
Prob_matrix::~Prob_matrix(){
	delete [] p_matrix; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	delete [] p_prior; // 8(ref ACTGNNNN) * diploid(4x4)
	delete [] base_freq; // 4 base
	delete [] type_likely; //The 17th element rate_t[16] will be used in comparison
	delete [] type_prob;
	if (ran_sum_mode == true)
	{
		delete [] p_rank; // 6bit: N; 5bit: n1; 11bit; T1
		delete [] p_binom; // Total * case;
	}
}

/**
 * DATE:  
 * FUNCTION:    soap file generate correction matrix 
 * PARAMETER:   alignment :soap file ; para ; genome :reference and dbsnp inforamtion
 * RETURN :     
 */
int Prob_matrix::matrix_gen(igzstream & alignment, Parameter * para, Genome * genome) {
	// Read Alignment files
	ubit64_t * count_matrix = new ubit64_t [256*256*4*4];
	memset(count_matrix, 0, sizeof(ubit64_t)*256*256*4*4);
	Soap_format soap;
	map<Chr_name, Chr_info*>::iterator current_chr;
	current_chr = genome->chromosomes.end();
	std::string::size_type coord;
	for(std::string line; getline(alignment, line);) {
		std::istringstream ss(line);
		if (ss >> soap)
			deal_reads(count_matrix, genome, soap, current_chr); //get count matrix 
	}
	count_qual(count_matrix, para);
	delete [] count_matrix;

	//memset(&p_matrix[((ubit64_t)para->q_max-para->q_min+1)<<12], 0, 256*256*4*4 - (para->q_max-para->q_min+1)*256*4*4);
	// Note: from now on, the first 8 bit of p_matrix is its quality score, not the FASTQ char
	return 1;
}

int Prob_matrix::matrix_read(std::fstream &mat_in, Parameter * para) 
{
	int q_char, type;
	std::string::size_type coord;
	for(std::string line; getline(mat_in, line);) 
	{
		std::istringstream s(line);
		s>>q_char>>coord;
		for(type=0;type!=16;++type)
		{
			s>>p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type];
			//cerr<<q_char<<"|"<<coord<<"|"<<p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type]<<endl;
			//exit(0);
		}
	}
	return 1;
}

int Prob_matrix::matrix_write(std::fstream &mat_out, Parameter * para) 
{
	for( char q_char = para->q_min; q_char <= para->q_max; q_char++ )
	{
		for( std::string::size_type coord=0; coord != para->read_length; coord++)
		{
			mat_out<<((ubit64_t)q_char-para->q_min)<<'\t'<<coord;
			for(char type=0;type!=16;type++) 
			{
				mat_out<<'\t'<<scientific<<showpoint<<setprecision(16)<<p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | type];
			}
			mat_out<<endl;
		}
	}
	return 1;
}

/* program by Bill */
int Prob_matrix::matrix_gen(SamCtrl &alignment, Parameter * para, Genome * genome) {
	// Read Alignment files
	Soap_format soap;
	ubit64_t * count_matrix = new ubit64_t [256*256*4*4];
	memset(count_matrix, 0, sizeof(ubit64_t)*256*256*4*4);
	map<Chr_name, Chr_info*>::iterator current_chr;
	current_chr = genome->chromosomes.end();
	int r;
	std::string line;

	while((r = alignment.readline(line)) >=0) {
		line = alignment_format(line);
		if (line == NOUSE_ALIGNMENT)
			continue;
		std::istringstream ss(line);
		if (ss >> soap)
			deal_reads(count_matrix, genome, soap, current_chr);
	}

	count_qual(count_matrix, para);
	delete [] count_matrix;

	//memset(&p_matrix[((ubit64_t)para->q_max-para->q_min+1)<<12], 0, 256*256*4*4 - (para->q_max-para->q_min+1)*256*4*4);
	// Note: from now on, the first 8 bit of p_matrix is its quality score, not the FASTQ char
	return 1;
}

/**
 * DATE: 2010-8-5
 * FUNCTION: deal with the reads, get count matrix
 * PARAMETER: count_matrix: is the point to the matrix. genome: is the point to the Genome.
				soap: the soapformat to be deal with. current_chr: the current chr.
 * RETURN: void
 */
void Prob_matrix::deal_reads(ubit64_t *count_matrix, Genome *genome, Soap_format &soap, map<Chr_name, Chr_info*>::iterator &current_chr) {
	std::string::size_type coord;
	ubit64_t ref(0);
	if(soap.get_pos() < 0) {
		return;
	}
	//cerr<<soap<<endl;
	// In the overloaded "+" above, soap.position will be substracted by 1 so that coordiates start from 0
	//current _chr is a new chromose 
	cerr << __FUNCTION__ << __LINE__ << "\tpos\t" << soap.get_pos() << "\tsoap.get_read_len()\t" << soap.get_read_len() << endl;
	if (current_chr == genome->chromosomes.end() || current_chr->first != soap.get_chr_name()) {
		current_chr = genome->chromosomes.find(soap.get_chr_name()); // get the current chromose name
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
	//the end of this soap position is out of this chromose length   
	if (soap.get_pos()+soap.get_read_len() >= current_chr->second->length()) {
		return;
	}
	//only process the soap which best hit==1 .
	if (soap.is_unique()) {
		for(coord = 0; coord != soap.get_read_len(); coord++) {
			if (soap.is_N(coord)) { 
				;             //the corrd soap is N ,do nothing
			} 
			else {
				if(! (soap.get_pos()+coord<current_chr->second->length())) {        //the position out of chromose length,it's error
					cerr<<soap<<endl;
					cerr<<"The program found the above read has exceed the reference length:\n";
					cerr<<"The read is aligned to postion: "<<soap.get_pos()<<" with read length: "<<soap.get_read_len()<<endl;
					cerr<<"Reference: "<<current_chr->first<<" FASTA Length: "<<current_chr->second->length()<<endl;
					exit(255);
				}
				ref = current_chr->second->get_bin_base(soap.get_pos()+coord); //get this position in the referece
				if ( (ref&12) !=0 ) {
					// This is an N on reference or a dbSNP which should be excluded from calibration
					;
				}
				else {
					if(soap.is_fwd()) { //"+"
						// forward strand, soap.get_qual : quality , coord :position in the reads,ref : reference genotype,  soap.get_base : genotype in the reads
						//allele type, quality score, coordinates on the read, obeserve genotype. 
		cerr << __FUNCTION__ << __LINE__ << "\tcoord\t" << coord << endl;
						count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | (coord<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
		cerr << __FUNCTION__ << __LINE__ << endl;
					}
					else {//"-"
						// reverse strand
						count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | ((soap.get_read_len()-1-coord)<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
					}
				}
			}
		}

	}
}
/**
 * DATE: 2010-8-5
 * FUNCTION: count the quality.
 * PARAMETER: count_matrix: is the point to the matrix. para: is the point to the Parameter.
 * RETURN: void
 */
void Prob_matrix::count_qual(ubit64_t *count_matrix, Parameter *para) {
	std::string::size_type coord;
	ubit64_t o_base/*o_based base*/, t_base/*theorecical(supposed) base*/, type, sum[4], same_qual_count_by_type[16], same_qual_count_by_t_base[4], same_qual_count_total, same_qual_count_mismatch;
	char q_char/*fastq quality char*/;

	const ubit64_t sta_pow=10; // minimum number to say statistically powerful

	//the q_min default is 64, q_max is 64+40
	for(q_char=para->q_min; q_char<=para->q_max ;q_char++) 
	{
		memset(same_qual_count_by_type, 0, sizeof(ubit64_t)*16);
		memset(same_qual_count_by_t_base, 0, sizeof(ubit64_t)*4);
		same_qual_count_total = 0;
		same_qual_count_mismatch = 0;
		for(coord=0; coord != para->read_length ; coord++) { //para->read_length is get from "-L"
			for(type=0;type!=16;type++) {//type is the 2 genotype from ACGT
				// If the sample is small, then we will not consider the effect of read cycle.
				//quality == q_char , coordinates on the read is the same,  
				same_qual_count_by_type[type] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				//same_qual_count_by_t_base[0]~[3],(type>>2)&3 is referce genotype
				same_qual_count_by_t_base[(type>>2)&3] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				same_qual_count_total += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				//cerr<<(int)type<<'\t'<<same_qual_count_by_type[type]<<'\t'<<same_qual_count_by_t_base[0]<<endl;
				if(type % 5 != 0) {
					// Mismatches
					same_qual_count_mismatch += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				}
			}
		}
		for(coord=0; coord != para->read_length ; coord++) 
		{
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
}
