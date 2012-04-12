#include "soap_snp.h"
#include <cassert>

int Call_win::initialize(ubit64_t start) {
	std::string::size_type i;
	for(i=0; i != read_len + win_size ; i++) {
		sites[i].pos = i+start;
	}
	return 1;
}
//update 2010-12-20 by guyue
int Call_win::deep_init(ubit64_t start) {
	memset(sites, 0, pos_size * (read_len + win_size));
	int i;
	for(i=0; i != read_len + win_size ; i++) {
		sites[i].pos = i+start;
		sites[i].ori = 0xFF;
		//sites[i].depth = 0;
		//sites[i].repeat_time = 0;
		//sites[i].dep_uni = 0;
		//memset(sites[i].count_uni,0,sizeof(int)*4);
		//memset(sites[i].q_sum,0,sizeof(int)*4);
		//// add by Bill
		//memset(sites[i].count_sfs,0,sizeof(int)*4);
		//memset(sites[i].base_info,0,sizeof(small_int)*4*2*64*256);
		//memset(sites[i].count_all,0,sizeof(int)*4);
	}
	return 1;
}
//update 2010-12-20 by guyue
int Call_win::recycle() 
{
	std::string::size_type i;
	//for(i=0; i != read_len ; i++) 
	//{
	//	sites[i].pos = sites[i+win_size].pos;
	//	sites[i].ori = sites[i+win_size].ori;
	//	sites[i].depth = sites[i+win_size].depth;
	//	sites[i].repeat_time = sites[i+win_size].repeat_time;
	//	sites[i].dep_uni = sites[i+win_size].dep_uni;
	//	memcpy(sites[i].base_info, sites[i+win_size].base_info, sizeof(small_int)*4*2*64*256); // 4 types of bases, 2 strands, max quality score is 64, and max read length 256
	//	memcpy(sites[i].count_uni, sites[i+win_size].count_uni, sizeof(int)*4);
	//	memcpy(sites[i].q_sum, sites[i+win_size].q_sum, sizeof(int)*4);
	//	memcpy(sites[i].count_all, sites[i+win_size].count_all, sizeof(int)*4);
	//	// add by Bill
	//	memcpy(sites[i].count_sfs, sites[i+win_size].count_sfs, sizeof(int)*4);
	//}
	 
	memcpy(&sites[0], &sites[win_size], pos_size*read_len);
	memset(&sites[read_len], 0, pos_size * (win_size));
	for(i=read_len; i != read_len+win_size; i++) 
	{
		//sites[i].ori = 0xFF;
		//sites[i].pos = sites[i-1].pos+1;
		//sites[i].depth = 0;
		//sites[i].repeat_time = 0;
		//sites[i].dep_uni = 0;
		//memset(sites[i].count_uni,0,sizeof(int)*4);
		//memset(sites[i].q_sum,0,sizeof(int)*4);
		//memset(sites[i].base_info,0,sizeof(small_int)*4*2*64*256);
		//memset(sites[i].count_all,0,sizeof(int)*4);
		//// add by Bill
		//memset(sites[i].count_sfs,0,sizeof(int)*4);
		sites[i].ori = 0xFF;
		sites[i].pos = sites[i-1].pos+1;
	}
	return 1;
}


int Call_win::call_cns(Chr_name call_name, Chr_info* call_chr, ubit64_t call_length, Prob_matrix * mat, Parameter * para, gzoutstream * consensus, my_ofstream & baseinfo, SfsMethod &sfsMethod, const int id) 
{
	std::string::size_type coord;
	small_int k;
	ubit64_t o_base, strand;
	char allele1, allele2, genotype, type, type1/*best genotype*/, type2/*suboptimal genotype*/, base1, base2, base3;
	int i, q_score, q_adjusted, qual1, qual2, qual3, q_cns, all_count1, all_count2, all_count3;
	//int *pcr_dep_count;
	//pcr_dep_count = new int [para->read_length*2];
	double  rank_sum_test_value, binomial_test_value;
	bool is_out;
//	double * real_p_prior = new double [16];
//	double * likelihoods = new double [10];
//	memset(likelihoods, 0, sizeof(double)*10);
	//int n_pcr_dep_count;//update by zhukai on 2010-12-27

	//std::cerr<<"Call length="<<call_length<<endl;
	for(std::string::size_type j=0; j != call_length; ++j)  //modify the j++ to ++j    update by zhukai on 2010-12-21
	{
		//cerr<<sites[j].pos<<endl;	
		if(para->region_only) 
		{
			if (NULL== call_chr->get_region() ) 
			{
				break;
			}
			else if (! call_chr->is_in_region(sites[j].pos)) 
			{
				continue;
		    }
		    else 
			{
				;
			}
		}
		sites[j].ori = (call_chr->get_bin_base(sites[j].pos))&0xF;
		if( ((sites[j].ori&4) !=0)/*an N*/ && sites[j].depth == 0) 
		{			// CNS text format:
			// ChrID\tPos\tRef\tCns\tQual\tBase1\tAvgQ1\tCountUni1\tCountAll1\tBase2\tAvgQ2\tCountUni2\tCountAll2\tDepth\tRank_sum\tCopyNum\tSNPstauts\n"
			if(!para->glf_format && ! para->is_snp_only)
			{
				// the output file was not opened,
				if (!consensus->is_open())
				{
					continue;
				}
				(*consensus)<<call_name<<'\t'<<(sites[j].pos+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.00000\t255.00000\t0"<<"\n";
			}
			else if (para->glf_format) 
			{
				memset(likelihoods, 0, sizeof(double)*10);
				consensus->write(reinterpret_cast<char*>(likelihoods), sizeof(double)*10);
				consensus->flush();
				if(!consensus->good()) 
				{
					cerr<<"Broken ofstream after writting Position "<<(sites[j].pos+1)<<" at "<<call_name<<endl;
					exit(255);
				}
				baseinfo<<call_name<<'\t'<<(sites[j].pos+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.00000\t255.00000\t0"<<"\n";
			}
			else
			{
				;
			}
			continue;
		}
		base1 = 0, base2 = 0, base3 = 0;
		qual1 = -1, qual2= -2, qual3 = -3;
		all_count1 = 0, all_count2 =0, all_count3 =0;
		if(sites[j].dep_uni) 
		{
			// This position is uniquely covered by at least one nucleotide
			for(i=0;i!=4; ++i) 
			{
				// i is four kind of alleles
				if(sites[j].q_sum[i] >= qual1) 
				{
					base3 = base2;
					qual3 = qual2;
					base2 = base1;
					qual2 = qual1;
					base1 = i;
					qual1 = sites[j].q_sum[i];
				}
				else if (sites[j].q_sum[i]>=qual2) 
				{
					base3 = base2;
					qual3 = qual2;
					base2 = i;
					qual2  = sites[j].q_sum[i];
				}
				else if (sites[j].q_sum[i]>=qual3) 
				{
					base3 = i;
					qual3  = sites[j].q_sum[i];
				}
				else 
				{
					;
				}
			}
			if(qual1 == 0) 
			{
				// Adjust the best base so that things won't look ugly if the pos is not covered
				base1 = (sites[j].ori&7);
			}
			else if(qual2 ==0 && base1 != (sites[j].ori&7)) 
			{
				base2 = (sites[j].ori&7);
			}
			else 
			{
				;
			}
		}
		else 
		{
			// This position is covered by all repeats
			for(i=0;i!=4; ++i)
			{
				if(sites[j].count_all[i] >= all_count1) 
				{
					base3 = base2;
					all_count3 = all_count2;
					base2 = base1;
					all_count2 = all_count1;
					base1 = i;
					all_count1 = sites[j].count_all[i];
				}
				else if (sites[j].count_all[i]>=all_count2) 
				{
					base3 = base2;
					all_count3 = all_count2;
					base2 = i;
					all_count2  = sites[j].count_all[i];
				}
				else if (sites[j].count_all[i]>=all_count3)
				{
					base3 = i;
					all_count3  = sites[j].count_all[i];
				}
				else 
				{
					;
				}
			}
			if(all_count1 ==0) 
			{
				base1 = (sites[j].ori&7);
			}
			else if(all_count2 ==0 && base1 != (sites[j].ori&7))
			{
				base2 = (sites[j].ori&7);
			}
			else 
			{
				;
			}
		}
		// Calculate likelihood
		for(genotype=0;genotype!=16; ++genotype)
		{
			mat->type_likely[genotype] = 0.0;
		}
		for(o_base=0;o_base!=4; ++o_base) 
		{
			//cerr<<o_base<<endl;
			if(sites[j].count_uni[o_base]==0) 
			{
				continue;
			}
			global_dep_count = -1;
			memset(pcr_dep_count, 0, sizeof(int)*2*para->read_length);
		
			for(q_score=para->q_max-para->q_min; q_score != -1; --q_score) 
			{
				for(coord=0; coord != para->read_length; ++coord)
				{
					for(strand=0; strand!=2; ++strand) 
					{
						base_i = sites[j].base_info[o_base<<15|strand<<14|q_score<<8|coord]; //update by zhukai on 2010-12-28
						for(k=0; k != base_i; ++k)
						{							
							if(pcr_dep_count[strand*para->read_length+coord]==0)
							{
								global_dep_count += 1;
							}
							pcr_dep_count[strand*para->read_length+coord] += 1;
							q_adjusted = int( pow(10, (log10(q_score) +(pcr_dep_count[strand*para->read_length+coord]-1)*para->pcr_dependency +global_dep_count*para->global_dependency)) +0.5);
							if(q_adjusted < 1)
							{
								q_adjusted = 1;
							}
							for(allele1 = 0;allele1 != 4; ++allele1)
							{
								for(allele2 = allele1; allele2 != 4; ++allele2)
								{
									mat->type_likely[allele1<<2|allele2] += log10(0.5*mat->p_matrix[ ((ubit64_t)q_adjusted<<12) | (coord <<4) | (allele1<<2) | o_base] +0.5*mat->p_matrix[((ubit64_t)q_adjusted<<12) | (coord <<4) | (allele2<<2) | o_base]);
								}
							}
						}
					}
				}
			}
		}

		/* it is used to sfs, add by Bill*/
		if (para->sfs != NULL)
		{
			sfsMethod.getMapData(sites[j], mat, call_name, id);
		}

		// the output file was not opened,
		if (!consensus->is_open())
		{
			continue;
		}

		if(para->glf_format) 
		{
			for(int i=0; i!=10; ++i) 
			{
				consensus->write(reinterpret_cast<char*>(&mat->type_likely[glf_type_code[i]]), sizeof(double));
			}
			consensus->flush();
			if(!consensus->good())
			{
				cerr<<"Broken ofstream after writting Position "<<(sites[j].pos+1)<<" at "<<call_name<<endl;
				exit(255);
			}
			//continue;
		}
		// Calculate Posterior Probability
		memcpy(real_p_prior, &mat->p_prior[((ubit64_t)sites[j].ori&0x7)<<4], sizeof(double)*16);
		if ( (sites[j].ori & 0x8) && para->refine_mode) 
		{
			// a dbSNP
			snp_p_prior_gen(real_p_prior, call_chr->find_snp(sites[j].pos), para, sites[j].ori);
		}
		memset(mat->type_prob,0,sizeof(rate_t)*17);
		type2=type1=16;
		for (allele1=0; allele1!=4; ++allele1) 
		{
			for (allele2=allele1; allele2!=4; ++allele2) 
			{
				genotype = allele1<<2|allele2;
				if (para->is_monoploid && allele1 != allele2) 
				{
					continue;
				}
				mat->type_prob[genotype] = mat->type_likely[genotype] + log10(real_p_prior[genotype]) ;

				if (mat->type_prob[genotype] >= mat->type_prob[type1] || type1 == 16) 
				{
					// find the genotype which type_prob is biggest.
					type2 = type1;
					type1 = genotype;
				}
				else if (mat->type_prob[genotype] >= mat->type_prob[type2] || type2 ==16) 
				{
					// find the genotype which type_prob is second biggest.
					type2 = genotype;
				}
				else
				{
					;
				}
			}
		}
		is_out = true; // Check if the position needs to be output, useful in snp-only mode

		if (para->rank_sum_mode)
		{
			rank_sum_test_value = rank_test(sites[j], type1, mat->p_rank, para);
		}
		else 
		{
			rank_sum_test_value = 1.0;
		}
		if(rank_sum_test_value==0.0) 
		{
			// avoid double genotype overflow
			q_cns = 0;
		}
		else 
		{
			q_cns = (int)(10*(mat->type_prob[type1]-mat->type_prob[type2])+10*log10(rank_sum_test_value));
		}

		if ( (type1&3) == ((type1>>2)&3)) 
		{ // Called Homozygous
			if (qual1>0 && base1 != (type1&3))
			{
				//Wired: best base is not the consensus!
				q_cns = 0;
			}
			else if (/*qual2>0 &&*/ q_cns > qual1-qual2) 
			{
				// Should not bigger than this
				q_cns = qual1-qual2;
			}
			else 
			{
				;
			}
		}
		else 
		{	// Called Heterozygous
			if(sites[j].q_sum[base1]>0 && sites[j].q_sum[base2]>0 && type1 == (base1<base2 ? (base1<<2|base2):(base2<<2|base1))) 
			{
				// The best bases are in the heterozygote
				if (q_cns > qual2-qual3) 
				{
					q_cns = qual2-qual3;
				}
			}
			else 
			{	// Ok, wired things happened
				q_cns = 0;
			}
		}
		if(q_cns>99) 
		{
			q_cns = 99;
		}
		if (q_cns<0) 
		{
			q_cns = 0;
		}
		if(! para->glf_format) 
		{
			// ChrID\tPos\tRef\tCns\tQual\tBase1\tAvgQ1\tCountUni1\tCountAll1\tBase2\tAvgQ2\tCountUni2\tCountAll2\tDepth\tRank_sum\tCopyNum\tSNPstauts\n"
			if(! para->is_snp_only || (abbv[type1] != "ACTGNNNN"[(sites[j].ori&0x7)] && sites[j].depth > 0)) 
			{
				if(base1<4 && base2<4)
				{
					(*consensus)<<call_name<<'\t'<<(sites[j].pos+1)<<'\t'<<("ACTGNNNN"[(sites[j].ori&0x7)])<<'\t'<<abbv[type1]<<'\t'<<q_cns<<'\t'
					<<("ACTGNNNN"[base1])<<'\t'<<(sites[j].q_sum[base1]==0?0:sites[j].q_sum[base1]/sites[j].count_uni[base1])<<'\t'<<sites[j].count_uni[base1]<<'\t'<<sites[j].count_all[base1]<<'\t'
					<<("ACTGNNNN"[base2])<<'\t'<<(sites[j].q_sum[base2]==0?0:sites[j].q_sum[base2]/sites[j].count_uni[base2])<<'\t'<<sites[j].count_uni[base2]<<'\t'<<sites[j].count_all[base2]<<'\t'
					<<sites[j].depth<<'\t'<<rank_sum_test_value<<'\t'<<(sites[j].depth==0?255:(double)(sites[j].repeat_time)/sites[j].depth)<<'\t'<<((sites[j].ori&8)?1:0)<<"\n";
				}
				else if(base1<4) 
				{
					(*consensus)<<call_name<<'\t'<<(sites[j].pos+1)<<'\t'<<("ACTGNNNN"[(sites[j].ori&0x7)])<<'\t'<<abbv[type1]<<'\t'<<q_cns<<'\t'
					<<("ACTGNNNN"[base1])<<'\t'<<(sites[j].q_sum[base1]==0?0:sites[j].q_sum[base1]/sites[j].count_uni[base1])<<'\t'<<sites[j].count_uni[base1]<<'\t'<<sites[j].count_all[base1]<<'\t'
					<<"N\t0\t0\t0\t"
					<<sites[j].depth<<'\t'<<rank_sum_test_value<<'\t'<<(sites[j].depth==0?255:(double)(sites[j].repeat_time)/sites[j].depth)<<'\t'<<((sites[j].ori&8)?1:0)<<"\n";
				}
				else
				{
					(*consensus)<<call_name<<'\t'<<(sites[j].pos+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<"\n";
				}
			}
		}
		else 
		{
			if(base1<4 && base2<4)
			{
				baseinfo<<call_name<<'\t'<<(sites[j].pos+1)<<'\t'<<("ACTGNNNN"[(sites[j].ori&0x7)])<<'\t'<<abbv[type1]<<'\t'<<q_cns<<'\t'
					<<("ACTGNNNN"[base1])<<'\t'<<(sites[j].q_sum[base1]==0?0:sites[j].q_sum[base1]/sites[j].count_uni[base1])<<'\t'<<sites[j].count_uni[base1]<<'\t'<<sites[j].count_all[base1]<<'\t'
					<<("ACTGNNNN"[base2])<<'\t'<<(sites[j].q_sum[base2]==0?0:sites[j].q_sum[base2]/sites[j].count_uni[base2])<<'\t'<<sites[j].count_uni[base2]<<'\t'<<sites[j].count_all[base2]<<'\t'
					<<sites[j].depth<<'\t'<<showpoint<<rank_sum_test_value<<'\t'<<(sites[j].depth==0?255:(double)(sites[j].repeat_time)/sites[j].depth)<<'\t'<<((sites[j].ori&8)?1:0)<<"\n";
			}
			else if(base1<4) 
			{
				baseinfo<<call_name<<'\t'<<(sites[j].pos+1)<<'\t'<<("ACTGNNNN"[(sites[j].ori&0x7)])<<'\t'<<abbv[type1]<<'\t'<<q_cns<<'\t'
					<<("ACTGNNNN"[base1])<<'\t'<<(sites[j].q_sum[base1]==0?0:sites[j].q_sum[base1]/sites[j].count_uni[base1])<<'\t'<<sites[j].count_uni[base1]<<'\t'<<sites[j].count_all[base1]<<'\t'
					<<"N\t0\t0\t0\t"
					<<sites[j].depth<<'\t'<<showpoint<<rank_sum_test_value<<'\t'<<(sites[j].depth==0?255:(double)(sites[j].repeat_time)/sites[j].depth)<<'\t'<<((sites[j].ori&8)?1:0)<<"\n";
			}
			else 
			{
				baseinfo<<call_name<<'\t'<<(sites[j].pos+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<"\n";
			}
		}
		//if(sites[j].pos == 250949) {
		//	exit(1);
		//}
	}
//	delete [] real_p_prior;
//	delete [] pcr_dep_count;
//	delete [] likelihoods;
	return 1;
}

//update 2010-11-08 guyue
/**
 * DATE:  
 * FUNCTION:    call the consensus from soap file . 
 * PARAMETER:   
 * RETURN :     
 */
int Call_win::soap2cns(igzstream & alignment, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id)
{
	Soap_format soap;
	current_chr = genome->chromosomes.end();
	for(std::string line; getline(alignment, line);) 
	{
		std::istringstream ss(line);
		if( ss >> soap ) 
			deal_read(soap, consensus, baseinfo, genome, mat, para, sfsMethod, id);
	}

	deal_tail(consensus, baseinfo, genome, mat, para, sfsMethod, id);

	alignment.close();
	consensus->close();
	baseinfo.close();
	return 1;
}

/* program by Bill */
int Call_win::soap2cns(SamCtrl &alignment, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id) {
	Soap_format soap;
	current_chr = genome->chromosomes.end();
	int r;
	std::string line;
	while((r = alignment.readline(line)) >=0) {
		line = alignment_format(line);   //from sam/bam to soap file 

		if (line == NOUSE_ALIGNMENT)
			continue;
		std::istringstream ss(line);
		if( ss >> soap ) 
			deal_read(soap, consensus, baseinfo, genome, mat, para, sfsMethod, id);
	}

	deal_tail(consensus, baseinfo, genome, mat, para, sfsMethod, id);

	consensus->close();
	baseinfo.close();
	return 1;
}

/**
 * DATE: 2010-8-5
 * FUNCTION: deal with the win.
 * PARAMETER: current_chr: is the chr to be processed. last_start: the last position that be processed
 *				consensus: output file handle. genome: is the point to the Genome.
 *				mat: the matrix. para: the parameter.
 *				sfsMethod: a SfsMethod object, use to sfs.
 *				id: the individual's id, use to sfs
 * RETURN: void
 */
void Call_win::pro_win(gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter *para, SfsMethod &sfsMethod, const int id)
{
	//done_pro_win = true; // 2011-1-20
	
	if (current_chr == genome->chromosomes.end())
	{
		return;
	}

	if (!para->region_only || current_chr->second->is_in_region_win(last_start)) 
	{
		//cerr<<"at"<<soap.get_pos()<<"Called"<<last_start<<endl;
		if(last_start > sites[win_size-1].pos) 
		{
			initialize((last_start/win_size)*win_size);
		}
		call_cns(current_chr->first, current_chr->second, win_size, mat, para, consensus, baseinfo, sfsMethod, id);
		last_start = sites[win_size-1].pos;  //last start is the end position of the process window 
		recycled = false;
	}

	if (!para->region_only) 
	{
		//cerr<<"recycled"<<last_start<<endl;
		recycle();
		recycled = true;
		last_start = sites[win_size-1].pos;
	}
	else if (!recycled) 
	{
		//cerr<<"recycled"<<" "<<last_start<<endl;
		/*	if (last_start>(sites[win_size-1].pos)) {
		cerr<<"Unexpected "<<last_start<<">"<<(sites[win_size-1].pos)<<endl;
		exit(1);
		if ((last_start+1)%win_size!=0) {
		cerr<<"Assertion Error!"<<last_start<<">"<<sites[win_size-1].pos<<endl;
		exit(1); 
		//cerr << last_satrt <<" "<<sites[win_size-1].pos<<endl;
		}
		initialize(last_start-win_size+1);
		}
		else {*/
		if (current_chr->second->is_in_region_win(sites[win_size].pos)) 
		{
			recycle();
		}
		else
		{
			deep_init(sites[win_size].pos);
		}
		//}
		recycled = true;
		last_start = sites[win_size-1].pos;
	}
	else 
	{
		assert((last_start+1)%win_size==0);
		last_start += win_size;
	}
}

/**
 * DATE: 2010-8-5
 * FUNCTION: deal with one read.
 * PARAMETER: recycled: the flag. 
 *				soap: the soapformat to be processed. soap: SOAP format read from the file.
 *				consensus: output file handle. genome: is the point to the Genome.
 *				mat: the matrix. para: the parameter.
 *				sfsMethod: a SfsMethod object, use to sfs.
 *				id: the individual's id, use to sfs
 * RETURN: void
 */
void Call_win::deal_read(Soap_format &soap, gzoutstream * consensus, my_ofstream & baseinfo, Genome * genome, Prob_matrix * mat, Parameter * para, SfsMethod &sfsMethod, const int id) 
{
	int coord, sub;
	//cerr<<"A"<<soap.get_pos()<<endl;
	//exit(1);
	if(soap.get_pos() < 0)
	{
		return;
	}	
	// first time deal_read or come to a new chromose
	if (current_chr == genome->chromosomes.end() || current_chr->first != soap.get_chr_name()) 
	{
		//come to a new chromose
		if(current_chr != genome->chromosomes.end()) 
		{
			recycled = false;
			while(current_chr->second->length() > sites[win_size-1].pos && current_chr->second->length() > last_start && (current_chr->second->get_region() != NULL || current_chr->second->m_region_len > last_start)) 
			{
				pro_win(consensus, baseinfo, genome, mat, para, sfsMethod, id);
			}
	
			// The last window
			if(last_start > sites[win_size-1].pos) 
			{
				initialize((last_start/win_size)*win_size);
			}
			call_cns(current_chr->first, current_chr->second, current_chr->second->length()%win_size, mat, para, consensus, baseinfo, sfsMethod, id);
	
		}
		
		current_chr = genome->chromosomes.find(soap.get_chr_name());
		//initialize(0);
		/*if (last_start < current_chr->second->getStartPos())
			last_start = current_chr->second->getStartPos();*/
		last_start = current_chr->second->getStartPos();
		//initialize((last_start/win_size)*win_size);
		// use deep_init to make sites to be pure. 2010-12-16 Bill.
		deep_init((last_start/win_size)*win_size);
		cerr<<"Processing "<<current_chr->first<<endl;
	}
	else {
		;
	}

	//the end of this read is out of the chromorese 
	if (soap.get_pos()+soap.get_read_len() >= current_chr->second->length()) {
		return;
	}
	
	//is out of targ region
	if (para->region_only && (current_chr->second->get_region() == NULL || !current_chr->second->is_in_region(soap.get_pos()))) {
		return;
	}
	
	if(soap.get_pos() < last_start) 
	{
		//cerr << soap.get_pos() << "<" <<last_start<<endl;
		return;
		/*cerr<<"Errors in sorting:"<<soap.get_pos()<<"<"<<last_start<<endl;
		exit(255);*/
	}
	
	recycled = false;
	//some window before the start of this soap  isn't process
	while (soap.get_pos()/win_size > last_start/win_size ) 
	{ 
		//if ( soap.get_pos()/win_size - last_start/win_size >1)
		//cerr << soap.get_pos()/win_size <<" "<<last_start/win_size<<endl;
		pro_win(consensus, baseinfo, genome, mat, para, sfsMethod, id);  //deal the pre-window
	}
	// 2011-1-20
	done_pro_win = true;
	//cerr<<"die"<<endl;
	//exit(1);
	last_start=soap.get_pos();
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	
	for(coord=0; coord<soap.get_read_len(); coord++){
		if( (soap.get_pos()+coord)/win_size == soap.get_pos()/win_size ) {
			// In the same sliding window
			sub = (soap.get_pos()+coord) % win_size;
		}
		else {
			sub = (soap.get_pos()+coord) % win_size + win_size; // Use the tail to store the info so that it won't intervene the uncalled bases
		}

		if (sub >= (win_size + para->read_length))
		{
			// update at 2010-10-12
			cerr << "\tWARNING: the parameter '-L' was wrong, you should set the longest read's length!!!" 
				 << "  At position : " << soap.get_pos()
				 << endl;
			break;	// we don't need to exit the program. juet beak;
			// exit(0);
		}

		sites[sub].depth += 1;
		sites[sub].repeat_time += soap.get_hit();
		if((soap.is_N(coord)) || soap.get_qual(coord)<para->q_min || sites[sub].dep_uni >= 0xFF) {
			// An N, low quality or meaningless huge depth
			continue;
		}
		//just deal with the position best hit equal to 1
		if(soap.get_hit() == 1) {
			//exit((fprintf(stderr, "Wo Cao!\n")));
			sites[sub].dep_uni += 1;
			// Update the covering info: 4x2x64x64 matrix, base x strand x q_score x read_pos, 2-1-6-6 bits for each
			if(soap.is_fwd()) {
				// Binary strand: 0 for plus and 1 for minus
				sites[sub].base_info[(((ubit64_t)(soap.get_base(coord)&0x6)|0))<<14 | ((ubit64_t)(soap.get_qual(coord)-para->q_min))<<8 | coord ] += 1;
			}
			else {
				sites[sub].base_info[(((ubit64_t)(soap.get_base(coord)&0x6)|1))<<14 | ((ubit64_t)(soap.get_qual(coord)-para->q_min))<<8 | (soap.get_read_len()-1-coord) ] += 1;

			}
			sites[sub].count_uni[(soap.get_base(coord)>>1)&3] += 1;
			sites[sub].q_sum[(soap.get_base(coord)>>1)&3] += (soap.get_qual(coord)-para->q_min);

			// add by Bill at 2010-10-15
			if (para->sfs != NULL && soap.get_qual(coord) > para->sfs->qs)
			{
				sites[sub].count_sfs[(soap.get_base(coord)>>1)&3] += 1;
			}
		}
		else {
			;// Repeats
		}
		sites[sub].count_all[(soap.get_base(coord)>>1)&3] += 1;
		
	}
}

// processed tail of chromosome
/**
 * DATE: 2010-8-18
 * FUNCTION: processed tail of chromosome.
 * PARAMETER: recycled: the flag. consensus: output file handle. baseinfo: output file handle.
 *				genome: is the point to the Genome.
 *				mat: the matrix. para: the parameter.
 *				sfsMethod: a SfsMethod object, use to sfs.
 *				id: the individual's id, use to sfs
 * RETURN: void
 */
void Call_win::deal_tail(gzoutstream * consensus, my_ofstream& baseinfo, Genome* genome, Prob_matrix* mat, Parameter* para, SfsMethod &sfsMethod, const int id)
{
	//The unprocessed tail of chromosome
	recycled = false;
	while(current_chr->second->length() > sites[win_size-1].pos && current_chr->second->length() > last_start && (current_chr->second->get_region() != NULL || current_chr->second->m_region_len > last_start)) {

		pro_win(consensus, baseinfo, genome, mat, para, sfsMethod, id);

	}
	// The last window
	if(last_start > sites[win_size-1].pos) {
		initialize((last_start/win_size)*win_size);
	}  
	call_cns(current_chr->first, current_chr->second, current_chr->second->length()%win_size, mat, para, consensus, baseinfo, sfsMethod, id);
}
