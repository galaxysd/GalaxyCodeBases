#include "soap_snp.h"

// add chromosome to varible chromosomes.
bool Genome::add_chr(Chr_name & name) {
	Chr_info * new_chr = new Chr_info;
	pair<map<Chr_name, Chr_info*>::iterator, bool> insert_pair;
	insert_pair=chromosomes.insert(pair<Chr_name, Chr_info*>(name,new_chr));
	return insert_pair.second;
}

Genome::~Genome(){
	for( map<Chr_name, Chr_info*>::iterator iter=chromosomes.begin(); iter!= chromosomes.end(); iter++ ){
		;
	}
}
Chr_info::Chr_info(const Chr_info & other) {
	dbsnp = other.dbsnp;
	len = other.len;
	if (len%capacity==0) {
		bin_seq = new ubit64_t [len/capacity];
		memcpy(bin_seq, other.bin_seq, sizeof(ubit64_t)*len/capacity);
	}
	else {
		bin_seq = new ubit64_t [1+len/capacity];
		memcpy(bin_seq, other.bin_seq, sizeof(ubit64_t)*len/capacity);
	}
}


// allocate every sequence memory£¬ load sequence
int Chr_info::binarize(std::string & seq)
{
	len = seq.length();
	//update 11-26	

	//cerr<<len<<endl;
	// 4bit for each base
	// Allocate memory
	
	if (len%capacity==0) 
	{
		bin_seq = new ubit64_t [len/capacity];
		memset(bin_seq,0,sizeof(ubit64_t)*len/capacity);
	}
	else 
	{
		bin_seq = new ubit64_t [1+len/capacity];
		memset(bin_seq,0,sizeof(ubit64_t)*(1+len/capacity));
	}
	bool is_set_start = false;

	// Add each base, 7 is 0b111
	for(std::string::size_type i=0;i!=seq.length();i++) {
		//update 11-26
		
		if (!is_set_start && (seq[i] != 'N'))
		{
			//update 11-29
			if (i < 1000)
			{
				m_start_position = 0;
			}
			else
			{
				m_start_position = i - global_win_size;  //get the position is not N 
				m_start_position = m_start_position / global_win_size * global_win_size;
			}
			is_set_start = true;
		}
		bin_seq[i/capacity] |= ((((ubit64_t)seq[i]>>1)&7)<<(i%capacity*4));
	}
	return 1;
}

int Chr_info::insert_snp(std::string::size_type pos, Snp_info & snp_form) 
{
	Snp_info * new_snp = new Snp_info;
	*new_snp = snp_form;
	pair<map<ubit64_t, Snp_info*>::iterator, bool> insert_pair;
	insert_pair = dbsnp.insert(pair<ubit64_t, Snp_info*>(pos,new_snp));	
	if(insert_pair.second) 
	{
		// Successful insertion
		// Modify the binary sequence! Mark SNPs	
		
		bin_seq[pos/capacity] |= (1ULL<<(pos%capacity*4+3));
		
	}
	else 
	{
		cerr<<"Warning: Snp insertion failed\t"<<pos<<endl;  //this position is a duplication
		return 0;
	}
	return 1;
}
//set region ,if in the region , the postion in the region mask is 1
int Chr_info::set_region(int start, int end) {
	if(start<0) {
		start = 0;
	}
	else if (start >= len) { //start out of the chromose length
		start = len;
	}

	// add by Bill.
	if (end > m_region_len)
	{
		m_region_len = end;
	}

	if(end<0) {
		end = 0;
	}
	else if (end >= len) {
		end = len;
	}
	if (start > end) {
		cerr<<"Invalid region: "<<start<<"-"<<end<<endl;
		exit(255);
	}
	// Specific mask
	if(start/64 == end/64) {
		region_mask[start/64] |= ((~((~(0ULL))<<(end-start+1)))<<(63-end%64));
	}
	else {
		if(start % 64) {
			region_mask[start/64] |= (~((~(0ULL))<<(64-start%64)));
		}
		else {
			region_mask[start/64] = ~(0ULL);
		}
		region_mask[end/64] |= ((~(0ULL))<<(63-end%64));
		if(end/64-start/64>1) {
			memset(region_mask+start/64+1, 0xFF, sizeof(ubit64_t)*(end/64-start/64-1));
		}
	}
	// Window mask
	start /= global_win_size; end /= global_win_size;
	//cerr<<start<<"\t"<<end<<endl;
	if(start/64 == end/64) {
		region_win_mask[start/64] |= ((~((~(0ULL))<<(end-start+1)))<<(63-end%64));
	}
	else {
		if(start % 64) {
			region_win_mask[start/64] |= (~((~(0ULL))<<(64-start%64)));
		}
		else {
			region_win_mask[start/64] = ~(0ULL);
		}
		region_win_mask[end/64] |= ((~(0ULL))<<(63-end%64));
		if(end/64-start/64>1) {
			memset(region_win_mask+start/64+1, 0xFF, sizeof(ubit64_t)*(end/64-start/64-1));
		}
	}
	
	return 1;
}

/**
 * DATE:  
 * FUNCTION:    initialize region mask and ragion window,new array [chromorese len /64] or [chromorese len /64 +1]
 * PARAMETER:   
 * RETURN:   
 */
int Chr_info::region_mask_ini(){
	//Specific mask
	if(len%64==0) {
		region_mask = new ubit64_t [len/64];
		memset(region_mask, 0, sizeof(ubit64_t)*(len/64));
	}
	else {
		region_mask = new ubit64_t [len/64+1];
		memset(region_mask, 0, sizeof(ubit64_t)*(len/64+1));
	}
	//Window mask
	int win_len = len/global_win_size +1;
	if(win_len%64==0) {
		region_win_mask = new ubit64_t [win_len/64];
		memset(region_win_mask, 0, sizeof(ubit64_t)*(win_len/64));
	}
	else {
		region_win_mask = new ubit64_t [win_len/64+1];
		memset(region_win_mask, 0, sizeof(ubit64_t)*(win_len/64+1));
	}
	return 1;
}

// read every reads' starting and ending positions
int Genome::read_region(my_ifstream & region, Parameter * para) 
{
	Chr_name current_name(""), prev_name("");
	int start, end;
	map<Chr_name, Chr_info*>::iterator chr_iter;
	for(std::string buff;getline(region,buff);) 
	{
		std::istringstream s(buff);
		if(s>>current_name>>start>>end) 
		{
			if(current_name != prev_name) 
			{
				chr_iter = chromosomes.find(current_name);
				if(chr_iter == chromosomes.end()) 
				{
					cerr<<"Unexpected Chromosome:"<<current_name<<endl;
					continue;
				}
				//initialize the region pointer 
				if(NULL == chr_iter->second->get_region()) 
				{
					chr_iter->second->region_mask_ini();
				}
			}
			chr_iter->second->set_region(start-para->read_length, end-1);
			prev_name = current_name;
		}
		else 
		{
			cerr<<"Wrong format in target region file"<<endl;
			return 0;
		}
	}
	return 1;
}

Genome::Genome(my_ifstream &fasta,my_ifstream & known_snp) 
{
	std::string seq("");
	Chr_name current_name("");
	map<Chr_name, Chr_info*>::iterator chr_iter;
	for(std::string buff;getline(fasta,buff);) 
	{
		if('>' == buff[0]) 
		{	// Fasta id
			// Deal with previous chromosome
			if( chromosomes.find(current_name) != chromosomes.end()) 
			{
				chr_iter = chromosomes.find(current_name);
				chr_iter->second->binarize(seq); // initialize sequence
			}
			// Insert new chromosome
			std::string::size_type i;
			for(i=1;!isspace(buff[i]) && i != buff.length();i++) 
			{
				;
			}
			Chr_name new_chr_name(buff,1,i-1);
			if(! add_chr(new_chr_name)) 
			{
				std::cerr<<"Insert Chromosome "<<new_chr_name<<" Failed!\n";
			}
			current_name = new_chr_name;
			seq = "";
		}
		else 
		{
			seq += buff;
		}
	}
	if(seq.length() != 0 && chromosomes.find(current_name) != chromosomes.end()) 
	{
		chr_iter = chromosomes.find(current_name);
		chr_iter->second->binarize(seq);
	}
	
	// initialize known dbSNP
	if( known_snp ) 
	{
		Chr_name current_name;
		Snp_info snp_form;
		std::string::size_type pos;
		for(std::string buff;getline(known_snp, buff);) 
		{
			// Format: Chr\tPos\thapmap?\tvalidated?\tis_indel?\tA\tC\tT\tG\trsID\n
			std::istringstream s(buff);
			s>>current_name>>pos;
			s>>snp_form;
			if( chromosomes.find(current_name) != chromosomes.end()) 
			{
				// The SNP is located on an valid chromosme
				pos -= 1; // Coordinates starts from 0
				if( pos < chromosomes.find(current_name)->second->length())//update by zhukai on 2010-12-28
				{
					(chromosomes.find(current_name)->second)->insert_snp(pos, snp_form);
				}
				else
				{
					cerr<<"Warning :The current position is "<<pos + 1<<" out of the region of the reference sequence "<<endl;//update by zhukai on 2010-12-28
					
				}
			}
		}
	}
}
