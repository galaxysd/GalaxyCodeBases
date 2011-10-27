#include "simulate.h"

//from ASCII of A C G T to 0 1 2 3, auto dealing with upper or lower case.
//8bit char type, A=a=N=n=0, C=c=1, G=g=2, T=t=3, others as 4.
char alphabet[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//from 0 1 2 3 4 to A C G T N
char Bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//form 0 1 2 3 4 to T G C A N
char c_bases[5] ={
		'T', 'G', 'C', 'A', 'N'
};


//GC bias relative abundant distribution, GC rate 0%~100%, GC bias relative abundant 0 ~ 100
double GC_bias_abundant[101]=
{
	0.19,0.19,0.19,0.19,0.19,0.19,0.19,0.19,0.19,0.19,
	0.73,0.73,0.73,0.73,0.73,0.73,0.73,0.73,0.73,0.73,
	0.94,0.94,0.94,0.94,0.94,0.94,0.94,0.94,0.94,0.94,
	0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,
	0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,
	0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,
	0.63,0.63,0.63,0.63,0.63,0.63,0.63,0.63,0.63,0.63,
	0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,
	0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,
	0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,
};

//for split string
void split(string strLine, vector<string>& tokens, const char* delim)
{
  for(;;) {
  	//erase delimiter
  	int i = strLine.find_first_not_of(delim);
  	if(i == -1)
    	break;
		strLine.erase(0, i);
   	i = strLine.find_first_of(delim);
  	if(i == -1) {
     	tokens.push_back(strLine);
    	break;
   	} 
   	else {
    	string token = strLine.substr(0, i);
    	strLine.erase(0, i);
     	tokens.push_back(token);
    }
  }
}

//Rrealization of snp
char get_match(char base){
	char bases[4][3]={{'T','G','C'},
			 {'A','G','C'},
			 {'A','T','G'},
			 {'A','T','C'}};
	int a;
	char n='N';
	switch (base)
	{
	case 'A': a=0;break;
	case 'T': a=1;break;
	case 'C': a=2;break;
	case 'G': a=3;break;
	default: return n;
	}
	int num=int(rand()%3);
	return bases[a][num];
}

//Produce heterozygous SNPs in multiploid
string Get_snp(string &seq,ofstream &snp,string id, float hetersnp_rate){
	int seq_len=seq.size();
	int snp_num=int(seq_len*hetersnp_rate);
	for (int i=0;i<snp_num;i++)
	{
		int index=int(rand()%seq_len);
		snp<<id<<"\t"<<index+1<<"\t"<<seq[index]<<"\t";
		seq[index]=get_match(seq[index]);
		snp<<seq[index]<<endl;
	}
	return seq;
}

//Getting insert sequence when heterozygous indel rate is bigger than 0
string get_insertion(int num){
	char base[]={'A','T','G','C'};
	string s;
	for (int a=0;a<num;a++)
	{
		int index=int(rand()%4);
		s+=base[index];
	}
	return s;
}

//Produce heterozygous indels in multiploid
string Get_indel(string &seq,ofstream &indel,string id1,float heterindel_rate){
	int seq_len=seq.size();
	int indel_num=int(seq_len*heterindel_rate);
	int array[3]={2,3,6};
	int p=1;
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<indel_num/2/array[i];j++)
		{
			int num=int(rand()%seq_len);
			if (num+p>seq_len)
			{
				j--;
			}else{
				indel<<id1<<"\t"<<"-"<<"\t"<<num+1<<"\t"<<p<<"\t";
				for (int k=0;k<p;k++)
				{
					indel<<seq[num+k];
					seq[num+k]='N';
				}
				indel<<endl;
			}	
		}
		p++;
	}
	p=1;
	vector<uint8_t> insert(seq_len);
	string s;
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<indel_num/2/array[i];j++)
		{
			int num=int(rand()%seq_len);
			insert[num]=p;
		}
		p++;
	}
	for (int i=0;i<seq_len;i++)
	{
		if (insert[i]>=1)
		{
			string temp;
			temp=get_insertion(insert[i]);
			s+=temp;
			indel<<id1<<"\t+\t"<<i+1<<"\t"<<int(insert[i])<<"\t"<<temp<<endl;
			if (seq[i]!='N')
			{
				s+=seq[i];
			}
		}else{
			if (seq[i]!='N')
			{
				s+=seq[i];
			}
		}
	}
	return s;
}

//Rrealization of error sequencing
char get_error_match(char base)
{
	char bases[4][3]={{'C','G','T'},
			 {'A','G','T'},
			 {'A','C','T'},
			 {'A','C','G'}};
			 	
	//simulate error bias
	double error_bias[4][3]={{0.518336873,0.812283803,1},
		{0.575332872,0.721712965,1},
		{0.282213738,0.430879522,1},
		{0.188680942,0.490612013,1}};
			 	
	int a;
	char n='N';
	switch (base)
	{
	case 'A': a=0;break;
	case 'C': a=1;break;
	case 'G': a=2;break;
	case 'T': a=3;break;
	default: return n;
	}
	double num = double(rand())/double(RAND_MAX);
	int i = 0;
	for(i=0; i<3; i++)
	{
		double p = error_bias[a][i];
		if(num <= p){break;}
	}
	return bases[a][i];
}


//get the reverse and complement sequence
string reversecomplementary (string read)
{	
	string rc_read;
	for (int64_t i=read.size()-1; i>=0; i--)
	{	
		rc_read.push_back(c_bases[alphabet[read[i]]]);
	}
	return rc_read;
}

/*
//simulate the insertsize distribution with the model of normal distribution function
//The insertsize range is limited in (¦Ì-5¦Ò£¬¦Ì+5¦Ò), which covers almost all the data.
vector <int> insert_distribution(int reads_pair){
	double pi=3.1415926535;
	vector <double> insert;
	vector <int> insert_num;
	double total,temp1;
	int temp2=0,total2=0,num=0;
	for (int i=insertsize_mean-5*insertsize_sd;i<=insertsize_mean+5*insertsize_sd;i++)
	{
		temp1=1/sqrt(2*pi)/insertsize_sd/exp(pow((i-insertsize_mean),2)/(2*pow(insertsize_sd,2)));
		insert.push_back(temp1);
		total+=temp1;
	}
	for (int i=insertsize_mean-5*insertsize_sd;i<=insertsize_mean+5*insertsize_sd;i++){
		temp2=int(0.5+insert[num++]/total*reads_pair);
		insert_num.push_back(temp2);
		total2+=temp2;
	}
	insert_num[5*insertsize_sd]+=reads_pair-total2;
	insert.clear();
	return insert_num;
}
*/

//simulate normal distribution insertsize by Box-muller method
int simulate_insertsize(int mean, int sd)
{
	int insertsize = 0;
	double x1,x2,radius_pow2,y1 = 0.0;
	
	do
	{
		x1 = 2.0*double(rand())/double(RAND_MAX)-1.0;
		x2 = 2.0*double(rand())/double(RAND_MAX)-1.0;
		radius_pow2 = x1*x1+x2*x2;
	}while(radius_pow2>=1.0 || radius_pow2 == 0);
	
	radius_pow2 = sqrt(-2.0*log(radius_pow2)/radius_pow2);
	y1 = x1*radius_pow2;
	insertsize = int(mean+y1*sd);
	
	return insertsize;
}
/*
//Simulate illumina error distribution on different cycles,
//with the model function f(x)=0.00001*x**4
vector <double> error_distribution(int rd_pair, float error_rate, int read_length){
	double basic_error_rate=0.001;
	double total_error=(error_rate-basic_error_rate)*read_length;
	vector <double> error_dist;
	double total;
	for (int i=1;i<=read_length;i++)
	{
		double temp=0.00001*pow(i,4);
		error_dist.push_back(temp);
		total+=temp;
	}
	for (int i=1;i<=read_length;i++)
	{
		double temp=(basic_error_rate+error_dist[i-1]/total*total_error)*rd_pair;
		error_dist[i-1]=temp;
	}
	return error_dist;
}
*/

//Simulate illumina error distribution on different cycles,
//with the model function f(x)=0.00001*x**4
vector <double> error_distribution(float error_rate, int read_length){
	double basic_error_rate=0.001;
	double total_error=(error_rate-basic_error_rate)*read_length;
	vector <double> error_dist;
	double total;
	for (int i=1;i<=read_length;i++)
	{
		double temp=0.00001*pow(i,4);
		error_dist.push_back(temp);
		total+=temp;
	}
	for (int i=1;i<=read_length;i++)
	{
		double temp=(basic_error_rate+error_dist[i-1]/total*total_error);
		error_dist[i-1]=temp;
	}
	return error_dist;
}


//simulate GC bias
int simulate_GC_bias(string insert_str){
	int is_ignore = 0;
	int GC_num = 0;
	int insert_size = insert_str.size();
	//get GC rate
	for (int i=0; i<insert_size; i++)
	{	
//		if(alphabet[insert_str[i]] == 1 || alphabet[insert_str[i]] == 2){GC_num++;} //auto dealing with upper or lower case.
		if(insert_str[i] == 'G' || insert_str[i] == 'C'){GC_num++;}
	}
	
	double GC_rate = double(GC_num)/double(insert_size);

	//get relative abundance
	double bias_abund = GC_bias_abundant[int(GC_rate*100)];
	
	//decide whether ignore this insert string.
	double num = double(rand())/double(RAND_MAX);
	if(num > bias_abund){is_ignore = 1;}
		
	return is_ignore;
}



//simulate quality value
string simulate_quality(vector<int> error_pos, int read_len, double** correct_base_quality, double** error_base_quality)
{
	string quality_str;
	//simulate correct bases quality line
	for(int i=0; i<read_len; i++)
	{
		double num = double(rand())/double(RAND_MAX);
		for(int j=2; j<41; j++)
		{
			if(num<=correct_base_quality[i][j])
			{
				int ascii_value=j+64;
				quality_str.push_back(char(ascii_value));
				break;
			}
		}
	}
	//simulate error bases quality value
	int error_num = error_pos.size();
	if(error_num > 0){
  	for(int i=0; i<error_num; i++)
  	{
  		double num = double(rand())/double(RAND_MAX);
  		for(int j=2; j<41; j++)
  		{
  			if(num<=error_base_quality[error_pos[i]][j])
  			{
  				int ascii_value=j+64;
  				quality_str[error_pos[i]] = char(ascii_value);
  				break;
  			}
  		}
  	}
	}
	
	return quality_str;
}
