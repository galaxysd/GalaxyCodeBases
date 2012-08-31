//本程序根据参考序列上mapped的reads来计算测序样本的genotype, 读取bwa mpielup格式输入文件
//依据bayes theory，以sample出现纯合和杂合SNP的概率和测序错误分值模型为先验概率，计算观测数据下每一种genotype的后验概率。
//测序错误分值模型，与模拟reads所用的模型相同，由分析重测序数据而得。它是一个四维矩阵：ref_base x cycle x seq_base x quality。
//解释：我们实际想得到的是P(D|genotype), 此处ref_base即为geotype, 而D包含三个属性(cycle,seq_base,quality)。
//我们想得到是P(cycle,seq_base,quality|ref_base), 而它又等价于P(cycle|ref_base) * P(seq_base,quality|cycle,ref_base)。
//因为P(cycle|ref_base)可以认为是一个常数，在整个bayes公式的计算中不起作用，因此我们把它忽略。
//这样可以将P(seq_base,quality|cycle,ref_base)直接用于bayes公式，也就形成了ref_base*cycle作为列，seq_base*quality作为行的二维matrix文件。
//本程序采用long double (128bit浮点数)，来高精度的表示概率。
//coverage depth最大允许值为400,以防计算过程中数值指数部分超过long double界限。超过400的部分将被忽略，程序不报错不中止。
//选取概率最大的genotype作为最终结果，并将其错误率转换为phred quality score, 并设上限为100。

//pirs和soapsnp的matrix是等效的，为保持一致，此程序采用soapsnp格式的matrix，  P(base|ref,qual,cycle)

//Author: Fan Wei, fanweisz09@gmail.com
//Date: 2012/6/4


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"

using namespace std;

typedef long double rate_t;
rate_t PTi[4][10]; //store the prior probability for each substitutions
rate_t *PdTi;   //store the distribution of error and quality
int Ref_Base_num;
int Cycle_num;
int Seq_Base_num;
int Quality_num;
int Qual_add_val;
int Cycle_add_val;
int Ref_add_val;
int Ref_Cyc_base_qual;
string Error_profile_file;
string Out_prefix = "output";
float Trans_Tranv_ratio = 2.0;
float Hom_SNP_rate = 0.0005;
float Het_SNP_rate = 0.0005;
int SNP_prior_mode = 1;
int Ascii_Qual_Start = 33;
int Only_output_snp = 1;

//由０１２３ 4到ＡＣＧＴ N
char Bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//用0，1，2，3，分别代表A,C,G,T
int DipGeno[10][2] = {0,0, 1,1, 2,2, 3,3, 0,1, 0,2, 0,3, 1,2, 1,3, 2,3};

map<int,char> TwoBases;
map<int, int> TransSubs;


void define_TransSubs()
{
	TransSubs[0] = 2;  TransSubs[2] = 0;
	TransSubs[1] = 3;  TransSubs[3] = 1;
}


void define_two_bases_markers()
{

	TwoBases[0] = 'A';	TwoBases[1] = 'C';	TwoBases[2] = 'G';	TwoBases[3] = 'T';
	TwoBases[4] = 'M';	TwoBases[5] = 'R';  TwoBases[6] = 'W';
	TwoBases[7] = 'S'; TwoBases[8] = 'Y'; TwoBases[9] = 'K';
}


void usage()
{	cout << "calculate_diploid_genotype  <*.mpileup>" << endl;
	cout << "   version 1.0" << endl;
	cout << "   -r <float>  homozygous SNP prior rate, default=" << Hom_SNP_rate << endl;
	cout << "   -e <float>  heterozygous SNP prior rate, default=" << Het_SNP_rate << endl;
	cout << "   -t <float>  transition vs transversion prior ratio, default=" << Trans_Tranv_ratio << endl;
	cout << "   -d <int>    set the prior mode, 0: no prior; 1: with prior;  default="<< SNP_prior_mode << endl;
	cout << "   -m <file>   matrix file of error and quality distributions, default=exe_path/XieSunney.blood.matrix.soapsnp" << endl;
	cout << "   -q <int>    ASCII chracter standing for quality==0, default=" << Ascii_Qual_Start << endl;
	cout << "   -o <string> prefix for the output files" << endl;
	cout << "   -s <int>    two output mode, 1, only snp; 0, cns; default=" << Only_output_snp << endl;
	cout << "   -h          get help information" << endl;
	exit(0);
}


//由ＡＣＧＴ N得到０１２３ 4
int get_base (char refBase)
{
	int refBaseId;
	switch (refBase)
	{
		case 'A': refBaseId = 0; break;
		case 'C': refBaseId = 1; break;
		case 'G': refBaseId = 2; break;
		case 'T': refBaseId = 3; break;
		case 'a': refBaseId = 0; break;
		case 'c': refBaseId = 1; break;
		case 'g': refBaseId = 2; break;
		case 't': refBaseId = 3; break;
		default:  refBaseId = 4;
	}
	return refBaseId;
}


//限定freq_rate在1到1e-9之间，统计样本量级1G
void load_soapsnp_error_profile (string &matrix_file)
{
	igzstream infile;
	infile.open(matrix_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << matrix_file << endl;
		exit(0);
	}

	int pi = 0; //PdTi的循环变量
	int line_field_number = 18;
	vector<string> textlines;

	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{	if (lineStr[0] == '#')
		{	continue;
		}

		vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
		if (lineVec.size() != line_field_number)
		{	continue;
		}

		Quality_num = boost::lexical_cast<int>(lineVec[0]) + 1;
		Cycle_num = boost::lexical_cast<int>(lineVec[1]) + 1;

		textlines.push_back(lineStr);
	}

	Ref_Base_num = 4;
	Seq_Base_num = 4;
	Ref_Cyc_base_qual = Ref_Base_num * Seq_Base_num * Quality_num * Cycle_num;
	Qual_add_val = Ref_Base_num * Seq_Base_num * Cycle_num;
	Cycle_add_val = Ref_Base_num * Seq_Base_num;
	Ref_add_val = Seq_Base_num;

	PdTi = new rate_t[Ref_Cyc_base_qual];
	for (int j=0; j<textlines.size(); j++)
	{
		vector<string> lineVec;
		boost::split(lineVec,textlines[j], boost::is_any_of(" \t\n"), boost::token_compress_on);
		int this_qual = boost::lexical_cast<int>(lineVec[0]);
		int this_cycl = boost::lexical_cast<int>(lineVec[1]);

		for (int i=2; i<line_field_number; i++)
		{	int k = i - 2;
			int this_ref = k / Seq_Base_num;
			int this_base = k % Seq_Base_num;
			pi = this_qual*Qual_add_val + this_cycl*Cycle_add_val + this_ref*Ref_add_val  + this_base;
			rate_t freq_rate = boost::lexical_cast<rate_t>(lineVec[i]);

			if (pi < Ref_Cyc_base_qual)
			{	PdTi[pi] = (freq_rate > 1e-9L) ? freq_rate : 1e-9L;
			}else{
				cerr << "pi exceeds its range " << Ref_Cyc_base_qual <<  endl;
				exit(0);
			}
		}
	}

}


//根据Trans_Tranv_ratio,Hom_SNP_rate,Het_SNP_rate,SNP_prior_mode来计算genotype的prior probability
void assign_PTi_rates ()
{
	double TsTvR = Trans_Tranv_ratio;
	double NonSnpRate = 1 - Hom_SNP_rate - Het_SNP_rate;
	double hom_trans_rate = Hom_SNP_rate * TsTvR*2 / (TsTvR*2 + 1 + 1);
	double hom_tranv_rate = Hom_SNP_rate * 1 / (TsTvR*2 + 1 + 1);

	double HapSnpRate = Hom_SNP_rate + Het_SNP_rate;  //近似，非准确
	double HapNonSnpRate = 1 - HapSnpRate;
	double HapTransSnpRate = HapSnpRate * TsTvR*2 / (TsTvR*2 + 1 + 1);
	double HapTranvSnpRate = HapSnpRate * 1 / (TsTvR*2 + 1 + 1);

	double hetNonTransRate;
	double hetNonTranvRate;
	double hetTransTranvRate;
	double hetTranvTranvRate;
	double relativeTotal;

	relativeTotal = HapNonSnpRate*HapTransSnpRate * 1 + HapNonSnpRate*HapTranvSnpRate *2 + HapTransSnpRate*HapTranvSnpRate * 2 + HapTranvSnpRate*HapTranvSnpRate  * 1;
	hetNonTransRate = (relativeTotal > 0) ? HapNonSnpRate*HapTransSnpRate / relativeTotal * Het_SNP_rate : 0;
	hetNonTranvRate = (relativeTotal > 0) ? HapNonSnpRate*HapTranvSnpRate / relativeTotal * Het_SNP_rate : 0;
	hetTransTranvRate = (relativeTotal > 0) ? HapTransSnpRate*HapTranvSnpRate / relativeTotal * Het_SNP_rate : 0;
	hetTranvTranvRate = (relativeTotal > 0) ? HapTranvSnpRate*HapTranvSnpRate / relativeTotal * Het_SNP_rate : 0;

	if (SNP_prior_mode == 0)
	{	NonSnpRate = hom_trans_rate = hom_tranv_rate = hetNonTransRate = hetNonTranvRate = hetTransTranvRate = hetTranvTranvRate = 0.1;
	}

	for (int i=0; i<4; i++)
	{	for (int j=0; j<10; j++)
		{	if (i == DipGeno[j][0]  && i == DipGeno[j][1])   // non-SNP
			{
				PTi[i][j] = NonSnpRate;
			}else if(i != DipGeno[j][0]  && DipGeno[j][0] == DipGeno[j][1])  //homo SNP
			{
				PTi[i][j] = (TransSubs[i] == DipGeno[j][0]) ? hom_trans_rate : hom_tranv_rate ;
			}else //het SNP
			{	if (i == DipGeno[j][0]  || i == DipGeno[j][1]) //het-one
				{	PTi[i][j] = (TransSubs[i] == DipGeno[j][0] || TransSubs[i] == DipGeno[j][1]) ? hetNonTransRate : hetNonTranvRate;
				}else //het-two
				{	PTi[i][j] = (TransSubs[i] != DipGeno[j][0] && TransSubs[i] != DipGeno[j][1]) ? hetTranvTranvRate : hetTransTranvRate;
				}
			}
		}
	}

}

//从pileup格式的base_str中得到干净的bases信息
void get_clean_bases(int refBase, string &raw_base_str,vector<int> &clean_base_vec)
{
	int pos = 0;
	int len = raw_base_str.size();
	while (pos < len)
	{
		if (raw_base_str[pos] == '^')
		{	pos += 2;
		}
		else if (raw_base_str[pos] == '+' || raw_base_str[pos] == '-')
		{	pos += 1;
			string indel_num;
			while (pos < len)
			{	if (raw_base_str[pos] >= 48 && raw_base_str[pos] <= 57 )  //ASCII: 0-9
				{	indel_num.push_back(raw_base_str[pos]);
					pos ++;
				}
				else
				{	break;
				}
			}
			pos += atoi(indel_num.c_str());
		}
		else
		{
			switch (raw_base_str[pos])
			{
				case '.': clean_base_vec.push_back(refBase); break;
				case ',': clean_base_vec.push_back(refBase); break;
				case 'A': clean_base_vec.push_back(0); break;
				case 'a': clean_base_vec.push_back(0); break;
				case 'C': clean_base_vec.push_back(1); break;
				case 'c': clean_base_vec.push_back(1); break;
				case 'G': clean_base_vec.push_back(2); break;
				case 'g': clean_base_vec.push_back(2); break;
				case 'T': clean_base_vec.push_back(3); break;
				case 't': clean_base_vec.push_back(3); break;
				case 'N': clean_base_vec.push_back(4); break;
				case 'n': clean_base_vec.push_back(4); break;
				case '*': clean_base_vec.push_back(4); break;
				case '>': clean_base_vec.push_back(4); break;
				case '<': clean_base_vec.push_back(4); break;
				//自动省略$符号的判断
			}
			pos += 1;
		}

	}

}

//从pileup格式的quality_str中得到quality信息，并转换成phred-scale, 33开始
void get_quality(vector<int> &qual_vec, string &qual_str)
{
	for (int i=0; i<qual_str.size(); i++)
	{	int qual_val = qual_str[i] - Ascii_Qual_Start;
		qual_val = (qual_val < Quality_num) ? qual_val : Quality_num - 1;  //质量值大于40的全部当成40
		qual_vec.push_back(qual_val);
	}
}

//从pileup格式的cycle_str中得到cycle信息
void get_cycle (vector<int> &cycl_vec, string &cycl_str)
{
	vector<string> temp_vec;
	boost::split(temp_vec,cycl_str, boost::is_any_of(","), boost::token_compress_on);
	for (int i=0; i<temp_vec.size(); i++)
	{	//int cycle_val = boost::lexical_cast<int>(temp_vec[i]);  //runs much slow
		int cycle_val = atoi(temp_vec[i].c_str()); //runs fast
		if (cycle_val > Cycle_num)
		{	cycle_val = Cycle_num;
		}
		cycl_vec.push_back(cycle_val - 1);  //cycle 0-99代表100个cycles
	}
}


//根据指定的参考序列碱基，计算每一种genotype的后验概率P(Ti|D)
//coverage depth最大允许值为400,以防计算过程中数值指数部分超过long double界限
void calculte_genotype_probability(int refBase, vector<rate_t> &pTiDvec, vector<int> &baseVec, vector<int> &qualVec, vector<int> &cyclVec)
{
	//cerr << "refbase: " << refBase << "  " << Bases[refBase] << endl;
	rate_t pD = 0.0;
	int len = (baseVec.size() < 400) ? baseVec.size() : 400;

	//vector<rate_t> pDTiVec;
	for (int genotype=0; genotype<10; genotype++)
	{
		rate_t pTi = PTi[refBase][genotype];  //
		rate_t pDTi = 1.0;

		if (pTi != 0)
		{	for (int j=0; j<len; j++)
			{
				int num1 = qualVec[j]*Qual_add_val + cyclVec[j]*Cycle_add_val + DipGeno[genotype][0]*Ref_add_val +  + baseVec[j];
				int num2 = qualVec[j]*Qual_add_val + cyclVec[j]*Cycle_add_val + DipGeno[genotype][1]*Ref_add_val +  + baseVec[j];
				if (num1 >= Ref_Cyc_base_qual || num2 >= Ref_Cyc_base_qual)
				{	cerr << "num1 or num2 exceeds its range " << num1 << "\t" << num2 << endl;
					exit(0);
				}
				pDTi *=    (PdTi[num1] + PdTi[num2] ) / 2;
				//cerr << "# " << DipGeno[genotype][0] << "\t" << DipGeno[genotype][1] << "\t" << cyclVec[j] << "\t" << baseVec[j] << "\t" << qualVec[j] << "\t"<< PdTi[num1] << "\t" << PdTi[num2]  << "\t" << num1 << "\t" << num2 << endl;

			}
		}
		//pDTiVec.push_back(pDTi);
		pTiDvec.push_back( pTi*pDTi );
		pD += pTi*pDTi;
		//cerr << "calculate: " << genotype << "\t" << pTi << "\t" << pDTi << "\t" << pTi*pDTi << "\t" << pD << endl;
	}

	for (int genotype=0; genotype<10; genotype++)
	{	pTiDvec[genotype] = pTiDvec[genotype] / pD;
	}

	//output middle results
	/*for (int genotype=0; genotype<10; genotype++)
	{	rate_t pTi = PTi[refBase][genotype];
		rate_t pDTi = pDTiVec[genotype];
		cout << genotype << "\t" << pTi << "\t" << pDTi << "\t" << pTi*pDTi << "\t" << pD << "\t" << pTiDvec[genotype] << endl;
	}*/
}

//find the genotype with maximum posterior likelihood and calcualte the phred-scale score (up to 100)
void find_genotype_score(vector<rate_t> &pTiDvec, int &sample_genotype, int &phred_score)
{
	int max_id = 0;
	rate_t max_likelihood = 0.0;
	for (int i=0; i<pTiDvec.size(); i++)
	{	if (pTiDvec[i] > max_likelihood)
		{	max_id = i;
			max_likelihood = pTiDvec[i];

		}
		//cerr << "find max: " << i << "  " << pTiDvec[i] << endl;
	}

	sample_genotype = max_id;

	rate_t error_rate = 1 - pTiDvec[max_id];
	if (error_rate < 1e-10L)
	{	error_rate = 1e-10L;
	}
	phred_score = int(-10 * log10(error_rate));

}

void output_prior_file(string &prior_file)
{
	ofstream PRIOR (prior_file.c_str());
	if (! PRIOR)
	{	cerr << "fail create prior file" << prior_file << endl;
	}

	PRIOR << "The 10 genotypes: \n";
	for (int i=0; i<10; i++)
	{	PRIOR << i << "(" << Bases[DipGeno[i][0]] << Bases[DipGeno[i][1]] << ")  ";
	}
	PRIOR << endl << endl;


	PRIOR << "The compressed representaiton for each genotype:\n";
	map<int,char>::const_iterator map_it = TwoBases.begin();
	while (map_it != TwoBases.end())
	{	PRIOR << map_it->first << ":" << map_it->second << ";  ";
		++map_it;
	}
	PRIOR << "\n\n";

	PRIOR << "The transison substitution matrix:\n";
	map<int, int>::const_iterator map_it2 = TransSubs.begin();
	while (map_it2 != TransSubs.end())
	{	PRIOR << Bases[map_it2->first] << ":" << Bases[map_it2->second] << "; ";
		++map_it2;
	}
	PRIOR << endl  << endl;


	PRIOR << "Hom_SNP_rate: " << Hom_SNP_rate << endl;
	PRIOR << "Het_SNP_rate: " << Het_SNP_rate << endl;
	PRIOR << "All_SNP_rate: " << Hom_SNP_rate + Het_SNP_rate << endl;
	PRIOR << "Trans_Tranv_ratio: " << Trans_Tranv_ratio << endl;

	PRIOR << "\nThe prior probabilites of genotypes for each reference base:\n\n";
	PRIOR << "R\\G";
	for (int j=0; j<10; j++)
	{	PRIOR << "\t" << Bases[DipGeno[j][0]] << Bases[DipGeno[j][1]];
	}
	PRIOR << endl;
	for (int i=0; i<4; i++)
	{	PRIOR << Bases[i] ;
		for (int j=0; j<10; j++)
		{
			PRIOR << "\t" << PTi[i][j];
		}
		PRIOR << endl;

	}

	PRIOR << "\n	Ref_Base_num " << Ref_Base_num << endl;
	PRIOR << "	Cycle_num " << Cycle_num << endl;
	PRIOR << "	Seq_Base_num " << Seq_Base_num << endl;
	PRIOR << "	Quality_num " << Quality_num << endl;
	PRIOR << "	PdTi array size: " << Ref_Cyc_base_qual << endl;
	PRIOR << "	Qual_add_val " << Qual_add_val << endl;
	PRIOR << "	Cycle_add_val " << Cycle_add_val << endl;
	PRIOR << "	Ref_add_val " << Ref_add_val << endl;

}

int main(int argc, char *argv[])
{
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "t:r:e:d:m:q:o:s:h")) !=-1) {
		switch(c) {
			case 't': Trans_Tranv_ratio=atof(optarg); break;
			case 'r': Hom_SNP_rate=atof(optarg); break;
			case 'e': Het_SNP_rate=atof(optarg); break;
			case 'd': SNP_prior_mode=atoi(optarg); break;
			case 'm': Error_profile_file=optarg; break;
			case 'q': Ascii_Qual_Start=atoi(optarg); break;
			case 'o': Out_prefix=optarg; break;
			case 's': Only_output_snp=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 2) usage();

	if (! Error_profile_file.size())
	{	string exepath = argv[0];
		int numh = exepath.rfind("/");
		Error_profile_file = exepath.substr(0,numh+1) + "XieSunney.blood.matrix.soapsnp";
	}
	cerr << "\nUsed error profile matrix file: " << Error_profile_file << endl << endl;

	string pileup_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数

	clock_t time_start, time_end;
	time_start = clock();

	//caluate the prior probability of each genotype
	define_two_bases_markers();
	define_TransSubs();
	assign_PTi_rates();
	cerr << "Assign prior probability to genotype done:\n\n";

	//load the error and quality profile matrix into memory
	load_soapsnp_error_profile(Error_profile_file);
	cerr << "Dimensions of error profile:\n";

	string prior_file = Out_prefix + ".prior";
	output_prior_file(prior_file);

	time_end = clock();
	cerr << "\nLoad the error profile done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	//parse the pileup file, and calculate the likelihood of each genotypes
	igzstream infile;
	infile.open(pileup_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << pileup_file << endl;
	}

	string gentoytpe_file;
	if (Only_output_snp == 1)
	{	gentoytpe_file = Out_prefix + ".snp";
	}else
	{	gentoytpe_file = Out_prefix + ".cns";
	}
	ofstream outfile (gentoytpe_file.c_str());
	if ( ! outfile )
	{	cerr << "fail to open output file" << gentoytpe_file << endl;
	}

	long int loop = 0;
	outfile << "#Chr\tPos\tRef\tGeno\tScore\tDepth\tBases\tQuals\tCycles\t[all mapped bases]\n";
	string lineStr;
	while (getline( infile, lineStr, '\n' ))
	{
		vector<string> lineVec;
		boost::split(lineVec,lineStr, boost::is_any_of("\t"), boost::token_compress_on);
		int refBase = get_base(lineVec[2][0]);

		vector<int> baseVec;
		get_clean_bases(refBase,lineVec[4],baseVec);

		vector<int> qualVec;
		get_quality(qualVec, lineVec[5]);

		vector<int> cyclVec;
		get_cycle(cyclVec, lineVec[7]);

		////////////////debug///////////////
		//cerr << "refBase: " << Bases[refBase] << endl;
		//for (int i=0; i<baseVec.size(); i++)
		//{	cerr << i+1 << "\t" << Bases[baseVec[i]] << "\t" << qualVec[i] << "\t" << cyclVec[i] << endl;
		//}

		////////////////debug///////////////

		//check the result of parsing pileup file
		if (baseVec.size() != qualVec.size() || baseVec.size() != cyclVec.size())
		{	cerr << "Error happens at : " << lineVec[4] << "\t" << lineVec[5] << "\t" << lineVec[7] << endl << endl;
			exit(0);
		}

		vector<int> FbaseVec;
		vector<int> FqualVec;
		vector<int> FcyclVec;

		//filter the observed data, base not N, and qual > 5;
		for (int i=0; i<baseVec.size(); i++)
		{	if (baseVec[i] != 4 && qualVec[i] >= 5) //当质量值小于5的时候，P(D|G)几乎是相等的，因此P(G|D) = P(G), 数据没有发挥作用，结果将只反映先验概率，造成错误。
			{	FbaseVec.push_back(baseVec[i]);
				FqualVec.push_back(qualVec[i]);
				FcyclVec.push_back(cyclVec[i]);
			}
		}
		if (FbaseVec.size() == 0)
		{	continue;
		}

		vector<rate_t> pTiDvec;
		calculte_genotype_probability(refBase,pTiDvec, FbaseVec, FqualVec, FcyclVec);
		int sample_genotype = 0;
		int phred_score = 0;
		find_genotype_score(pTiDvec,sample_genotype,phred_score);
		if (Only_output_snp == 0 || sample_genotype != refBase)
		{

			outfile << lineVec[0] << "\t" << lineVec[1] << "\t" << lineVec[2] << "\t";
			if (sample_genotype < 4)
			{	outfile << Bases[sample_genotype];
			}else{
				outfile << TwoBases[sample_genotype] << "(" << Bases[DipGeno[sample_genotype][0]] << Bases[DipGeno[sample_genotype][1]] << ")";
			}

			outfile << "\t"	<< phred_score << "\t" << FbaseVec.size() << "\t" ;
			for (int i=0; i<FbaseVec.size(); i++)
			{	outfile << Bases[FbaseVec[i]];
			}
			outfile << "\t";
			for (int i=0; i<FbaseVec.size(); i++)
			{	outfile << FqualVec[i] << ",";
			}
			outfile << "\t";
			for (int i=0; i<FbaseVec.size(); i++)
			{	outfile << FcyclVec[i] << ",";
			}
			outfile << "\t" << baseVec.size() << "\t"<< lineVec[4] << "\t" << lineVec[5] << "\t" << lineVec[7] << endl;

		}

		//loop ++;
		//if (loop > 10000)
		//{	break;
		//}
	}

	infile.close();
	outfile.close();
	time_end = clock();
	cerr << "\nParse the pileup file and genotype calculation done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}

