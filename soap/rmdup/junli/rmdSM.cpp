#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <functional>

#include <unistd.h>
#include <cstring>
#define DEF_ONE  10000000
using namespace std;

#define USAGE \
	"--------------------------------------------------------------------------------------------------------\n"\
	" Program: rmdSM   remove duplication form soap result,then sort and merge by chr \n"\
	" Vision:  v1.0.0 low memory use, and correct the mismatch numbers of soap2 \n"\
	" Contact:         lijun3@genomics.org.cn\n\n"\
	" Usage:[option]\n\n"\
	"       -in_list  <str>   input list of soap result files, support read gzip file \"*.soap,*.single\"\n"\
	"                         only reads pairs in *.soap files used to remove duplication\n"\
	"                         list file Format: soap_result_file_name [lane_id] \n"\
	"                         if only give the soap_result_file_name, do not change the soap reasult,\n"\
	"                         just sort,merge and remove duplication reads \n"\
	"       -chr      <str>   input list of chr name  \n"\
	"       -out_dir  <str>   the directory of reslut \n"\
	"       -dp               out put the duplication reads which have been removed in *.dp file\n"\
	"       -prefix   <str>   prefix of output file,such as \"sample1\"\n\n"\
	" memory prediction:      1. (reads number(75b,of all files in list) / chr number)*0.35 G/M\n"\
	"                         2. 1M reads(75b,of the best covered chromosome) needs 0.35g\n"\
	"---------------------------------------------------------------------------------------------------------\n"\


typedef unsigned long long _uint64;

typedef vector<string> Vstring;

typedef struct PosInfo
{
	_uint64 rOne_pos;
	_uint64 rTwo_pos;
	int total_len;
	int used;
	string chrname;
	_uint64 readsIndex;

}posInfo;

typedef struct Reads
{
	_uint64 pos;
	int len;
	string new_id;
	bool stat;
	string info;

}reads;

typedef struct pe{
	string id;
	string ser;
	int length;
	string chrname;
	_uint64 pos;
	int hit_type;
	
}PE;

typedef struct ReadsCount{
	_uint64 reads;
	_uint64 bases;
}readsCount;




bool less_read(const reads & m1, const reads & m2) {
        return m1.pos <  m2.pos;
}


bool less_two(const posInfo & m1, const posInfo & m2) {
	return (m1.rOne_pos <  m2.rOne_pos)|| ((m1.rOne_pos ==  m2.rOne_pos) && (m1.rTwo_pos <  m2.rTwo_pos));
}

////////////

string int2str(int i)
{
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}

string uint_642str(_uint64 i)
{
	string s;
	  stringstream ss(s);
	 ss << i;
	 return ss.str();

}

/////////////////


typedef struct LIST
{
        string file;
	int nu;

}list;
typedef vector<list> Vlist;

typedef vector<posInfo> Vposinfo;

typedef map<string, Vposinfo> MVposinfo;

typedef map<string,string> MCMD;

typedef vector<string> Vchr;
typedef map<string,_uint64> MCHRSIZE;


typedef vector<reads> VChrReads;

MCMD mCmd;
MCHRSIZE mChrSize;

string outDir;
string LibName;
bool control = true;
Vchr vChr;
string suff = "";

bool outdp = false;


/////////////////////////////////////////////////

bool readList(Vlist &vs,string listFile)
{

	ifstream ifs(listFile.c_str(),ifstream::in);

        if(!ifs.good())
        {
                cerr << "list of soap files is wrong,please check : " << listFile<< endl;
                return false;
        }
	cerr << "soap result files: " << endl;
	while(!ifs.eof())
        {
			string line;
			getline(ifs,line);

		if(line.length() > 0)
		{

			istringstream isone (line,istringstream::in);

			string file;
			int nu=-1; 
			list ls;

			isone >> file >> nu; 
			ls.file = file;
			ls.nu = nu;

			vs.push_back(ls);

			cerr << ls.file << " " << ls.nu << endl; 
		}
	
	}

	ifs.close();

	return true;	


}



bool readChr(Vstring &vs,string listFile)
{

	ifstream ifs(listFile.c_str(),ifstream::in);

        if(!ifs.good())
        {
                cerr << "list of chr is wrong,please check : " << listFile<< endl;
                return false;

        }

	while(!ifs.eof())
        {
		string line;
		ifs >> line;

		if(line.length() > 0)
		{
			vs.push_back(line);
			
			cerr << "chr : " << line << endl;
					
		}
	
	}

	ifs.close();

	return true;	


}

_uint64  readpeNew(readsCount &sum,string file , MVposinfo & mvpi,string chr,VChrReads &vChrReads,_uint64 &Index,bool filter,string nu,bool isrmd)
{
	
	ifstream ifs(file.c_str(),ifstream::in);
	cerr << "read: "<<file << " IS rm duplication : " << isrmd << endl;
	if(!ifs.good())
	{
		cerr << "open soap file: " << file  << " error!" << endl;
		return false;

	}
	
	_uint64 totalreads=0;
	_uint64 totalbases=0;
	PE beforReads;
	beforReads.id="first";
	reads line;

	while(!ifs.eof())
	{
		string readsOne;

		getline(ifs,readsOne);

		 string id,ser,quality,ab,fr,chrname,end;
                 _uint64 pos,s_pos;
                 int hit_type,length;
                 int nhits;

                istringstream isone (readsOne,istringstream::in);
 
               	isone >>id >> ser >> quality >> nhits >> ab >> length >> fr >> chrname >> pos >> hit_type >> end;
		int trim= 0;
		string s_end = end;

		if(chrname == chr)

		{

			size_t found = readsOne.find_last_of("\t");
	
			if(found == string::npos)
       			 {

                        // error line ,not split by table !

       			 }
			string last_seq= readsOne.substr(found+1);
			size_t last_pos = last_seq.find_first_of("ATGC");
			int mismatch =0;

		  	while (last_pos!=string::npos) // count mismatch number
  			{
    			
				mismatch++;
	
 		   		last_pos =last_seq.find_first_of("ATGC",last_pos+1);
	  		}

			string setion;


			s_pos = pos;

			line.pos = pos;
			line.stat= false;
			line.len = length;

			string OneId = id.substr(0,id.find_last_of("/"));
			string TwoId = beforReads.id.substr(0,beforReads.id.find_last_of("/"));


			if(s_end.find('S') != string::npos) {   // to define is trim or not
                                        trim =1;
			}

                        size_t m = s_end.find_first_of('M');

                        if(fr == "+") {         // get reads position in 5' before trim

                                                size_t s = s_end.find_first_of('S');

                                                if(m != string::npos && s != string::npos&&s < m) {	// trim 5'

                                                        string dis = s_end.substr(0,s);
                                                        int d = atoi(dis.c_str());
                                                        s_pos = pos -d;

                                                 }


                         }
                         else {                  // get reads position in 5' before trim

                                                size_t s = s_end.find_last_of('S');
						
                                                if(m != string::npos && s != string::npos && s>m)	// trim 5'
                                                {

                                                        string dis = s_end.substr(m+1,s-m-1);

                                                        int d = atoi(dis.c_str());
                                                        s_pos = pos + length + d -1;

                                                }
						else {	// trim 3'

							s_pos = pos + length - 1;
		
						}

                        }

			if(filter)
			{
		
				end = int2str(trim);
		
				string new_id = nu;
				new_id +="-";
				new_id +=OneId;

				line.new_id = new_id;
				line.new_id.reserve(new_id.size());

				setion += ser;	setion += "\t";
				setion += quality;	setion += "\t";
				setion += int2str(nhits); setion += "\t";
				setion += ab; setion += "\t";
				setion += int2str(length); setion += "\t";
				setion += fr;	setion += "\t";
				setion += chrname; setion += "\t";

				setion += uint_642str(pos); 
				setion += "\t";
				setion += int2str(mismatch);
				setion += "\t";
				setion += end;

				line.info = setion;

				line.info.reserve(setion.size());

			}
			else
                        {

                                line.new_id = id;
                                line.new_id.reserve(id.size());

				string setion2 = readsOne.substr(readsOne.find_first_of("\t")+1);
			
				line.info = setion2;	
				line.info.reserve(setion.size());

                        }    
			

			vChrReads.push_back(line);

	 		if((isrmd && beforReads.chrname == chrname && OneId == TwoId ))
			{
				posInfo pinfo;
				pinfo.chrname=chrname;
				pinfo.chrname.reserve(chrname.size());
				
				if(s_pos < beforReads.pos)
				{
					pinfo.rOne_pos = s_pos;
					pinfo.rTwo_pos = beforReads.pos;

					pinfo.total_len = length + beforReads.length;

				}

				else
				{
					pinfo.rOne_pos = beforReads.pos;
                 		        pinfo.rTwo_pos = s_pos;

					pinfo.total_len = length + beforReads.length;

				}
					pinfo.used = 1;
					pinfo.readsIndex = Index;

				MVposinfo::iterator it = mvpi.find(pinfo.chrname);

				if(mvpi.empty() || it == mvpi.end())
				{
					Vposinfo vpi;
				 	vpi.push_back(pinfo);	
					mvpi.insert(pair<string,Vposinfo>(pinfo.chrname,vpi));	
					vChr.push_back(pinfo.chrname);

				}
				else
				{
			
					it->second.push_back(pinfo);

				}

			}

		     	Index++;	
			totalreads++;
			totalbases += length;
			beforReads.id = id;
			beforReads.pos = s_pos;
			beforReads.length = length;
			beforReads.chrname = chrname;

		}

	}

	ifs.close();
	cerr << "index: " << Index << " " << totalreads << endl;
	sum.reads = totalreads;

	sum.bases = totalbases;
	return totalreads;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
///	out put result

bool  outPutMerge(string file,VChrReads &vchr)
{

	string outfile = file;
	string outdbfile = file;
	outdbfile += ".dp";

	_uint64 num_reads = 0;
	_uint64 num_bases = 0;
	_uint64 nu_duplication=0;

cerr << "Output: " << outfile << endl;

	ofstream out_fs(outfile.c_str());
	ofstream out_dp;

	if(outdp) {

		out_dp.open(outdbfile.c_str());
		if(!out_dp.good()) {
			cerr << " open outdpput file error ! " << outfile << endl;
			out_dp.close();
			return false;

		}

	}
	if(!out_fs.good())	{
		cerr << " open output file error ! " << outfile << endl;	
	
		out_fs.close();
		return false;

	}

	sort(vchr.begin(),vchr.end(),less_read);
cerr << "Size " << vchr.size() << endl;

	for(int i=0;i < vchr.size();i++)
	{
		if(! vchr.at(i).stat)
		{
			out_fs <<vchr.at(i).new_id << "\t" << vchr.at(i).info << endl;

			num_reads++;
			num_bases += vchr.at(i).len;

		}
		else {

			out_dp << vchr.at(i).new_id << "\t" << vchr.at(i).info << endl;

		}

	}
	
	cout << num_reads << "\t" << num_bases << endl;
	out_fs.close();

	if(outdp) {
		out_dp.close();
	}
	return true;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
////	delete duplication

_uint64 reset(MVposinfo& vpi,VChrReads &vChrReads)
{


	_uint64 totalline = 0;

	for(MVposinfo::iterator it =vpi.begin();it!=vpi.end();++it)
	{

	cerr << "size " << it->first << "\t"<< it->second.size() << "\t" << it->second.capacity() << "\t" << sizeof(it->second.at(1)) << endl;
	cerr << "vector" << "\t" << vChrReads.size() << "\t" << vChrReads.capacity() << "\t" << sizeof(vChrReads.at(1)) << endl;

		
		sort(it->second.begin(),it->second.end(),less_two);

		_uint64 rOpos = 0;
		_uint64 rTpos = 0;


		totalline += it->second.size();
	
		long long index =0;

		for(unsigned long long  i = 1; i<it->second.size(); ++i)
		{
		

			if(it->second.at(i).rOne_pos == it->second.at(i-1).rOne_pos && it->second.at(i).rTwo_pos == it->second.at(i-1).rTwo_pos)
			{
			
			
				if(it->second.at(i).total_len <= it->second.at(index).total_len)
                         	{
					it->second.at(i).used =0 ;
					vChrReads.at(it->second.at(i).readsIndex).stat=true;
					vChrReads.at(it->second.at(i).readsIndex-1).stat=true;
			 	}
                         	else
				{
					it->second.at(index).used = 0;
					vChrReads.at(it->second.at(i).readsIndex).stat=true;
					vChrReads.at(it->second.at(i).readsIndex-1).stat=true;
					index = i;

				}
	
			}

			else
			{
				index = i;


			}
		
		}

	}

	return totalline;

}


void doMerge(Vlist& soapList,string chr,string outdir,bool control)
{


	MVposinfo vpi;
	VChrReads vReads;

	string outFile = outdir;
	if(suff.length()!=0)
	{
		
		outFile += suff;
		outFile += ".";

	}
	outFile += chr;

        _uint64 Index =0;	
	_uint64 t_reads= 0;
	_uint64 t_bases = 0;

	_uint64 lastreads = 0;

cerr << "do chromosome " << chr << endl;

	for(int i=0; i < soapList.size();i++)
   	{


		bool filt =true;
		readsCount sum;
		bool isrmd = true;

		
		string nu = int2str(soapList.at(i).nu);
		if(soapList.at(i).file.find("single") != string::npos )
		{
			isrmd= false;

		}

		if(soapList.at(i).nu == -1 )
		{

			nu = int2str(i);
			filt= false;
			cerr << "IS filter: " << "no" << endl;
		}
		else
		{
	
			cerr << "IS filter: " <<  "yes " << endl; 
		}

		readpeNew(sum,soapList.at(i).file,vpi,chr,vReads,Index,filt,nu,isrmd);
		
		if(sum.reads == 0) {
			cerr << " may donot have chromosome " << chr << " in file: " << soapList.at(i).file << endl;

		}
		t_reads += sum.reads;
		t_bases += sum.bases;

	}

	reset(vpi,vReads);

cout << chr << "\t" << t_reads << "\t" << t_bases << "\t";

	{
		outPutMerge(outFile,vReads);
	}	
	

}




//////////////////////////////////////////////


int main(int argc ,char** argv)
{

	Vstring chrList;
	Vlist soapList;
	string fn;

	string soapListFile,chrListFile,cn;
	string cmd;

	--argc ; ++argv ;
	if(argc ==0) {
		cerr << USAGE;
		exit(0);

	}
	
	while(argc) {
		
		if( !strcmp(argv[0], "-in_list") && argc >=2) {
			soapListFile= argv[1];
			argc -=2; argv +=2;

			cerr << "soap_results: "<<  soapListFile << endl;

		}
		else if( !strcmp(argv[0],"-chr") && argc >=2) {
			chrListFile = argv[1];
                        argc -=2; argv +=2;
			cerr << "chr: " << chrListFile << endl;
			
		}
		else if( !strcmp(argv[0],"-out_dir") && argc >=2) {
                        outDir = argv[1];
                        argc -=2; argv +=2;
			cerr << "Out direcotry: " << outDir << endl;

                        
                }
		else if(!strcmp(argv[0],"-prefix") && argc >=2) {
                        suff = argv[1];
                        argc -=2; argv +=2;
                    	
			cerr << "prefix: " << suff << endl;    
                }
		else if(!strcmp(argv[0],"-dp") && argc >=1) {

			outdp = true;
			argc -=1; argv +=1;

		}
		else {

			cerr << USAGE;
			exit(0);
		}


	}
	
	if(!readList(soapList,soapListFile))
        {
                cerr << "read list error !" << endl;
                exit(0);
        }
        if(!readChr(chrList,chrListFile))
        {
                cerr << "read chr error !" << endl;
                exit(0);
        }

	outDir += "/";
	control = false;

	for(Vstring::iterator it = chrList.begin();it!=chrList.end();it++)
	{
					
		doMerge(soapList,*it,outDir,control);
				
	}

	

}
