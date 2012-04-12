
#include "soap_snp.h"
#include <getopt.h>
#include <ctime>
#include "FileListManager.h"
#include "MatrixManager.h"
#include "Readwin.h"
#include "Call_winManager.h"
#include "accessControl.h"
#include "SfsMethod.h"
using namespace std;
#include <istream>
#include <sys/mman.h>

#define VERSIONSTR "v2.0"

int usage() {

	cerr<<"		****************************************************************			"<<endl;
	cerr<<"		*		  Software name :SOAPsnp                       *			"<<endl;
	cerr<<"		*	     	        Version : v 2.0                        *			"<<endl;
	cerr<<"		*		    Last Update : 2010.10.15		       *			"<<endl;
	cerr<<"		*		         Author : Bill Tang                    *			"<<endl;
	cerr<<"		*		         E_mail : tangzhoubiao@genomics.org.cn *			"<<endl;
	cerr<<"		*				  zhouguyue@genomics.org.cn    *			"<<endl;
	cerr<<"		*		      Copyright : BGI. All Rights Reserved.    *			"<<endl;
	cerr<<"		****************************************************************			"<<endl;
	cerr<<"`Compulsory Parameters:"<<endl;
	cerr<<"		-i <FILE> Input SORTED Soap Result"<<endl;
	cerr<<"		-S <FILE> Input SORTED SAM Result"<<endl;
	cerr<<"		-B <FILE> Input SORTED BAM Result"<<endl;
	cerr<<"		-l if the Soap is filelist input"<<endl;
	cerr<<"		-d <FILE> Reference Sequence in fasta format"<<endl;
	cerr<<"		-o <FILE> (<DIR> if input is file list, if want to call sfs, you don't need to set this option, but you want call both, you can set this option) Output consensus file Optional Parameters:(Default in [])"<<endl;
	cerr<<"		-z <Char> ASCII chracter standing for quality==0 [@]"<<endl;
	cerr<<"		-g <Double> Global Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.9]"<<endl;
	cerr<<"		-p <Double> PCR Error Dependency Coefficient, 0.0(complete  dependent)~1.0(complete independent)[0.5]"<<endl;
	cerr<<"		-r <Double> novel altHOM prior probability [0.0005]"<<endl;
	cerr<<"		-e <Double> novel HET prior probability [0.0010]"<<endl;
	cerr<<"		-t set transition/transversion ratio to 2:1 in prior probability"<<endl;
	cerr<<"		-s <FILE> Pre-formated dbSNP information"<<endl;
	cerr<<"		-2 specify this option will REFINE SNPs using dbSNPs information [Off]"<<endl;
	cerr<<"		-a <Double> Validated HET prior, if no allele frequency known [0.1]"<<endl;
	cerr<<"		-b <Double> Validated altHOM prior, if no allele frequency known[0.05]"<<endl;
	cerr<<"		-j <Double> Unvalidated HET prior, if no allele frequency known [0.02]"<<endl;
	cerr<<"		-k <Double> Unvalidated altHOM rate, if no allele frequency known[0.01]"<<endl;
	cerr<<"		-u Enable rank sum test to give HET further penalty for better accuracy.[Off]"<<endl;
	//cerr<<"-n Enable binomial probability calculation to give HET for better accuracy. [Off]"<<endl;
	cerr<<"		-m Enable monoploid calling mode, this will ensure all consensus as HOM  and you probably should SPECIFY"<<endl<<
	      "		   higher altHOM rate. [Off]"<<endl;	
	cerr<<"		-q Only output potential SNPs. Useful in Text output mode. [Off]"<<endl;
	cerr<<"		-M <FILE>(<DIR> if set -l parameter) Output the quality calibration matrix;the matrix can be reused with -I "<<endl;
	cerr<<"	           if you rerun the program"<<endl;
	cerr<<"		-I <FILE>(<FILELIST> if set -l parameter) Input previous quality calibration matrix. It cannot be used "<<endl<<
	      "		   simutaneously with -M ."<<endl;
	cerr<<"		-L <short> maximum length of read [45] ."<<endl;
	cerr<<"		-Q <short> maximum FASTQ quality score [40] ."<<endl;
//	cerr<<"		-F <int> Output format. 0: Text; 1: GLFv2; 2: GPFv2.[0]"<<endl;
	cerr<<"		-E <String> Extra headers EXCEPT CHROMOSOME FIELD specified in GLFv2 output. "<<endl<<
	      "	           Format is \"TypeName1:DataName1:TypeName2:DataName2\"[""] ."<<endl;
	cerr<<"		-T <FILE> Only call consensus on regions specified in FILE. Format: ChrName\\tStart\\tEnd."<<endl;
	cerr <<"		-C <int> The CPU number which you want to set(It only work with -l).[4]" << endl;
	cerr << "		-f <FILE> the sfs parameter file. The file format: (must be used with -l parameter)" << endl << endl;
	cerr << "			doBay	0/1	" << endl;
	cerr << "			doJoint	0/1		" << endl;
	cerr << "			writeFr 0/1			" << endl;
	cerr << "			outfiles	<FILE PATH>	" << endl;
//	cerr << "			start	position	" << endl;
//	cerr << "			stop	position	" << endl;
	cerr << "			underFlowProtect	0/1	" << endl;
	cerr << "			allowPhatZero	0/1	" << endl;
	cerr << "			bias	[0,1]	" << endl;
	cerr << "			pvar	[0,1]	" << endl;
	cerr << "			eps	[0,1]	" << endl;
	cerr << "			pThres 	[0,1]	" << endl;
	cerr << "			sigleMinorJoint	0/1" << endl;
	cerr << "			sigleMinorBay	0/1" << endl << endl;

	//cerr<<"-S <FILE> Output summary of consensus"<<endl;
	cerr << "		-h Display this help" << endl;
	exit(1);
	return 0;
}

int readme() {
	return usage();
}



int main ( int argc, char * argv[]) {
	AccessControl("/share/backup/luoruibang/accessControl/SOAPsnp-v2", VERSIONSTR);

	// This part is the default values of all parameters
	Parameter *para = new Parameter();
	std::string alignment_name, consensus_name(""),soapresult_name;
	std::string matrix_file("");
	std::string sfs_path("");
	bool is_matrix_in = false; // Generate the matrix or just read it?
	bool list = false;
	FileListManager fileListManager;
	SfsMethod sfsMethod; 
	int c;
	Files *files = new Files();
	int ret;
	char in_mode[5] = {0};
	int CPU = 4;

	while((c=getopt(argc,argv,"i:d:o:z:g:p:r:e:ts:2a:b:j:k:unmqM:I:L:Q:S:B:F:E:T:h:f:C:l")) != -1) {
		switch(c) {
			case 'l':
			{
				/*alignment is a file list*/
				list = true;
				break;
			}	
			case 'i':
			{
				// Soap Alignment Result
				/*files->soap_result.clear();
				//files->soap_result.open(optarg);
				//if( ! files->soap_result) {
				cerr<<"No such file or directory:"<<optarg<<endl;
				exit(1);
				}*/
				alignment_name = optarg;
				break;
			}
			case 'd':
			{
				// The reference genome in fasta format
				files->ref_seq.clear();
				files->ref_seq.open(optarg);
				if( ! files->ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files->ref_seq.clear();
				break;
			}
			case 'o':
			{
				consensus_name = optarg;
				break;
			}
			case 'z':
			{
				// The char stands for quality==0 in fastq format
				para->q_min = optarg[0];
				if(para->q_min == 33) {
					clog<<"Standard Fastq System Set"<<endl;
				}
				else if(para->q_min == 64) {
					clog<<"Illumina Fastq System Set"<<endl;
				}
				else {
					clog<<"Other types of Fastq files?? Are you sure?"<<endl;
				}
				para->q_max = para->q_min + 40;
				break;
			}
			case 'g':
			{
				para->global_dependency= log10(atof(optarg));
				break;
			}
			case 'p':
			{
				para->pcr_dependency= log10(atof(optarg));
				break;
			}
			case 'r':
			{
				para->althom_novel_r = atof(optarg);
				break;
			}
			case 'e':
			{
				para->het_novel_r=atof(optarg);
				break;
			}
			case 't':
			{
				para->transition_dominant=true;
				break;
			}
			case 's':
			{
				// Optional: A pre-formated dbSNP table
				files->dbsnp.clear();
				files->dbsnp.open(optarg);
				if( ! files->ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files->dbsnp.clear();
				break;
			}
			case '2':
			{
				// Refine prior probability based on dbSNP information
				para->refine_mode = true;
				break;
			}
			case 'a':
			{
				para->althom_val_r=atof(optarg);
				break;
			}
			case 'b':
			{
				para->het_val_r=atof(optarg);
				break;
			}
			case 'j':
			{
				para->althom_unval_r=atof(optarg);
				break;
			}
			case 'k':
			{
				para->het_unval_r=atof(optarg);
				break;
			}
			case 'u':
			{
				para->rank_sum_mode = true;
				break;
			}
			case 'n':
			{
				para->binom_mode = true;
				break;
			}
		 	case 'm':
			{
				para->is_monoploid=1;
				break;
			}
			case 'q':
			{
				para->is_snp_only=1;
				break;
			}
			case 'M':
			{
				matrix_file = optarg;
				
				break;
			}
			case 'I':
			{
				matrix_file = optarg;
				
				is_matrix_in = true;
				break;
			}
			case 'S':
			{
				//files->summary.open(optarg);
				//// Output the summary of consensus
				//if( ! files->summary ) {
				//	cerr<<"No such file or directory: "<<optarg<<endl;
				//	exit(1);
				//}
				alignment_name = optarg;
				in_mode[0] = 'r';
				break;
			}
			case 'B':
			{
				alignment_name = optarg;
				in_mode[0] = 'r';
				in_mode[1] = 'b';
				break;
			}
			case 'L':
			{
				para->read_length = atoi(optarg);
				break;
			}
			case 'Q':
			{
				para->q_max = optarg[0];
				if(para->q_max < para->q_min) {
					cerr<< "FASTQ quality character error: Q_MAX > Q_MIN" <<endl;
				}
				break;
			}
			case 'F': {
				para->glf_format = atoi(optarg);
				break;
			}
			case 'E': {
				para->glf_header = optarg;
				break;
			}
			case 'T': {
				files->region.clear();
				files->region.open(optarg);
				files->region.clear();
				para->region_only = true;
				break;
			}
			case 'f': {
				sfs_path = optarg;
				break;
					  }
			case 'C': {
				CPU = atoi(optarg);
				break;
					  }
			case 'h':readme();break;
			case '?':usage();break;
			default: cerr<<"Unknown error in command line parameters"<<endl;
		}
	}
	
	if (alignment_name == "")
	{
		usage();
		exit(1);
	}
	if (sfs_path != "")
	{
		if (list != true)
		{
			cerr << "ERROR: -f must be used with -l" << endl;
			exit(0);
		}
	}
	
	if (CPU < 1)
	{
		cerr << "\tERROR: You have set the wrong CPU number " << CPU << endl
			 << "\tUse -h for more information." << endl;
		exit(0);
	}

	/*if is file list,open file list*/
	if(list == true)                            
	{
		ret = fileListManager.openAli(alignment_name,in_mode);
		if (BASEFILE_ERROR == ret)
		{
			cerr << "ERROR: File list " << alignment_name << "cannot open" << endl;
			exit(1);
		}
		else if (SOAPFILE_ERROR == ret)
		{
			cerr << "ERROR: FILE OPEN" << endl;
			exit(1);
		}
		
		// if need to call sfs, then init the para->sfs
		if (sfs_path != "")
		{
			para->sfs = sfs_init();
			int ret = init_sfs_para(sfs_path, para->sfs);
			if (ret == SFS_PARA_ERROR)
			{
				exit(1);
			}
			// init the SfsMethod object.
			sfsMethod.init(fileListManager.getFileNum());
			
			ret = files->OpenSfsfile(para->sfs->outfiles, para->sfs->writeFr, para->sfs->doBay, para->sfs->doJoint);
			if (ret == -1)
			{
				cerr << "ERROR: FILE OPEN" << para->sfs->outfiles << endl;
				exit(0);
			}
		}
		
	}
	else 
	{	
		/*if it is sam\bam open sam_result*/
		if(in_mode[0] == 'r')
		{
			files->sam_result.open(alignment_name.c_str(), in_mode);
		}
		/*open a soap file*/
		else 
		{
			files->soap_result.clear();
			files->soap_result.open(alignment_name.c_str());
			if( !files->soap_result) 
			{
				cerr<<"No such file or directory:"<<alignment_name.c_str()<<endl;
				exit(1);
			}
		}

		if (matrix_file != "")
		{
			if (is_matrix_in)
			{
				// Input the calibration matrix
				files->matrix_file.open(matrix_file.c_str(), fstream::in);
				if( ! files->matrix_file) {
					cerr<<"No such file or directory:"<<matrix_file<<endl;
					exit(1);
				}
				files->matrix_file.clear();
			}
			else 
			{
				// Output the calibration matrix
				files->matrix_file.open(matrix_file.c_str(), fstream::out);
				if( ! files->matrix_file) {
					cerr<<"Cannot creat file :"<<matrix_file<<endl;
					exit(1);
				}
				files->matrix_file.clear();
			}
		}
	}

	if((consensus_name == "" && sfs_path == "") || !files->ref_seq || (!files->sam_result.isOpened() && !files->soap_result && !fileListManager.getFileNum())) 
	{
		// These are compulsory parameters
		usage();
	}

	/*if it is noraml soapsnp format open,else change to binary and open */
	if ( ! para->glf_format ) {
		if(list == true)
		{
			ret = fileListManager.openCnsFile(alignment_name, consensus_name, "");
			if (CNSFILE_ERROR == ret)
			{
				cerr << "ERROR: cannot create consensus file" << endl;
				exit(1);
			}
		}
		else
		{
			// Normal SOAPsnp tab-delimited text format
			files->consensus.clear();
			files->consensus.open(consensus_name.c_str());
			if( ! files->consensus ) {
				cerr<<"Cannot creat file:" <<consensus_name<<endl;
				exit(1);
			}
			files->consensus.clear();
		}
	}
	else {
		// SOAPsnp-defined GLF and baseinfo format
		files->consensus.clear();
		files->consensus.open(consensus_name.c_str(), ios::binary);
		if(!files->consensus) {
			cerr<<"Cannot creat file:" <<consensus_name<<endl;
			exit(1);
		}
		files->consensus.clear();

		/*ouput consensus_name.baseinfo*/
		files->baseinfo.clear();
		string baseinfo_name = consensus_name + ".baseinfo";
		files->baseinfo.open(baseinfo_name.c_str());
		if(!files->baseinfo) 
		{
			cerr<<"Cannot creat file:" <<baseinfo_name<<endl;
			exit(1);
		}
		files->baseinfo.clear();

		/*ouput consensus_name.index*/
		files->o_region.clear();
		string o_region_name = consensus_name + ".index";
		files->o_region.open(o_region_name.c_str());
		if(!files->o_region) 
		{
			cerr<<"Cannot creat file:" <<o_region_name<<endl;
			exit(1);
		}
		files->o_region.clear();
	}

	time_t timep;
    time(&timep);
    clog << "genome start " << ctime(&timep);
	//Read the chromosomes into memory
	Genome * genome = new Genome(files->ref_seq, files->dbsnp);
	files->ref_seq.close();
	files->dbsnp.close();
	clog<<"Reading Chromosome and dbSNP information Done."<<endl;
	if(para->region_only && files->region) {
		genome->read_region(files->region, para);
		clog<<"Read target region done."<<endl;
	}

	time(&timep);
    clog << "genome end " << ctime(&timep);
	
	int file_num = fileListManager.getFileNum();

	if (list == true)
	{
		time(&timep);
		clog << ctime(&timep);

		MatrixManager matrixManager;
		Call_winManager call_winManager(para->read_length, file_num, genome);
		vector<Readwin> readwin_vec;
		ThreadManager threadManager(CPU);

		if (is_matrix_in)
		{
			// open matrix files to input.
			fileListManager.openMatrixFile(matrix_file, ios::in);
			for (int i = 0; i < fileListManager.getFileNum(); ++i)
			{
				if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager.getMatrixFile(i)), para))
				{
					cerr << "ERROR: add matrix failed!" << endl;
					exit(0);
				}
				Readwin new_readwin;
				readwin_vec.push_back(new_readwin);
			}
		}
		else 
		{
			matrixManager.setMatrixNum(file_num);
			if (matrix_file != "")
			{
				// open matrix file to output.
				fileListManager.openMatrixFile(matrix_file, ios::out, alignment_name);
			}

			if (in_mode[0] == 'r')
			{
				for (int i = 0; i < fileListManager.getFileNum(); ++i)
				{
					time(&timep);
					clog << ctime(&timep);

					MATRIX_ARGS *tmp_args = new MATRIX_ARGS(&matrixManager, para, genome, 0,
															fileListManager.getMatrixFile(i),
															fileListManager.getSamFile(i),i);
					// add the data to ThreadManager.
					threadManager.AddThread(_matrixManager_addMatrix, tmp_args);
					
					/*fstream *mat_out_p = fileListManager.getMatrixFile(i);

					if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager.getSamFile(i)), para, genome, mat_out_p))
					{
						cerr << "ERROR: add matrix failed!" << endl;
						exit(0);
					}*/

					Readwin new_readwin;
					readwin_vec.push_back(new_readwin);
				}
			}
			else
			{
				for (int i = 0; i < fileListManager.getFileNum(); ++i)
				{
					MATRIX_ARGS *tmp_args = new MATRIX_ARGS(&matrixManager, para, genome,
							fileListManager.getSoapFile(i),
							fileListManager.getMatrixFile(i),
							0,i);
					// add the data to ThreadManager.
					threadManager.AddThread(_matrixManager_addMatrix, tmp_args);

					/*fstream *mat_out_p = fileListManager.getMatrixFile(i);
					if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager.getSoapFile(i)), para, genome, mat_out_p))
					{
						cerr << "ERROR: add matrix failed!" << endl;
						exit(0);
					}*/

					Readwin new_readwin;
					readwin_vec.push_back(new_readwin);
				}
			}

			threadManager.Run();
			threadManager.Reset();
		}
		clog << "Correction Matrix Done!" << endl;

		fileListManager.closeAliFiles();
		fileListManager.openAli(alignment_name, in_mode);

		int count = 0;
		int sub_count = 0;
		time(&timep);
		clog << "begin to call !!!!!  time : " <<  ctime(&timep);
		CThreadPool threadpool(CPU);
		CALL_WIN_ARGS * call_win_args = NULL;
		Call_win_Task * cw_Task;
		// two vector used to release resource.
		vector<CALL_WIN_ARGS*> call_win_args_vec;
		vector<Call_win_Task*> cw_Task_vec;
	
		do 
		{
			ret = fileListManager.readWin(readwin_vec);
			if (COUNT_ERROR == ret)
			{
				cerr << "ERROR: something wrong with file number!" << endl;
				exit(0);
			}
			else if (INPUT_ERROR == ret)
			{
				cerr << "ERROR: something wrong with input files!" << endl;
				exit(0);
			}
			else if (COME_NEW_CHR == ret)
			{
				for (int i = 0; i < fileListManager.getFileNum(); ++i)
				{
					call_winManager.soap2cns(readwin_vec[i].getReadwin(), *(fileListManager.getCnsFile(i)), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), sfsMethod);
					readwin_vec[i].winChange();
					call_winManager.soap2cns(readwin_vec[i].getReadwin(), *(fileListManager.getCnsFile(i)), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), sfsMethod);
					readwin_vec[i].reset();
				}
			}
			for (int i = 0; i < fileListManager.getFileNum(); ++i)
			{
				call_win_args = new CALL_WIN_ARGS(&call_winManager, &(readwin_vec[i].getReadwin()), fileListManager.getCnsFile(i), &(files->baseinfo), genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), &sfsMethod);
				call_win_args_vec.push_back(call_win_args);
				cw_Task = new Call_win_Task();
				cw_Task->SetData(call_win_args);
				cw_Task_vec.push_back(cw_Task);
				threadpool.AddTask(cw_Task);
				//call_winManager.soap2cns(readwin_vec[i].getReadwin(), *(fileListManager.getCnsFile(i)), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), sfsMethod);
			}

			while (threadpool.getTaskSize() != 0)
			{
				usleep(1000);
			}

			for (int j = 0; j < cw_Task_vec.size(); ++j)
			{
				delete call_win_args_vec[j];
				delete cw_Task_vec[j];
			}
			call_win_args_vec.clear();
			cw_Task_vec.clear();

			if (para->sfs != NULL)
			{
				// do sfs
				sfsMethod.call_SFS(para, files);
			}
		} while(ret != FILE_END);// end of do

		time(&timep);
		clog << "Main:" << "end of call!!!" << ctime(&timep);
		
		for (int i = 0; i < fileListManager.getFileNum(); ++i)
		{
			call_winManager.dealTail(*(fileListManager.getCnsFile(i)), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, sfsMethod);
		}

		if (para->sfs != NULL)
		{
			// do sfs
			sfsMethod.call_SFS(para, files);
		}

		threadpool.StopAll();
		fileListManager.closeAliFiles();
		fileListManager.closeCnsFiles();
	}
	else
	{
		/***********************************/
		Prob_matrix * mat = new Prob_matrix;
		if( ! is_matrix_in) {
			//Read the soap result and give the calibration matrix
			if(in_mode[0] == 'r') {
				mat->matrix_gen(files->sam_result, para, genome);
			} else {
				mat->matrix_gen(files->soap_result, para, genome);
			}
			if (files->matrix_file) {
				mat->matrix_write(files->matrix_file, para);
			}
		}
		else {
			mat->matrix_read(files->matrix_file, para);
		}
		files->matrix_file.close();
		clog<<"Correction Matrix Done!"<<endl;
		mat->prior_gen(para);
		mat->rank_table_gen();
		Call_win info(para->read_length);
		info.initialize(0);

		//Call the consensus
		if (in_mode[0] != 'r')
		{
			files->soap_result.close();
			files->soap_result.clear();
			files->soap_result.open(alignment_name.c_str());
			files->soap_result.clear();
			info.soap2cns(files->soap_result, files->consensus, files->baseinfo, genome, mat, para, sfsMethod, 0);
			files->soap_result.close();
			files->consensus.close();
		} 
		else 
		{
			files->sam_result.close();
			files->sam_result.open(alignment_name.c_str(), in_mode);
			info.soap2cns(files->sam_result, files->consensus, files->baseinfo, genome, mat, para, sfsMethod, 0);
			files->sam_result.close();
			files->consensus.close();
		}
		/***********************************/
	}
	delete files;
	delete para;
	cerr<<"Consensus Done!"<<endl;
	return 0;
}

