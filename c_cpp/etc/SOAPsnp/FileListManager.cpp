/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-9
 *CLASS NAME: FileListManager
 *FUNCTION : a class that can control file list. support open, readline, close etc
 *			file functions.
 *FILE NAME : FileListManager.cpp
 *UPDATE DATE : 2010-8-19
 *UPDATE BY : Bill Tang
 *UPDATE DATE : 2010-9-3
 *UPDATE BY : Bill Tang
 *
 *UPDATE: 2010-9-3 add class member m_endfile_count, change function readWin, 
 *			readWinFromSoap and readWinFromSam. 
 *			Make it can judge the files number that reach tail.
 *UPDATE: 2010-10-4 change fucntion closeAliFiles, make it can close matrix files and delete all object!
 *UPDATE: 2010-10-8 change function openCnsFile: add the function that if don't need output consensus, create empty objects.
 *UPDATE: 2010-11-2 change the m_soap_file_vec's type, to vector<*igzstream>.
					change the interface which is related to soap file.
 *UPDATA BY :BIll Tang
 *UPDATE: 2010-11-15 add the sem's control function.
 *******************************************************************************
 */

#include "FileListManager.h"
#include "SamCtrl.h"
#include <sys/resource.h>


FileListManager::FileListManager()
	: m_stop_pos(0)
	, m_endfile_count(0)
{
	m_soap_file_vec.clear();
	m_sam_file_vec.clear();
	m_cns_file_vec.clear();
	m_file_number = 0;
	m_base_file.clear();
}



FileListManager::~FileListManager()
{
	/*delete f_output stream*/
	for (vector<gzoutstream*>::iterator iter = m_cns_file_vec.begin();
		 iter != m_cns_file_vec.end(); 
		 ++iter)
	{
		delete *iter;
	}

	/*delete f_input stream*/
	for (vector<igzstream*>::iterator iter = m_soap_file_vec.begin();
		 iter != m_soap_file_vec.end(); 
		 ++iter)
	{
		delete *iter;
	}


	for (vector<SamCtrl*>::iterator iter = m_sam_file_vec.begin();
		 iter != m_sam_file_vec.end(); 
		 ++iter)
	{
		delete *iter;
	}

	for (vector<fstream*>::iterator iter = m_matrix_file_vec.begin();
		 iter != m_matrix_file_vec.end(); 
		 ++iter)
	{
		delete *iter;
	}
}


/**
 * DATE: 2010-8-9
 * FUNCTION: open output file function.
 * PARAMETER: path: the output consensus files path. mode: the open mode.
 * RETURN: the number of consensus files that be opened.
 */
int FileListManager::openCnsFile(const string& path, const string& outdir, const char* mode, Parameter * para)
{
	string cns_Path, cns_name;         //the one consensus file name
	gzoutstream* f_output;

	// if don't need to out put consensus, then create empty objects.
	if (outdir == "")
	{
		for (int i = 0; i < m_file_number; ++i)
		{
			// create new ofstream object.update by zhukai on 2010-12-09
			if(para->is_cns_gz)
			{
				f_output = new myogzstream();
			}
			else
			{
				f_output = new myofstream();
			}
			m_cns_file_vec.push_back(f_output);
		}
		return m_cns_file_vec.size();
	}

	m_base_file.open(path.c_str());  //open filelist
	if (!m_base_file)         //can not open filelist
	{
		cerr << "Cannot open file:" << path << endl;
		return BASEFILE_ERROR;
	}

	string cmd = "mkdir -p " + outdir;
	system(cmd.c_str());

	/*open each cnsfile and push in cns vector*/
	while (getline(m_base_file, cns_Path))          
	{
		//f_output = new gzoutstream();   

		/*find cons name from cns path*/
		int pos = cns_Path.rfind('/');   
		
		//update by zhukai on 2010-12-09
		if(para->is_cns_gz)
		{
			f_output = new myogzstream();
			cns_name = outdir + cns_Path.substr(pos, string::npos) + ".consus.gz";
		}
		else
		{
			f_output = new myofstream();
			cns_name = outdir + cns_Path.substr(pos, string::npos) + ".consus";
		}

		
		
		/*put cns in to cns vector*/
		f_output->open(cns_name.c_str(),std::ios::out);
		
		if (!f_output->is_open()) 
		{
			cerr << "Cannot open file:" << cns_name << endl;
			return CNSFILE_ERROR;
		}

		// set the output format. 2010-12-15 Bill
		f_output->set_f(std::ios::showpoint);
		m_cns_file_vec.push_back(f_output);
	}

	m_base_file.close();                        //close consensus files list 

	return m_cns_file_vec.size();
}


/**
 * DATE: 2010-8-9
 * FUNCTION: open all the file in the file list
 * PARAMETER: listfilename: the path of file that contain list of infile. mode: the open mode.
 * RETURN: 
 */
int FileListManager::openAli(const string& listfilename,const char* mode)
{
	if (mode[0] == 'r')                       //is sam/bam file	
	{
		return openSamFile(listfilename,mode);            //open sam/bam file ,return file number
	}
	else 
	{
		//cerr << "open soap file" <<endl;
		return openSoapFile(listfilename);                   //open soap file , return file number
	}
}


/**
 * DATE: 2010-8-30
 * FUNCTION: open sam/bam file
 * PARAMETER: mode: the open mode.
 * RETURN: BASEFILE_ERROR when cannot open file list. SAMFILE_ERROR when cannot open sam file.
 */
int FileListManager::openSamFile(const string& listfilename, const char* mode)
{
	string filepath;
	SamCtrl *sctrl;

	m_base_file.open(listfilename.c_str());              //open sam/bam files list
	if (!m_base_file) 
	{
		cerr << "Cannot open file:" << listfilename << endl;
		return BASEFILE_ERROR;
	}
	
	/*open each sam/bam file in the filelist and input sam vector*/
	while (m_base_file >> filepath)
	{
		/*set max limit */
	//	if(setlimit(m_sam_file_vec.size()))
		{

			sctrl = new SamCtrl();
			if (!(sctrl->open(filepath.c_str(), mode)))
			{
				cerr << "Cannot open file:" << filepath << endl;
				return SAMFILE_ERROR;
			}
			m_sam_file_vec.push_back(sctrl);
		}
	//	else
		{
	//		return SETLIMIT_ERROR;
		}
	}
	m_file_number = m_sam_file_vec.size();   // bam/sam file number in file list
	m_base_file.close();                        //close sam/bam files list
	
	return m_file_number;
}


/**
 * DATE: 2010-8-30
 * FUNCTION: open soap format file.
 * PARAMETER: void.
 * RETURN: BASEFILE_ERROR when cannot open file list. SAMFILE_ERROR when cannot open soap file.
 */
int FileListManager::openSoapFile(const string& listfilename)
{
	string filepath;
	igzstream* f_input; 

	m_base_file.open(listfilename.c_str());              //open soap filelist
	if (!m_base_file) 
	{
		cerr<<"Cannot open file:" <<listfilename<<endl;
		return BASEFILE_ERROR;
	}
	/*open each soap file in the filelist and input soap vector*/
	while (m_base_file >> filepath)
	{
		/*set max limit */
		//if ( setlimit( m_soap_file_vec.size() ) )
		{
			f_input = new igzstream();
			f_input->open(filepath.c_str()) ;
			if (!(*f_input)) 
			{
				cerr << "Cannot open file:" << filepath << endl;
				return SOAPFILE_ERROR;
			}
			m_soap_file_vec.push_back(f_input);
			f_input->clear();
		}
		//else
		{
		//	return SETLIMIT_ERROR;
		}
	}

	m_file_number = m_soap_file_vec.size();   // soap file number in file list	
	m_base_file.close();                    //close soap listfile
		
	return m_file_number;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: read a line from the *vec[index]
 * PARAMETER: line: the read buffer. index: the index of file handles vec.
 * RETURN: the length of string line.
 */
int FileListManager::readLine(string& line, const int index)
{
	return 0;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: close all the files
 * PARAMETER: void.
 * RETURN: void
 * UPDATE ON 2010-10-4
 *			delete the vector object and add delete m_matrix_file_vec.
 */
void FileListManager::closeAliFiles(void)
{
	int vec_size = 0;

	/*close soap files*/
	while (vec_size < m_soap_file_vec.size())
	{
		m_soap_file_vec[vec_size]->close();
		delete m_soap_file_vec[vec_size];
		vec_size++;
	}
	
	vec_size = 0;
	/*close sam/bam files*/
	while(vec_size < m_sam_file_vec.size())
	{
		m_sam_file_vec[vec_size]->close();
		delete m_sam_file_vec[vec_size];
		vec_size++;
	}

	vec_size = 0;
	/*close sam/bam files*/
	while(vec_size < m_matrix_file_vec.size())
	{
		m_matrix_file_vec[vec_size]->close();
		delete m_matrix_file_vec[vec_size];
		vec_size++;
	}

	m_soap_file_vec.clear();
	m_sam_file_vec.clear();
	m_matrix_file_vec.clear();
}

/**
 * DATE: 2010-8-9
 * FUNCTION: close all the files
 * PARAMETER: void.
 * RETURN: void
 */
void FileListManager::closeCnsFiles(void)
{
	int vec_size=0;
	/*close cns file*/
	while (vec_size < m_cns_file_vec.size())
	{
		m_cns_file_vec[vec_size]->close();	
		vec_size++;
	}
}



/**
 * DATE: 2010-8-9
 * FUNCTION: read lines to the readwin_vec. The lines' start position must inside a win size. 
 *				When the lines' end position exceed the win,put it to the next win.
 * PARAMETER: readwin_vec: the vector of Readwin objs. storing the info of every readwin.
 * RETURN: COUNT_ERROR when readwin_vec.size() not equal to the files' number.
 */
int FileListManager::readWin(vector<Readwin>& readwin_vec, Parameter * para)
{
	//CThreadPool threadpool(para->CPU);
	if (readwin_vec.size() != m_cns_file_vec.size())
	{
		// file size not right.
		return COUNT_ERROR;
	}

	// get signal pointers.
	sem_t * sem_read_p = &(para->sem_read);
	sem_t * sem_call_cns_p = &(para->sem_call_cns);
	sem_t * sem_readwin_return_p = &(para->sem_readwin_return);
	
	////CThreadPool threadpool(cpu);
	//READ_WIN_ARGS * read_win_args = NULL;
	//Read_win_Task * rw_Task = NULL;
	//// two vector used to release resource.
	//vector<READ_WIN_ARGS*> read_win_args_vec ;
	//vector<Read_win_Task*> rw_Task_vec ;
	int ret = 0;
	while(1)
	{
		// waiting for the signal.
		sem_wait(sem_read_p);
		sem_wait(sem_readwin_return_p);

		if (m_soap_file_vec.size() == readwin_vec.size())
		{
			for (int i = 0; i < readwin_vec.size(); ++i)
			{
				// read records from soap file to the win.
				ret = readWinFromSoap(readwin_vec[i], *m_soap_file_vec[i]);

				if (COME_NEW_CHR == ret)
				{
					para->ret = COME_NEW_CHR;
					goto next_cycle;
					//return COME_NEW_CHR;
				}
				else if (FILE_END == ret)
				{
					m_endfile_count ++;
				}
			}

	//		int *ret = new int[m_file_number];
	//		memset(ret, 0, sizeof(int) * m_file_number);
	//		for (int i = 0; i < readwin_vec.size(); ++i)
	//		{
	//			// read records from soap file to the win.
	//			read_win_args = new READ_WIN_ARGS(&readwin_vec[i], m_soap_file_vec[i],NULL, ret, i);
	//			read_win_args_vec.push_back(read_win_args);
	//			rw_Task = new Read_win_Task();
	//			rw_Task->SetData(read_win_args);
	//			rw_Task_vec.push_back(rw_Task);
	//			threadpool.AddTask(rw_Task);
	//		}
	//		while (threadpool.getTaskSize() != 0)
	//		{
	//			usleep(1);
	//		}
	//		//clear vector
	//		for (int j = 0; j < rw_Task_vec.size(); ++j)
	//		{
	//			delete read_win_args_vec[j];
	//			delete rw_Task_vec[j];
	//		}
	//		read_win_args_vec.clear();
	//		rw_Task_vec.clear();
	////cerr << __FUNCTION__ << __LINE__ << endl;

	//		for (int i = 0; i < readwin_vec.size(); ++i)
	//		{
	//			if (COME_NEW_CHR == ret[i])
	//			{
	//				para->ret = COME_NEW_CHR;
	//				//return COME_NEW_CHR;
	//			} 
	//			else if (FILE_END == ret[i])
	//			{
	//				m_endfile_count ++;
	//			}
	//		}
	////cerr << __FUNCTION__ << __LINE__ << endl;
	//		delete [] ret;
		}
		else if (m_sam_file_vec.size() == readwin_vec.size())
		{
			for (int i = 0; i < readwin_vec.size(); ++i)
			{
				// read records from soap file to the win.
				ret = readWinFromSam(readwin_vec[i], *m_sam_file_vec[i]);
				if (COME_NEW_CHR == ret)
				{
					para->ret = COME_NEW_CHR;
					goto next_cycle;
					//return COME_NEW_CHR;
				}
				else if (FILE_END == ret)
				{
					m_endfile_count ++;
				}
			}
			//int *ret = new int[m_file_number];
			//memset(ret, 0, sizeof(int) * m_file_number);
			//for (int i = 0; i < readwin_vec.size(); ++i)
			//{
			//	read_win_args = new READ_WIN_ARGS(&readwin_vec[i], NULL, m_sam_file_vec[i], ret, i);
			//	read_win_args_vec.push_back(read_win_args);
			//	rw_Task = new Read_win_Task();
			//	rw_Task->SetData(read_win_args);
			//	rw_Task_vec.push_back(rw_Task);
			//	threadpool.AddTask(rw_Task);
			//}
			//while (threadpool.getTaskSize() != 0)
			//{
			//	usleep(10);
			//}
			////clear vector
			//for (int j = 0; j < rw_Task_vec.size(); ++j)
			//{
			//	delete read_win_args_vec[j];
			//	delete rw_Task_vec[j];
			//}
			//read_win_args_vec.clear();
			//rw_Task_vec.clear();

			//for (int i = 0; i < readwin_vec.size(); ++i)
			//{
			//	if (COME_NEW_CHR == ret[i])
			//	{
			//		para->ret = COME_NEW_CHR;
			//		// return COME_NEW_CHR;
			//	} 
			//	else if (FILE_END == ret[i])
			//	{
			//		m_endfile_count ++;
			//	}
			//}
			//delete [] ret;
		} 
		else
		{
			// something wrong with input files.
			para->ret = INPUT_ERROR;
			//return INPUT_ERROR;
		}

		if (m_endfile_count != m_file_number)
		{
			para->ret = 0;
			//return 0;
		}
		else
		{
			// all file reach tail.
			para->ret = FILE_END;
			sem_post(sem_call_cns_p);
			break;
			//return FILE_END;
		}

		next_cycle:
		sem_post(sem_call_cns_p);
	}
	return 1;
}


// get the soap file handle pointer m_soap_file_vec[index]
/**
 * DATE: 2010-8-9
 * FUNCTION: get the soap file handle pointer m_soap_file_vec[index]
 * PARAMETER: index: the vector index.
 * RETURN: return soap file pointer m_soap_file_vec[index].failed return INDEX_OVER_FLOW
 */
igzstream * FileListManager::getSoapFile(const int index)
{
	if (index < 0 || index >= m_soap_file_vec.size())
	{
		// index overflow
		return INDEX_OVER_FLOW;
	}

	return m_soap_file_vec[index];
}


/**
 * DATE: 2010-8-9
 * FUNCTION: get the soap file handle pointer m_soap_file_vec[index]
 * PARAMETER: index: the vector index.
 * RETURN: return sam/bam file pointer m_soap_file_vec[index].failed return INDEX_OVER_FLOW
 */
SamCtrl* FileListManager::getSamFile(const int index)
{
	if (index < 0 || index >= m_sam_file_vec.size())
	{
		// index overflow
		return INDEX_OVER_FLOW;
	}

	return m_sam_file_vec[index];
}


/**
 * DATE: 2010-8-9
 * FUNCTION: read record from soapfile and send to the readwin.
 * PARAMETER:	readwin: the Readwin object that would be processed.
 *				soapfile: the input file handls.
 * RETURN: FILE_END when file is over else 0.
 */
int FileListManager::readWinFromSoap(Readwin& readwin, igzstream& soapfile)
{
	readwin.winChange();	// start from next win.
	std::string line;	// a temp string object.
	int ret;	// store the flag from function Readwin::addRead.

	
	// judge if the win can be added.
	if (readwin.isAbleToAdd())
	{
		// read record form file.
		while (getline(soapfile,line))
		{
			ret = readwin.addRead(line);

			if (ADD_SUCCESSFUL == ret)
			{
				// successful.
				continue;
			}
			else if (READ_EXCEED_WIN == ret)
			{
				// record exceed.
				break;
			}
			else if (READ_POS_ERROR == ret)
			{
				// something error with the record's position.
				cerr << "ERROR : Something wrong with the read's position with line:" << endl;
				cerr << line << endl;
				//break;
				continue;
			}
			else if (COME_NEW_CHR == ret)
			{
				if (readwin.getLastCount() == m_file_number)
				{
					// all the samples are reach a new chromosome.
					return COME_NEW_CHR;
				}
				break;
			}
		}	
		if (soapfile.eof())
		{
			// file is over.
			readwin.setLastPos(FILE_END);
			return FILE_END;
		}
	}
	return 0;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: read record from samCtrl and send to the readwin.
 * PARAMETER:	readwin: the Readwin object that would be processed.
 *				samCtrl: the input file handls.
 * RETURN: FILE_END when file is over else 0.
 */
int FileListManager::readWinFromSam(Readwin& readwin, SamCtrl& samCtrl)
{
	readwin.winChange();	// start from next win.
	std::string line;	// a temp string object.
	int ret;	// store the flag from function Readwin::addRead.
	int size = 0;	// store the readline flag.

	// judge if the win can be added.
	if (readwin.isAbleToAdd())
	{
		// read record form file.
		while ((size = samCtrl.readline(line)) != -1)
		{
			line = alignment_format(line);
			// add record to the Readwin object.
			ret = readwin.addRead(line);
			
			if (ADD_SUCCESSFUL == ret)
			{
				// successful.
				continue;
			}
			else if (READ_EXCEED_WIN == ret)
			{
				// record exceed.
				break;
			}
			else if (READ_POS_ERROR == ret)
			{
				// something error with the record's position.
				cerr << "ERROR : Something wrong with the read's position with line:" << endl;
				cerr << line << endl;
				//break;
				continue;
			}
			else if (COME_NEW_CHR == ret)
			{
				if (readwin.getLastCount() == m_file_number)
				{
					// all the samples are reach a new chromosome.
					return COME_NEW_CHR;
				}
				break;
			}
		}

		if (size == -1)
		{
			// file is over.
			readwin.setLastPos(FILE_END);
			return FILE_END;
		}
	}
	

	return 0;
}


// get the file number
/**
 * DATE: 2010-8-9
 * FUNCTION: get the file number
 * PARAMETER:	
 * RETURN: m_file_number
 */
int FileListManager::getFileNum(void)
{
	return m_file_number;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: get the consensus file handle pointer m_cns_file_vec[index]
 * PARAMETER: index: the vector index.
 * RETURN: return consensus file pointer m_cns_file_vec[index].failed return INDEX_OVER_FLOW
 */
gzoutstream* FileListManager::getCnsFile(const int index)
{
	if (index < 0 || index >= m_cns_file_vec.size())
	{
		// index overflow
		return INDEX_OVER_FLOW;
	}

	return m_cns_file_vec[index];
}


/**
 * DATE: 2010-8-25
 * FUNCTION: open matrix files with open_mode mode
 * PARAMETER: matrix_list: the output or input matrix files path. mode: the open mode.
 *				alignment_list: alignment file list.
 * RETURN: the number of matrix files that be opened.
 */
int FileListManager::openMatrixFile(const std::string matrix_list, std::ios_base::open_mode mode, const std::string alignment_list)
{
	fstream *fsp;
	if (mode == std::ios::in)       //matrix file is existed , read it
	{
		string filepath;
		m_base_file.open(matrix_list.c_str());              //open matrix filelist
		if (!m_base_file) 
		{
			cerr<<"Cannot open file:" << matrix_list << endl;
			return BASEFILE_ERROR;
		}
		//cerr <<"open matrix file 642"<<endl;
		/*open each matrix file in the filelist and input soap vector*/
		while (m_base_file >> filepath)
		{
			//cerr <<"open matrix file 646"<<endl;
			fsp = new fstream();
			//cerr <<"open matrix file 647"<<endl;
			fsp->open(filepath.c_str(), fstream::in) ;
			//cerr <<"open matrix file 649"<<endl;
			if (!fsp->is_open()) 
			{
				cerr << "Cannot open file:" << filepath << endl;
				return SOAPFILE_ERROR;
			}
			m_matrix_file_vec.push_back(fsp);
			fsp->clear();
		}

		m_base_file.close();    
	}
	else if (alignment_list != "")
	{
		string mat_Path, mat_name;         //the one matrix file name

		m_base_file.open(alignment_list.c_str());  //open matrix file path

		if (!m_base_file)        
		{
			cerr << "Cannot open file:" << alignment_list << endl;
			return BASEFILE_ERROR;
		}

		string cmd = "mkdir -p " + matrix_list;
		system(cmd.c_str());
		string list_path = matrix_list + "/matrix_list";
		ofstream list_out(list_path.c_str());

		/*open each matrix file and push in cns vector*/
		while (getline(m_base_file, mat_Path))          
		{
			fsp = new fstream();   
			/*find cons name from cns path*/
			int pos = mat_Path.rfind('/');
			mat_name = matrix_list + mat_Path.substr(pos, string::npos) + ".matrix";
			/*put cns in to cns vector*/
			fsp->open(mat_name.c_str(), fstream::out);

			if (!fsp->is_open()) 
			{
				cerr << "Cannot open file:" << mat_name << endl;
				return CNSFILE_ERROR;
			}

			list_out << mat_name << endl;
			m_matrix_file_vec.push_back(fsp);
		}

		list_out.close();
		m_base_file.close();
	}
	return m_matrix_file_vec.size();
}


/**
 * DATE: 2010-8-25
 * FUNCTION: get the matrix file handle pointer m_matrix_file_vec[index]
 * PARAMETER: index: the vector index.
 * RETURN: return matrix file pointer m_matrix_file_vec[index].failed return INDEX_OVER_FLOW
 */
fstream* FileListManager::getMatrixFile(const int index)
{
	if (index >= m_matrix_file_vec.size())
	{
		// index overflow
		return INDEX_OVER_FLOW;
	}

	return m_matrix_file_vec[index];
}

/**
 * DATE: 2010-8-30
 * FUNCTION: set max limit file number 
 * PARAMETER: filenumber : open filenumber now
 * RETURN: return SETLIMIT_ERROR or return SETLIMIT_OK
 */
//int FileListManager::setlimit(int filenumber)
//{
//	struct rlimit r;
//
//	/*if open file number up to limit */
//	if (filenumber >= r.rlim_max)
//	{
//		/*set limit*/ 
//		r.rlim_cur = 2000;        
//		r.rlim_max = 2000;
//
//		/*if set limit error*/
//		if (setrlimit(RLIMIT_NOFILE, &r) < 0)
//		{
//			cerr << "setrlimit error\n";
//			return SETLIMIT_ERROR;
//		}
//	}
//	return SETLIMIT_OK;
//}

 
  

//2010-11-08
int Read_win_Task::Run()
{
	READ_WIN_ARGS *read_win_args = (READ_WIN_ARGS*)this->m_ptrData;
	
	Readwin & readwin = *(read_win_args->readwin_p);
	igzstream * m_soap_file = read_win_args->m_soap_file_p ;
	SamCtrl * m_sam_file = read_win_args->m_sam_file_p;
	int * ret = read_win_args->ret_p;
	int  index =  read_win_args->index;
	if (read_win_args->m_sam_file_p == NULL)
	{
		ret[index] = fileListManager.readWinFromSoap(readwin, *m_soap_file);
	}
	else
	{
		ret[index] = fileListManager.readWinFromSam(readwin,* m_sam_file);
	}
	return 1;
}

/**
 * DATE: 2010-11-12
 * FUNCTION: a thread function used to run FlieListManager::readWin function.
 * PARAMETER: __Args the parameter structure.
 * RETURN: 
 */
void *_flieListManager_readWin(void * __Args)
{
	// get args.
	BIG_READ_WIN_ARGS * args = (BIG_READ_WIN_ARGS*)__Args;
	// run the function.
	args->fileListManager->readWin((*args->readwin_vec), args->para);
	return NULL;
}
