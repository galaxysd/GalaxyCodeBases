/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-10
 *CLASS NAME: Call_winManager
 *FUNCTION : a class that can control Call_win vector.
 *FILE NAME : Call_winManager.cpp
 *UPDATE DATE : 2010-9-15
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */


#include "Call_winManager.h"


Call_winManager::Call_winManager(void)
{
	m_cw_vec.clear();

}

/**
 * DATE: 2010-10-12
 * FUNCTION: a constructor that can initialize the m_cw_vec vector.
 * PARAMETER: read_len: the longest read's length. file_name
 * RETURN: the flag of addCallwin result.
 */
Call_winManager::Call_winManager(ubit64_t read_len, int file_num, Genome* genome)
{
	m_cw_vec.clear();
	for (int i = 0; i != file_num; ++i)
	{
		Call_win *call_win = new Call_win(read_len);
		if (call_win == NULL)
		{
			cerr << "\tERROR: Can't create Call_win object!" << endl;
			exit(0);
		}

		// add it to the vector
		m_cw_vec.push_back(call_win);
		// initializion.
		call_win->initialize(0);
		call_win->current_chr = genome->chromosomes.end();
	}
}

Call_winManager::~Call_winManager(void)
{
	for (int i = 0; i < m_cw_vec.size(); ++i)
	{
		// release the memory.
		delete m_cw_vec[i];
	}
}


/**
 * DATE: 2010-8-9
 * FUNCTION: add a Call_win object to the Call_win vector.
 * PARAMETER: read_len: the longest read's length.
 * RETURN: the flag of addCallwin result.
 */
int Call_winManager::addCallwin(ubit64_t read_len)
{
	Call_win *call_win = new Call_win(read_len);
	if (call_win == NULL)
	{
		return CRE_WIN_FAILED;
	}

	m_cw_vec.push_back(call_win);
	return m_cw_vec.size();
}


/**
 * DATE: 2010-8-9
 * FUNCTION: soap to consensus.
 * PARAMETER: alignment_vec: alignments to be processed. consensus: the out file handle.
 *				baseinfo: the out file handle. para: the Parameter pointer.
 *				genome: pointer of object that contian reference information.
 *				index: index of the Call_win which to be processed.
 *				ali_index: the index of m_ali_vec that can be processed in the beginning of the 'soap2cns'
 *				sfsMethod: a SfsMethod object, use to sfs.
 * RETURN: CRE_WIN_FAILED when create Call_win failed. INDEX_EXCEED when index exceed
 *			SOAP2CNS_SUCCESS when successful.
 * UPDATE: 2010-10-12, reduce the addCallwin part.
 * UPDATE : 2010-11-22 process the win before new reads outsid the window
 */
int Call_winManager::soap2cns(vector<Soap_format>& alignment_vec, \
							  gzoutstream * consensus,\
							  my_ofstream& baseinfo, \
							  Genome* genome, \
							  Prob_matrix* mat,\
							  Parameter* para, \
							  int index, \
							  int ali_index, \
							  SfsMethod &sfsMethod)
{
/*	if (index == m_cw_vec.size())
	{
		// create a new Call_win object.
		if (addCallwin(para->read_length) != (index + 1))
		{
			return CRE_WIN_FAILED;
		}

		// initialize.
		m_cw_vec[index]->initialize(0);
		m_cw_vec[index]->current_chr = genome->chromosomes.end();
	}
	*/

	if (index > m_cw_vec.size())
	{
		// index is over flow.
		return INDEX_EXCEED;
	}
	
	//process the win before new reads outsid the window
	for (int i = ali_index; i < alignment_vec.size(); ++i)
	{
		m_cw_vec[index]->deal_read(alignment_vec[i], consensus, baseinfo, genome, mat, para, sfsMethod, index);
	}
	// add by guyue 2010-11-25
	if (!m_cw_vec[index]->done_pro_win && m_cw_vec[index]->last_start > 999)
	{
		m_cw_vec[index]->recycled = false;
		m_cw_vec[index]->pro_win(consensus, baseinfo, genome, mat, para, sfsMethod, index);
	}
	m_cw_vec[index]->done_pro_win = false;
	// clear up the alignments.
	alignment_vec.clear();
	return SOAP2CNS_SUCCESS;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: deal the tail of sample.
 * PARAMETER:  consensus: the out file handle.
 *				baseinfo: the out file handle. para: the Parameter pointer.
 *				genome: pointer of object that contian reference information.
 *				index: index of the Call_win which to be processed.
 *				sfsMethod: a SfsMethod object, use to sfs.
 * RETURN:  INDEX_EXCEED when index exceed
 *			SOAP2CNS_SUCCESS when successful.
 */
int Call_winManager::dealTail(gzoutstream * consensus, my_ofstream& baseinfo, Genome* genome, Prob_matrix* mat, Parameter* para, int index, SfsMethod &sfsMethod)
{
	if (index >= m_cw_vec.size())
	{
		// index is over flow.
		return INDEX_EXCEED;
	}
	
	m_cw_vec[index]->deal_tail(consensus, baseinfo, genome, mat, para, sfsMethod, index);

	return SOAP2CNS_SUCCESS;
}

/**
 * DATE: 2010-10-14
 * FUNCTION: the call back function that used to mutli_thread.
 * PARAMETER:  __Args is the parameter list.
 * RETURN: 
 */
//void* _call_winManager_soap2cns_norm(void * __Args)
//{
//	_call_win_args * _args=(_call_win_args*) __Args;
//	_args->call_winManager->soap2cns(_args->readwin->getReadwin(), *(_args->consensus),*(_args->baseinfo),_args->genome,_args->mat, _args->para, _args->index, _args->readwin->getSoap2cnsIdx());
//
//}

/**
 * DATE: 2010-10-13
 * FUNCTION: overload the Run function. Do the Call_winManager's soap2cns function.
 * PARAMETER: 
 * RETURN: 
 */
int Call_win_Task::Run()
{
	CALL_WIN_ARGS *call_win_args = (CALL_WIN_ARGS*)this->m_ptrData;
	// set all the parameters.
	Call_winManager & call_winManager = *(call_win_args->call_winManager_p);
	vector<Soap_format> & alignment_vec_p = *(call_win_args->alignment_vec_p);
	gzoutstream * consensus = call_win_args->consensus_p;
	my_ofstream & baseinfo = *(call_win_args->baseinfo_p);
	Genome * genome = call_win_args->genome;
	Prob_matrix * mat = call_win_args->mat;
	Parameter * para = call_win_args->para;
	int index = call_win_args->index;
	int ali_index = call_win_args->ali_index;
	SfsMethod & sfsMethod = *(call_win_args->sfsMethod_p);

	// do the soap2cns.
	call_winManager.soap2cns(alignment_vec_p, consensus, baseinfo, genome, mat, para, index, ali_index, sfsMethod);
	return 1;
}

void Call_winManager::setLast_start(int start_pos)
{
	if (start_pos <= 0)
	{
		return;
	}
	for (int i = 0; i < m_cw_vec.size(); ++i)
	{
		m_cw_vec[i]->last_start = start_pos;
		m_cw_vec[i]->m_is_set_ls = true;
	}
}
