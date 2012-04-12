/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-9
 *CLASS NAME: FileListManager
 *FUNCTION : a class that can control Prob_matrix vector. support add, change, get 
 *			 Matrix ect functions.
 *FILE NAME : MatrixManager.cpp
 *UPDATE DATE : 2010-8-25
 *UPDATE BY : Bill Tang
 *UPDATE: 2010-11-2 change the interface which is related to soap file.
 *UPDATA BY :BIll Tang
 *******************************************************************************
 */
#include "MatrixManager.h"


MatrixManager::MatrixManager(void)
{
}


MatrixManager::~MatrixManager(void)
{
	// delete the Prob_matrix objects.
	for (int i = 0; i < m_mat_vec.size(); i++)
	{
		delete m_mat_vec[i];
	}
}

/**
 * DATE: 2010-8-9
 * FUNCTION: add a Prob_matrix object to the m_mat_vec.
 * PARAMETER: soap_vec: the alignment vector will be processed. para: the Parameter pointer.
 *				genome: pointer of boject that contian reference information. index: the index of the m_mat_vec
 * RETURN: if add failed then return CREATE_MAT_FAILED. successful return m_mat_vec.size().
 * UPDATE: 2010-10-12 add a parameter index
 */
int MatrixManager::addMatrix(SamCtrl &alignment, Parameter* para, Genome* genome, fstream *outFile, int index)
{
	Prob_matrix *prob_matrix = new Prob_matrix(para->rank_sum_mode);

	if (prob_matrix == NULL)
	{
		// create Prob_matrix object failed
		return CREATE_MAT_FAILED;
	}

	// generate matrix.
	prob_matrix->matrix_gen(alignment, para, genome);
	// write matrix to file.
	if (outFile != NULL)
	{
		prob_matrix->matrix_write(*outFile, para);
	}
	// count prior
	prob_matrix->prior_gen(para);
	if (para->rank_sum_mode == true)
	{
		prob_matrix->rank_table_gen();
	}
	m_mat_vec[index] = prob_matrix;

	return m_mat_vec.size();
}

/**
 * DATE: 2010-8-9
 * FUNCTION: add a Prob_matrix object to the m_mat_vec.
 * PARAMETER: soap_vec: the alignment vector will be processed. para: the Parameter pointer.
 *				genome: pointer of boject that contian reference information. index: the index of the m_mat_vec
 * RETURN: if add failed then return CREATE_MAT_FAILED. successful return m_mat_vec.size().
 * UPDATE: 2010-10-12 add a parameter index
 */
int MatrixManager::addMatrix(igzstream &alignment, Parameter* para, Genome* genome, fstream *outFile, int index)
{
	Prob_matrix *prob_matrix = new Prob_matrix(para->rank_sum_mode); //initialize array will be used in the prob matrix

	if (prob_matrix == NULL)
	{
		// create Prob_matrix object failed
		return CREATE_MAT_FAILED;
	}

	// generate correction matrix.
	prob_matrix->matrix_gen(alignment, para, genome);

	// write matrix to file.
	if (outFile != NULL)
	{
		prob_matrix->matrix_write(*outFile, para);
	}

	// generate the prior matrix 
	prob_matrix->prior_gen(para);

	if (para->rank_sum_mode == true)
	{
		prob_matrix->rank_table_gen();
	}
	m_mat_vec[index] = prob_matrix;

	return m_mat_vec.size();
}


/**
 * DATE: 2010-8-9
 * FUNCTION: change the Prob_matrix with new soap_vec.
 * PARAMETER: soap_vec: the alignment vector will be processed. para: the Parameter pointer.
 *				genome: pointer of object that contian reference information.
 *				index: index of the Prob_matrix point which to be change.
 * RETURN: the flag of change result.
 */
int MatrixManager::changeMatrix(vector<Soap_format>& soap_vec, Parameter* para, Genome* genome, int index)
{
	return 0;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: get a Prob_matrix point that point to the m_mat_vec[index].
 * PARAMETER: index: the index of Prob_matrix pointer which to be get.
 * RETURN: a Prob_matrix pointer. If index over flow then return INDEX_OVER_FLOW.
 */
Prob_matrix* MatrixManager::getMatrix(int index)
{
	// judge if the index exceed the region of the vector.
	if (index < 0 || index > m_mat_vec.size())
	{
		return INDEX_OVER_FLOW;
	}

	return m_mat_vec[index];
}


 
/**
 * DATE: 2010-8-25
 * FUNCTION: add a matrix from the prior matrixfile.
 * PARAMETER: soap_vec: the alignment vector will be processed. para: the Parameter pointer.
 * RETURN: if add failed then return CREATE_MAT_FAILED. successful return m_mat_vec.size().
 */
int MatrixManager::addMatrix(std::fstream& matrix_file, Parameter *para)
{
	//update 11-11 for ran_sum_mode flag
	Prob_matrix *prob_matrix = new Prob_matrix(para->rank_sum_mode);

	if (prob_matrix == NULL)
	{
		// create Prob_matrix object failed
		return CREATE_MAT_FAILED;
	}

	prob_matrix->matrix_read(matrix_file, para);
	// count prior
	prob_matrix->prior_gen(para);
	//2010-11-11
	if (para->rank_sum_mode == true)
	{
		prob_matrix->rank_table_gen();
	}
	m_mat_vec.push_back(prob_matrix);

	return m_mat_vec.size();
}


/**
 * DATE: 2010-10-12
 * FUNCTION: set the m_mat_vec's size.
 * PARAMETER: number is the size which m_mat_vec should has.
 * RETURN: if add failed then return CREATE_MAT_FAILED. successful return m_mat_vec.size().
 */
void MatrixManager::setMatrixNum(const int number)
{
	int tmp_num = number;
	if (number < 0)
		tmp_num = 1;
	Prob_matrix *prob_matrix = NULL;
	while (m_mat_vec.size() < tmp_num)
	{
		m_mat_vec.push_back(prob_matrix);
	}
}


/**
 * DATE: 2010-10-12
 * FUNCTION: the funtion be used to do mutli_thread.
 * PARAMETER: __Args is the parameter array.
 * RETURN: 
 */
void *_matrixManager_addMatrix(void * __Args)
{
	MATRIX_ARGS * _args = (MATRIX_ARGS*)__Args;

	if(_args->sam_alignment != 0)
	{
		if (CREATE_MAT_FAILED == _args->matrixManager->addMatrix(*(_args->sam_alignment),_args->para,_args->genome,_args->outFile,_args->index))
		{
			cerr << "ERROR: add matrix failed!" << endl;
			exit(0);
		}
	}
	else
	{
		if (CREATE_MAT_FAILED ==_args->matrixManager->addMatrix(*(_args->soap_alignment),_args->para,_args->genome,_args->outFile,_args->index))
		{
			cerr << "ERROR: add matrix failed!" << endl;
			exit(0);
		}
	}
	return NULL;
}

