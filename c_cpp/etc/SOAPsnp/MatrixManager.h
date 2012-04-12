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
 *FILE NAME : MatrixManager.h
 *UPDATE DATE : 2010-8-25
 *UPDATE BY : Bill Tang
 *UPDATE: 2010-11-2 change the interface which is related to soap file.
 *UPDATA BY :BIll Tang
 *******************************************************************************
 */

#pragma once
#include <fstream>
#include <vector>
#include "soap_snp.h"

#define INDEX_OVER_FLOW	NULL
#define CREATE_MAT_FAILED -1

using namespace std;
class MatrixManager
{
public:
	MatrixManager(void);
	virtual ~MatrixManager(void);
private:
	// a Prob_matrix vectior point
	vector<Prob_matrix*> m_mat_vec;
public:
	// add Matrix to the vector
	virtual int addMatrix(SamCtrl &alignment, Parameter* para, Genome* genome, fstream *outFile, int index);
	// add a matrix from the prior matrixfile.
	virtual int addMatrix(igzstream &alignment, Parameter* para, Genome* genome, fstream *outFile, int index);
	// change the Prob_matrix with new soap_vec
	virtual int changeMatrix(vector<Soap_format>& soap_vec, Parameter* para, Genome* genome, int index);
	// get a Prob_matrix point that point to the m_mat_vec[index]
	virtual Prob_matrix* getMatrix(int index);
	// add a matrix from the prior matrixfile.
	virtual int addMatrix(std::fstream& matrix_file, Parameter *para);
	// set the m_mat_vec's size
	virtual void setMatrixNum(const int number);
};


void *_matrixManager_addMatrix(void * __args);

/**
 * DATE: 2010-10-12
 * FUNCTION: a parameter structure that about matrix.
 * PARAMETER: 
 * RETURN: 
 */
typedef struct _matrix_args
{
	MatrixManager * matrixManager;
	//FileListManager * fileListManager;
	Parameter * para;
	Genome* genome;
	igzstream *soap_alignment;
	std::fstream *outFile;
	SamCtrl *sam_alignment;
	int index;

	/**
	 * DATE: 2010-10-12
	 * FUNCTION: initialize the struct.
	 * PARAMETER: 
	 * RETURN: 
	 */
	inline _matrix_args(MatrixManager * a, Parameter * b, Genome* c,
			igzstream * d, std::fstream * e, SamCtrl * f, int g)
	{
		matrixManager = a;
		para = b;
		genome = c;
		soap_alignment = d;
		outFile = e;
		sam_alignment = f;
		index = g;
	}
	;
	/**
	 * DATE: 2010-10-12
	 * FUNCTION: initialize the struct.
	 * PARAMETER: 
	 * RETURN: 
	 */
	inline _matrix_args()
	{
		matrixManager = 0;
		para = 0;
		genome = 0;
		soap_alignment = 0;
		outFile = 0;
		sam_alignment = 0;
	}
	;

}MATRIX_ARGS;

