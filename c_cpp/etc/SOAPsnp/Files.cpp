/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-9-15
 *UPDATE DATE : 2010-9-15
 *UPDATE BY :	Zhou Guyue
 *******************************************************************************
 */



#include "soap_snp.h"
#include <cassert>

/**
 * DATE: 2010-9-15
 * FUNCTION: open sfs file
 * PARAMETER: outfiles : outfile name. 
			  writeFr : whether writer frequence file.
			  doBay : whether do Bay.
			  doJoint : whether do Joint
 * RETURN:	 OPENSFS_ERROR:if can not open file  ,OPENSFS_SUCC : if can open file.
 */
int Files::OpenSfsfile(const string outfiles, const int writeFr, const int doBay, const int doJoint)
{
	 //generete the output filenames
	if(!outfiles.c_str())
	{
		printf("Must supply -outfiles (-fai)\n");
		return OPENSFS_ERROR;
	}
	 
	fSFSall = outfiles + ".sfs";
	fFreq = outfiles + ".frq";
	fJoint = outfiles + ".bjoint";
 
    //open the persite FILE streams
	if(writeFr)
	{
		//freqfile.clear();
		//freqfile.open(fFreq.c_str());
		freqfile= getFILE(fFreq.c_str(),"w");
	}
	if(doBay)
	{
		//sfsfile.clear();
		//sfsfile.open(fSFSall.c_str());
		sfsfile = getFILE(fSFSall.c_str(),"w");
		//sfsfile-open(fSFSall.c_str());
	}
	if(doJoint)
	{
		//jointSfsfile.clear();
		//jointSfsfile.open(fJoint.c_str(), ios::binary);
		jointSfsfile = getFILE(fJoint.c_str(),"w");
	}
	if ((writeFr && !freqfile) || (doBay && !sfsfile) || (doJoint && !jointSfsfile))
	{
		cerr << "\topen sfs file failed!" << endl;
		return OPENSFS_ERROR;
	}
	return OPENSFS_SUCC;

}

 /**
 * DATE: 2010-9-15
 * FUNCTION: get file in mode
 * PARAMETER: fname : file name , mode : write. 			   
 * RETURN:	 FILE .
 */
FILE* Files::getFILE(const char*fname,const char* mode)
{
	FILE *fp;
	if(NULL == (fp = fopen(fname, mode)))
	{
		fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
		exit(0);
	}
	return fp;
}