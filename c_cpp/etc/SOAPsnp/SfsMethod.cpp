/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-9-13
 *CLASS NAME: SfsMethod
 *FUNCTION : a class that contain method about sfs.
 *FILE NAME : SfsMethod.cpp
 *UPDATE DATE : 2010-9-13
 *UPDATE BY : Bill Tang
  *UPDATE DATE : 2010-9-27
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */

#include "SfsMethod.h"
#include "soap_snp.h"

#define BE_PROCESSED -1

bool operator <= (const loci& first,const  loci& second) {
	if(first.chromo==second.chromo)
		return first.position<=second.position;
	else
		return first.chromo<=second.chromo;
}


bool operator > (const loci& first,const  loci& second) {
	return (!(first <= second));
}

bool operator == (const loci& first,const  loci& second) {
	if(first.chromo==second.chromo)
		return first.position==second.position;
	else
		return first.chromo==second.chromo;
}

SfsMethod::SfsMethod()
{
	binomLookup = NULL;
	likes = NULL;
	bayes = NULL;
	asso[0].clear();
	asso[1].clear();
	asso[2].clear();
	//end of update
	m_map_idx = 0;
	m_map_idx_process = BE_PROCESSED;
	sem_init(&sem_call_sfs, 0, 0);
	sem_init(&sem_call_sfs_return, 0, 1);
	//sem_init(&sem_map_number, 0, 2);
	file_end_flag = 0;
}

pthread_mutex_t SfsMethod::m_pthreadMutex = PTHREAD_MUTEX_INITIALIZER;

SfsMethod::SfsMethod(int numInds)
{
	/*
	depending on major/minor we will only use 3 genotypes,
	so we will need a small lookup table to decide which columns to extract from the glf file
	  A C G T
	A 0 1 2 3
	C   4 5 6
	G     7 8
	T       9

	*/
	int n = 0;
	for(int i = 0; i < 4; i++)
	{
		for(int j = i; j < 4; j++)
		{
			glfLookup[i][j] = n;
			glfLookup[j][i] = n;
			n++;
		}
	}

	generate_likeLookup();
	if (numInds < 0)
	{
		generate_binom(0);
		likes = allocDoubleArray(1);
		bayes = allocDoubleArray(1);
		m_numInds = numInds;
	}
	else
	{
		generate_binom(numInds);
		likes = allocDoubleArray(2 * numInds + 1);
		bayes = allocDoubleArray(2 * numInds + 1);
		m_numInds = numInds;
	}
}

/**
 * DATE: 2010-9-16
 * FUNCTION: initialization :likes , bays, binomLookup matrix
 * PARAMETER:  numInds: the individual number.	
 * RETURN:	 0
 */
int SfsMethod::init(const int numInds)
{
	/*
	depending on major/minor we will only use 3 genotypes,
	so we will need a small lookup table to decide which columns to extract from the glf file
	  A C G T
	A 0 1 2 3
	C   4 5 6
	G     7 8
	T       9

	*/
	int n = 0;
	for(int i = 0; i < 4; i++)
	{
		for(int j = i; j < 4; j++)
		{
			glfLookup[i][j] = n;
			glfLookup[j][i] = n;
			n++;
		}
	}

	generate_likeLookup();
	if (numInds < 0)
	{
		generate_binom(0);
		likes = allocDoubleArray(1);
		bayes = allocDoubleArray(1);
		m_numInds = numInds;
	}
	else
	{
		generate_binom(numInds);
		likes = allocDoubleArray(2 * numInds + 1);
		bayes = allocDoubleArray(2 * numInds + 1);
		m_numInds = numInds;
	}
	//allocVec();
	return 0;
}

SfsMethod::~SfsMethod(void)
{
	if (binomLookup != NULL)
	{
		delete [] binomLookup;
		binomLookup = NULL;
	}
	if (likes != NULL)
	{
		delete [] likes;
		likes = NULL;
	}
	if (bayes != NULL)
	{
		delete [] bayes;
		bayes = NULL;
	}
	//delVector(asso);
	cleanUpMap(asso[0]);
	cleanUpMap(asso[1]);
	cleanUpMap(asso[2]);
	//clear sem
	sem_destroy(&sem_call_sfs);
	sem_destroy(&sem_call_sfs_return);
}

void outMap(aMap::iterator &itr, int numInds)
{
	cerr << "chr : " << itr->first.chromo << " pos : " << itr->first.position;
	datum p = itr->second;
	cerr << " ref :" << p.ref << "  major : " << p.major << " minor: " << p.minor << endl;
	int *lk = p.lk;
	cerr << "lk :" << endl;
	for (int j = 0; j < numInds; j++)
	{
		for (int i = 0; i < 10; i++)
			cerr << lk[4 * j + i] << "\t";
		cerr << endl;
	}
	cerr << endl;
	int *locus = p.locus;
	cerr << "locus: " << endl;
	for (int j = 0; j <= numInds; j++)
	{
		for (int i = 0; i < 4; i++)
			cerr << locus[4 * j + i] << "\t";
		cerr << endl;
	}
	cerr << "phat : " << p.phat;
	cerr << endl << endl;
}

/**
 * DATE: 2010-9-13
 * FUNCTION: the function for algo
 * PARAMETER:	asso : the maps of the datums
 *				numInds : the number of individuals
 *				eps : the errorate
 *				pvar : the probability of being variable
 *				B :the bias correction coefficient
 *				sfsfile : sfs file pionter
 *				normalize : 
 *				alternative : 
 *				singleMinor : we  should loop over 3 possible minor or one.
 *				pThres : select the out put threshold of snp £¬default is 0.01
 *				allowPhatZero
 * RETURN: void
 * UPDATE:2010-10-21, add a set phat condition.
 * UPDATE:2010-11-29, control the sfs file column 
 */
void SfsMethod::algo(aMap & asso, int numInds, double eps, double pvar, double B, FILE * sfsfile, 
					 int normalize, int alternative, int singleMinor, double pThres, int allowPhatZero, int firstlast)
{
	normalize=0;

	//algorithm goes on by a site on site basis
	double *pis = new double[numInds];
	double *wis = new double[numInds];
	double *cf = new double[numInds];//the bias coefiecients
	double *sumMinors = allocDoubleArray(2*numInds+1);//hte sum of the 3 different minors
	double *hj = allocDoubleArray(2*numInds+1);

	for(aMap::iterator it=asso.begin(); it!=asso.end(); it++) {

		//outMap(it, numInds);
		datum p = it->second;
		//part one
		if (p.major > 3 || p.minor==p.major){
			//printf("never here\n");
			continue;
		}

		//set the resultarray to zeros
		for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
			sumMinors[sm] = 0;

		//loop through the 3 different minors
		for(int minor_offset=0;minor_offset<4;minor_offset++) {

			if(minor_offset == p.major)
				continue;

			//hook for only calculating one minor
			if(singleMinor)
				minor_offset = p.minor;

			int major_offset = p.major;
			int Aa_offset = glfLookup[minor_offset][major_offset];//0-9
			int AA_offset = glfLookup[minor_offset][minor_offset];//0-9
			int aa_offset = glfLookup[major_offset][major_offset];//0-9
			

			for(int i=0;i<numInds;i++) {
				int ni = p.locus[i*4+minor_offset];
				int nt = p.locus[i*4+minor_offset] + p.locus[i*4+major_offset];

				if(B!=0)//bias version
					cf[i] = pow((0.5-B),ni)*pow((0.5+B),nt-ni)*pow(0.5,-nt);
				else
					cf[i] =1;

				//FIX if we have a very large number ( paralogue sites), the cf[i] will overflow, and it should be zero
				if(isinf_local(cf[i])||isnan_local(cf[i])){
					// fprintf(stderr,"POSSIBLE ERROR: Problems calculating bias correction coefficient, site might be paralogue?\t%d\t%d\n",it->first.chromo,it->first.position);
					cerr << "POSSIBLE ERROR: Problems calculating bias correction coefficient, site might be paralogue?\t" << it->first.chromo << "\t" << it->first.position << endl;
					cf[i] = 1;
				}
				if(nt==0){//if we dont have any reads for individual 'i'
					pis[i] = 0;
					wis[i] = 0;
					continue;
				}
				pis[i] = (ni-eps*nt)/(nt*(1-2*eps));
				wis[i] = 2.0*nt/(nt+1.0);
			}
			
			double tmp=0;
			for(int i=0;i<numInds;i++)
				tmp+= (pis[i]*wis[i]);

			double phat = std::max(0.0,tmp/sum(wis,numInds));//original

			//    phat = std::max(pvar/(4*numInds),tmp/sum(wis,numInds)); //ramus old
			//      
			if(allowPhatZero==0){
				int depth = 0;
				for (int pl=0;pl<4;pl++) depth += p.locus[4*numInds+pl];
				phat = std::max(pvar/std::min(10,depth),phat); //thorfinn new
			}
			//
			//add bias stuff
			double newc = (0.5-B)/(0.5+B);
			phat = std::min(1.0,phat/(newc+phat*(1.0-newc)));
			//end bias

			// change at 2010-10-21
			if (minor_offset == p.minor)
			{
				it->second.phat = phat;
			}

			//part two
			// change on 2010-9-18
			memset(hj, 0, sizeof(double) * (2*numInds+1));
			/*for(int i = 0; i < 2*numInds+1; ++i)
			{
				hj[i] = 0.0;
			}*/

			if(0)
				for(int index=0;index<(2*numInds+1);index++)
					hj[index]=log(hj[index]);
			double PAA,PAa,Paa;
			//    phat =  0.224518;

			double totmax = 0.0;
			double gaa_prod = 1.0; //product of the major genotyp likelihood
			for(int i=0 ; i<numInds ;i++) {
				double GAA,GAa,Gaa;
#ifndef RASMUS_INPUT
				GAA = likeLookup[p.lk[i*10+AA_offset]];
				GAa = likeLookup[p.lk[i*10+Aa_offset]];
				Gaa = likeLookup[p.lk[i*10+aa_offset]];
#else
				GAA = exp(p.lk[i*10+AA_offset]);
				GAa = exp(p.lk[i*10+Aa_offset]);
				Gaa = exp(p.lk[i*10+aa_offset]);

				//      printf("GAA=%.15f\tGAa=%.15f\tGaa=%.15f\n",GAA,GAa,Gaa);
				//printf("sums g,=%.15f\n",GAa+GAA+Gaa);
#endif
				if(1){
					double sums= GAA+GAa+Gaa;//this is essential for floating point precision
					GAA = GAA/(sums);
					GAa = GAa/(sums);
					Gaa = Gaa/(sums);

				}

				gaa_prod *= Gaa;
				double MAA,MAa,Maa;
				if(phat==0||1){
					MAA = GAA * phat*phat;
					MAa = GAa * 2*(1-phat)*phat*cf[i]; //DRAGON remember to add cf[i]
					Maa = Gaa * (1-phat)*(1-phat);
				}else{
					//	  fprintf(stderr,"hereee\n");
					MAA = GAA ;
					MAa = GAa * 2;
					Maa = Gaa ;
				}
				//printf("mAA=%.15f\tmAa=%.15f\tMaa=%.15f\n",MAA,MAa,Maa);

				if(1 &&(isnan_local(Maa) || isnan_local(MAa) || isnan_local(MAA))){
					//fprintf(stderr,"Possible underflow at: \t%d\t%d\n",it->first.chromo,it->first.position);
					cerr << "Possible underflow at: \t" << it->first.chromo << "\t" << it->first.position << endl;
					goto next_allele;// change 2010.9.18
					//exit(0);
				}
				if(0){
					PAA =log(MAA);///(MAA+MAa+Maa);
					PAa =log(MAa);///(MAA+MAa+Maa);
					Paa =log(Maa);///(MAA+MAa+Maa);
					//      printf("sum af P's=%f\n",PAA+PAa+Paa);
				}else{
					PAA =(MAA);///(MAA+MAa+Maa);
					PAa =(MAa);///(MAA+MAa+Maa);
					Paa =(Maa);///(MAA+MAa+Maa);
				}

				if(normalize){
					PAA =log(MAA/(MAA+MAa+Maa));
					PAa =log(MAa/(MAA+MAa+Maa));
					Paa =log(Maa/(MAA+MAa+Maa));
				}

				//check for underflow error, this should only occur once in a blue moon
				if(isnan_local(Paa)||isnan_local(PAa)||isnan_local(PAA)){
					printf("fixing nan error at : \t%s\t%d\n",it->first.chromo.c_str(),it->first.position); 
					goto next_allele;// change 2010.9.18 
					//exit(0);
				}
				if(isnan_local(Paa)||isnan_local(PAa)||isnan_local(PAA)){
					//fprintf(stderr,"Possible underflow at: \t%d\t%d\n",it->first.chromo,it->first.position);
					cerr << "Possible underflow at: \t" << it->first.chromo << "\t" << it->first.position << endl;
					printf("PAA=%f\tPAa=%f\tPaa=%f\tMAA=%f\tMAa=%f\tMaa=%f\n",PAA,PAa,Paa,MAA,MAa,Maa);
					
				}

				double mymax;
				if (Paa > PAa && Paa > PAA) mymax = Paa;
				else if (PAa > PAA) mymax = PAa;
				else mymax = PAA;


				Paa=Paa/mymax;
				PAa=PAa/mymax;
				PAA=PAA/mymax;
				totmax = totmax + log(mymax);

				if(i==0){
					hj[0] =Paa;
					hj[1] =PAa;
					hj[2] =PAA;
				}else{

					for(int j=2*(i+1); j>1;j--){
						double tmp;
						if(0)
							tmp =log(exp(PAA+hj[j-2])+exp(PAa+hj[j-1])+exp(Paa+hj[j]));
						else
							tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];

						if(isnan_local(cf[i]))
							//fprintf(stderr,"cfi is inf\n");
							cerr << "cfi is inf" <<endl; 
						if(isnan_local(tmp)){
							//fprintf(stderr,"jis nan in algorithm:%d\n",j );
							cerr << "jis nan in algorithm: " << j <<endl;
							//fprintf(stderr,"%f\t%f\t%f\t%f\t%f\n",PAA,PAa,Paa,mymax,cf[i]);
							cerr << PAA << "\t" << PAa << "\t" << Paa << "\t" << mymax << "\t" << cf[i] <<endl;
							//fprintf(stderr,"at : \t%d\t%d\n",it->first.chromo,it->first.position);
							cerr << "at : " << it->first.chromo << "\t" << it->first.position << endl;
							hj[j] = 0;
							goto next_allele;// change 2010.9.18
							//exit(0);
						}else
							hj[j]  =tmp;
					}
					if(0){
						hj[1] = log(exp(Paa+hj[1])+exp(PAa+hj[0]));
						hj[0] = Paa+hj[0];
					}
					else{
						hj[1] = Paa*hj[1] + PAa*hj[0];
						hj[0] = Paa*hj[0];
					}
				}
			}
			//recursion done, now do the normalization, and extrastuff

			if(alternative){
				//populate the log(M) vector
				double *m = allocDoubleArray(2*numInds+1);
				for(int mi=0;mi<(2*numInds+1);mi++)
					m[mi] = binomLookup[mi]*pow(phat,mi)*pow(1-phat,2*numInds-mi);

				//now add the mi to the hi,
				for(int i=0;i<(2*numInds+1) ;i++)
					hj[i] = hj[i] * m[i];

				// change on 2010-9-25
				delete [] m;
			}

			if(0) // this is a artifact from the no-underflow version
				for (int i=0;i<(2*numInds+1);i++)
					hj[i] = exp(hj[i]);

			for (int i=0;i<(2*numInds+1);i++)
				hj[i] = exp(log(hj[i])+totmax);

			double bottom  = pvar*sum(hj,2*numInds+1)+(1-pvar)*gaa_prod;


			if(bottom==0||isnan_local(bottom)||isinf_local(bottom)){
				//fprintf(stderr,"bottom equal zero should not happen: %f\n",bottom);
				cerr << "bottom equal zero should not happen: " << bottom <<endl;
				//fprintf(stderr,"at : \t%d\t%d\n",it->first.chromo,it->first.position);
				cerr << "at : " << "\t" <<it->first.chromo << "\t" << it->first.position <<endl;
				//fprintf(skipped_sites_fp,"%s\t%d\n",it->first.chromo.c_str(),it->first.position);
				//sfsfile << it->second.chromo << "\t" << it->second.position << "\t-1" << endl;
				fprintf(sfsfile,"%s\t%d\t-1\n",it->first.chromo.c_str(),it->first.position);
				
				minor_offset=4;
				goto next_allele; // change on 2010.9.18
				break;

				//	exit(0);
			}


			//now input the normalized hj array
			hj[0] = ((pvar*hj[0])+(1-pvar)*gaa_prod)/bottom;//prod(p.lk,numInds,2))/bottom;
			for(int i=1;i<(2*numInds+1) ; i++)
				hj[i] = (hj[i]*pvar/bottom);

			//validate that everything looks fine
			if(sum(hj,2*numInds+1)>1.0000001){
				printf("This shouldnt occur, it means that the sfs for a loci has an improper prob space\n");
				printf("sums is:%f\n",sum(hj,2*numInds+1));
				printf("at : \t%s\t%d\n",it->first.chromo.c_str(),it->first.position);
				goto next_allele; // change 2010.9.18
				//exit(0);
			}
			//update the sumSFS of the site
			for(int ss=0;ss<(2*numInds+1);ss++)
				sumMinors[ss] += hj[ss];
			likes[getMaxId(hj,2*numInds+1)]++;
			
			if(singleMinor)
				break;
		}


		for(int i=0;i<(2*numInds+1);i++)
			if(singleMinor==0)
				sumMinors[i] = sumMinors[i]/3.0;

		//update the global bayes/like estimates
		for(int i=0;i<(2*numInds+1);i++)
			bayes[i] += sumMinors[i];

		//printit
		if((1-sumMinors[0]-sumMinors[2*numInds])>pThres){
			//sfsfile << it->second.chromo << "\t" << it->second.position << "\t";
			fprintf(sfsfile,"%s\t%d\t",it->first.chromo.c_str(),it->first.position);
			//sfsfile << getChar(it->second.major) << "\t" << getChar(it->second.minor) << "\t";
			fprintf(sfsfile,"%c\t%c\t",getChar(it->second.major),getChar(it->second.minor));
		//update 11-29 for control the sfs file column
			//fprintf(sfsfile,"%f\t",sumMinors[i]);
			for(int i=0;i<(2*numInds);i++)
			{
				if ( firstlast == 1 )  //only output the first and last column to the sfs file
				{
					fprintf(sfsfile,"%f\t",sumMinors[i]);
					break; 
				}
				else
				{
					fprintf(sfsfile,"%f\t",sumMinors[i]);
				}
			}
				//sfsfile << sumMinors[i] << "\t";
			//sfsfile << sumMinors[2*numInds] << endl;
			fprintf(sfsfile,"%f\n",sumMinors[2*numInds]);
		}

next_allele:
		;
		
		//cleanup


	}
	delete []  hj;
	delete [] pis;
	delete [] wis;
	delete [] sumMinors;
	delete [] cf;
}


/**
 * DATE: 2010-8-9
 * FUNCTION:  do algo joint function
 * PARAMETER:	asso : the maps of the datums
 *				numInds : the number of individuals
 *				eps : the errorate
 *				pvar : the probability of being variable
 *				B : the bias correction coefficient
 *				sfsfile : sfs file pionter
 *				underFlowProtect: 
 *				alternative : 
 *				singleMinor : should we loop over 3 possible minor. MAYBE
 *				fold : 
 * RETURN: void
 */
void SfsMethod::algoJoint(aMap & asso, int numInds, double eps, double pvar, double B, FILE * sfsfile, 
						  int underFlowProtect, int alternative, int singleMinor, int fold)
{
	
	//algorithm goes on by a site on site basis //pretty much the same as 'algo' without the prior phat

	double sumMinors[2*numInds+1];  //the sum of the 3 different minors

	for(aMap::iterator it=asso.begin(); it!=asso.end(); it++) {//loop over sites
		datum p = it->second;

		//part one

		if (p.major > 3 || p.minor==p.major){
		//	printf("never runs\n");
			continue;
		}

		//set the resultarray to zeros
		for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
			sumMinors[sm] = 0;

		//loop through the 3 different minors
		for(int minor_offset=0;minor_offset<4;minor_offset++) {

			if(minor_offset == p.major)
				continue;

			//hook for only calculating one minor
			if(singleMinor)
				minor_offset = p.minor;

			int major_offset = p.major;
			int Aa_offset = glfLookup[minor_offset][major_offset];//0-9
			int AA_offset = glfLookup[minor_offset][minor_offset];//0-9
			int aa_offset = glfLookup[major_offset][major_offset];//0-9

			//part two
			double hj[2*numInds+1];
			for(int index=0;index<(2*numInds+1);index++)
				if(underFlowProtect==0)
					hj[index]=0;
				else
					hj[index]=log(0);
			double PAA,PAa,Paa;

			for(int i=0 ; i<numInds ;i++) {
				//printf("AA=%d\tAa=%d\taa=%d\n",p.lk[i*3+AA],p.lk[i*3+Aa],p.lk[i*3+aa]);
				double GAA,GAa,Gaa;
#ifndef RASMUS_INPUT //abi input
				GAA = log(likeLookup[p.lk[i*10+AA_offset]]);
				if(p.lk[i*10+AA_offset]==0 && p.lk[i*10+Aa_offset]==0 && p.lk[i*10+aa_offset]==0)
					GAa = log(likeLookup[p.lk[i*10+Aa_offset]]);
				else
					GAa = log(2.0*likeLookup[p.lk[i*10+Aa_offset]]);
				Gaa = log(likeLookup[p.lk[i*10+aa_offset]]);
#else
				GAA = p.lk[i*10+AA_offset];
				GAa = log(2.0)+p.lk[i*10+Aa_offset];
				Gaa = p.lk[i*10+aa_offset];

#endif
				if(underFlowProtect==0){
					GAA=exp(GAA);
					GAa=exp(GAa);
					Gaa=exp(Gaa);
				}


				PAA =(GAA);///(MAA+MAa+Maa);
				PAa =(GAa);///(MAA+MAa+Maa);
				Paa =(Gaa);///(MAA+MAa+Maa);


				//check for underflow error, this should only occur once in a blue moon
				if(isnan_local(Paa)||isnan_local(PAa)||isnan_local(Paa)){
					//fprintf(stderr,"Possible underflow at: \t%d\t%d\n",it->first.chromo,it->first.position);
					cerr << "Possible underflow at: \t" << it->first.chromo  << "\t" << it->first.position << endl;
					printf("PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
				}

				if(i==0){
					hj[0] =Paa;
					hj[1] =PAa;
					hj[2] =PAA;
				}else{

					for(int j=2*(i+1); j>1;j--){
						//print_array(hj,2*numInds+1);
						double tmp;
						if(underFlowProtect==1)
							tmp =addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
						else
							tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];

						if(isnan_local(tmp)){
							printf("jis nan:%d\n",j );
							print_array(hj,2*numInds+2);
							hj[j] = 0;
							//	    exit(0);
							break;
						}else
							hj[j]  =tmp;
					}
					if(underFlowProtect==1){
						hj[1] = addProtect2(Paa+hj[1],PAa+hj[0]);
						hj[0] = Paa+hj[0];
					}
					else{
						hj[1] = Paa*hj[1] + PAa*hj[0];
						hj[0] = Paa*hj[0];
					}
				}


			}
			for(int i=0;i<(2*numInds+1);i++)
				if(underFlowProtect==0)
					sumMinors[i] +=  exp(log(hj[i])-log(bico(2*numInds,i)));
				else
					sumMinors[i] +=  exp(hj[i]-log(bico(2*numInds,i)));
			if(singleMinor)
				break;


		}
		//printit by rasmus 0.98 modified by thorfinn 0.981

		for(int i=0;i<(2*numInds+1);i++)
			if(singleMinor)
				sumMinors[i] = log(sumMinors[i]);
			else
				sumMinors[i] = log(sumMinors[i]/3);

		if(fold){
			int newDim = numInds+1;
			for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
				sumMinors[i] = log(exp(sumMinors[i]-log(2.0)) + exp(sumMinors[2*numInds-i]-log(2.0)));//THORFINN NEW
			//sfsfile << sumMinors;
			fwrite(sumMinors,sizeof(double),newDim,sfsfile);
		}else
			//sfsfile << sumMinors;
			fwrite(sumMinors,sizeof(double),2*numInds+1,sfsfile);

	}

}

/**
 * DATE: 2010-9-13
 * FUNCTION: write frequnce file
 * PARAMETER:	
 *				pFile : a file pointer 
 *				numInds : the number of individuals			
 *				asso : the maps of the datums
				
 * RETURN: void
 */
void SfsMethod::writeFreq(FILE * pFile, int numInds, aMap & asso){
	FILE *db =NULL;
	int printCounts=0;
	if(printCounts)
		db=fopen("indcounts.txt","w");
	for(aMap::iterator it=asso.begin(); it!=asso.end(); it++){
		//loci l = it->first;
		//datum d = it->second;
		//outMap(it, numInds);

		datum d = it->second;
		if (d.major > 3 || d.minor > 3 || it->first.position == 0)
		{
			continue;
		}

		//pFile << d.chromo << "\t" << d.position << "\t";
		fprintf(pFile,"%s\t%d\t",it->first.chromo.c_str(),it->first.position) ;
		//pFile << getChar(d.major) << "\t" << getChar(d.minor) << "\t";
		fprintf(pFile,"%c\t%c\t",getChar(d.major),getChar(d.minor)) ;
		int *bases = d.locus;
		for(int i=0;i<4;i++)
			fprintf(pFile,"%d\t",bases[4*numInds+i]);
			//pFile << bases[4*numInds+i] << "\t";
		if(d.major!=d.minor)
			fprintf(pFile,"%d\t",bases[4*numInds+d.major]+bases[4*numInds+d.minor]);
			//pFile << bases[4*numInds+d.major]+bases[4*numInds+d.minor] << "\t";
		else
			fprintf(pFile,"%d\t",bases[4*numInds+d.major]);
			//pFile << bases[4*numInds+d.major] << "\t";
		//pFile << d.phat << endl;
		fprintf(pFile,"%f\n",d.phat);

		if(printCounts){
			for(int i=0;i<4*numInds-1;i++)
				fprintf(db,"%d ",bases[i]);
			fprintf(db,"%d\n",bases[4*numInds-1]);
		}  
	}
	if(printCounts)
		fclose(db);
	fflush(pFile);
}

/**
 * DATE: 2010-9-13
 * FUNCTION: allocation functions, this shouldn't look magical for anyone
 * PARAMETER: len is the requested length of double array.
 * RETURN: the pointer to the new double array.
 */
double* SfsMethod::allocDoubleArray(int len)
{
	double *ret= new double[len];
	memset(ret, 0, sizeof(double) * len);
	/*for(int i=0;i<len;i++)
		ret[i]=0;*/
	return ret;
}

/**
 * DATE: 2010-9-13
 * FUNCTION: test for infinity
 * PARAMETER: x is the double to be tested.
 * RETURN: 
 */
int SfsMethod::isinf_local(double x)
{
#ifdef __INTEL_COMPILER
	return isinf(x);
#else
	return std::isinf(x);
#endif
}


/**
 * DATE: 2010-9-13
 * FUNCTION: test for a NaN
 * PARAMETER: x is the double to be tested.
 * RETURN: 
 */
int SfsMethod::isnan_local(double x)
{
#ifdef __INTEL_COMPILER
	return isnan(x);
#else
	return std::isnan(x);
#endif
}

/**
 * DATE: 2010-9-13
 * FUNCTION: add up all the value form the ary
 * PARAMETER: ary is the double array pointer, len is the array's length.
 * RETURN: the count.
 */
double SfsMethod::sum(const double* ary, int len)
{
	double s =0;
	for(int i=0;i<len ; i++)
		s+=ary[i];
	return s;
}

/**
 * DATE: 2010-9-13
 * FUNCTION: generate like look up table.
 * PARAMETER: void
 * RETURN: void
 */
void SfsMethod::generate_likeLookup()
{
	for (short int i=0;i<=255;i++)
		likeLookup[i] =  calcLikeRatio(i); //(pow(10.0,-i/10.0));
}


/**
 * DATE: 2010-9-13
 * FUNCTION: generate binom. input should be number of individuals
 * PARAMETER: k is the individuals' number.
 * RETURN: void
 */
void SfsMethod::generate_binom(int k)
{
	binomLookup = new double[k*2+1];
	for(int i=0;i<(k*2+1);i++)
		binomLookup[i] = bico(k*2,i);
}

/**
 * DATE: 2010-9-13
 * FUNCTION :  To calculate the nCk ( n combination k)
 * PARAMETER : n is the number of individulals.  
 *			   k is the individuals' number.
 * RETURN: double
 */
double SfsMethod::bico(int n, int k)
{
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

/**
 * DATE: 2010-9-13
 * FUNCTION: generate factln
 * PARAMETER: n : number of individuals
 * RETURN: double
 */
double SfsMethod::factln(int n)
{
	static double a[101]; //A static array is automatically initialized to zero.
	if (n < 0) printf("Negative factorial in routine factln: %d \n", n);
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); 
	else return gammln(n+1.0);
}

/**
 * DATE: 2010-9-13
 * FUNCTION: generate gammln
 * PARAMETER: xx 
 * RETURN: double
 */
double SfsMethod::gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

/**
 * DATE: 2010-9-13
 * FUNCTION: get max id
 * PARAMETER: d : the double array pointer
 *			  len  £º the array's length.
 * RETURN: int.
 */
int SfsMethod::getMaxId(double *d,int len){
		int m=0;
		for(int i=1; i<len; i++)
			if(d[i]>d[m])
				m=i;
		return m;
}

/**
 * DATE: 2010-9-13
 * FUNCTION: clean map
 * PARAMETER: asso : the maps of the datums
 * RETURN: void.
 */
void SfsMethod::cleanUpMap(aMap& asso)
{
	for(aMap::iterator it=asso.begin();it!=asso.end();it++)
	{
		delete [] it->second.lk;
		delete [] it->second.locus;
	}
}

/**
 * DATE: 2010-9-13
 * FUNCTION: get char from int
 * PARAMETER: i : genetype number 
 * RETURN: genetype.
 */
char SfsMethod::getChar(int i){
	if(i==0)
		return 'a';
	else if(i==1)
		return 'c';
	else if(i==2)
		return 'g';
	else if(i==3)
		return 't';
	else
		puts("Error getting char value of base");
	return 0;
}

/**
 * DATE: 2010-9-13
 * FUNCTION:  add protect 
 * PARAMETER: ary : the double array pointer
 *			  len  £º the array's length.
 * RETURN: protect add result .
 */
double SfsMethod::addProtect(double *ary,int len){
	double maxVal = ary[0];
	for(int i=1;i<len;i++)
		if(ary[i]>maxVal)
			maxVal=ary[i];

	double sumVal = 0;
	for(int i=0;i<len;i++)
		sumVal += exp(ary[i]-maxVal);
	return log(sumVal) + maxVal;
}


/**
 * DATE: 2010-9-13
 * FUNCTION:  add three protect 
 * PARAMETER: a,b,c : double 
 *			   
 * RETURN: three protect add result .
 */
double SfsMethod::addProtect3(double a,double b, double c){
	double maxVal;// = std::max(a,std::max(b,c));
	if(a>b&&a>c)
		maxVal=a;
	else if(b>c)
		maxVal=b;
	else
		maxVal=c;
	double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
	return log(sumVal) + maxVal;
}

/**
 * DATE: 2010-9-13
 * FUNCTION:  add two protect 
 * PARAMETER: a,b : double 
 *			   
 * RETURN: two protect add result .
 */
double SfsMethod::addProtect2(double a,double b){
	double maxVal = std::max(a,b);

	double sumVal = exp(a-maxVal)+exp(b-maxVal);
	return log(sumVal) + maxVal;
}

/**
 * DATE: 2010-9-13
 * FUNCTION:  print array
 * PARAMETER: ary : the double array pointer
 *			  len  £º the array's length.
 * RETURN: void.
 */
void SfsMethod::print_array(double *ary,int len){
	for (int i=0;i<len;i++)
		printf("%f\t",(ary[i]));
	printf("\n");
}

/*
allocates the arrays needed for one site for all individuals
*/
/**
 * DATE: 2010-9-14
 * FUNCTION:  allocates the arrays needed for one site for all individuals
 * PARAMETER: numInds : the individuals' number.
 * RETURN: datum structure.
 */
datum SfsMethod::allocDatum(int numInds){
	datum data;
	//4 times the number of individuals +4 (the last 4elems will be the base sum of all inds)
	data.locus=allocIntArray(4*numInds+4);
	data.lk=allocIntArray(10*numInds);//allocUnsignedCharArray(3*numInds);
	data.phat =0.0;
//	data.position = 0;
//	data.chromo = "";
//	data.is_be_record = false;
	return data;
}

/**
 * DATE: 2010-9-14
 * FUNCTION:  allocates the arrays needed for one site for all individuals
 * PARAMETER: len : array length.
 * RETURN: array pointer.
 */
int* SfsMethod::allocIntArray(int len){
	//  printf("allocing a int array:%d\n",len);
	int *ret = NULL;
	ret = new int[len];
	if (ret == NULL)
	{
		cerr << "\tallocIntArray: new operation failed!" << endl;
		exit(0);
	}
	/*for(int i=0;i<len;i++)
		ret[i]=0;*/
	memset(ret, 0, sizeof(int) * len);
	return ret;
}

/**
 * DATE: 2010-9-14
 * FUNCTION:  will calculate the sums and input the values in the last 4 slots, 
 *				and set the minor, major accordingly
 * PARAMETER: asso : the maps of the datums
 *				numInds : the number of individuals
 * RETURN: void
 */
void SfsMethod::calcSumBias(aMap & asso,int numInds){
	// for(aMap::iterator it = asso.begin(); it != asso.end(); ++it){
	for(aMap::iterator it=asso.begin(); it!=asso.end(); it++){
		int *bases = it->second.locus; //4length array of basesums
		for(int i=0;i<numInds;i++){
			bases[4*numInds+A] += bases[i*4+A];
			bases[4*numInds+C] += bases[i*4+C];
			bases[4*numInds+G] += bases[i*4+G];
			bases[4*numInds+T] += bases[i*4+T];
		}
		//now get the major/minor

		//    int maj = it->second.major;

		int maj = 0;
		maj = it->second.ref;

		int temp=0;

		int min = maj;

		for(int i=0;i<4;i++)
		{
			if(maj==i) //we should check the against the major allele
				continue;
			else if (bases[numInds*4+i] > temp){
				min = i;
				temp = bases[numInds*4+i];
			}
		}

		if (maj==min){
			min=0;
			while (maj==min)	
				min++;
		}
		it->second.minor = min;
		it->second.major = maj;
	}
}


/**
 * DATE: 2010-9-15
 * FUNCTION:  soaplikelihood to sfslikelihood
 * PARAMETER:  type_likely:sopa likelihood. likelihood : sfs likelihood
 * RETURN:	 void
 */
void SfsMethod::soaplk_sfslk(int * likelihood, double * type_likely, int id)
{		 
	assert(type_likely != NULL && likelihood != NULL);

	char allele1;
	char allele2;
	char genotype;
	char type1 = 0;
	for (allele1 = 0; allele1 != 4; allele1++) 
	{
		// Find the largest likelihood
		for (allele2 = allele1; allele2 != 4; allele2++) 
		{
			genotype = allele1<<2|allele2;
			if (type_likely[genotype] > type_likely[type1]) 
			{
				type1 = genotype;
			}
		}
	}

	for(int i = 0; i != 10; i++) 
	{
		if(type_likely[type1] - type_likely[sfs_type[i]] > 25.5)
		{
			likelihood[10 * id + i] = 255;
		}
		else
		{
			likelihood[10 * id + i] = 10 * (type_likely[type1] - type_likely[sfs_type[i]]);// sfs_type[10]={0,1,3,2,5,7,6,15,11,10}; AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
		}
	}
}

/**
 * DATE: 2010-9-15
 * FUNCTION:  copy source's coverage to the dest
 * PARAMETER:  dest: the object array. source : the source array. id the individual's index.
 * RETURN:	 void
 */
void SfsMethod::cov_Cpy(int *dest, int const *source, int id)
{
	assert(source != NULL && dest != NULL);
	dest[4 * id + 0] = source[0];
	dest[4 * id + 1] = source[1];
	dest[4 * id + 2] = source[3];
	dest[4 * id + 3] = source[2];
}


/**
 * DATE: 2010-9-15
 * FUNCTION:  alloc vector memery
 * PARAMETER:  void
 * RETURN:	 void
 */
void SfsMethod::allocVec(void)
{
/*	//asso = new aMap();
	asso = NULL;
	asso = new aVector[global_win_size];
	if (asso == NULL)
	{
		cerr << "\tallocVec : new opreation failed!" << endl;
		exit(0);
	}
	for (int i = 0; i < global_win_size; ++i)
		asso[i] = allocDatum(m_numInds);
*/
}


/**
 * DATE: 2010-9-15
 * FUNCTION:  free the map memery
 * PARAMETER:  amap the pointer to be deleted.
 * RETURN:	 void
 */
void SfsMethod::delMap(aMap &amap)
{
	// clean up the map's node
	cleanUpMap(amap);
	amap.clear();
}


/**
 * DATE: 2010-9-15
 * FUNCTION:  do SFS
 * PARAMETER:  para: the parameters  files: Files pointer.
 * RETURN:	 void
 */
int SfsMethod::call_SFS(Parameter* para, Files* files)
{
	if (files == NULL)
	{
		return POINTER_NULL;
	}
		//update 11-`6
	sem_t * sem_call_sfs_p = &(sem_call_sfs);
	sem_t * sem_call_sfs_return_p = &(sem_call_sfs_return);
	///////////////////////time //////////////////
	//static time_t sumtime = 0;
	//time_t start = 0, end = 0;
	//
	////////////////////////////////////////////// 

	//sem_t * sem_map_number_p = &(sem_map_number);
	while (1)
	{
		sem_wait(sem_call_sfs_p);
		/*time(&start);*/
		int tmp = getidxProcess();
		mapChange();
		if (file_end_flag == 0)
		{
			sem_post(sem_call_sfs_return_p);
		}
		SFS_PARA *sfs = para->sfs;
		//int tmp = getidxProcess();
		calcSumBias(asso[tmp], m_numInds);//get major and minor
		if (sfs->doBay)
		{
			// do bay.
			algo(asso[tmp], m_numInds, sfs->eps, sfs->pvar, sfs->bias, files->sfsfile, sfs->under_FP, sfs->alternative, sfs->sigle_MB, sfs->pThres, sfs->allow_PathZ, sfs->sfs_first_last);
		}
		if (sfs->doJoint)
		{
			// do joint.
			algoJoint(asso[tmp], m_numInds, sfs->eps, sfs->pvar, sfs->bias, files->jointSfsfile, sfs->under_FP, sfs->alternative, sfs->sigle_MJ, sfs->jointFold);
		}
		if (sfs->writeFr)
		{
			// output frequence.
			writeFreq(files->freqfile, m_numInds, asso[tmp]);
		}
		// free the map.
		delMap(asso[tmp]);
		//////////////////////////////////////////
		//time(&end);
		//sumtime += end - start; 
		//clog << "sfs  tmp sumtime : "<< sumtime <<endl;
		//////////////////////////////////////////
		//sem_post(sem_map_number_p);
		if (file_end_flag == 1)
			return 0;
	}
}


// get the data what sfs want
/**
 * DATE: 2010-9-15
 * FUNCTION:  get the data what sfs want, include site's pos, chr name, depth, and the 10 likelihoods.
 * PARAMETER:  site: the allete information.
 *				mat: Prob_matrix point, include likelihoods information.
 *				chr: chromosome name.
 *				id: the individual's index.
 * RETURN:	 void
 * UPDATE: 2010-10-14 Add a pthread_mutex_lock.
 */
int SfsMethod::getMapData(const Pos_info& site, const Prob_matrix* mat, const std::string & chr, const int id)
{
	assert(id < m_numInds);
	loci tmp = {chr, site.pos + 1};
	// lock.
	pthread_mutex_lock(&m_pthreadMutex);
	aMap::iterator itr = asso[m_map_idx].find(tmp);
	int *the_line = NULL;
	int *lk = NULL;
	
	if (itr == asso[m_map_idx].end())
	{
		// new site.
		datum dats = allocDatum(m_numInds);
		// because soapsnp not the same with realsfs in base sort, snp is ACTG, sfs is ACGT.
		if (site.ori == 2)
			dats.ref = 3;
		else if (site.ori == 3)
			dats.ref = 2;
		else
			dats.ref = site.ori;
		the_line = dats.locus;
		lk = dats.lk;
		asso[m_map_idx].insert(std::make_pair(tmp, dats));
		//	asso[tmp] = dats;
	}
	else
	{
		the_line = itr->second.locus;
		lk = itr->second.lk;
		assert(itr->first.position == tmp.position);
	}
	// unlock
	pthread_mutex_unlock(&m_pthreadMutex);

	// get coverage.
	cov_Cpy(the_line, site.count_sfs, id);
	// get likelihoods.
	soaplk_sfslk(lk, mat->type_likely, id);
	
	return 0;
}





// This builds a map from: char* -> int
/**
 * DATE: 2010-9-15
 * FUNCTION:  This builds a map from: char* -> int
 *				such that each chromosome is represented as an integer
 *				input is a file looking like
 *				chr1
 *				chr2
 *				or a .fai generated file from samtools
 * PARAMETER: fname : file name
 * RETURN:	 void
 */
void SfsMethod::buildMap(const char* fname)
{
	my_ifstream pfile;
	pfile.open(fname,std::ios::in);
	char buffer[1024];
	int id=0;
	//fprintf(stderr,"\t-> Building chromosome index\n");
	cerr << "\t-> Building chromosome index" <<endl;
	while(pfile.getline(buffer,1024) ) {
		char* chr = strdup(strtok(buffer,"\t"));
		m_faiIndex.insert(std::make_pair(chr,id++));
		//fprintf(stderr,"\t-> %s->%d\n",chr,id-1);
		cerr << "\t->" << chr << "->" << id-1 <<endl;
	}
	//fprintf(stderr,"\t-> Done building chromosome index\n");
	cerr << "\t-> Done building chromosome index" << endl;
}


// delete the vector's data
/**
 * DATE: 2010-9-26
 * FUNCTION:  free the vector's data memery
 * PARAMETER:  avec the aVector to be deleted.
 * RETURN:	 void
 */
void SfsMethod::delVector(aVector * avec)
{
	for (int i = 0; i < global_win_size; i++)
	{
		if (avec[i].lk != NULL)
			delete [] avec[i].lk;
		if (avec[i].locus != NULL)
			delete [] avec[i].locus;
	}
	delete [] avec;
}


/**
 * DATE: 2010-9-26
 * FUNCTION:  get the data what sfs want, include site's pos, chr name, depth, and the 10 likelihoods.
 * PARAMETER:  site: the allete information.
 *				mat: Prob_matrix point, include likelihoods information.
 *				chr: chromosome name.
 *				id: the individual's index.
 * RETURN:	 void
 */
int SfsMethod::getSFSData(const Pos_info & site, const Prob_matrix * mat, const std::string & chr, const int id)
{
/*	assert(id < m_numInds);
	int index = site.pos % global_win_size;
	assert(index < global_win_size);

	if (!asso[index].is_be_record)
	{
		asso[index].chromo = chr;
		asso[index].is_be_record = true;
		if (site.ori == 2)
			asso[index].ref = 3;
		else if (site.ori == 3)
			asso[index].ref = 2;
		else
			asso[index].ref = site.ori;
		asso[index].position = site.pos + 1;
	}

	assert(asso[index].position == (site.pos + 1));
	// get coverage.
	cov_Cpy(asso[index].locus, site.count_all, id);
	// get likelihoods.
	soaplk_sfslk(asso[index].lk, mat->type_likely, id);
	*/
	return 0;
}


// clean up the vector's member.
/**
 * DATE: 2010-9-30
 * FUNCTION:  clean up the vector's member.
 * PARAMETER:  void
 * RETURN:	 void
 */
void SfsMethod::cleanVec(void)
{
/*	for (int i = 0; i < global_win_size; ++i)
	{
		asso[i].chromo = "";
		asso[i].is_be_record = false;
		asso[i].major = 0;
		asso[i].minor = 0;
		asso[i].ref = 0;
		asso[i].phat = 0;
		asso[i].position = 0;
		memset(asso[i].lk, 0, sizeof(int) * (10 * m_numInds));
		memset(asso[i].locus, 0, sizeof(int) * (4 * m_numInds + 4));
	}*/
}
/**
 * DATE: 2010-11-17
 * FUNCTION:  change to the next map
 * PARAMETER:  void
 * RETURN:	 void
 */
void SfsMethod::mapChange(void)
{
	// move behind another map
	m_map_idx = 1 - m_map_idx;
	
}

//get the map can be processed
/**
 * DATE: 2010-11-16
 * FUNCTION: get the index can be process
 * PARAMETER:  void
 * RETURN: the index of map can be process
 */
int SfsMethod::getidxProcess(void)
{
	while(m_map_idx_process == BE_PROCESSED) //wait the get index .
	{
		usleep(1);
	}
	int tmp;
	(m_map_idx_process == BE_PROCESSED) ? (tmp = 2) : (tmp = m_map_idx_process);
	m_map_idx_process = BE_PROCESSED;
	return tmp;
}

//set the map index
/**
 * DATE: 2010-11-16
 * FUNCTION: get the index can be process
 * PARAMETER:  void
 * RETURN: the index of map can be process
 */
void SfsMethod::setidxProcess(void)
{
	while(m_map_idx_process != BE_PROCESSED) //wait the get index .
	{
		usleep(1);
	}
	m_map_idx_process = m_map_idx ;
}

/**
 * DATE: 2010-11-16
 * FUNCTION: a thread function used to run SfsMethod::call_sfs() function.
 * PARAMETER: __Args the parameter structure.
 * RETURN: 
 */
void *_sfsMethod_callsfs(void * __Args)
{
	// get args.
	BIG_CALL_SFS_ARGS * args = (BIG_CALL_SFS_ARGS*)__Args;
	// run the function.
	args->sfsMethod->call_SFS(args->para, args->files);
	return NULL;
}
