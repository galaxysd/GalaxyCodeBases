#include "pairs.h"

using namespace std;

extern Param param;

PairAlign::PairAlign()
{
	n_filtered_pairs=n_filtered_a=n_filtered_b=0;
	SetFlag();
}
void PairAlign::SetFlag()
{
	_sa.SetFlag('a');
	_sb.SetFlag('b');
}
void PairAlign::ImportFileFormat(int format1, int format2)
{
	_sa.ImportFileFormat(format1);
	_sb.ImportFileFormat(format2);
}
void PairAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2)
{
	_sa.ImportBatchReads(n, a1);
	_sb.ImportBatchReads(n, a2);
	num_reads=n;
}
int PairAlign::GetExactPairs()
{
	int i, j, k;
	PairHit pp;
	pp.na=pp.nb=0;
	//a+, b-
	pp.chain=0;
	if(_sa._cur_n_hit[0] && _sb._cur_n_chit[0]) {
		i=j=0;
		while(i<_sa._cur_n_hit[0]) {
			while((j<_sb._cur_n_chit[0]) &&((_sb.chits[0][j].chr<_sa.hits[0][i].chr) 
				||((_sb.chits[0][j].chr==_sa.hits[0][i].chr) &&(_sb.chits[0][j].loc+_sb._pread->seq.size()<_sa.hits[0][i].loc+param.min_insert))))
				j++;
			k=j;
			while((k<_sb._cur_n_chit[0]) &&(_sb.chits[0][k].chr==_sa.hits[0][i].chr) 
				&&(_sb.chits[0][k].loc+_sb._pread->seq.size()<=_sa.hits[0][i].loc+param.max_insert)) {
				if(_cur_n_hits[0]>=MAXHITS)
					return 0;
				pp.a=_sa.hits[0][i];
				pp.b=_sb.chits[0][k];
				pairhits[0][_cur_n_hits[0]++]=pp;
				k++;
			}
			i++;
		}
	}
	//a-, b+
	pp.chain=1;
	if(_sa._cur_n_chit[0] && _sb._cur_n_hit[0]) {
		i=j=0;
		while(i<_sa._cur_n_chit[0]) {
			while((j<_sb._cur_n_hit[0]) &&((_sb.hits[0][j].chr<_sa.chits[0][i].chr) 
				||((_sb.hits[0][j].chr==_sa.chits[0][i].chr) &&(_sb.hits[0][j].loc<_sa.chits[0][i].loc+_sa._pread->seq.size()-param.max_insert))))
				j++;
			k=j;
			while((k<_sb._cur_n_hit[0]) &&(_sb.hits[0][k].chr==_sa.chits[0][i].chr) 
				&&(_sb.hits[0][k].loc<=_sa.chits[0][i].loc+_sa._pread->seq.size()-param.min_insert)) {
				if(_cur_n_hits[0]>=MAXHITS)
					return 0;				
				pp.a=_sa.chits[0][i];
				pp.b=_sb.hits[0][k];
				pairhits[0][_cur_n_hits[0]++]=pp;
				k++;
			}
			i++;
		}
	}
	if(_cur_n_hits[0]>0)
		return 0;
	return -1;	
}
int PairAlign::GetExact2SnpPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,nsnp;
	//a+ vs b-
	if(_sa._cur_n_hit[0]) {
		pp.chain=0;
		pp.na=0;
		for(i=0; i<_sa._cur_n_hit[0]; i++) {
			nsnp=_sb.SnpAlign_range(1, _sa.hits[0][i].chr, _sa.hits[0][i].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[0][i].loc+param.max_insert-_sb._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.a=_sa.hits[0][i];
				for(j=0; j<_sb._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=MAXHITS)
						break;
					pp.b=_sb.bound_hits[nsnp][j];
					pp.nb=nsnp;					
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}
		}
	}
	//a- vs b+
	if(_sa._cur_n_chit[0]) {
		pp.chain=1;
		pp.na=0;
		for(i=0; i<_sa._cur_n_chit[0]; i++) {
			nsnp=_sb.SnpAlign_range(0, _sa.chits[0][i].chr, _sa.chits[0][i].loc-param.max_insert+_sa._pread->seq.size(), _sa.chits[0][i].loc-param.min_insert+_sa._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.a=_sa.chits[0][i];
				for(j=0; j<_sb._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=MAXHITS)
						break;		
					pp.b=_sb.bound_hits[nsnp][j];
					pp.nb=nsnp;			
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}	
	}
	//b+ vs a-
	if(_sb._cur_n_hit[0]) {
		pp.chain=1;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_hit[0]; i++) {
			nsnp=_sa.SnpAlign_range(1, _sb.hits[0][i].chr, _sb.hits[0][i].loc+param.min_insert-_sa._pread->seq.size(), _sb.hits[0][i].loc+param.max_insert-_sa._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.b=_sb.hits[0][i];
				for(j=0; j<_sa._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=MAXHITS)
						break;	
					pp.a=_sa.bound_hits[nsnp][j];
					pp.na=nsnp;				
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}		
	}
	//b- vs a+
	if(_sb._cur_n_chit[0]) {
		pp.chain=0;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_chit[0]; i++) {
			nsnp=_sa.SnpAlign_range(0, _sb.chits[0][i].chr, _sb.chits[0][i].loc-param.max_insert+_sb._pread->seq.size(), _sb.chits[0][i].loc-param.min_insert+_sb._pread->seq.size(), ref);
			if(-1 !=nsnp) {
				pp.b=_sb.chits[0][i];
				for(j=0; j<_sa._cur_n_boundhit[nsnp]; j++) {
					if(_cur_n_hits[nsnp]>=MAXHITS)
						break;
					pp.a=_sa.bound_hits[nsnp][j];
					pp.na=nsnp;					
					pairhits[nsnp][_cur_n_hits[nsnp]++]=pp;
				}
			}			
		}		
	}
	for(i=0; i<=param.max_snp_num; i++) {
		if(_cur_n_hits[i]>0)
			return i;
	}
	return -1;
}
int PairAlign::GetSnp2SnpPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,k,h,nsnp;
	for(i=1; i<=param.max_snp_num; i++) {
		//a+ vs b-
		if(_sa._cur_n_hit[i]) {
  		pp.chain=0;
  		pp.na=i;			
  		for(j=0; j<_sa._cur_n_hit[i]; j++) {
  			nsnp=_sb.SnpAlign_range(1, _sa.hits[i][j].chr, _sa.hits[i][j].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[i][j].loc+param.max_insert-_sb._pread->seq.size(), ref);
  			if(-1 !=nsnp) {
  				pp.a=_sa.hits[i][j];
  				for(k=0; k<_sb._cur_n_boundhit[nsnp]; k++) {
  					if(_cur_n_hits[nsnp+i]>=MAXHITS)
  						break;					
  					pp.b=_sb.bound_hits[nsnp][k];
  					pp.nb=nsnp;
  					pairhits[nsnp+i][_cur_n_hits[nsnp+i]++]=pp;
  				}
  			}		
  		}
  	}
		//a- vs b+
		if(_sa._cur_n_chit[i]) {
  		pp.chain=1;
  		pp.na=i;			
  		for(j=0; j<_sa._cur_n_chit[i]; j++) {
  			nsnp=_sb.SnpAlign_range(0, _sa.chits[i][j].chr, _sa.chits[i][j].loc-param.max_insert+_sa._pread->seq.size(), _sa.chits[i][j].loc-param.min_insert+_sa._pread->seq.size(), ref);
  			if(-1 !=nsnp) {
  				pp.a=_sa.chits[i][j];
  				for(k=0; k<_sb._cur_n_boundhit[nsnp]; k++) {
  					if(_cur_n_hits[nsnp+i]>=MAXHITS)
  						break;						
  					pp.b=_sb.bound_hits[nsnp][k];
  					pp.nb=nsnp;
  					pairhits[nsnp+i][_cur_n_hits[nsnp+i]++]=pp;
  				}
  			}			
  		}
  	}
		for(h=2; h<=i*2; h++) {
			if(_cur_n_hits[h]>0)
				return h;
		}	
	}
	return -1;
}
int PairAlign::GetExact2GapPairs(RefSeq &ref)
{
	PairHit pp;
	int i,j,k,h,gapsize;
	//a+ vs b-
	if(_sa._cur_n_hit[0]) {
		pp.chain=0;
		pp.na=0;
		for(i=0; i<_sa._cur_n_hit[0]; i++) {
			gapsize=_sb.GapAlign_range(1, _sa.hits[0][i].chr, _sa.hits[0][i].loc+param.min_insert-_sb._pread->seq.size(), _sa.hits[0][i].loc+param.max_insert-_sb._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.a=_sa.hits[0][i];
				for(j=0; j<_sb._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=MAXHITS)
						break;					
					pp.b=_sb.bound_gaphits[j];
					pp.nb=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}
	}
	//a- vs b+
	if(_sa._cur_n_chit[0]) {
		pp.chain=1;
		pp.na=0;
		for(i=0; i<_sa._cur_n_chit[0]; i++) {
			gapsize=_sb.GapAlign_range(0, _sa.chits[0][i].chr, _sa.chits[0][i].loc-param.max_insert+_sa._pread->seq.size(), _sa.chits[0][i].loc-param.min_insert+_sa._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.a=_sa.chits[0][i];
				for(j=0; j<_sb._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=MAXHITS)
						break;						
					pp.b=_sb.bound_gaphits[j];
					pp.nb=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}		
	}
	//b+ vs a-
	if(_sb._cur_n_hit[0]) {
		pp.chain=1;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_hit[0]; i++) {
			gapsize=_sa.GapAlign_range(1, _sb.hits[0][i].chr, _sb.hits[0][i].loc+param.min_insert-_sa._pread->seq.size(), _sb.hits[0][i].loc+param.max_insert-_sa._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.b=_sb.hits[0][i];
				for(j=0; j<_sa._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=MAXHITS)
						break;						
					pp.a=_sa.bound_gaphits[j];
					pp.na=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}
	}
	//b- vs a+
	if(_sb._cur_n_chit[0]) {
		pp.chain=0;
		pp.nb=0;
		for(i=0; i<_sb._cur_n_chit[0]; i++) {
			gapsize=_sa.GapAlign_range(0, _sb.chits[0][i].chr, _sb.chits[0][i].loc-param.max_insert+_sb._pread->seq.size(), _sb.chits[0][i].loc-param.min_insert+_sb._pread->seq.size(), ref);
			if(-1!=gapsize) {
				pp.b=_sb.chits[0][i];
				for(j=0; j<_sa._cur_n_boundgaphit; j++) {
					if(_cur_n_gaphits[gapsize]>=MAXHITS)
						break;						
					pp.a=_sa.bound_gaphits[j];
					pp.na=gapsize;
					pairgaphits[gapsize][_cur_n_gaphits[gapsize]++]=pp;
				}
			}
		}		
	}
	for(i=1; i<=MAXGAP; i++) {
		if(_cur_n_gaphits[i])
			return i;
	}
	return -1;
}
int PairAlign::RunAlign(RefSeq &ref)
{
	int i;
	for(i=0; i<=param.max_snp_num*2; i++)
		_cur_n_hits[i]=0;
	for(i=0; i<=MAXGAP; i++)
		_cur_n_gaphits[i]=0;
	_sa.ClearHits();
	_sb.ClearHits();
	_sa.ConvertBinaySeq();
	_sb.ConvertBinaySeq();
	//get exact+exact pairs
	_sa.GenerateSeeds_1(0);
	_sb.GenerateSeeds_1(0);
	_sa.SnpAlign_0(ref);
	_sb.SnpAlign_0(ref);
	if(((_sa._cur_n_hit[0]&&_sb._cur_n_chit[0]) ||(_sa._cur_n_chit[0]&&_sb._cur_n_hit[0])) &&(-1!=GetExactPairs()))
		return 1;
	//get exact+snp pairs
	_sa.GenerateSeeds_1(1);
	_sa.GenerateSeeds_1(2);
	_sa.GenerateSeeds_2(3);
	_sa.GenerateSeeds_2(5);
	_sa.GenerateSeeds_3(4);
	_sb.GenerateSeeds_1(1);
	_sb.GenerateSeeds_1(2);
	_sb.GenerateSeeds_2(3);
	_sb.GenerateSeeds_2(5);
	_sb.GenerateSeeds_3(4);	
	if(-1!=GetExact2SnpPairs(ref))
		return 2;
	//snp alignment for a, then do snp align for b in the flanking region, get pairs
	_sa.SnpAlign_1(ref);
	_sa.SnpAlign_2(ref);
	if(-1!=GetSnp2SnpPairs(ref))
	{
		return 3;
	}	
	//gap align in flanking region and get exact+gap pairs
	if(param.max_gap_size &&(_sa._cur_n_hit[0] ||_sb._cur_n_chit[0] ||_sa._cur_n_chit[0] ||_sb._cur_n_hit[0]) &&(-1!=GetExact2GapPairs(ref)))
	{
		return 4;
	}
	return 0;		
}
void PairAlign::Do_Batch(RefSeq &ref)
{
	_str_align.clear();
	int tt=0;
	int filter1, filter2;
	for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
		filter1=_sa.FilterReads();
		filter2=_sb.FilterReads();
		if(filter1 && filter2) {
			n_filtered_pairs++;
			continue;
		}
		else if(filter1) {
			n_filtered_a++;
			if(_sb.RunAlign(ref))
				_sb.StringAlign(ref, _str_align);
		}
		else if(filter2) {
			n_filtered_b++;
			if(_sa.RunAlign(ref))
				_sa.StringAlign(ref, _str_align);
		}
		else
			if(RunAlign(ref))
				StringAlign(ref, _str_align);
	}
//	cout<<_str_align<<endl;
}
void PairAlign::StringAlign(RefSeq &ref, string &os)
{
	_sa.Reverse_Seq();
	_sa.Reverse_Qual();
	_sb.Reverse_Seq();
	_sb.Reverse_Qual();	
	int i, j;
	//snp hits
	for(i=0; i<=param.max_snp_num*2; i++) {
		if(0==_cur_n_hits[i])
			continue;
		if(1==_cur_n_hits[i]) {
			if(pairhits[i][0].na==0)
				_sa.s_OutHit(pairhits[i][0].chain, 1, pairhits[i][0].na, &pairhits[i][0].a, 1, ref, os);
			else
				_sa.s_OutHit(pairhits[i][0].chain, 1, pairhits[i][0].na, &pairhits[i][0].a, 1, ref, os);
			if(pairhits[i][0].nb==0)
				_sb.s_OutHit(!pairhits[i][0].chain, 1, pairhits[i][0].nb, &pairhits[i][0].b, 1, ref, os);
			else
				_sb.s_OutHit(!pairhits[i][0].chain, 1, pairhits[i][0].nb, &pairhits[i][0].b, 1, ref, os);			
		}
		else if(1==param.report_repeat_hits) {   //randomly pick up one
			j=rand()%_cur_n_hits[i];
			_sa.s_OutHit(pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].na, &pairhits[i][j].a, 1, ref, os);
			_sb.s_OutHit(!pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].nb, &pairhits[i][j].b, 1, ref, os);
		}
		else if(2==param.report_repeat_hits) {   //output all repeat hits
			for(j=0; j<_cur_n_hits[i]; j++) {
				_sa.s_OutHit(pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].na, &pairhits[i][j].a, 1, ref, os);
				_sb.s_OutHit(!pairhits[i][j].chain, _cur_n_hits[i], pairhits[i][j].nb, &pairhits[i][j].b, 1, ref, os);				
			}
		}
		return;
	}
	//gap hits
	for(i=1; i<=MAXGAP; i++) {
		if(0==_cur_n_gaphits[i])
			continue;
		if(1==_cur_n_gaphits[i]) {
			_sa.s_OutGapHit(pairgaphits[i][0].chain, _cur_n_gaphits[i], i, &pairgaphits[i][0].a, ref, os);
			_sb.s_OutGapHit(!pairgaphits[i][0].chain, _cur_n_gaphits[i], i, &pairgaphits[i][0].b, ref, os);
		}
		else if(1==param.report_repeat_hits) {
			j=rand()%_cur_n_gaphits[i];
			_sa.s_OutGapHit(pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].a, ref, os);
			_sb.s_OutGapHit(!pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].b, ref, os);
		}
		else if(2==param.report_repeat_hits) {
			for(j=0; j<_cur_n_gaphits[i]; j++) {						
				_sa.s_OutGapHit(pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].a, ref, os);
				_sb.s_OutGapHit(!pairgaphits[i][j].chain, _cur_n_gaphits[i], i, &pairgaphits[i][j].b, ref, os);
			}
		}
		return;
	}
}

