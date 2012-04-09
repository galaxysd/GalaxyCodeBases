//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22
//function: 在指定insert size的条件下，将一条序列与adapter进行比较
//默认的比对标准：最小比对长度为5bp, 最大错配率为10%,即10个碱基最多允许1个错配。


#include "locate_adapter.h"

using namespace std;


//将所给序列与adapter序列进行比对，如果能比上则返回1，比不上反回0。
int locate_adapter ( string &seq, string adapter, int insert_size, int &align_len, int &mis_match )
{	
	int min_align_len = 5;
	float max_mismatch_rate = 0.1;

	//如果insert_size为负数，则截取部分后一部分adapter，并将insert_size重置为0
	int abs_insert_size = abs(insert_size);
	if (insert_size < 0 && abs_insert_size < adapter.size())
	{	adapter = adapter.substr(abs_insert_size,adapter.size()-abs_insert_size);
		insert_size = 0;
	}
	
	//insert_size统一作正数处理
	int is_adapter = 0;
	mis_match = 0;
	align_len = ( insert_size + adapter.size() <= seq.size() ) ? adapter.size() : ( seq.size() - insert_size );
	int max_mismatch_num = align_len * max_mismatch_rate; //最大允许的错配碱基数
	
	if (align_len >= min_align_len) //当符合最小比对长度时才操作
	{	for (int i=0; i<align_len; i++)
		{	if ( seq[insert_size+i] != adapter[i] )
			{	mis_match++;
				if (mis_match > max_mismatch_num)
				{	break;
				}
			}
		}
	}
	
	if (mis_match <= max_mismatch_num)
	{	is_adapter = 1;
	}
	
	return is_adapter;

}

