#include <iostream>
#include <map>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <boost/progress.hpp>
#include <boost/progress.hpp>
#include "gzstream.h"
#include "self_util.h"
#include "stat_soap_coverage.h"
using namespace std;

stat_soap_coverage::stat_soap_coverage(string str_ref_file_name, 
        string str_output_prefix, vector<string> vec_soap_file_name, 
        vector<string> vec_width, bool b_gcdump, bool b_depwindump)
{
    this->str_ref_file_name = str_ref_file_name;
    this->str_output_prefix = str_output_prefix;
    this->vec_soap_file_name = vec_soap_file_name;
    this->vec_width = vec_width;
    this->b_gcdump = b_gcdump;
    this->b_depwindump = b_depwindump;

    Run();
}

stat_soap_coverage::~stat_soap_coverage()
{

}

void stat_soap_coverage::DealReference()
{
    boost::progress_timer timer;
    igzstream in(str_ref_file_name.c_str());
    string line;
    string keyname = "";
    string sequence = "";
    //uint64_t countLine = 0;


    while(getline(in, line))
    {
        TrimLeft(line);
        TrimRight(line);
        if(line[0] == '>')
        {
            if(sequence.length() != 0)
            {
                map_reference_base[keyname] = sequence;
                cerr << "map_reference_base:" << keyname << ":" << sequence.length() << endl;
                vec_chr_keyname.push_back(keyname);
            }

            int index;
            if(((index = line.find(" ")) != string::npos) || ((index = line.find("\t")) != string::npos) || ((index = line.find("\n")) != string::npos))
            {
                keyname = line.substr(1, index-1);
            }
            else
            {
                keyname = line.substr(1);
            }

            sequence.clear();
            continue;
        }

        sequence += line;
    }

    if(sequence.length() != 0)
    {
        map_reference_base[keyname] = sequence;
        vec_chr_keyname.push_back(keyname);
    }
    
    in.close();
    cout << "deal reference time: " ;
}

void stat_soap_coverage::DealSoapCoverage()
{
  
    boost::progress_timer timer;
    boost::progress_display display(vec_soap_file_name.size());

    for(vector<string>::iterator it = vec_soap_file_name.begin(); it != vec_soap_file_name.end(); ++it)
    {
        igzstream in((*it).c_str());
        string line;
        string keyname;
        string strTemp = "";
        vector<unsigned int> chr_soap_coverage;
        //uint64_t countLine = 0;
        int in_temp;



        while(getline(in, line))
        {
            TrimLeft(line);
            TrimRight(line);
            
            if(line[0] == '>')
            {
                if(chr_soap_coverage.size() == 0)
                {
                    int index;
                    if(((index = line.find(" ")) != string::npos) || ((index = line. find("\t")) != string::npos) || ((index = line.find("\n")) != string::npos))
                    {
                        keyname = line.substr(1, index-1);    
                    }
                    else
                    {
                        keyname = line.substr(1);
                    }
                    for(unsigned int i=0; i<map_reference_base[keyname].size(); ++i)
                    {
                        
                        in >> in_temp;
                        if(in_temp == 65535)
                        {
                            in_temp = 0;
                        }
                        chr_soap_coverage.push_back((unsigned int)in_temp);
                    }
                    if(map_soap_coverage.count(keyname) == 0)
                    {
                        map_soap_coverage[keyname] = chr_soap_coverage;
                    }
                    else
                    {
                        for(unsigned int i=0; i<map_soap_coverage[keyname].size(); ++i)
                        {
							if ( uint64_t(map_soap_coverage[keyname][i]) + chr_soap_coverage[i] > uint64_t(UINT32_MAX) ) {
								map_soap_coverage[keyname][i] = UINT32_MAX;
							} else {
								map_soap_coverage[keyname][i] += chr_soap_coverage[i];
							}
/*
                            uint64_t temp = map_soap_coverage[keyname][i] + chr_soap_coverage[i];
                            if(temp > 65534)
                            {
                                map_soap_coverage[keyname][i] = 65534;
                            }
                            else
                            {
                                map_soap_coverage[keyname][i] += chr_soap_coverage[i];
                            }
*/
                        }
                    }

                  
                }

                keyname = "";
                chr_soap_coverage.clear();

            }

        }
        in.close();
        ++display;
    }
    cout << "deal soapcoverage time: ";
}

void stat_soap_coverage::StatGC()
{
    boost::progress_timer timer;
    boost::progress_display display(vec_width.size());
    bool b_depth = true;
    for(unsigned int i=0; i<vec_width.size(); i++)
    {
        unsigned int width = toInt(vec_width[i]);
        map<double, vector<double> > map_soap_gc_depth;
        vector<double> gc_keyname;
        map<double, uint64_t> map_temp_wincount;

        ogzstream gzgcdump, gzdepwindump;
        winCountN[width] = 0;
        if(b_gcdump)
        {
            gzgcdump.open((str_output_prefix+"_"+toStr(width)+".refgc.gz").c_str());
        }

        if(b_depwindump)
        {
            gzdepwindump.open((str_output_prefix+"_"+toStr(width)+".windep.gz").c_str());
        }


        for(map<string, string>::iterator it = map_reference_base.begin(); it != map_reference_base.end(); ++it)
        {

            string keyname = it->first;
            if(b_gcdump)
            {
                gzgcdump << ">" << keyname <<endl;
            }

            if(b_depwindump)
            {
                gzdepwindump << ">" << keyname << endl;
            }
            uint64_t countPos = 0;
            int countGC = 0;
            int countN = 0;
            uint64_t sumBase = 0;

            for(unsigned int j=0; j<it->second.length(); ++j)
            {
                countPos++;

                if(b_depth)
                {
                    if(map_stat_depth.count(map_soap_coverage[keyname][j]) == 0)
                    {
                        map<string, uint64_t> temp_depth;
                        for(unsigned int chr_key = 0; chr_key < vec_chr_keyname.size(); ++chr_key)
                        {
                            temp_depth[vec_chr_keyname[chr_key]] = 0;
                        }
                        temp_depth[keyname] = 1;
                        map_stat_depth[map_soap_coverage[keyname][j]] = temp_depth;
                    }
                    else
                    {
                        map_stat_depth[map_soap_coverage[keyname][j]][keyname] += 1;
                    }
                }



                if(countPos < width)
                {
                    if((it->second[j] == 'N') || (it->second[j] == 'n'))
                    {
                        countN++;
                    }
                    else
                    {
                        sumBase += map_soap_coverage[keyname][j];
                    }

                    if((it->second[j] == 'G') || (it->second[j] == 'C') || (it->second[j] == 'g') || (it->second[j] == 'c'))
                    {
                        countGC++;
                    }
                }
                else
                {
                    if(map_sumwincount.count(width) == 0)
                    {
                        map_sumwincount[width] = 1;
                    }
                    else
                    {
                        map_sumwincount[width] += 1;
                    }

                    if((width-countN) != 0)
                    {
                        double rate = (double)countGC/(width-countN) * 100;
                        double key = int(rate) + 0.5;
                        if(map_temp_wincount.count(key) == 0)
                        {
                            map_temp_wincount[key] = 1;
                        }
                        else
                        {
                            map_temp_wincount[key] += 1;
                        }
                    }
                    else
                    {
                        winCountN[width] += 1;
                    }

                    if(((double)countN/width >= 0.75) || ((width - countN) <= 30))
                    {
                        countGC = 0;
                        sumBase = 0;
                        countN = 0;
                        if(b_gcdump)
                        {
                            gzgcdump << -1 << endl;
                        }

                        if(b_depwindump)
                        {
                            gzdepwindump << -1 << endl;
                        }

                    }
                    else
                    {
                        double gc_rate = (double)countGC/(width-countN) * 100;
                        if(b_gcdump)
                        {
                            gzgcdump << gc_rate << endl;
                        }

                        double gc_key = int(gc_rate) + 0.5;
                        
                        if(map_sumdepthcount.count(width) == 0)
                        {
                            map_sumdepthcount[width] = 1;
                        }
                        else
                        {
                            map_sumdepthcount[width] += 1;
                        }

                        if(map_soap_gc_depth.count(gc_key) == 0)
                        {
                            vector<double> temp;
                            temp.push_back((double)sumBase/(width-countN));
                            map_soap_gc_depth[gc_key] = temp;
                            gc_keyname.push_back(gc_key);
                        }
                        else
                        {
                            map_soap_gc_depth[gc_key].push_back((double)sumBase/(width-countN));
                        }

                        if(b_depwindump)
                        {
                            gzdepwindump << ((double)sumBase/(width-countN)) << endl;
                        }

                        countGC = 0;
                        sumBase = 0;
                        countN = 0;
                    }

                    if((it->second[j] == 'N') || (it->second[j] == 'n'))
                    {
                        countN++;
                    }
                    else
                    {
                        sumBase += map_soap_coverage[keyname][j];
                    }
                    
                    if((it->second[j] == 'G') || (it->second[j] == 'C') || (it->second[j] == 'g') || (it->second[j] == 'c'))
                    {
                        countGC++;
                    }
                    countPos = 1;
                }
            }
        }

        if(b_gcdump)
        {
            gzgcdump.close();
        }
        
        if(b_depwindump)
        {
            gzdepwindump.close();
        }

        map_wincount[width] = map_temp_wincount;
        b_depth = false;
        if(map_width_soap_gc_depth.count(width) == 0)
        {
            map_width_soap_gc_depth[width] = map_soap_gc_depth;
            map_gc_keyname[width] = gc_keyname;
        }
        else
        {
            cerr << "error !" << __FILE__ << ", " << __LINE__ <<endl;
            exit(EXIT_FAILURE);
        }

        ++display;
    }
    cout << "stat GC time: ";
}

void stat_soap_coverage::StatCoverage()
{
    boost::progress_timer timer;
    boost::progress_display display(map_reference_base.size());
    uint64_t all_countN = 0;
    uint64_t all_coverageNum = 0;
    uint64_t all_sumBase = 0;
    uint64_t all_sum = 0;
    for(map<string, string>::iterator it = map_reference_base.begin(); it != map_reference_base.end(); ++it)
    {
        string keyname = it->first;
        uint64_t countN = 0;
        uint64_t coverageNum = 0;
        uint64_t sumBase = 0;

        for(unsigned int i=0; i<it->second.length(); ++i)
        {
            all_sum++;
            if((it->second[i] == 'N') || (it->second[i] == 'n'))
            {
                all_countN++;
                countN++;
                continue;
            }
            
            if(map_soap_coverage[keyname][i] != 0)
            {
                all_coverageNum++;
                coverageNum++;
                sumBase += map_soap_coverage[keyname][i];
                all_sumBase += map_soap_coverage[keyname][i];
            }
        }

        if(map_stat_coverage.count(keyname) == 0)
        {
            vector<double> temp;
            temp.push_back((double)sumBase/(it->second.length()-countN));
            temp.push_back((double)coverageNum);
            temp.push_back((double)coverageNum/(it->second.length()-countN));
            temp.push_back((double)(it->second.length()-countN));
            temp.push_back((double)countN);
            temp.push_back((double)it->second.length());
            map_stat_coverage[keyname] = temp;
        }
        else
        {
            cerr << "error !" << __FILE__ << ", " << __LINE__<< endl;
            exit(EXIT_FAILURE);
        }
        ++display;
    }

    if(map_stat_coverage.count("_All_") == 0)
    {
        vector<double> temp;
        temp.push_back((double)all_sumBase/(all_sum - all_countN));
        temp.push_back((double)all_coverageNum);
        temp.push_back((double)all_coverageNum/(all_sum - all_countN));
        temp.push_back((double)(all_sum - all_countN));
        temp.push_back((double)all_countN);
        temp.push_back((double)all_sum);
        map_stat_coverage["_All_"] = temp;
    }
    cout << "stat coverage time: ";
}

void stat_soap_coverage::DealStat()
{
    cout << "statGC begin!" << endl;
    StatGC();
    cout << "statGC end!" << endl;
    cout << "statCoverage begin!" << endl;
    StatCoverage();
    cout << "statCoverage end!" << endl;
    
    boost::progress_timer timer;
    boost::progress_display display(map_width_soap_gc_depth.size());
    for(map<int, map<double, vector<double> > >::iterator it = map_width_soap_gc_depth.begin(); it != map_width_soap_gc_depth.end(); ++it)
    {
        int width = it->first;
        ofstream out((str_output_prefix + "_" + toStr(width) + ".dat").c_str());
        
        if(!out)
        {
            cerr << "can't open the file " << (str_output_prefix + "_" + toStr(width) + ".dat") << endl;
            exit(EXIT_FAILURE);
        }
        vector<double> gc_keyname;
        gc_keyname = map_gc_keyname[width];
        sort(gc_keyname.begin(), gc_keyname.end());
        map<double, vector<double> > temp_gc_output;

        for(unsigned int i=0; i<gc_keyname.size(); ++i)
        {
            vector<double> temp = map_width_soap_gc_depth[width][gc_keyname[i]];
            //double gc_rate = gc_keyname[i];
            uint64_t ref_count =  map_width_soap_gc_depth[width][gc_keyname[i]].size();
            double sum_coverage = 0;

            for(unsigned int j=0; j<temp.size(); ++j)
            {
                sum_coverage += temp[j];
            }

            double avg_value = sum_coverage/temp.size();
            double min_value = *min_element(temp.begin(), temp.end());
            double max_value = *max_element(temp.begin(), temp.end());
            sort(temp.begin(), temp.end());
            double Q1, Q2, Q3;
            if(temp.size() % 2 == 0)
            {
                if(temp.size() < 4)
                {
                    Q1 = 0;
                    Q2 = 0;
                    Q3 = 0;
                }
                else
                {
                    Q1 = temp[int((temp.size()+1)/4)-1] + (temp[int((temp.size()+1)/4)] - temp[int((temp.size()+1)/4)-1]) * (((double)(temp.size()+1)/4)-(int((temp.size()+1)/4)));
                    Q2 = temp[int((temp.size()+1)/2) - 1] + (temp[int((temp.size()+1)/2)] - temp[int((temp.size()+1)/2) - 1]) * (((double)(temp.size()+1)/2)-(int((temp.size()+1)/2)));
                    Q3 = temp[int((temp.size()+1)/4*3) - 1] + (temp[int((temp.size()+1)/4*3)] - temp[int((temp.size()+1)/4*3)]) * (((double)(temp.size()+1)/4*3)-(int((temp.size()+1)/4*3)));
            
                }
            }
            else
            {
                if(temp.size() < 3)
                {
                    Q1 = 0;
                    Q2 = 0;
                    Q3 = 0;
                }
                else
                {
                    Q1 = temp[int(temp.size() + 1)/4 - 1];
                    Q2 = temp[int(temp.size() + 1)/2 - 1];
                    Q3 = temp[int(temp.size() + 1)/4 * 3 - 1];
                }
            }
            double small_value, big_value;
            double small_temp_value = Q1-1.5*(Q3-Q1);
            double big_temp_value = Q3+1.5*(Q3-Q1);

            for(unsigned int small_i = 0; small_i < temp.size(); small_i++)
            {
                if(small_temp_value < temp[small_i])
                {
                    small_value = temp[small_i];
                    break;
                }
                small_value = temp[small_i];
            }

            for(unsigned int big_i=0; big_i < temp.size();big_i++)
            {
                if(big_temp_value < temp[big_i])
                {
                    if(big_i == 0)
                        big_value = temp[big_i];
                    else
                        big_value = temp[big_i - 1];
                    break;
                }
                big_value = temp[big_i];
            }

            vector<double> temp_output;
            temp_output.push_back(double(ref_count));
            temp_output.push_back(double(avg_value));
            temp_output.push_back(small_value);
            temp_output.push_back(double(Q1));
            temp_output.push_back(double(Q2));
            temp_output.push_back(double(Q3));
            temp_output.push_back(big_value);
            temp_output.push_back(double(min_value));
            temp_output.push_back(double(max_value));
            
            temp_gc_output[gc_keyname[i]] = temp_output;
        }

        double sum_avg = 0;
        double sum_ref_count = 0;
        for(unsigned int i=0; i<gc_keyname.size(); ++i)
        {
            sum_avg += temp_gc_output[gc_keyname[i]][1];
            sum_ref_count += temp_gc_output[gc_keyname[i]][0];
        }

        double k = sum_avg/sum_ref_count;
        out << "#WinSize=" << width << "\tWinCount=" << map_sumwincount[width] << "\tDepthCount=" << map_sumdepthcount[width] << endl
            << "#All-N windows count: " << winCountN[width] << endl
            << "#GC%\tRefCnt\tDepthCnt\tMean\tSmall\tQ1\tMid\tQ3\tBig\tMin\tMax\tRefcntcal"
            << endl;
        for(unsigned int i=0; i<gc_keyname.size(); ++i)
        {
            out << gc_keyname[i] << "\t" << map_wincount[width][gc_keyname[i]] 
                << "\t" << uint64_t(temp_gc_output[gc_keyname[i]][0])
                << "\t" << temp_gc_output[gc_keyname[i]][1]
                << "\t" << temp_gc_output[gc_keyname[i]][2]
                << "\t" << temp_gc_output[gc_keyname[i]][3]
                << "\t" << temp_gc_output[gc_keyname[i]][4]
                << "\t" << temp_gc_output[gc_keyname[i]][5]
                << "\t" << temp_gc_output[gc_keyname[i]][6]
                << "\t" << temp_gc_output[gc_keyname[i]][7]
                << "\t" << temp_gc_output[gc_keyname[i]][8]
                << "\t" << temp_gc_output[gc_keyname[i]][0]*k

                << endl;
        }

        out.close();
        ++display;
    }

    ofstream log((str_output_prefix + "_stat" + ".dat").c_str());
    if(!log)
    {
        cerr << "can't open file " << (str_output_prefix + "_stat" + ".dat") << __FILE__ << ", " << __LINE__ << endl;
        exit(EXIT_FAILURE);
    }
    
    log << "#chrid\tdepth\tcovered\tcvgratio\tchrlen_no_N\tNzone\tchrlen" << endl;
    for(map<string, vector<double> >::iterator it=map_stat_coverage.begin(); it!=map_stat_coverage.end(); ++it)
    {
        log << it->first << "\t" << it->second[0] << "\t" << uint64_t(it->second[1]) << "\t" << it->second[2] << "\t" << uint64_t(it->second[3]) << "\t" << uint64_t(it->second[4]) << "\t" << uint64_t(it->second[5]) << endl;
    }

    log.close();
    cout << "deal stat time: ";
}

void stat_soap_coverage::Run()
{
    cout << "dealReference begin!" << endl;
    DealReference();
    cout << "dealReference end!" << endl;
    cout << "dealDealSoapCoverage begin!" << endl;
    DealSoapCoverage();
    cout << "dealDealSoapCoverage end!" << endl;
    DealStat();
    statDepth();
}


void stat_soap_coverage::statDepth()
{
    cout << "stat depth time: ";
    boost::progress_timer timer;
    ofstream out((str_output_prefix + "_" + "stat.depth").c_str());

    out << "#Depth\t_All_";
    for(unsigned int j=0; j<vec_chr_keyname.size(); j++)
    {
        out << "\t" << vec_chr_keyname[j];
    }

    out << endl;
    vector<double> temp;
    for(map<double, map<string, uint64_t> >::iterator it2 = map_stat_depth.begin(); it2 != map_stat_depth.end(); it2++)
    {
        temp.push_back(it2->first);
    }
    sort(temp.begin(), temp.end());
    for(unsigned int i=0; i<temp.size(); ++i)
    {
        uint64_t sum = 0;
        for(unsigned int j=0; j<vec_chr_keyname.size(); j++)
        {
            sum += map_stat_depth[temp[i]][vec_chr_keyname[j]];
        }

        out << temp[i] << "\t" << sum << "\t";
        for(unsigned int j=0; j<vec_chr_keyname.size(); j++)
        {
            out << "\t" << uint64_t(map_stat_depth[temp[i]][vec_chr_keyname[j]]);
        }
        out << endl;
    }
    out.close();
}
