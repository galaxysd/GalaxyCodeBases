// Project     : M-Vicuna (Modularized-Vicuna)
// Name        : main.cpp
// Author      : Xiao Yang
// Created on  : May 13, 2013
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2013)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or
//				 guaranteed support whatsoever. The Broad Institute cannot
//				 be responsible for its use,	misuse, or functionality.
// Description :

#include <iostream>

#include "xutil.h"
#include "xny/seq_cmp.hpp"

int main (int argc, char** argv){

	std::string s0, s1, s2;
	std::cin >> s0 >> s1;// >> s2;

	double timing = get_time();
	double start_time = timing;
	
	bio::global_alignment galnF (2, -3, -5, -2, 1, true);
	//bio::global_alignment galnR (2, -3, -5, -2, 1, false);
	galnF.set_alignment_type(1);
	//galnR.set_alignment_type(1);

	galnF (s0, s1);
	std::cout << "Path1: " << galnF.path() << "\n";

	//galnR (s0, s2);
	//std::cout << "Path2: " << galnR.path() << "\n";

	std::cout << "s0: " << s0 << "\ns1: " << s1 << "\n";//s2: " << s2 << "\n";

	print_time("Whole program takes \t", start_time);
	std::cout << "DONE!\n";

	return (EXIT_SUCCESS);
}






