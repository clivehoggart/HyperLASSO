// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   InputData.h 
 *   header file for InputData class
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef INPUT_DATA_H
#define INPUT_DATA_H 1
#include <iostream>
#include <fstream>
#include <vector>
//#include <string>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <string.h>
#include "StringSplitter.h"
#include "StringConvertor.h"
#include "common.h"

bool isWhiteLine(const char *p);
void readFile(const char *fname, Matrix_s& data);
void LoadInts( const char *fname,
	       std::vector< std::vector<short> >& data,
	       std::vector<std::string>& labels,
	       std::vector< std::vector<bool> >& missing );
void LoadDoubles( const char *fname,
	         std::vector< std::vector<double> >& data,
	         std::vector<std::string>& labels,
		 std::vector< std::vector<bool> >& missing );

#endif /* !defined INPUT_DATA_H */
