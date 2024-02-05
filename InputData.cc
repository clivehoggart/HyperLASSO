/**
 *   ADMIXMAP
 *   InputData.cc 
 *   Class to read and check all input data files
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
#include "InputData.h"

using namespace std;

/**
 *  Auxilary functions to read data from file.
 */
bool isWhiteLine(const char *p)
{
    while (*p) {
        if (!isspace(*p++)) {
            return false;
        }
    }

    return true;
}

void readFile(const char *fname, Matrix_s& data)
{
    if (0 == fname || 0 == strlen(fname)) return;

    ifstream in(fname);
    if (!in.is_open()) {
        cout << "Cannot open file for reading: " << fname << endl;
        exit(0);
    }
    else {
       cout << "Loading " << fname << endl;
    }

    data.clear();
    try {
        StringSplitter splitter;

        string line;        

        while (getline(in, line)) {
            if (!isWhiteLine(line.c_str())) {
	      data.push_back(splitter.split(line.c_str()));
            }
        }
    } catch (...) {
        in.close();
        throw;
    }
}

void LoadDoubles( const char *fname,
		  vector< vector<double> >& data, vector<string>& labels, vector< vector<bool> >& missing )
{
  if (0 == fname || 0 == strlen(fname)) return;
  
  ifstream in(fname);
  if (!in.is_open()) {
    cout << "Cannot open file for reading: " << fname << endl;
    exit(0);
  }
  else {
    cout << "Loading " << fname << endl;
  }
  
  try {
    StringSplitter splitter;
    string line;
    int j=0;
    
    while (getline(in, line)) {
      if (!isWhiteLine(line.c_str())) {
	Vector_s vline(splitter.split(line.c_str()));
	if( j == 0 ){
	  labels = vline;
	  data.resize( labels.size() );
	  missing.resize( labels.size() );
	}else{
	  vector<float> vec;
	  for( unsigned int i=0; i<vline.size(); i++ ){
             if( vline[i] != "NA" ){
                data[i].push_back(StringConvertor::toFloat(vline[i]));
                missing[i].push_back(false);
             }else{
                data[i].push_back(0);
                missing[i].push_back(true);
             }
	  }
	}
	j++;
      }
    }
  } catch (...) {
    in.close();
    throw;
  }
}

void LoadInts( const char *fname,
	       vector< vector<short> >& data, vector<string>& labels, vector< vector<bool> >& missing )
{
  if (0 == fname || 0 == strlen(fname)) return;
  
  ifstream in(fname);
  if (!in.is_open()) {
    cout << "Cannot open file for reading: " << fname << endl;
    exit(0);
  }
  else {
    cout << "Loading " << fname << endl;
  }
  
  try {
    StringSplitter splitter;
    string line;
    int j=0;
    
    while (getline(in, line)) {
      if (!isWhiteLine(line.c_str())) {
	Vector_s vline(splitter.split(line.c_str()));
	if( j == 0 ){
	  labels = vline;
	  data.resize( labels.size() );
	  missing.resize( labels.size() );
	}else{
	  vector<float> vec;
	  for( unsigned int i=0; i<vline.size(); i++ ){
             if( vline[i] != "NA" ){
                data[i].push_back(StringConvertor::toInt(vline[i]));
                missing[i].push_back(false);
             }else{
                data[i].push_back(0);
                missing[i].push_back(true);
             }
	  }
	}
	j++;
      }
    }
  } catch (...) {
    in.close();
    throw;
  }
}
