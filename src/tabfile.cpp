#include "tabfile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

void tabfile::load(char *fn, int skip)
{
    n_skipped_line=skip;
    std::ifstream in(fn);
    std::string buf;
    std::string tempstr;

    // read skipped line
    for(int i=0;i<n_skipped_line;i++)
    {
        std::getline(in,buf);
        skipped_line.push_back(buf);
    }

    // test read first raw and get column count
    std::streampos spos=in.tellg();
    std::getline(in,buf);
    int slen=buf.size();
    n_column=1;
    for(int i=0;i<slen-1;i++) if(buf[i]<=32&&buf[i+1]>32) n_column++;
    in.seekg(spos);

    // initialize data
    data.resize(n_column);
    
    // load data
    while(!in.eof())
    {
        for(int i=0;i<n_column;i++)
        {
            in>>tempstr;
            if(in.eof()) break;
            data[i].push_back(tempstr);
        }
    }
    in.close();
};

template <class T>
void tabfile::get_data(std::vector<T> &v, int column)
{
    int size=data[column].size();
    v.resize(size);
    for(int i=0;i<size;i++) std::stringstream(data[column][i])>>v[i];
};

template void tabfile::get_data(std::vector<int>&,int);
template void tabfile::get_data(std::vector<float>&,int);
template void tabfile::get_data(std::vector<double>&,int);
template void tabfile::get_data(std::vector<std::string>&,int);

template <class T>
void tabfile::get_map(std::map<std::string,T> &m, int key_column, int data_column)
{
    std::vector<T> d;
    get_data(d,data_column);
    int size=data[key_column].size();
    for(int i=0;i<size;i++) m[data[key_column][i]]=d[i];
};

template void tabfile::get_map(std::map<std::string,int>&,int,int);
template void tabfile::get_map(std::map<std::string,float>&,int,int);
template void tabfile::get_map(std::map<std::string,double>&,int,int);
template void tabfile::get_map(std::map<std::string,std::string>&,int,int);
