#ifndef _TABFILE_H
#define _TABFILE_H

#include <string>
#include <vector>
#include <map>

struct tabfile
{
    int n_skipped_line;
    int n_column;
    
    std::vector <std::string> skipped_line;
    std::vector <std::vector<std::string> > data;

    void load(char *fn, int skip=1);
    template<class T> void get_data(std::vector<T> &v, int column);
    template<class T> void get_map(std::map<std::string,T> &m, int key_column, int data_column);
};
#endif
