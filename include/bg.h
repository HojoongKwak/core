#ifndef _BEDGRAPH_H
#define _BEDGRAPH_H

#include <string>
#include <vector>

struct bgdata
{
    // Load data from file
    int load(char *fn);
    // Get bedgraph level from a range
    void getlevel(std::vector<float> &o, std::string &c, int s, int e, int r=1, bool isaverage=false);
    // Get bedgraph level of a chromosome
    void getChr(std::vector<float> &o, std::string &c, int r=1);

    std::vector<std::string> chr;
	std::vector<std::vector<int> > start;
    std::vector<std::vector<int> > end;
	std::vector<std::vector<float> > val;
    std::vector<float> v;

    int ch;
    int res;
    int len;
	bool set(std::string &c, int r=1);
    void getcoverage(std::vector<float> &o, std::string &c, int s, int e, int r=1);
    int getChrID(std::string &c);
    float &operator[](int i);
    float operator()(int c,int i);
    float operator()(std::string &c,int i);
};
#endif
