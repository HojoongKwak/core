#include "bed.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

bool bedTrack::read(std::ifstream &in)
{
    std::string buf,line;
    while(in>>buf,std::getline(in,line),buf=="track"||buf=="browser"||buf[0]=='#');
    if(in.eof()) return false;
    chrom=buf;
    std::stringstream ss(line);
    ss>>chromStart>>chromEnd;
    ss>>name>>score>>strand;
    ss>>thickStart>>thickEnd>>itemRGB>>blockCount>>blockSizes>>blockStarts;
    return true;
};

void bedTrack::write(std::ofstream &out, int ncol)
{
    out<<chrom<<"\t"<<chromStart<<"\t"<<chromEnd;
    if(ncol>3) out<<"\t"<<name<<"\t"<<score<<"\t"<<strand;
    if(ncol>6) out<<"\t"<<thickStart<<"\t"<<thickEnd<<"\t"<<itemRGB;
    if(ncol>9) out<<"\t"<<blockCount<<"\t"<<blockSizes<<"\t"<<blockStarts;
    out<<std::endl;
};

double bedTrack::pair_to_score(int a, int b)
{
    int lg_a = log10((double)a+1.0)*200.0+0.5;
    int lg_b = log10((double)b+1.0)*200.0+0.5;
    int lg_s = lg_a + lg_b;
    score = (double)(lg_s*(lg_s+1)/2+lg_b)/2000.0;
    return score;
};

void bedTrack::score_to_pair(int &a, int &b)
{
    double w = floor((sqrt(score*16000.0+1.0)-1.0)/2.0);
    double t = w*(w+1.0)/2.0;
    double lg_b = score*2000.0-t;
    double lg_a = w - lg_b;
    a = pow(10,lg_a/200.0)-0.5;
    b = pow(10,lg_b/200.0)-0.5;
};

