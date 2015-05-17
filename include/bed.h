#ifndef _BED_H
#define _BED_H

#include <string>
#include <iostream>
#include <fstream>

struct bedTrack
{
    std::string chrom;
    int chromStart, chromEnd;
    std::string name;
    double score;
    char strand;
    int thickStart, thickEnd;
    std::string itemRGB;
    int blockCount;
    std::string blockSizes, blockStarts;

    bool read(std::ifstream &in); 
    void write(std::ofstream &out, int ncol=3);
    double pair_to_score(int a, int b);
    void score_to_pair(int &a, int &b);
};
#endif
