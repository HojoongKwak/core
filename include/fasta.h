#ifndef _FASTA_H
#define _FASTA_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

struct fastaTrack
{
    std::string name, seq;
    bool read(std::ifstream &in);
    void write(std::ofstream &out);
};

struct fastaFile
{
	std::vector<fastaTrack> track;
	void load(char *fileName);
};
#endif
