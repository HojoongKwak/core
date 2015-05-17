#include "fasta.h"

#include <string>
#include <iostream>
#include <fstream>

bool fastaTrack::read(std::ifstream &in)
{
	std::string line;
	std::getline(in,line);
	if(line[0]!='>') return false;
	name=line.substr(1);
	while(!in.eof())
	{
		std::streampos pos=in.tellg();
		std::getline(in,line);
		if(line[0]=='>')
		{
			in.seekg(pos);
			break;
		}
		seq+=line;
	}
	return true;
};

void fastaTrack::write(std::ofstream &out)
{
	out<<">"<<name<<std::endl;
	int curPos=0;
	while(curPos<seq.size())
	{
		out<<seq.substr(curPos,80)<<std::endl;
		curPos+=80;
	}
};

void fastaFile::load(char *fileName)
{
	std::ifstream in(fileName);
	fastaTrack t;
	while(t.read(in)) track.push_back(t);
};
