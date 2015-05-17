#ifndef _SAMTRACK_H
#define _SAMTRACK_H

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

struct samtrack
{
	std::string id;
	std::string chr;
	int pos;
	int mapq;
	char strand;
	std::string seq,phred;
	std::string MDtag;
	std::string CIGAR;
	std::vector<int> mapping;
	std::string mismap;
	void parse_CIGAR()
	{
		int cigarN=0;
		std::string delim="MIDNSHP=X";
		std::size_t current=0;
		std::size_t next;
		std::vector<int> cigarBN;
		std::string cigarOP;
		while(next=CIGAR.find_first_of(delim,current),next!=std::string::npos)
		{
			cigarOP.push_back(CIGAR[next]);
			cigarBN.push_back(0);
			std::stringstream(CIGAR.substr(current,next-current))>>cigarBN[cigarN];
			current=next+1;
			++cigarN;
		}
		int cpos=0,it=0;
		mapping.resize(seq.size());
		mismap.clear();
		mismap.resize(seq.size(),' ');
		for(int i=0;i<cigarN;++i)
		{
			switch(cigarOP[i])
			{
				case 'M':
					for(int j=0;j<cigarBN[i];++j,++cpos,++it)
					{
						mapping[it]=pos+cpos;
						mismap[it]='.';
					}
					break;
				case 'I':
					for(int j=0;j<cigarBN[i];++j,++it) mapping[it]=pos+cpos;
					break;
				case 'N':
				case 'D':
					cpos+=cigarBN[i];
					break;
				case 'S':
					for(int j=0;j<cigarBN[i];++j,++it) mapping[it]=-1;
					break;
			}
		}
	};

	void parse_MDtag()
	{
		std::string m=MDtag.substr(5);
		
		std::size_t cpos,ppos=0;
		std::size_t it=mismap.find('.');
		while(cpos=m.find_first_of("ACGT^",ppos),cpos!=std::string::npos)
		{
			if(cpos>ppos)
			{
				int match;
				std::stringstream(m.substr(ppos,cpos-ppos))>>match;
				for(int i=0;i<match;++i) it=mismap.find('.',it+1);
				ppos=cpos;
			}
			else
			{
				if(m[cpos]!='^')
				{
					mismap[it]=m[cpos];
					it=mismap.find('.',it+1);
					ppos++;
				}
				else
				{
					ppos++;
					while(ppos<m.size())
					{
						char c=m[ppos];
						if(c=='A'||c=='C'||c=='G'||c=='T') ppos++;
						else break;
					}
				}
			}
		}
	};

	bool read(std::ifstream &in)
	{
		std::string buf;
		in>>id;
		if(in.eof()) return false;
		int flag;
		in>>flag;
		if(flag&0x16) strand='-'; else strand='+';
		in>>chr;
		in>>pos;
		in>>mapq;
		in>>CIGAR;	// CIGAR
		in>>buf>>buf>>buf;
		in>>seq>>phred;
		getline(in,buf);
		std::stringstream ss(buf);
		ss>>MDtag;	// MDtag
		return true;
	};
};
#endif
