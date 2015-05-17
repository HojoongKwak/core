#ifndef _PSLTRACK_H
#define _PSLTRACK_H

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

struct psltrack
{
	int match,mismatch,repmatch,Ns,Qgc,Qgb,Tgc,Tgb;
	char strand;
	std::string Qname;
	int Qsize,Qstart,Qend;
	std::string Tname;
	int Tsize,Tstart,Tend;
	int blockCount;
	std::string blockSizestring,Qstartstring,Tstartstring;
	std::vector<int> blockSizes,Qstarts,Tstarts;
	void parse_csv_string(std::string s,int c,std::vector<int> &v)
	{
		std::string str;
		int spacer;
		v.resize(0);
		for(int i=0;i<c;++i)
		{
			spacer=s.find_first_of(',');
			str=s.substr(0,spacer);
			v.push_back(atoi(str.c_str()));
			s=s.substr(spacer+1);
		}
	};
	int Tpos_from_Qpos(int Qpos, int &Bpos5, int &Bpos3)
	{
		int blockpos=0;
		if(Qpos<Qstart||Qpos>Qend) return -1;
		if(strand=='+')
		{
			while(blockpos<blockCount)
			{
				if(Qpos>=Qstarts[blockpos]) ++blockpos;
				else break;
			}
			--blockpos;
			if(Qpos-Qstarts[blockpos]>=blockSizes[blockpos])
			{
				Bpos5=-1;
				Bpos3=-1;
				return -1;
			}
			else
			{
				int Tpos=Tstarts[blockpos]+Qpos-Qstarts[blockpos]+1;
				Bpos5=Tpos-Tstarts[blockpos]+1;
				Bpos3=blockSizes[blockpos]+1-Bpos5;
				return Tpos;
			}
		}
		else
		{
			Qpos=Qsize-Qpos;
			while(blockpos<blockCount)
			{
				if(Qpos>Qstarts[blockpos]) ++blockpos;
				else break;
			}
			--blockpos;
			if(Qpos-Qstarts[blockpos]>blockSizes[blockpos])
			{
				Bpos5=-1;
				Bpos3=-1;
				return -1;
			}
			else
			{
				int Tpos=Tstarts[blockpos]+Qpos-Qstarts[blockpos];
				Bpos3=Tpos-Tstarts[blockpos]+1;
				Bpos5=blockSizes[blockpos]+1-Bpos3;
				return Tpos;
			}
		}
	};
			
	bool read(std::ifstream &in)
	{
		in>>match;
		if(in.eof()) return false;
		in>>mismatch>>repmatch>>Ns>>Qgc>>Qgb>>Tgc>>Tgb;
		in>>strand;
		in>>Qname>>Qsize>>Qstart>>Qend;
		in>>Tname>>Tsize>>Tstart>>Tend;
		in>>blockCount>>blockSizestring>>Qstartstring>>Tstartstring;
		parse_csv_string(Qstartstring,blockCount,Qstarts);
		parse_csv_string(Tstartstring,blockCount,Tstarts);
		parse_csv_string(blockSizestring,blockCount,blockSizes);
		std::string buf;
		std::getline(in,buf);
		return true;
	};
	
	void write(std::ofstream &out)
	{
		out<<match<<"\t"<<mismatch<<"\t"<<repmatch<<"\t";
		out<<Ns<<"\t"<<Qgc<<"\t"<<Qgb<<"\t"<<Tgc<<"\t"<<Tgb<<"\t";
		out<<strand<<"\t";
		out<<Qname<<"\t"<<Qsize<<"\t"<<Qstart<<"\t"<<Qend<<"\t";
		out<<Tname<<"\t"<<Tsize<<"\t"<<Tstart<<"\t"<<Tend<<"\t";
		out<<blockCount<<"\t"<<blockSizestring<<"\t"<<Qstartstring<<"\t"<<Tstartstring<<std::endl;
	};
};
#endif
