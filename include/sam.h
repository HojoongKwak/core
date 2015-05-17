#ifndef _SAM_H
#define _SAM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

struct SAMtrack
{
    std::string QNAME;
    int FLAG;
    std::string RNAME;
    int POS;
    int MAPQ;
    std::string CIGAR;
    std::string RNEXT;
    int PNEXT;
    int TLEN;
    std::string SEQ;
    std::string QUAL;
    std::string optional;
    bool read(std::ifstream &in)
    {
        in>>QNAME;
        if(in.eof()) return false;
        in>>FLAG>>RNAME>>POS>>MAPQ>>CIGAR>>RNEXT>>PNEXT>>TLEN>>SEQ>>QUAL;
        getline(in,optional);
        return true;
    };
    void write(std::ofstream &out)
    {
        out<<QNAME<<"\t"<<FLAG<<"\t"<<RNAME<<"\t"<<POS<<"\t"<<MAPQ<<"\t";
        out<<CIGAR<<"\t"<<RNEXT<<"\t"<<PNEXT<<"\t"<<TLEN<<"\t"<<SEQ<<"\t"<<QUAL<<std::endl;
    };

	std::vector<char> cigarOP;
	std::vector<int> cigarBN;
	int cigarN;

	SAMtrack()
	{
		cigarOP.resize(100);
		cigarBN.resize(100);
	};

	void parse_CIGAR()
	{
		std::string delim="MIDNSHP=X";
		std::size_t current=0;
		std::size_t next;
		cigarN=0;
		while(next=CIGAR.find_first_of(delim,current),next!=std::string::npos)
		{
			cigarOP[cigarN]=CIGAR[next];
			std::stringstream(CIGAR.substr(current,next-current))>>cigarBN[cigarN];
			current=next+1;
			++cigarN;
		}
	};

	int get_mapped_length()
	{
		int maplen=0;
		for(int i=0;i<cigarN;++i)
		{
			switch(cigarOP[i])
			{
				case 'M':
				case 'N':
				case 'D':
					maplen+=cigarBN[i];
					break;
			}
		}
		return maplen;
	};
};

struct SAMfile
{
    std::vector<std::string> SN;
    std::vector<int> LN;

    void read_header(std::ifstream &in)
    {
        std::streampos pos;
        std::string buf;
        do
        {
            pos=in.tellg();
            getline(in,buf);
            if(buf.substr(0,3)=="@SQ")
            {
                std::stringstream nameStream(buf.substr(buf.find("SN:")+3));
                std::string name;
                nameStream>>name;
                std::stringstream lenStream(buf.substr(buf.find("LN:")+3));
                int len;
                lenStream>>len;
                SN.push_back(name);
                LN.push_back(len);
            }
        } while(buf[0]=='@');
        in.seekg(pos);
    };
		
	void convert_to_bedgraph(std::ifstream &in, std::ofstream &outp, std::ofstream &outm, int cov, int minlen=16)
	{
		read_header(in);
		int nchr=LN.size();

		// allocate memory for bedgraph
		std::cout<<"Allocating memory"<<std::endl;
		std::vector< std::vector<int> > bg_p(LN.size()), bg_m(LN.size());
		long long usedmem=0LL;
		for(int i=0;i<nchr;++i)
		{
			try
			{
				bg_p[i].resize(LN[i]);
				bg_m[i].resize(LN[i]);
			}
			catch(std::bad_alloc const&)
			{
				std::cout<<"Memory allocation failure."<<std::endl;
				exit(0);
			}
			usedmem+=(long long)LN[i]*2*sizeof(int);
			// std::cout<<"... "<<usedmem/1048576<<"Mb allocated"<<std::endl;
		}

		// Read all the tracks
		std::cout<<"Reading tracks"<<std::endl;
		SAMtrack t;
		int readCount=0;
		while(t.read(in))
		{
			++readCount;
			if(readCount%500000==0) std::cout<<"... "<<(float)(readCount/100000)/10<<"M reads\r"<<std::flush;
			if(t.FLAG&0x04) continue;
			int chr = std::find(SN.begin(),SN.end(),t.RNAME)-SN.begin();
			t.parse_CIGAR();
			int maplen=t.get_mapped_length();
			if(maplen<minlen) continue;
			if(cov==5)		// coverage of 5' end
			{
				if(t.FLAG&0x10)	++bg_m[chr][t.POS+maplen-1];
				else ++bg_p[chr][t.POS];
			}
			else if(cov==3)	// coverage of 3' end
			{
				if(t.FLAG&0x10) ++bg_m[chr][t.POS];
				else ++bg_p[chr][t.POS+maplen-1];
			}
			else			// base coverage
			{
				int curpos=t.POS;
				for(int i=0;i<t.cigarN;++i)
				{
					switch(t.cigarOP[i])
					{
						case 'M':
							if(t.FLAG&0x10) for(int j=0;j<t.cigarBN[i];++j) ++bg_m[chr][curpos++];
							else for(int j=0;j<t.cigarBN[i];++j) ++bg_p[chr][curpos++];
							break;
						case 'D':
						case 'N':
							curpos+=t.cigarBN[i];
							break;
					}
				}
			}
		}
		std::cout<<readCount<<"reads loaded."<<std::endl;

		// write bedgraph file
		std::cout<<"Writing bedgraph"<<std::endl;
		for(int i=0;i<nchr;++i)
		{
			if(SN[i].find("chr")==std::string::npos) SN[i].insert(0,"chr");
			std::cout<<"... "<<SN[i]<<"\r"<<std::flush;
			int prev_pp=0,prev_mp=0;
			int prev_pv=bg_p[i][0],prev_mv=bg_m[i][0];

			for(int j=0;j<LN[i];++j)
			{
				int cur_pv=bg_p[i][j],cur_mv=bg_m[i][j];
				if(cur_pv!=prev_pv)
				{
					if(prev_pv>0) outp<<SN[i]<<"\t"<<prev_pp<<"\t"<<j<<"\t"<<prev_pv<<std::endl;
					prev_pp=j;
					prev_pv=cur_pv;
				}
				if(cur_mv!=prev_mv)
				{
					if(prev_mv>0) outm<<SN[i]<<"\t"<<prev_mp<<"\t"<<j<<"\t"<<-prev_mv<<std::endl;
					prev_mp=j;
					prev_mv=cur_mv;
				}
				if(j==LN[i]-1)
				{
					if(prev_pv>0) outp<<SN[i]<<"\t"<<prev_pp<<"\t"<<j<<"\t"<<prev_pv<<std::endl;
					if(prev_pv>0) outm<<SN[i]<<"\t"<<prev_mp<<"\t"<<j<<"\t"<<-prev_mv<<std::endl;
				}
			}
		}
		std::cout<<"... complete."<<std::endl;
	};

	void convert_to_bedgraph(std::ifstream &in, std::ofstream &out, int cov, int minlen=16)
	{
		std::cout<<"Reading header"<<std::endl;
		read_header(in);
		std::cout<<"... complete."<<std::endl;
		int nchr=LN.size();

		// allocate memory for bedgraph
		std::cout<<"Allocating memory"<<std::endl;
		std::vector< std::vector<int> > bg(nchr);
		long long usedmem=0LL;
		for(int i=0;i<nchr;++i)
		{
			try
			{
				bg[i].resize(LN[i]);
			}
			catch(std::bad_alloc const&)
			{
				std::cout<<"Memory allocation failure."<<std::endl;
				exit(0);
			}
			usedmem+=(long long)LN[i]*sizeof(int);
			std::cout<<"... "<<usedmem/1048576<<"Mb allocated"<<std::endl;
		}
		std::cout<<"... complete."<<std::endl;

		// Read all the tracks
		std::cout<<"Reading tracks"<<std::endl;
		SAMtrack t;
		int readCount=0;
		while(t.read(in))
		{
			++readCount;
			if(readCount%500000==0) std::cout<<"... "<<(float)(readCount/100000)/10<<"M reads\r"<<std::flush;
			if(t.FLAG&0x04) continue;
			int chr = std::find(SN.begin(),SN.end(),t.RNAME)-SN.begin();
			t.parse_CIGAR();
			int maplen=t.get_mapped_length();
			if(maplen<minlen) continue;
			if(cov==5)		// coverage of 5' end
			{
				if(t.FLAG&0x10)	++bg[chr][t.POS+maplen-1];
				else ++bg[chr][t.POS];
			}
			else if(cov==3)	// coverage of 3' end
			{
				if(t.FLAG&0x10) ++bg[chr][t.POS];
				else ++bg[chr][t.POS+maplen-1];
			}
			else			// base coverage
			{
				int curpos=t.POS;
				for(int i=0;i<t.cigarN;++i)
				{
					switch(t.cigarOP[i])
					{
						case 'M':
							for(int j=0;j<t.cigarBN[i];++j) ++bg[chr][curpos++];
							break;
						case 'D':
						case 'N':
							curpos+=t.cigarBN[i];
							break;
					}
				}
			}
		}
		std::cout<<readCount<<"reads loaded."<<std::endl;

		// write bedgraph file
		std::cout<<"Writing bedgraph"<<std::endl;
		for(int i=0;i<nchr;++i)
		{
			if(SN[i].find("chr")==std::string::npos) SN[i].insert(0,"chr");
			std::cout<<"... "<<SN[i]<<"\r"<<std::flush;
			int prev_p=0;
			int prev_v=bg[i][0];

			for(int j=0;j<LN[i];++j)
			{
				int cur_v=bg[i][j];
				if(cur_v!=prev_v)
				{
					if(prev_v>0) out<<SN[i]<<"\t"<<prev_p<<"\t"<<j<<"\t"<<prev_v<<std::endl;
					prev_p=j;
					prev_v=cur_v;
				}
				if(j==LN[i]-1)
				{
					if(prev_v>0) out<<SN[i]<<"\t"<<prev_p<<"\t"<<j<<"\t"<<prev_v<<std::endl;
				}
			}
		}
		std::cout<<"... complete."<<std::endl;
	};
};

#endif
