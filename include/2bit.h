#ifndef _2BIT_HEADER
#define _2BIT_HEADER

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>

struct _2bit_header
{
    int signature;
    int version;
    int sequenceCount;
    int reserved;
    void read(std::ifstream &in)
    {
        in.read((char *)&signature,4);
        in.read((char *)&version,4);
        in.read((char *)&sequenceCount,4);
        in.read((char *)&reserved,4);
    };
    void write(std::ofstream &out)
    {
        out.write((char *)&signature,4);
        out.write((char *)&version,4);
        out.write((char *)&sequenceCount,4);
        out.write((char *)&reserved,4);
    };
};

struct _2bit_index
{
    int nameSize;
    char *name;
    int offset;
    void read(std::ifstream &in)
    {
        nameSize=0;
        in.read((char *)&nameSize,1);
        name=new char[nameSize+1];
        in.read(name,nameSize);
        name[nameSize]='\0';
        in.read((char *)&offset,4);
    };
    void write(std::ofstream &out)
    {
        out.write((char *)&nameSize,1);
        out.write(name,nameSize);
        out.write((char *)&offset,4);
    };
};

struct _2bit_sequence_record
{
   int dnaSize;
   int nBlockCount;
   int *nBlockStarts;
   int *nBlockSizes;
   int maskBlockCount;
   int *maskBlockStarts;
   int *maskBlockSizes;
   char *packedDna;
   void read(std::ifstream &in)
   {
        in.read((char *)&dnaSize,4);
        in.read((char *)&nBlockCount,4);
        nBlockStarts=new int[nBlockCount];
        for(int i=0;i<nBlockCount;i++) in.read((char *)&(nBlockStarts[i]),4);
        nBlockSizes=new int[nBlockCount];
        for(int i=0;i<nBlockCount;i++) in.read((char *)&(nBlockSizes[i]),4);

        in.read((char *)&maskBlockCount,4);
        maskBlockStarts=new int[maskBlockCount];
        for(int i=0;i<maskBlockCount;i++) in.read((char *)&(maskBlockStarts[i]),4);
        maskBlockSizes=new int[maskBlockCount];
        for(int i=0;i<maskBlockCount;i++) in.read((char *)&(maskBlockSizes[i]),4);
        int buf;
        in.read((char *)&buf,4);
        int maplen=(dnaSize/16+1)*4;
        packedDna=new char[maplen];
        in.read(packedDna,maplen);
   };
   
   void write(std::ofstream &out)
   {
        out.write((char *)&dnaSize,4);
        int maplen=(dnaSize/16+1)*4;
        packedDna=new char[maplen];
        out.write(packedDna,maplen);
   };

    inline int X(int i)
    {
        return (packedDna[i>>2]>>(6-(i&3)-(i&3)))&3;
    };
   
    void setX(int a, int i)
    {
        int offset=i>>2;
        int segment=i&3;
        packedDna[offset]&=~(3<<(6-segment-segment));
        packedDna[offset]|=a<<(6-segment-segment);
    };

    void applynBlock()
    {
        srand(time(NULL));
        for(int i=0;i<nBlockCount;i++)
        {
            for(int j=nBlockStarts[i];j<nBlockStarts[i]+nBlockSizes[i];j++)
            {
                setX(rand()%4,j);
            }
        }
    };
};

struct _2bit_file
{
    _2bit_header header;
    _2bit_index *indices;
    _2bit_sequence_record *records;
    int seqCount;
    char *base, *rbase;
	_2bit_file()
	{
		base=new char[4];
		base[0]='T';
		base[1]='C';
		base[2]='A';
		base[3]='G';
		rbase=new char[4];
		rbase[0]='A';
		rbase[1]='G';
		rbase[2]='T';
		rbase[3]='C';
	};
	_2bit_file(char *fn)
	{
		base=new char[4];
		base[0]='T';
		base[1]='C';
		base[2]='A';
		base[3]='G';
		rbase=new char[4];
		rbase[0]='A';
		rbase[1]='G';
		rbase[2]='T';
		rbase[3]='C';
		load(fn);
	};

	void init(char *fn)
    {
        std::ifstream in(fn,std::ios::binary);
        header.read(in);
        seqCount=header.sequenceCount;
        indices=new _2bit_index[seqCount];
        for(int i=0;i<seqCount;i++) indices[i].read(in);
        records=new _2bit_sequence_record[seqCount];
        in.close();
    };

    void readRecord(char *fn,int index)
    {
        std::ifstream in(fn,std::ios::binary);
        in.seekg(indices[index].offset,std::ios::beg);
        records[index].read(in);
        in.close();
    };
    void freeRecord(int index)
    {
        delete[] records[index].packedDna;
    };
    void save(char *fn)
    {
        std::ofstream out(fn,std::ios::binary);
        header.write(out);
        for(int i=0;i<seqCount;i++) indices[i].write(out);
        for(int i=0;i<seqCount;i++) records[i].write(out);
        out.close();
    };

	std::vector<std::string> chrName;
	
	void load(char *fn)
	{
		init(fn);
		for(int i=0;i<seqCount;++i)
		{
			readRecord(fn,i);
			records[i].applynBlock();
			chrName.push_back(indices[i].name);
		}
	};
	int getIndex(std::string chr)
	{
		for(int i=0;i<seqCount;i++)
			if(chr==chrName[i]) return i;
		return -1;
	};
	char get(std::string chr,int pos)
	{
		int i=getIndex(chr);
		if(i>=0) return base[records[i].X(pos-1)];
		else return 'N';
	};
	char rget(std::string chr,int pos)
	{
		int i=getIndex(chr);
		if(i>=0) return rbase[records[i].X(pos-1)];
		else return 'N';
	};
	void get(std::string chr, int from, int to, std::string &seq)
	{
		int i=getIndex(chr);
		seq.resize(to-from+1);
		if(i>=0)
		{
			int start=from-1;
			int end=to-1;
			if(start<0) start=0;
			if(end>=records[i].dnaSize) end=records[i].dnaSize-1;
			for(int j=start;j<=end;++j)
				seq[j-from+1]=base[records[i].X(j)];
		}
	};
	void rget(std::string chr, int from, int to, std::string &seq)
	{
		int i=getIndex(chr);
		int n=to-from+1;
		seq.resize(n);
		if(i>=0)
		{
			int start=from-1;
			int end=to-1;
			if(start<0) start=0;
			if(end>=records[i].dnaSize) end=records[i].dnaSize-1;
			for(int j=start;j<=end;++j)
				seq[n-j+from-2]=rbase[records[i].X(j)];
		}
	};
};
#endif
