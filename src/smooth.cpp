#include "smooth.h"

#include <vector>
#include <cmath>
#include <algorithm>

double smoothFunction::polynomial_smooth_function(double x, double w)
{
    double r=x/w;
    r*=r;
    r=1-r;
    r*=r*r*36/35/w;
    return r;
};
	
void smoothFunction::set_smooth(double w)
{
    _gaussian_std=(double)w/2;
    _gaussian_table_range=3*_gaussian_std;
    _gaussian_distribution.resize(2*_gaussian_table_range+1);

    if(_gaussian_table_range==0)
    {
        _gaussian_distribution[_gaussian_table_range]=1;
        return;
    }

    double x,s=0;
    for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
    {
        x=(double)i;
        x=x*x/_gaussian_std/_gaussian_std/2;
        x=exp(-x);
        s+=x;
        _gaussian_distribution[i+_gaussian_table_range]=x;
    }

                                
    for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
        _gaussian_distribution[i+_gaussian_table_range]/=s;
};

smoothFunction::smoothFunction(double w)
{
    set_smooth((double)w);
};

template<class classType>
    void smoothFunction::apply(std::vector<classType> &source)
{
    if(_gaussian_table_range==0) return;
    int lim=source.size();
    std::vector<classType> dest(lim);
    for(int i=0;i<lim;i++)
    {
        if(source[i]!=0) for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
            if(p>=0&&p<lim) dest[p]+=_gaussian_distribution[j+_gaussian_table_range]*source[i];
    }
    if(source[0]!=0) for(int i=-_gaussian_table_range;i<0;i++)
    {
        for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
            if(p>=0&&p<lim) dest[p]+=_gaussian_distribution[j+_gaussian_table_range]*source[0];
    }
    if(source[lim-1]!=0) for(int i=lim;i<=lim+_gaussian_table_range;i++)
    {
        for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
            if(p>=0&&p<lim) dest[p]+=_gaussian_distribution[j+_gaussian_table_range]*source[lim-1];
    }
    copy(dest.begin(),dest.end(),source.begin());
};

template void smoothFunction::apply(std::vector<int>&);
template void smoothFunction::apply(std::vector<float>&);
template void smoothFunction::apply(std::vector<double>&);
