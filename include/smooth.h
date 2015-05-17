#ifndef _SMOOTH_H
#define _SMOOTH_H

#include <vector>

struct smoothFunction
{
    smoothFunction(double w);
    template<class classType>
        void apply(std::vector<classType> &source);
    
    void set_smooth(double w);
    std::vector<double> _gaussian_distribution;
    double _gaussian_std;
    int _gaussian_table_range;
    double polynomial_smooth_function(double x, double w);
};
#endif
