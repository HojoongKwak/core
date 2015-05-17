#include <cmath>
#include <algorithm>
#define _N_GAUSSIAN_TABLE 10

using namespace std;

namespace smooth
{
    double *_gaussian_distribution[_N_GAUSSIAN_TABLE];
    double _gaussian_std;
    int _gaussian_table_range;

    double polynomial_smooth_function(double x, double w)
	{
		double r;
		r=x/w;
		r=r*r;
		r=1-r;
		r=r*r*r;
		r*=36/35/w;
		return r;
	};
	
    template <class classType>
    void set_smooth(classType w, int table_id=0)
    {
        _gaussian_std=(double)w/2;
        _gaussian_table_range=3*_gaussian_std;
        _gaussian_distribution[table_id]=new double[2*_gaussian_table_range+1];

        if(_gaussian_table_range==0)
        {
            _gaussian_distribution[table_id][0]=1;
            return;
        }

        double x,s=0;
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
        {
            x=(double)i;
            x=x*x/_gaussian_std/_gaussian_std/2;
            x=exp(-x);
            s+=x;
            _gaussian_distribution[table_id][i+_gaussian_table_range]=x;
        }

                                    
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++) _gaussian_distribution[table_id][i+_gaussian_table_range]/=s;
    };

    template<class classType>
    void add_smooth(classType *array, classType value, int pos, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        for(int i=-_gaussian_table_range, p=pos-_gaussian_table_range;i<=_gaussian_table_range;i++,p++)
            if(p>=0&&p<lim) array[p]+=dstr[i+_gaussian_table_range]*value;
    };

    template<class classType1, class classType2>
    void smooth_to(classType1 *dest, classType2 *source, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        for(int i=0;i<lim;i++) dest[i]=0;
        for(int i=0;i<lim;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[i];
        }
        for(int i=-_gaussian_table_range;i<0;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[0];
        }
        for(int i=lim;i<=lim+_gaussian_table_range;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[lim-1];
        }
    };

    template<class classType>
    void smooth_it(classType *source, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        classType *dest;
        dest=new classType[lim];
        fill_n(dest,lim,0);
        for(int i=0;i<lim;i++)
        {
            if(source[i]!=0) for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[i];
        }
        if(source[0]!=0) for(int i=-_gaussian_table_range;i<0;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[0];
        }
        if(source[lim-1]!=0) for(int i=lim;i<=lim+_gaussian_table_range;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[lim-1];
        }
        copy(dest,dest+lim,source);
        delete[] dest;
    };
}

namespace stat
{
    template<class classType>
    void ntiles(classType *output, classType *array, int n, int lim)
    {
        classType *temp;
        temp = new classType[lim];
        for(int i=0;i<lim;i++) temp[i]=array[i];
        sort(temp,temp+lim);
        output[0]=temp[0];
        for(int i=1;i<=n;i++) output[i]=temp[i*lim/n-1];
        delete[] temp;
    };
    
    template<class classType>
    double mean(classType *array, int lim)
    {
        double sum=0;
        for(int i=0;i<lim;i++) sum+=array[i];
        return sum/lim;
    };

    template<class classType>
    double stdev(classType *array, int lim)
    {
        double average,sum_square=0;
        average=mean(array,lim);
        for(int i=0;i<lim;i++) sum_square+=(array[i]-average)*(array[i]-average);
        return sqrt(sum_square/lim);
    };

    struct histo
    {
        double *f;
        double *x_start;
        double *x_end;
        int lim;

        template<class classType>
        void set(classType x_begin, classType x_stop, classType width)
        {
            double range=(double)(x_stop-x_begin);
            lim = range/width+1;
            f=new double[lim];
            x_start=new double[lim];
            x_end=new double[lim];
            for(int i=0;i<lim;i++)
            {
                x_start[i]=(double)x_begin+width*i;
                x_end[i]=x_start[i]+width;
            }
        };
        
        template<class classType>
        void add(classType *array, int n)
        {
            for(int i=0;i<lim;i++) f[i]=0;
            int pos;
            for(int i=0;i<n;i++)
            {
                pos=0;
                while((double)array[i]>=x_end[pos]) pos++;
                f[pos]+=(double)1/n;
            }
        };
        
        void write(std::ofstream &out)
        {
            for(int i=0;i<lim;i++)
            {
                out<<x_start[i]<<"\t"<<f[i]<<endl;
                out<<x_end[i]<<"\t"<<f[i]<<endl;
            }
        };
    };



    using namespace smooth;

    struct density_map
    {
        double *f;
        double *x;
        double res;
        int lim;

        template<class classType>
        void set(classType x_start, classType x_end, classType width, classType resolution)
        {
            res=(double)resolution;
            double range=(double)(x_end-x_start);
            lim=range/res+1;
            f=new double[lim];
            x=new double[lim];
            for(int i=0;i<lim;i++) x[i]=(double)x_start+res*i;
            set_smooth(width/res);
        };

        template<class classType>
        void add(classType *array, int n)
        {
            for(int i=0;i<lim;i++) f[i]=0;
            for(int i=0;i<n;i++)
            {
                add_smooth(f,(classType)1/n/res,((double)array[i]-x[0])/res,lim);
            }
        };
        
        void write(std::ofstream &out)
        {
            for(int i=0;i<lim;i++)
            {
                out<<x[i]<<"\t"<<f[i]<<endl;
            }
        };
    };
}


