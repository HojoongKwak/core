#ifndef _ARGUMENTS_H
#define _ARGUMENTS_H

namespace arg
{
    template <class Type>
        void push(const char *id, const char *output, Type &vp, bool isoptional=false);
    bool get(int argc, char *argv[]);
    struct argInfo;
};
#endif
