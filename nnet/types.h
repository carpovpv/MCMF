#ifndef __TYPE_NET
#define __TYPE_NET

typedef unsigned char NetType;
typedef struct {
    int target;
    double value;
} output;

#define BACKPROP 0
#define PNNET 	 1
#define OCCK     2
#define OCCDIABOLO 3

#define EMPTY 666
#endif
