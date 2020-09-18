#ifndef TYPES_H_
#define TYPES_H_

#include <time.h>

typedef double* vector;
typedef unsigned int num;
typedef num* int_vector;
typedef int_vector Subgroup;
typedef int boolean;
typedef double DATA;
typedef double scalar;
typedef unsigned long long_num;
#define TRUE 1
#define FALSE 0

#define EPSILON 0.00001
#define IS_POSITIVE(X) ((X) > (EPSILON))
#define NON_ZERO(X) (IS_POSITIVE(X) || ((X) < (-(EPSILON))))

#define SHIFT TRUE
#define NO_SHIFT FALSE

#endif /* TYPES_H_ */