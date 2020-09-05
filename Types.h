#ifndef TYPES_H_
#define TYPES_H_

typedef double* vector;
typedef unsigned int num;
typedef num* int_vector;
typedef int_vector Subgroup;
typedef int boolean;
typedef double DATA;
typedef double scalar;
#define TRUE 1
#define FALSE 0

#define SHIFT TRUE
#define NO_SHIFT FALSE

#undef DEBUG_SPMAT
#undef DEBUG_MODMAT
#undef DEBUG_EIGEN
#undef DEBUG_DIV_TWO
#undef DEBUG

#endif /* TYPES_H_ */