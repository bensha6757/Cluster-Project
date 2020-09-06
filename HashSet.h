#ifndef HASHSET_H_
#define HASHSET_H_

#include "io_mem_errors.h"
#include "Types.h"
#include <stdlib.h>
#include <math.h>
#define SIZE (sizeof(long_num) * 8)

long_num* allocate_hash_set(num size);
long_num* allocate_hash_set_unmoved(num size);

void setFlag(long_num* set, num size, num i);

void resetFlag(long_num* set, num size, num i);

long_num getFlag(long_num* set, num size, num i);

int getNextSetFlag(long_num* set, num size, num lastFlag);


#endif /* HASHSET_H_ */