#ifndef HASHSET_H_
#define HASHSET_H_

#include "io_mem_errors.h"
#include "Types.h"
#include <stdlib.h>

num* allocate_hash_set(num size);

void setFlag(num *set, num size, num i);

num getFlag(num* set, num size, num i);


#endif /* HASHSET_H_ */