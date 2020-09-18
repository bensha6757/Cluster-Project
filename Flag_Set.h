#ifndef Flag_Set_H_
#define Flag_Set_H_

#include "IO_Mem_Errors.h"
#include "Types.h"
#include <stdlib.h>
#include <math.h>
#define BYTES_IN_CELL (sizeof(long_num))
#define BITS_IN_CELL (BYTES_IN_CELL*8)

long_num* allocate_flag_set(num size);

void set_flag(long_num* set, num size, num i);

void reset_flag(long_num* set, num size, num i);

num get_flag(long_num* set, num size, num i);

int get_next_set_flag(long_num* set, num size, num lastFlag, boolean first);


#endif /* Flag_Set_H_ */