#ifndef Flag_Set_H_
#define Flag_Set_H_

#include "IO_Mem_Errors.h"
#include "Types.h"
#include <stdlib.h>
#include <math.h>
#define BYTES_IN_CELL (sizeof(long_num))
#define BITS_IN_CELL (BYTES_IN_CELL*8)

/************************************************************************************************************
 * a data structure storing Unmoved vertices for the Maximization algorithm                                 *
 * allowing a quick get_next_set_flag in order to get the next unmoved vertex for the algorithm's purposes  *
 ************************************************************************************************************/

/* allocating a new flag set, with an amount of @param size set bits */
long_num* allocate_flag_set(num size);

/* setting the i'th bit */
void set_flag(long_num* set, num size, num i);

/* resetting the i'th bit */
void reset_flag(long_num* set, num size, num i);

/* check if the i'th bit is set */
num get_flag(long_num* set, num size, num i);

/* based on the last flag location, return the next set bit */
int get_next_set_flag(long_num* set, num size, num lastFlag, boolean first);


#endif /* Flag_Set_H_ */