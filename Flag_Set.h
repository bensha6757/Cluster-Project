#ifndef Flag_Set_H_
#define Flag_Set_H_

#include "IO_Mem_Errors.h"
#include "Types.h"
#include <stdlib.h>
#include <math.h>
#define BYTES_PER_CELL (sizeof(long_num))
#define BITS_PER_CELL (BYTES_PER_CELL*8)

/************************************************************************************************************
 * a data structure storing Unmoved vertices for the Maximization algorithm.                                *
 * allowing a quick get_next_set_flag in order to get the next unmoved vertex for the algorithm's purposes  *
 ************************************************************************************************************/

/* allocates a new flag set, with an amount of @param size set bits */
long_num* allocate_flag_set(num size);

/** Inset element i to set.
 *  @param set - the set to be probed.
 *  @param size - original size of set. Not to be changed after insertion/removal of element.
 *  @param i - the elemnt to be inserted to set.
 */
void set_flag(long_num* set, num size, num i);

/** Remove element i from set.
 *  @param set - the set to be probed.
 *  @param size - original size of set. Not to be changed after insertion/removal of element.
 *  @param i - the elemnt to be removed from set.
 */
void reset_flag(long_num* set, num size, num i);

/** Returns if element i is in set.
 *  @param set - the set to be probed.
 *  @param size - original size of set. Not to be changed after insertion/removal of element.
 *  @param i - the elemnt to be checked if is in set.
 */
num get_flag(long_num* set, num size, num i);

/** Returns an element in set, which is nearest to the last element removed from set.
 *  @param set - the set to be probed.
 *  @param size - original size of set. Not to be changed after insertion/removal of element.
 *  @param lastFlag - last element removed from set.
 *  @param first - a boolean indicating if any element was previously removed from set.
 */
int get_next_set_flag(long_num* set, num size, num lastFlag, boolean first);


#endif /* Flag_Set_H_ */