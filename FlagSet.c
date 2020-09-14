#include "FlagSet.h"

long_num* allocate_flag_set(num size){
    long_num *set, *p;
    num cells = size / BITS_IN_CELL;
    num r = size % BITS_IN_CELL > 0 ? 1 : 0; 
    set = (long_num*)malloc((cells+r) * BYTES_IN_CELL);
    for (p = set ; p < set + cells+r ; p++){
        *p = ~(*p & 0);
    }
    *(--p) >>= (BITS_IN_CELL - (size % BITS_IN_CELL));
    return set;
}

void set_flag(long_num* set, num size, num i){
    num locate = i / BITS_IN_CELL,  r = i % BITS_IN_CELL;
    long_num mask = (long_num)1 << r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    set[locate] |= mask;
}

void reset_flag(long_num* set, num size, num i){
    num locate = i / BITS_IN_CELL,  r = i % BITS_IN_CELL;
    long_num mask = ~((long_num)1 << r);
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    set[locate] &= mask;
}

num get_flag(long_num* set, num size, num i){
    num locate = i/BITS_IN_CELL,  r = i % BITS_IN_CELL;
    long_num cell = set[locate];
    long_num mask =  1 << r;
    long_num masked_cell = cell & mask;
    long_num i_bit = masked_cell >> r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    return (num) (i_bit % 2);
}

int get_next_set_flag(long_num* set, num size, num lastFlag, boolean first){
    num lastLoc = lastFlag / BITS_IN_CELL,  r = lastFlag % BITS_IN_CELL;
    num cap = (size / BITS_IN_CELL) + (size % BITS_IN_CELL > 0 ? 1 : 0);
    long_num *p = set + lastLoc, mask, cell = *p;
    num res;
    if (r + 1 != BITS_IN_CELL){
        if (first == FALSE){
            mask = ~(cell & 0) << (r + 1);
            cell &= mask;
        }
        res = (num)(log(cell & -cell)/log(2)) + (lastLoc * BITS_IN_CELL);
        if (res < size && cell != 0)
            return res;
        if (res >= size)
            return -1;
    }
    
    p++;
    while (p-set < cap && *p == 0){
        p++;
    }
    if (p-set >= cap)
        return -1;
    cell = *p;
    res = (num)(log(cell & -cell)/log(2)) + ((p-set) * BITS_IN_CELL);
    if (res < size)
        return res;
    return -1;
}