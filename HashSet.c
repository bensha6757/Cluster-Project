#include "HashSet.h"

long_num* allocate_hash_set(num size){
    long_num *set;
    num cap = size / SIZE;
    num r = size % SIZE > 0 ? 1 : 0; 
    set = (long_num*)calloc(cap+r,SIZE);
    VERIFY(set!=NULL,MEM_ALLOC_ERROR)
    return set;
}

long_num* allocate_hash_set_unmoved(num size){
    long_num *set, *p;
    num cap = size / SIZE;
    num r = size % SIZE > 0 ? 1 : 0; 
    set = (long_num*)calloc(cap+r,SIZE);
    VERIFY(set!=NULL,MEM_ALLOC_ERROR)
    for (p = set ; p < set + cap+r ; p++){
        *p = ~(*p);
    }
    return set;
}

void setFlag(long_num* set, num size, num i){
    num locate = i/SIZE,  r = i % SIZE;
    num mask = 1 << r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    set[locate] = set[locate] | mask;
}

void resetFlag(long_num* set, num size, num i){
    num locate = i/SIZE,  r = i % SIZE;
    num mask = ~(1 << r);
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    set[locate] &= mask;
}

long_num getFlag(long_num* set, num size, num i){
    num locate = i/SIZE,  r = i % SIZE;
    long_num cell = set[locate];
    long_num mask =  1 << r;
    long_num masked_cell = cell & mask;
    long_num i_bit = masked_cell >> r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    return i_bit % 2;
}

int getNextSetFlag(long_num* set, num size, num lastFlag){
    num lastLoc = lastFlag / SIZE,  r = lastFlag % SIZE, cap = (size / SIZE) + (size % SIZE > 0 ? 1 : 0);
    long_num *p = set + lastLoc;
    long_num mask = 0xffffffffffffffff << (r + 1);
    long_num cell = *p;
    num res;
    cell &= mask;
    res = (num)(log(cell & -cell)/log(2))  + (lastLoc * SIZE);
    #ifdef DEBUG_HASHSET
    printf("cap: %d \n", cap);
    printf("res: %d \n", res);
    #endif
    if (res < size && cell != 0)
        return res;
    if (res >= size)
        return -1;
    p++;
    while (p-set < cap && *p == 0){
        p++;
    }
    if (p-set >= cap)
        return -1;
    cell = *p;
    res = (num)(log(cell & -cell)/log(2)) + ((p-set) * SIZE);
    if (res < size)
        return res;
    return -1;
}