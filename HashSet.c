#include "HashSet.h"

num* allocate_hash_set(num size){
    num* set;
    num cap=size/sizeof(num);
    num r=size%sizeof(num) > 0 ? 1 : 0; 
    set=(num*)calloc(cap+r,sizeof(num));
    VERIFY(set!=NULL,MEM_ALLOC_ERROR)
    return set;

}

void setFlag(num* set, num size, num i){
    num locate = i/sizeof(num),  r=i%sizeof(num);
    num mask = 1 << r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    set[locate] = set[locate] | mask;
}

num getFlag(num* set, num size, num i){
    num locate = i/sizeof(num),  r=i%sizeof(num);
    num cell = set[locate];
    num mask =  1 << r;
    num masked_cell = cell & mask;
    num i_bit = masked_cell >> r;
    VERIFY(i<size, OUT_OF_BOUNDS_ERROR)
    return i_bit % 2;
}