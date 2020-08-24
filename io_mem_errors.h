/*
 * io_mem_errors.h
 *
 *  Created on: 22 באוג׳ 2020
 *      Author: גל
 */

#ifndef IO_MEM_ERRORS_H_
#define IO_MEM_ERRORS_H_

typedef enum error_t {
	FILE_READ_ERROR,
	FILE_WRITE_ERROR,
	MEM_ALLOC_ERROR,
	INFINITE_LOOP_ERROR
} ERROR;


#endif /* IO_MEM_ERRORS_H_ */
