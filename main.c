/*
 * main.c
 *
 *  Created on: 14 באוג׳ 2020
 *      Author: גל
 */

#include "modMat.h"



void fileToMatrix(char* filename, modMat *mat){
	size_t N;
	fopen(filename,"r");
	readVnumFromFile(filename,&N);
	mat=*allocateModMat(N);
	loadModMatrixFromFile(filename, mat);
	fclose(filename);
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	modMat *mat;

	fileToMatrix(inputFileName,mat);
	initPartition(divO);
	initPartition(divP);
	TrivialPartition(divP);
	while (!divP->isEmpty()){
		g=pickSubset(divP);
		DivideIntoTwo(g1,g2,g);
		if (g1->length()==0 or g2->length()==0)
			addToPartition(divO,g);
		else{

		}
	}
	PartitionsToFile(outputFileName,divO);
	return SUCCESS;
}
