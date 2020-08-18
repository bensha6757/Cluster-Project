/*
 * generateInput.c
 *
 *  Created on: 18 ����� 2020
 *      Author: ��
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>

int parseInt(char* str){
	int res=0;
	if (*str=='\0')
		exit(1);
	else
		res+=(int)((*str)-'0');
	str++;
	while (*str!='\0'){
		res*=10;
		res+=(int)((*str)-'0');
		str++;
	}
	return res;
}

int randInRangeExcluded(int t, int ex){
	int res=t+1;
	while (res>t || res==ex)
		res=rand();
	return res;
}

void writeRandomAdjList(int curr, int k, int n, FILE *file){
	int v=n+1, *adj,*p, *hash;
	adj=(int*)malloc(k*sizeof(int));
	if (adj==NULL)
		exit(1);
	for (p=adj; p<adj+k; p++){
		v=randInRangeExcluded(n,curr);
		*p=v;
	}
	assert(fwrite(adj,sizeof(int),k,file)==k);
	free(adj);
}

void generateGraphAdjencyFile(int num, FILE *file){
	int i, k;
	assert(fwrite(&num,sizeof(int),1,file)==1);	/* write n=|V| */
	for (i=0; i<num; i++){
		assert(fwrite(&i,sizeof(int),1,file)==1); /* write k_i - degree of v_i */
		k=randInRangeExcluded(num,num+1);
		writeRandomAdjList(i, k, num, file);
	}
}

int main(int argc, char* argv[]){
	int num;
	FILE *inputFile;
	SP_BUFF_SET();
	srand(time(NULL));

	assert(argc-1==2);
	num=parseInt(argv[1]);
	inputFile=fopen(argv[2],"w");
	assert(inputFile!=NULL);
	generateGraphAdjencyFile(num,inputFile);
	fclose(inputFile);
	return 0;
}
