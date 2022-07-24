#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "string.h"
#include "myProto.h"
#include "Result.h"
#include "mpiFunctions.h"



void resetVariables(int *firstRun, int *foundMaxScore, int *seq2MutantLen,
		int *maxOffset, float *maxScore, float coefficent, char *seq1,
		char *seq2) {
	*firstRun = TRUE;
	*foundMaxScore = FALSE;
	*seq2MutantLen = strlen(seq2) - 2;
	*maxScore = coefficent * (*seq2MutantLen);
	*maxOffset = strlen(seq1) - (*seq2MutantLen);
}

void createResultlMpiType(MPI_Datatype *ResultMPIType) { //creates Result MPI data type
	MPI_Datatype type[3] = { MPI_INT, MPI_INT, MPI_FLOAT };
	MPI_Aint disp[3];
	int blocklen[3] = { 2, 1, 1 };
	Result *result = (Result*) calloc(1, sizeof(Result));

	disp[0] = (char*) &result->ms - (char*) result;
	disp[1] = (char*) &result->offsetPlace - (char*) result;
	disp[2] = (char*) &result->score - (char*) result;
	MPI_Type_create_struct(3, blocklen, disp, type, ResultMPIType);
	MPI_Type_commit(ResultMPIType);
}

void compareResults(Result *result, int *firstRun, int offset, int n, int k,
		float scoreTemp) {
	if (*firstRun == TRUE) {
		setInfoToResult(result, offset, n, k, scoreTemp);
		*firstRun = FALSE;
	}
	if ((scoreTemp > result->score)) {
		setInfoToResult(result, offset, n, k, scoreTemp);
	}
}

void setInfoToResult(Result *result, int offset, int n, int k, float scoreTemp) {
	result->score = scoreTemp;
	result->offsetPlace = offset;
	result->ms[0] = n + 1;
	result->ms[1] = k + 1;
}

void printResult(Result *results, int numOfSeq2) {
	for (int i = 0; i < numOfSeq2; i++) {
		printf("line %d: offset = %d, score = %1.2f, MS(%d,%d)\n", i,
				results[i].offsetPlace, results[i].score, results[i].ms[0],
				results[i].ms[1]);
	}
}

void getLargestResult(Result *h_Result, Result *finalResult, int NUM_THREADS) {
	int maxResultIndex = 0;

	for (int i = 1; i < NUM_THREADS; i++) {
		if (h_Result[maxResultIndex].score < h_Result[i].score) {
			maxResultIndex = i;
		}
	}
	copyResult(finalResult, h_Result[maxResultIndex]);
}

void copyResult(Result *result, Result tempResult) {

	result->score = tempResult.score;
	result->offsetPlace = tempResult.offsetPlace;
	result->ms[0] = tempResult.ms[0];
	result->ms[1] = tempResult.ms[1];
}



