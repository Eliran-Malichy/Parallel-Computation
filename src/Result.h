/*
 * Results.h
 *
 *  Created on: 8 Feb 2022
 *      Author: eliran
 */

#ifndef RESULT_H_
#define RESULT_H_

#define TRUE 1
#define FALSE 0

typedef struct {
	int ms[2];
	int offsetPlace;
	float score;
} Result;

void compareResults(Result *result, int *firstRun, int offset, int n, int k,
		float scoreTemp);
void setInfoToResult(Result *result, int offset, int n, int k, float scoreTemp);
void printResult(Result *results, int numOfSeq2);
void getLargestResult(Result *h_Result, Result *finalResult, int NUM_THREADS);
void copyResult(Result *result, Result tempResult);
void resetVariables(int *firstRun, int *foundMaxScore, int *seq2MutantLen,
		int *maxOffset, float *maxScore, float coefficent, char *seq1,
		char *seq2);

#endif /* RESULT_H_ */
