#pragma once
#define MAX_NUMBERS 256
#define MAX_SEQ1_SIZE 5000
#define MAX_SEQ2_SIZE 3000
#define MAX_CONSERTIVE_GROUPS 9
#define MAX_SEMI_CONSERTIVE_GROUPS 11
#define MAX_COEFFICENTS 4
#define MAX_LETTERS 26
#define BUFFER 6
#include "Result.h"

typedef struct {
	int letter[MAX_LETTERS];
} LettersArray;

void sequenceRun(LettersArray *orderConservativeGroups,
		LettersArray *orderSemiConservativeGroups, char **seq2Array, char *seq1,
		Result *results, float coefficents[], int numOfSeq2);
int parallelRun(LettersArray *ConservativeGroups,
		LettersArray *SemiConservativeGroups, char *seq2, char *seq1,
		Result *result, float coefficents[]);
char** readFromFile(float *coefficents, char *seq1, int *numOfSeq2);
int comparison(LettersArray *conservativeGroups,
		LettersArray *semiConservativeGroups, int *counter, char charSeq1,
		char charSeq2);
char compareGroups(LettersArray *groups, int groupSize, char charSeq1,
		char charSeq2);
char* createMutant(char *seq2, int n, int k);
float calculateScore(float *coefficents, int *counter);
void createGroups(LettersArray *orderConservativeGroups,
		LettersArray *orderSemiConservativeGroups);
void markLettersInArray(LettersArray *newGroup,char const *group[], int size);
void resetCounter(int *counter);
void copyResultFromGpu(Result *h_Result, Result *d_Result);
char* allocateStringOnGpu(char *h_String);
LettersArray* allocateArraysOfLettersArray(LettersArray *array, int arraySize);
float* allocateCoefficentsOnGpu(float coefficents[]);
float* allocateFloatArrayZeroOnGpu(int arraySize);
Result* allocateResultOnGpu();
float** allocateArraysOfFloat(int numOfArrays, int arraySize);
void freeCharAllocatedMemory(char *d_char);
void freeFloatAllocatedMemory(float *d_float);
void freeAllThreadsMemoryFromGpu(float **d_BlocksScores,float *d_LinesTotalScores,Result *d_result);
void freeAllMemoryFromGpu(char *d_Seq1, char *d_Seq2, float *d_Coefficents,
		LettersArray *d_Conservative, LettersArray *d_SemiConservative);
void freeResultAllocatedMemory(Result *d_result);
void freeFloatArrays(float **array);

