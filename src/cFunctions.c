#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "string.h"
#include "myProto.h"
#include "Result.h"
#include "mpiFunctions.h"

void sequenceRun(LettersArray *orderConservativeGroups,
		LettersArray *orderSemiConservativeGroups, char **seq2Array, char *seq1,
		Result *results, float coefficents[], int numOfSeq2) {
	char *seq2Mutant;
	float scoreTemp = 0;
	float maxScore;
	int foundMaxScore;
	int firstRun = TRUE;
	int counter[MAX_COEFFICENTS] = { 0 };
	int maxOffset;
	int seq2MutantLen;

	for (int seq2Line = 0; seq2Line < numOfSeq2; seq2Line++) { // run on each seq2
		resetVariables(&firstRun, &foundMaxScore, &seq2MutantLen, &maxOffset,
				&maxScore, coefficents[0], seq1, seq2Array[seq2Line]);
		if (maxOffset < 0) { // if seq1 is shorter than seq2 mutant, return result with -1
			setInfoToResult(&results[seq2Line], -1, -2, -2, -1);
		} else {
			for (int n = 0;
					(n < strlen(seq2Array[seq2Line])) & (foundMaxScore == FALSE);
					n++) { // runs on n from 0 to seq2 len-1
				for (int k = n + 1;
						(k < strlen(seq2Array[seq2Line]))
								& (foundMaxScore == FALSE); k++) { // runs on k from n+1 to seq2 len-1
					seq2Mutant = createMutant(seq2Array[seq2Line], n, k);
					for (int offset = 0;
							(offset <= maxOffset) & (foundMaxScore == FALSE);
							offset++) { // runs on offset
						for (int i = 0; i < strlen(seq2Mutant); i++) { //compare chars of seq1 and seq2 mutant
							comparison(orderConservativeGroups,
									orderSemiConservativeGroups, counter,
									seq1[i + offset], seq2Mutant[i]);
						}
						scoreTemp = calculateScore(coefficents, counter);
						compareResults(&results[seq2Line], &firstRun, offset, n,
								k, scoreTemp);

						if (results[seq2Line].score == maxScore) {
							foundMaxScore = TRUE;
						}
					}
					free(seq2Mutant);
				}
			}
		}
	}
	printResult(results, numOfSeq2);
}

void createGroups(LettersArray *orderConservativeGroups,
		LettersArray *orderSemiConservativeGroups) {

	char const *conservativeGroups[MAX_CONSERTIVE_GROUPS] = { "NDEQ", "NEQK",
			"STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF" };
	char const *semiConservativeGroups[MAX_SEMI_CONSERTIVE_GROUPS] = { "SAG",
			"ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK",
			"HFY", "FVLIM" };

	markLettersInArray(orderConservativeGroups,conservativeGroups,MAX_CONSERTIVE_GROUPS);
	markLettersInArray(orderSemiConservativeGroups,semiConservativeGroups,MAX_SEMI_CONSERTIVE_GROUPS);
}

void markLettersInArray(LettersArray *newGroup,char const *group[], int size){//mark a letter in the arrays
	int index;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < strlen(group[i]); j++) {
			index = group[i][j] - 'A';
			newGroup[i].letter[index] = TRUE;
		}
	}
}
void resetCounter(int *counter) { // resets coefficents counter
	for (int i = 0; i < MAX_COEFFICENTS; i++) {
		counter[i] = 0;
	}
}

float calculateScore(float *coefficents, int *counter) {
	float score = 0;
	score = coefficents[0] * counter[0] - coefficents[1] * counter[1]
			- coefficents[2] * counter[2] - coefficents[3] * counter[3];
	resetCounter(counter);
	return score;
}

char* createMutant(char *seq2, int n, int k) {
	char *seq2Mutent = (char*) calloc(MAX_SEQ2_SIZE, sizeof(char));
	int j = 0;
	for (int i = 0; i < strlen(seq2) - 2; i++) {
		if ((j == n) | (j == k)) {
			j++;
		}
		seq2Mutent[i] = seq2[j];
		j++;
	}

	return seq2Mutent;
}
int comparison(LettersArray *conservativeGroups,
		LettersArray *semiConservativeGroups, int *counter, char charSeq1,
		char charSeq2) {
	if (charSeq1 == charSeq2) {	//Letters are equals
		counter[0]++;
		return 0;
	}

	if (compareGroups(conservativeGroups, MAX_CONSERTIVE_GROUPS, charSeq1,
			charSeq2) != ' ') {	//check if letters are in conservative Groups
		counter[1]++;
		return 0;
	}

	if (compareGroups(semiConservativeGroups,
	MAX_SEMI_CONSERTIVE_GROUPS, charSeq1, charSeq2) != ' ') {//check if letters are in semi Conservative Groups
		counter[2]++;
		return 0;
	}

	counter[3]++; //Letters in the pair are not equal, do not present both not in Conservative nor in Semi-Conservative groups
	return 0;
}

char compareGroups(LettersArray *groups, int groupSize, char charSeq1,
		char charSeq2) {
	int indexSeq1 = charSeq1 - 'A';
	int indexSeq2 = charSeq2 - 'A';
	for (int i = 0; i < groupSize; i++) {
		if (groups[i].letter[indexSeq1] == TRUE) {
			if (groups[i].letter[indexSeq2] == TRUE)
				return '$';
		}
	}
	return ' ';
}

char** readFromFile(float *coefficents, char *seq1, int *numOfSeq2) {
	char **seq2Array;

	scanf("%f  %f  %f  %f", &coefficents[0], &coefficents[1], &coefficents[2],
			&coefficents[3]);
	scanf("%s", seq1);
	scanf("%d", numOfSeq2);
	seq2Array = (char**) calloc(*numOfSeq2, sizeof(char*));
	for (int i = 0; i < *numOfSeq2; i++) {
		seq2Array[i] = (char*) calloc(MAX_SEQ2_SIZE, sizeof(char));
		scanf("%s", seq2Array[i]);
	}
	return seq2Array;
}

