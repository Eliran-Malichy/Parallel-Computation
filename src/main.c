/*
 ============================================================================
 Name        : FinalProject.c
 Author      : Eliran Malichy
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include "myProto.h"
#include "Result.h"
#include "mpiFunctions.h"
#define MASTER 0
#define MAX_SLAVES 2
#define COMPUTER_LOCAL 2
#define COMPUTER_OTHER 1
#define TAG_FINISH 0
#define TAG_RUN 1
#define TAG_READY 2
#define TAG_SENDING 3

int main(int argc, char *argv[]) {
	LettersArray *orderConservativeGroups = (LettersArray*) calloc(MAX_CONSERTIVE_GROUPS,
			sizeof(LettersArray));
	LettersArray *orderSemiConservativeGroups = (LettersArray*) calloc(
	MAX_SEMI_CONSERTIVE_GROUPS, sizeof(LettersArray));
	char **seq2Array;
	char *seq1 = (char*) calloc(MAX_SEQ1_SIZE, sizeof(char));
	char *seq2;
	Result *results;
	Result *ResultToSend;
	float coefficents[MAX_COEFFICENTS];
	int numOfSeq2Lines = 0;
	int processJob[MAX_SLAVES] = { 0 };
	int maxProcSize = MAX_SLAVES + 1;
	int procRank; /* rank of process */
	int procSize; /* number of processes */
	int source; /* rank of sender */
	int tag;
	int run = TRUE, countSend = 0, countRecive = 0, runningProcesses =
	MAX_SLAVES;
	int emptyPointer = 0;
	MPI_Datatype ResultMPIType;
	MPI_Status status; /* return status for receive */

	/* start up MPI */

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);/* find out process rank */
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);/* find out number of processes */

	if ((procSize != 1) & (procSize != maxProcSize)) {
		printf(
				"Run with one process for sequence run or three processes for parallel run only\n");
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	createResultlMpiType(&ResultMPIType); //creates mpi type

	if (procRank == 0) { //process 0 - read from file and create array of results
		seq2Array = readFromFile(coefficents, seq1, &numOfSeq2Lines);
		results = (Result*) calloc(numOfSeq2Lines, sizeof(Result)); //create array results
	}
	if (procSize == 1) { // sequence run
		createGroups(orderConservativeGroups,
				orderSemiConservativeGroups);
		sequenceRun(orderConservativeGroups, orderSemiConservativeGroups,
				seq2Array, seq1, results, coefficents, numOfSeq2Lines);
	}
	if (procSize == maxProcSize) { // parallel run
		MPI_Bcast(seq1, MAX_SEQ1_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD); //broadcast seq1 from 0 to all processes
		MPI_Bcast(coefficents, MAX_COEFFICENTS, MPI_FLOAT, 0, MPI_COMM_WORLD); //broadcast coefficents from 0 to all processes
		if (procRank == 0) { //process 0 job
			while (run == TRUE) {
				MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status); // Probe for an incoming message
				source = status.MPI_SOURCE;
				tag = status.MPI_TAG;
				if (tag == TAG_READY) {
					MPI_Recv(&emptyPointer, 1, MPI_INT, source, TAG_READY,
					MPI_COMM_WORLD, &status); //recive empty message
					if (countSend == numOfSeq2Lines) { //checks if all lines of seq2 were sent
						MPI_Send(NULL, 0, MPI_INT, source, TAG_FINISH,
						MPI_COMM_WORLD);
						runningProcesses--;
					} else {
						MPI_Send(seq2Array[countSend], MAX_SEQ2_SIZE, MPI_CHAR,
								source,
								TAG_RUN, MPI_COMM_WORLD); //send next seq2 line
						processJob[source - 1] = countSend; //saves line num for each process
						countSend++;
					}
				}
				if (tag == TAG_SENDING) { //recive result
					MPI_Recv(&results[processJob[source - 1]], 1, ResultMPIType,
							source,
							TAG_SENDING, MPI_COMM_WORLD, &status);
					countRecive++;
				}
				if ((runningProcesses == TAG_FINISH)
						& (countRecive == numOfSeq2Lines)) {

					run = FALSE;
					printResult(results, numOfSeq2Lines);
				}
			}
		} else { //processes 1 and 2 jobs
			createGroups(orderConservativeGroups,
					orderSemiConservativeGroups); //sort string arrays alphabetically

			while (run == TRUE) {
				seq2 = (char*) calloc(MAX_SEQ1_SIZE, sizeof(char));
				ResultToSend = (Result*) calloc(1, sizeof(Result));
				MPI_Send(NULL, 0, MPI_INT, 0, TAG_READY, MPI_COMM_WORLD); //notify process 0 that this process is ready to recive job
				MPI_Probe(0, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status); // Probe for an incoming message
				if (status.MPI_TAG == TAG_FINISH) {
					MPI_Recv(seq2, 0, MPI_CHAR, 0, TAG_FINISH, MPI_COMM_WORLD,
							&status);
					free(ResultToSend);
					free(seq2);
					break;
				}
				if (status.MPI_TAG == TAG_RUN) {
					MPI_Recv(seq2, MAX_SEQ2_SIZE, MPI_CHAR, 0, TAG_RUN,
					MPI_COMM_WORLD, &status); //recive seq 2 from process 0
					parallelRun(orderConservativeGroups,
							orderSemiConservativeGroups, seq2, seq1,
							ResultToSend, coefficents);
					MPI_Send(ResultToSend, 1, ResultMPIType, 0, TAG_SENDING,
					MPI_COMM_WORLD); //send result to process 0.
				}
				free(ResultToSend);
				free(seq2);
			}
		}
	}
	free(seq1);
	free(orderConservativeGroups);
	free(orderSemiConservativeGroups);
	fflush(stdout);
	/* shut down MPI */
	MPI_Finalize();

	return 0;
}
