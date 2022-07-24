#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"
#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include <omp.h>
#include "Result.h"
#define NUM_BLOCK_THREADS 64
#define FIND_MAX_BLOCK_THREADS 1024
#define NUM_THREADS 4


__device__ void resetCounters(int offset,float *LinesTotalScores,Result *resultsInBlock,int blockArrayLen,float **BlocksScores){
	int step=0;
	int numOfLoops = (offset/FIND_MAX_BLOCK_THREADS)+1;
	for(int i=0; i<numOfLoops;i++){
		if(threadIdx.x+step<offset){
			LinesTotalScores[threadIdx.x+step]=0;
		}
		step=step+FIND_MAX_BLOCK_THREADS;
	}
	step=0;
	for(int i=0; i<numOfLoops;i++){
		if(threadIdx.x+step<offset){
			for(int j=0;j<blockArrayLen;j++){
			BlocksScores[threadIdx.x+step][j]=0;
			}	
		}
		step=step+FIND_MAX_BLOCK_THREADS;
	}
	resultsInBlock[threadIdx.x].score=0;
	resultsInBlock[threadIdx.x].offsetPlace=0;
}


__device__ int getThreadIdDimX(){
	return  blockDim.x * blockIdx.x + threadIdx.x;
}


__device__ char compareGroupsOnGpu(LettersArray *groups, int groupSize, char charSeq1, char charSeq2) {
	int indexSeq1 = charSeq1-'A';
	int indexSeq2 = charSeq2-'A';
	for(int i=0;i<groupSize;i++){
		if(groups[i].letter[indexSeq1]==TRUE){
			if(groups[i].letter[indexSeq2]==TRUE)
				return '$';
		}
	}
	return ' ';
}


__device__ float comparisonOnGpu(LettersArray *conservativeGroups, LettersArray *semiConservativeGroups, char charSeq1, char charSeq2,float* coefficents) {
	if (charSeq1 == charSeq2) {	//Letters are equals
		return coefficents[0];
	}

	if (compareGroupsOnGpu(conservativeGroups, MAX_CONSERTIVE_GROUPS, charSeq1,charSeq2) != ' ') {	//check if letters are in conservative Groups
		return -coefficents[1];
	}

	if (compareGroupsOnGpu(semiConservativeGroups,
	MAX_SEMI_CONSERTIVE_GROUPS, charSeq1, charSeq2) != ' ') {//check if letters are in semi Conservative Groups
		return -coefficents[2];
	}

	return -coefficents[3];//Letters in the pair are not equal, do not present both not in Conservative nor in Semi-Conservative groups
}



__device__ void calculateBlockScore(float **BlocksScores,float *compareResultPerThread){//calculate each block score
		int threadId = getThreadIdDimX();
		int step=1;
		for(int i=0; (i<__log2f(NUM_BLOCK_THREADS));i++){
			if(threadIdx.x%(step*2)==0)
			{
				compareResultPerThread[threadIdx.x] = compareResultPerThread[threadIdx.x] + compareResultPerThread[threadIdx.x+step];
				step=step*2;		
			}
		__syncthreads();
		}
		if(threadIdx.x==0)
		{
			BlocksScores[blockIdx.y][blockIdx.x]=compareResultPerThread[threadIdx.x];
		}
}

__device__ char innerCheck(char *seq2,int k,int threadId){
		if(threadId+1==k){
			return seq2[threadId+2];
		}else{
			return seq2[threadId+1];
		}
}

__device__ char getCharFromSeq2(char *seq2, int n, int k,int threadId){//each thread get char from seq2, by that we "create" the seq2 mutant  

	if(threadId<n){
		return seq2[threadId];
	}
	if(threadId==n){
		return innerCheck(seq2, k, threadId);		
	}
	if(threadId<k){
		return innerCheck(seq2, k, threadId);
	}
 return seq2[threadId+2];
}



__device__ void calculateLinesScore(int numOfLoops,int threadId,int offset,int blockArrayLen,float *LinesTotalScores,float **BlocksScores){//calculate each offset total score 
		int step = 0;
		for(int i=0; i<numOfLoops;i++){
			if(threadId+step<offset){
				for(int i=0;i<blockArrayLen;i++){
					LinesTotalScores[threadId+step]=LinesTotalScores[threadId+step]+BlocksScores[threadId+step][i];
				}
			}
			step=step+FIND_MAX_BLOCK_THREADS;
		}

}

__device__ void findMaxScoreInlines(int numOfLoops,int threadId,int offset,float *LinesTotalScores,Result *resultsInBlock){//fineds the FIND_MAX_BLOCK_THREADS max scores 
	int step = 0;
	int l_offset=-1;
	float l_score=0;
	for(int i=0; i<numOfLoops;i++){
		if(threadId+step<offset){
			if(i==0){//first run
				l_offset=threadId;
				l_score=LinesTotalScores[threadId];
			}else{//second run and up
				if(LinesTotalScores[threadId+step]>l_score){
					l_offset=threadId+step;
					l_score=LinesTotalScores[threadId+step];
				}
			}
		}
		step=step+FIND_MAX_BLOCK_THREADS;
	}	
	resultsInBlock[threadId].score=l_score;
	resultsInBlock[threadId].offsetPlace=l_offset;
}

__device__ void findMaxScoreInBlock(int threadId,Result *resultsInBlock,Result *result){//findes the max score from scores in shared memory
	int jump=1;
	for(int i=0; i<__log2f(FIND_MAX_BLOCK_THREADS);i++){
		if(threadIdx.x%(jump*2)==0){
				if(resultsInBlock[threadId+jump].offsetPlace!=-1){
					if(resultsInBlock[threadId].score<resultsInBlock[threadId+jump].score){
						resultsInBlock[threadId].score=resultsInBlock[threadId+jump].score;
						resultsInBlock[threadId].offsetPlace=resultsInBlock[threadId+jump].offsetPlace;
					}
				}
			jump=jump*2;		
		}
			__syncthreads();
	}
	if(threadIdx.x==0){
		result->score=resultsInBlock[threadId].score;
		result->offsetPlace=resultsInBlock[threadId].offsetPlace;
	}
}
//each thread take a char from seq1 and seq2 meutant, compare them and save the result in shared memory.
//then the karnel calculate the score of each block by using the scores from the shared memory
//every line in the arrays represents an offset
__global__  void compare(int n, int k,char *seq1,char *seq2,int Seq2Len, float* coefficents,LettersArray *ConservativeGroups,LettersArray *SemiConservativeGroups,float **BlocksScores) {
	int threadId= getThreadIdDimX();
	char seq2Char;
	char seq1Char;
	__shared__ float compareResultPerThread[NUM_BLOCK_THREADS];

	
	//getGroupsToSharedMemory(orderConservativeGroups,orderSemiConservativeGroups,s_orderConservativeGroups,s_orderSemiConservativeGroups);//sends groups to shared memory for efficent run
	if(threadId<Seq2Len-2)
	{
		seq2Char = getCharFromSeq2( seq2, n, k, threadId);
		seq1Char = seq1[threadId+blockIdx.y];
		compareResultPerThread[threadIdx.x] = comparisonOnGpu(ConservativeGroups, SemiConservativeGroups, seq1Char, seq2Char,coefficents);
	}else{
		compareResultPerThread[threadIdx.x]=0;
	}
	__syncthreads();
	calculateBlockScore(BlocksScores,compareResultPerThread);	
}



__global__  void findMax(int offset,int blockArrayLen, Result *result, float *LinesTotalScores,float **BlocksScores){//findes the max score from array of offset scores
	int threadId= getThreadIdDimX();
	int numOfLoops = (offset/FIND_MAX_BLOCK_THREADS)+1;	
	__shared__ Result resultsInBlock[FIND_MAX_BLOCK_THREADS];
	result->score=0;
	result->offsetPlace=0;
	
	calculateLinesScore(numOfLoops,threadId,offset,blockArrayLen,LinesTotalScores,BlocksScores);//calculate the score of each line of blocks score
	__syncthreads();
	findMaxScoreInlines(numOfLoops,threadId,offset,LinesTotalScores,resultsInBlock);//findes FIND_MAX_BLOCK_THREADS max scores from the array of offsets scores
	__syncthreads();
	findMaxScoreInBlock(threadId,resultsInBlock,result);//findes the max score from the FIND_MAX_BLOCK_THREADS scores	
	__syncthreads();
	resetCounters(offset,LinesTotalScores,resultsInBlock,blockArrayLen,BlocksScores);//reset counters
	__syncthreads();
}
				

int parallelRun(LettersArray *ConservativeGroups,LettersArray *SemiConservativeGroups, char *seq2, char *seq1,Result *result, float coefficents[]) {
	
	volatile int run = TRUE;
	Result *resultTemp=(Result*) calloc(NUM_THREADS, sizeof(Result));
	Result *h_Result=(Result*) calloc(NUM_THREADS, sizeof(Result)); 
	float maxScore = coefficents[0] * (strlen(seq2) - 2);
	int firstRun = TRUE; //private
	int threadId; //private
  	int offset = (strlen(seq1)-(strlen(seq2) - 2))+1;
  	int sizeDimX=((strlen(seq2) - 2)/NUM_BLOCK_THREADS)+1;
  	dim3 gridShapeCompare(sizeDimX,offset);
  	dim3 gridShapeFindMax(1);
	cudaError_t err = cudaSuccess;
	cudaStream_t streams [NUM_THREADS];
	float *d_LinesTotalScores;
	float **d_BlocksScores;
	char* d_Seq1;
	char* d_Seq2;
	float *d_Coefficents;
	LettersArray *d_Conservative;
	LettersArray *d_SemiConservative;
	Result *d_result;
	
	
	if(offset<1){//checks if seq len is smaller than seq2 mutant len, if so returns result with -1 
	setInfoToResult(result,-1,-2,-2,-1);
	return 0;
	}
#pragma omp parallel num_threads(NUM_THREADS) firstprivate(firstRun) private(d_BlocksScores,d_result,d_LinesTotalScores,threadId)
{
		threadId = omp_get_thread_num();//get thread id
		cudaStreamCreate(& streams [ threadId ]);//init streams
		d_result=allocateResultOnGpu();//final result of one thread
		d_BlocksScores=allocateArraysOfFloat(offset,sizeDimX);//each place holds the score of one block
		d_LinesTotalScores=allocateFloatArrayZeroOnGpu(offset);//each place holds the total score of one line

#pragma omp single nowait
		d_Seq1=allocateStringOnGpu(seq1);
#pragma omp single nowait
		d_Coefficents = allocateCoefficentsOnGpu(coefficents);
#pragma omp single nowait
		d_Conservative=allocateArraysOfLettersArray(ConservativeGroups,MAX_CONSERTIVE_GROUPS);//	
#pragma omp single nowait
		d_SemiConservative=allocateArraysOfLettersArray(SemiConservativeGroups,MAX_SEMI_CONSERTIVE_GROUPS);
#pragma omp single nowait
		d_Seq2=allocateStringOnGpu(seq2);
#pragma omp barrier

#pragma omp for schedule(dynamic)
		for (int n = 0; n < strlen(seq2); n++) {
			for (int k = n + 1; k < strlen(seq2); k++) {
				if (run == TRUE) {
					// Launch the Kernel
  					compare<<<gridShapeCompare, NUM_BLOCK_THREADS,0,streams[threadId]>>>(n,k,d_Seq1,d_Seq2,strlen(seq2),d_Coefficents,d_Conservative,d_SemiConservative,d_BlocksScores);
 			 		err = cudaGetLastError();
 				   	if (err != cudaSuccess) 
 				   	{
 				       	fprintf(stderr, "Failed to launch compare kernel -  %s\n", cudaGetErrorString(err));
 				       	exit(EXIT_FAILURE);
 				   	}
					cudaStreamSynchronize(streams[threadId]);
					// Launch the Kernel
  					findMax<<<gridShapeFindMax, FIND_MAX_BLOCK_THREADS,0,streams[threadId]>>>(offset,sizeDimX,d_result,d_LinesTotalScores,d_BlocksScores);//findes the max score from array of offset scores
 			 		err = cudaGetLastError();
 				   	if (err != cudaSuccess) 
 				   	{
 				       	fprintf(stderr, "Failed to launch findMax kernel -  %s\n", cudaGetErrorString(err));
 				       	exit(EXIT_FAILURE);
 				   	}
					cudaStreamSynchronize(streams[threadId]);
					copyResultFromGpu(&resultTemp[threadId],d_result);
					compareResults(&h_Result[threadId],&firstRun,resultTemp[threadId].offsetPlace,n,k,resultTemp[threadId].score);
					if(h_Result[threadId].score==maxScore){
					# pragma omp atomic write
						run=FALSE;
					}								
				}
			}
		}
		#pragma omp barrier
		freeAllThreadsMemoryFromGpu(d_BlocksScores,d_LinesTotalScores,d_result);
		cudaStreamDestroy(streams[threadId]);	
	}
	getLargestResult(h_Result,result,NUM_THREADS);	
	freeAllMemoryFromGpu(d_Seq1,d_Seq2, d_Coefficents, d_Conservative, d_SemiConservative);
	free(resultTemp);
	free(h_Result);
	return 0;
}





////////////// other functions ////////////////////////////////




void copyResultFromGpu(Result *h_Result,Result *d_Result){
    // Copy the  result from GPU to the host memory.
   cudaError_t err = cudaSuccess;
    err = cudaMemcpy(h_Result, d_Result, 1* sizeof(Result), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy result from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

//////////// allocation functions//////////////////////////////


float** allocateArraysOfFloat(int numOfArrays, int arraySize){ //https://stackoverflow.com/questions/6561005/output-of-cuda-program-is-not-what-was-expected
	size_t bufferSize = arraySize*sizeof(float);
    float*  h_x[numOfArrays];
    float** d_x = 0;

    cudaMalloc( (void**)(&d_x), numOfArrays * sizeof(float*) );

    for ( int i=0; i < numOfArrays; i++ )
    {
        h_x[i] = NULL;
        cudaMalloc( (void**)(&h_x[i]), bufferSize );
    }

    cudaMemcpy( d_x, h_x, numOfArrays*sizeof(float*), cudaMemcpyHostToDevice);
  
return d_x;
}

LettersArray* allocateArraysOfLettersArray(LettersArray* array, int arraySize){ //https://stackoverflow.com/questions/6561005/output-of-cuda-program-is-not-what-was-expected
	
   LettersArray *d_LettersArray;
   size_t size = arraySize * sizeof(LettersArray);
   cudaError_t err = cudaSuccess;

    // Allocate memory on GPU to copy the data from the host
    err = cudaMalloc((void **)&d_LettersArray, size);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_LettersArray, array, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return d_LettersArray;
}

float* allocateCoefficentsOnGpu(float coefficents[]){

   float *d_Coefficents;
   size_t size = MAX_COEFFICENTS * sizeof(float);
   cudaError_t err = cudaSuccess;

    // Allocate memory on GPU to copy the data from the host
    err = cudaMalloc((void **)&d_Coefficents, size);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_Coefficents, coefficents, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return d_Coefficents;
}

float* allocateFloatArrayZeroOnGpu(int arraySize){
	float *arrayZero = (float*)calloc(arraySize,sizeof(float)); 
	float *d_array;
	size_t size = arraySize * sizeof(float);
	cudaError_t err = cudaSuccess;

	// Allocate memory on GPU to copy the data from the host
	err = cudaMalloc((void **)&d_array, size);
    if (err != cudaSuccess) {
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_array, arrayZero, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) 
    {
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    free(arrayZero);
    return d_array;
}

Result* allocateResultOnGpu(){
	Result *h_Result = (Result*)calloc(1,sizeof(Result));
	Result *d_Result;
	size_t size = 1 * sizeof(Result);
	cudaError_t err = cudaSuccess;

	// Allocate memory on GPU to copy the data from the host
	err = cudaMalloc((void **)&d_Result, size);
	if (err != cudaSuccess) 
	{
		fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	// Copy data from host to the GPU memory
	err = cudaMemcpy(d_Result, h_Result, size, cudaMemcpyHostToDevice);
	if (err != cudaSuccess) 
	{
		fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	free(h_Result);
    return d_Result;
}



char* allocateStringOnGpu(char* h_String){

    char *d_String;
   size_t size = strlen(h_String) * sizeof(char);
   cudaError_t err = cudaSuccess;

    // Allocate memory on GPU to copy the data from the host
    err = cudaMalloc((void **)&d_String, size);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_String, h_String, size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    return d_String;
}

//////////// free functions////////////////////////////////////

void freeAllThreadsMemoryFromGpu(float **d_BlocksScores,float *d_LinesTotalScores,Result *d_result){
		freeFloatArrays(d_BlocksScores);
		freeFloatAllocatedMemory(d_LinesTotalScores);
		freeResultAllocatedMemory(d_result);
}

void freeResultAllocatedMemory(Result* d_result){
	cudaError_t err = cudaSuccess;
	if (cudaFree(d_result) != cudaSuccess) {
        fprintf(stderr, "Failed to Result data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
} 

void freeCharAllocatedMemory(char* d_char){
	cudaError_t err = cudaSuccess;
	if (cudaFree(d_char) != cudaSuccess) {
        fprintf(stderr, "Failed to free char data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
} 



void freeFloatAllocatedMemory(float* d_float){
	cudaError_t err = cudaSuccess;
	if (cudaFree(d_float) != cudaSuccess) {
        fprintf(stderr, "Failed to free float data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
} 

void freeFloatArrays(float **d_float){
	cudaError_t err = cudaSuccess;
	if (cudaFree(d_float) != cudaSuccess) {
        fprintf(stderr, "Failed to free float data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

void freeLettersArrayAllocatedMemory(LettersArray* d_LettersArray){
	cudaError_t err = cudaSuccess;
	if (cudaFree(d_LettersArray) != cudaSuccess) {
        fprintf(stderr, "Failed to free float data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
} 

 

void freeAllMemoryFromGpu(char *d_Seq1,char *d_Seq2, float *d_Coefficents, LettersArray *d_Conservative, LettersArray *d_SemiConservative) {

	freeFloatAllocatedMemory(d_Coefficents);
 	freeCharAllocatedMemory(d_Seq1);
 	freeCharAllocatedMemory(d_Seq2);
	freeLettersArrayAllocatedMemory(d_Conservative);
	freeLettersArrayAllocatedMemory(d_SemiConservative);
}


