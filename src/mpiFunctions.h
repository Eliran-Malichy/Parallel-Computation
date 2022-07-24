/*
 * mpiFunctions.h
 *
 *  Created on: 18 Feb 2022
 *      Author: eliran
 */

#ifndef MPIFUNCTIONS_H_
#define MPIFUNCTIONS_H_
#include<mpi.h>
#define MAX_VALUE(a, b) ((a)>(b)? (a) : (b))
#define MIN_VALUE(a, b) ((a)>(b)? (b) : (a))

void createResultlMpiType(MPI_Datatype *ResultMPIType);

#endif /* MPIFUNCTIONS_H_ */
