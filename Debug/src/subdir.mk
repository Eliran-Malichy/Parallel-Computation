################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/Result.c \
../src/cFunctions.c \
../src/main.c 

O_SRCS += \
../src/Result.o \
../src/cFunctions.o \
../src/cudaFunctions.o \
../src/main.o 

C_DEPS += \
./src/Result.d \
./src/cFunctions.d \
./src/main.d 

OBJS += \
./src/Result.o \
./src/cFunctions.o \
./src/main.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	mpicc -fopenmp -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/Result.d ./src/Result.o ./src/cFunctions.d ./src/cFunctions.o ./src/main.d ./src/main.o

.PHONY: clean-src

