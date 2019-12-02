# Compiler:

PLAT =

CPP = /usr/local/bin/mpic++

# Available Compilation Flags for g++:
CPP_OPT = -O0
CPP_WARNINGS = 
#Wall
CPP_DEBUG = -g -Wall
CPP_VER = -std=c++11
CPP_STD = -lstdc++
OPEMMP = 
#-fopenmp

# Available Compilation Flags for MPI:

# Choosing the compilation Flags:

COMP_Flags = $(CPP_DEBUG) $(CPP_OPT) $(CPP_VER) $(OPEMMP)
EXEC_Flags = $(CPP_OPT) $(CPP_STD) $(OPEMMP)

# PathS to GMM library, source files and their headers: 
GMM = 
EIGEN = 
HDR = /Users/hamilton/Documents/GitHub/DHPC/include
SRC = /Users/hamilton/Documents/GitHub/DHPC/src

CXX_FLAGS =  -lpthread  -Wl, -Bdynamic -lm -lmpi

# Source file names:
CPP_SRCS += \
$(SRC)/main_dijkstra.cpp
#$(SRC)/main_A_Start.cpp  dijkstra

# Corresponding objects:
OBJS +=\
./main_dijkstra.o \

# Including the necessary libraries:
INCLUDES += \
-I$(EIGEN) \
-I$(GMM) \
-I$(HDR) \

# Generation of the executable:
%.o: $(SRC)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross C++ Compiler'
	$(CPP) $(COMP_Flags) $(INCLUDES) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" 
	@echo 'Finished building: $<'
	@echo ' '

all: aStart

aStart: $(OBJS)
	@echo 'Building Target: $@'
	$(CPP) $(EXEC_Flags) $(CXX_FLAGS) -o "aStart" $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

clean: 
	rm *.o *.d aStart
