CC=g++ -std=c++11 -O3 -D_GNU_SOURCE 
CFLAGS=-c
SOURCE=../src/
CMD=cmd/
all: build

build: vikNGS.o root math test parser Global.o vikNGScmd.o 
	$(CC) Log.o Request.o MemoryMapped.o \
VectorHelper.o GeneticsHelper.o RandomHelper.o StatisticsHelper.o \
StringTools.o VariantParser.o Filter.o SampleParser.o BEDParser.o  \
Test.o ScoreTestFunctions.o InputProcess.o vikNGS.o $(OUT)Global.o vikNGScmd.o \
-pthread -o vikNGS
	
root: Log.o Request.o MemoryMapped.o

Log.o:
	$(CC) $(CFLAGS) $(SOURCE)$(CMD)Log.cpp
Global.o:
	$(CC) $(CFLAGS) $(SOURCE)Global.cpp
Request.o:
	$(CC) $(CFLAGS) $(SOURCE)Request.cpp


math: VectorHelper.o GeneticsHelper.o RandomHelper.o StatisticsHelper.o

StatisticsHelper.o:
	$(CC) $(CFLAGS) -w $(SOURCE)Math/StatisticsHelper.cpp
RandomHelper.o:
	$(CC) $(CFLAGS) $(SOURCE)Math/RandomHelper.cpp
GeneticsHelper.o:
	$(CC) $(CFLAGS) $(SOURCE)Math/GeneticsHelper.cpp
VectorHelper.o: 
	$(CC) $(CFLAGS) $(SOURCE)Math/VectorHelper.cpp


test: Test.o ScoreTestFunctions.o
Test.o:
	$(CC) $(CFLAGS) $(SOURCE)Test/Test.cpp
ScoreTestFunctions.o: 
	$(CC) $(CFLAGS) $(SOURCE)Test/ScoreTestFunctions.cpp


parser: StringTools.o SampleParser.o VariantParser.o BEDParser.o Filter.o InputProcess.o  

MemoryMapped.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/MemoryMapped/MemoryMapped.cpp
InputProcess.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/InputProcess.cpp
StringTools.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/StringTools.cpp

Filter.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/Filter.cpp
VariantParser.o: 
	$(CC) $(CFLAGS) $(SOURCE)Parser/VariantParser.cpp
SampleParser.o:  
	$(CC) $(CFLAGS) $(SOURCE)Parser/SampleParser.cpp

BEDParser.o: 
	$(CC) $(CFLAGS) $(SOURCE)Parser/BEDParser.cpp



vikNGScmd.o: 
	$(CC) $(CFLAGS) $(SOURCE)$(CMD)vikNGScmd.cpp
vikNGS.o: 
	$(CC) $(CFLAGS) $(SOURCE)vikNGS.cpp
clean:
	rm -rf *o all
