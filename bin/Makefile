CC=g++ -std=c++11 -D_GNU_SOURCE
CFLAGS=-c
SOURCE=../src/

all: build

build: root parser test RVS.o 
	$(CC) RVS.o Log.o \
MemoryMapped.o InputParser.o ParserTools.o \
CompQuadForm.o RequestBuilder.o VectorHelper.o StatisticsHelper.o \
VariantFilterUtils.o VariantFilter.o SampleParserUtils.o SampleParser.o \
BEDParserUtils.o BEDParser.o VCFParserUtils.o VCFParser.o \
CommonTest.o RareTest.o RareTestObject.o \
-pthread -o vikNGS
	
parser: MemoryMapped.o VCFParser.o SampleParser.o BEDParser.o VariantFilter.o InputParser.o ParserTools.o 

MemoryMapped.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/MemoryMapped/MemoryMapped.cpp
InputParser.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/InputParser.cpp
ParserTools.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/ParserTools.cpp

VariantFilter.o: VariantFilterUtils.o
	$(CC) $(CFLAGS) $(SOURCE)Parser/Filter/VariantFilter.cpp
VariantFilterUtils.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/Filter/VariantFilterUtils.cpp

VCFParser.o: VCFParserUtils.o
	$(CC) $(CFLAGS) $(SOURCE)Parser/VCF/VCFParser.cpp
VCFParserUtils.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/VCF/VCFParserUtils.cpp

SampleParser.o: SampleParserUtils.o 
	$(CC) $(CFLAGS) $(SOURCE)Parser/Sample/SampleParser.cpp
SampleParserUtils.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/Sample/SampleParserUtils.cpp

BEDParser.o: BEDParserUtils.o
	$(CC) $(CFLAGS) $(SOURCE)Parser/BED/BEDParser.cpp
BEDParserUtils.o:
	$(CC) $(CFLAGS) $(SOURCE)Parser/BED/BEDParserUtils.cpp

root: Log.o StatisticsHelper.o CompQuadForm.o RequestBuilder.o VectorHelper.o 

Log.o:
	$(CC) $(CFLAGS) $(SOURCE)Log.cpp
StatisticsHelper.o:
	$(CC) $(CFLAGS) -w $(SOURCE)StatisticsHelper.cpp
CompQuadForm.o:
	$(CC) $(CFLAGS) $(SOURCE)CompQuadForm.cpp
RequestBuilder.o:
	$(CC) $(CFLAGS) $(SOURCE)RequestBuilder.cpp
VectorHelper.o: 
	$(CC) $(CFLAGS) -w $(SOURCE)VectorHelper.cpp

test: CommonTest.o RareTest.o RareTestObject.o
CommonTest.o:
	$(CC) $(CFLAGS) $(SOURCE)Test/CommonTest.cpp
RareTest.o:
	$(CC) $(CFLAGS) $(SOURCE)Test/RareTest.cpp
RareTestObject.o: 
	$(CC) $(CFLAGS) $(SOURCE)Test/RareTestObject.cpp


RVS.o: 
	$(CC) $(CFLAGS) $(SOURCE)RVS.cpp

clean:
	rm -rf *o all
