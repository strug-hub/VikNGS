CC=g++ -std=c++11 -D_GNU_SOURCE -O3
CFLAGS=-c
SOURCE=../src/

all: build

build: root math parser test RVS.o 
	$(CC) RVS.o Log.o \
MemoryMapped.o InputParser.o ParserTools.o \
CompQuadForm.o RequestBuilder.o VectorHelper.o StatisticsHelper.o \
VariantFilterUtils.o VariantFilter.o SampleParserUtils.o SampleParser.o \
BEDParserUtils.o BEDParser.o VCFParserUtils.o VCFParser.o \
CommonTest.o RareTest.o RareTestObject.o RareTestCollapseObject.o CommonTestObject.o \
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

root: Log.o RequestBuilder.o

Log.o:
	$(CC) $(CFLAGS) $(SOURCE)Log.cpp
RequestBuilder.o:
	$(CC) $(CFLAGS) $(SOURCE)RequestBuilder.cpp

math: StatisticsHelper.o CompQuadForm.o VectorHelper.o 

StatisticsHelper.o:
	$(CC) $(CFLAGS) -w $(SOURCE)Math/StatisticsHelper.cpp
CompQuadForm.o:
	$(CC) $(CFLAGS) $(SOURCE)Math/CompQuadForm.cpp
VectorHelper.o: 
	$(CC) $(CFLAGS) -w $(SOURCE)Math/VectorHelper.cpp

test: CommonTestObject.o RareTestObject.o RareTestCollapseObject.o CommonTest.o RareTest.o
CommonTest.o:
	$(CC) $(CFLAGS) $(SOURCE)Test/CommonTest.cpp
RareTest.o:
	$(CC) $(CFLAGS) $(SOURCE)Test/RareTest.cpp
RareTestObject.o: 
	$(CC) $(CFLAGS) $(SOURCE)Test/RareTestObject.cpp
CommonTestObject.o: 
	$(CC) $(CFLAGS) $(SOURCE)Test/CommonTestObject.cpp
RareTestCollapseObject.o: 
	$(CC) $(CFLAGS) $(SOURCE)Test/RareTestCollapseObject.cpp

RVS.o: 
	$(CC) $(CFLAGS) $(SOURCE)RVS.cpp

clean:
	rm -rf *o all
