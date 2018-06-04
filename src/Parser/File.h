#pragma once
#include "MemoryMapped/MemoryMapped.h"
#include <algorithm>

struct File {
    MemoryMapped mmap;
    uint64_t pos;

    uint64_t segmentSize;
    uint64_t pageSize;
    uint64_t pagesPerSegment;
    uint64_t currentPage;

    bool lastSegment;
    uint64_t lineNumber;
    uint64_t memory = 4e9; //4gb

    inline void open(std::string directory) {
        try {

            pageSize = mmap.getpagesize();
            pagesPerSegment = memory/pageSize;
            uint64_t one = 1;
            pagesPerSegment = std::max(one, pagesPerSegment);

            mmap.open(directory, pagesPerSegment * pageSize, MemoryMapped::CacheHint::SequentialScan);
            currentPage = pagesPerSegment;

            pos = 0;
            lineNumber = 0;
            segmentSize = mmap.mappedSize();

            lastSegment = false;
            if(segmentSize >= mmap.size())
                lastSegment = true;

        }
        catch (...) {
            throwError("file struct", "Cannot open file from provided directory.", directory);
        }
    }

    inline void close() {

        mmap.close();
    }

    inline std::string nextLine() {

        std::string line = "";

        while (true) {
            if(pos >= segmentSize){

                if(lastSegment)
                    break;
                else
                    remap();
            }

            if (mmap[pos] == '\n')
                break;

            line+= mmap[pos];
            pos++;
        }

        lineNumber++;
        pos++;

        return line;
    }

    inline int getLineNumber() {
        return lineNumber;
    }

    inline bool hasNext() {
        return !lastSegment || (pos < segmentSize);
    }

    void remap() {

        printInfo("Reading another 4GB");
        size_t offset = currentPage * pageSize;
        size_t mapSize = pagesPerSegment * pageSize;

        mmap.remap(offset, mapSize);
        segmentSize = mmap.mappedSize();

        lastSegment = (offset + mapSize) >= mmap.size();
        currentPage += pagesPerSegment;
        pos=0;

    }
};
