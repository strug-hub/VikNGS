#pragma once
#include "MemoryMapped/MemoryMapped.h"
#include <algorithm>
#include <stdexcept>

struct File {
    MemoryMapped mmap;
    uint64_t pos;

    uint64_t segmentSize;
    uint64_t pageSize;
    uint64_t pagesPerSegment;
    uint64_t currentPage;

    bool lastSegment;
    uint64_t lineNumber;
    uint64_t memory = static_cast<uint64_t>(4e9); //4gb

    inline void open(std::string directory) {
        try {

            pageSize = static_cast<uint64_t>(mmap.getpagesize());
            pagesPerSegment = memory/pageSize;
            uint64_t one = 1;
            pagesPerSegment = std::max(one, pagesPerSegment);

            mmap.open(directory, static_cast<size_t>(pagesPerSegment * pageSize), MemoryMapped::CacheHint::SequentialScan);
            currentPage = pagesPerSegment;

            pos = 0;
            lineNumber = 0;
            segmentSize = mmap.mappedSize();

            lastSegment = false;
            if(segmentSize >= mmap.size())
                lastSegment = true;

        }
        catch (...) {
            throw std::runtime_error("Cannot open file from provided directory: " + directory);
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

            line += static_cast<char>(mmap[pos]);
            pos++;
        }

        lineNumber++;
        pos++;

        return line;
    }

    inline uint64_t getLineNumber() {
        return lineNumber;
    }

    inline bool hasNext() {
        return !lastSegment || (pos < segmentSize);
    }

    void remap() {

        uint64_t offset = currentPage * pageSize;
        uint64_t mapSize = pagesPerSegment * pageSize;

        mmap.remap(offset, mapSize);
        segmentSize = mmap.mappedSize();

        lastSegment = (offset + mapSize) >= mmap.size();
        currentPage += pagesPerSegment;
        pos=0;

    }
};
