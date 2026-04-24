#pragma once
#include "IByteSource.h"
#include "MemoryMapped/MemoryMapped.h"

#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

// ----------------------------------------------------------------------------
// Mmap-backed line source — the original behavior. Reads in 4GB segments to
// tolerate files larger than addressable memory on 32-bit targets.
// ----------------------------------------------------------------------------
class MmapLineSource : public LineSource {
    MemoryMapped mmap;
    uint64_t pos = 0;
    uint64_t segmentSize = 0;
    uint64_t pageSize = 0;
    uint64_t pagesPerSegment = 0;
    uint64_t currentPage = 0;
    bool lastSegment = false;
    static constexpr uint64_t MEMORY = 4ULL * 1024 * 1024 * 1024; // 4 GB

public:
    void open(const std::string& directory) {
        pageSize = mmap.getpagesize();
        pagesPerSegment = std::max<uint64_t>(1, MEMORY / pageSize);
        mmap.open(directory, pagesPerSegment * pageSize, MemoryMapped::CacheHint::SequentialScan);
        currentPage = pagesPerSegment;
        pos = 0;
        segmentSize = mmap.mappedSize();
        lastSegment = (segmentSize >= mmap.size());
    }

    bool hasNext() override {
        return !lastSegment || (pos < segmentSize);
    }

    std::string nextLine() override {
        std::string line;
        while (true) {
            if (pos >= segmentSize) {
                if (lastSegment) break;
                remap();
            }
            if (mmap[pos] == '\n') break;
            line += mmap[pos];
            pos++;
        }
        pos++;
        return line;
    }

private:
    void remap() {
        size_t offset = currentPage * pageSize;
        size_t mapSize = pagesPerSegment * pageSize;
        mmap.remap(offset, mapSize);
        segmentSize = mmap.mappedSize();
        lastSegment = (offset + mapSize) >= mmap.size();
        currentPage += pagesPerSegment;
        pos = 0;
    }
};

// ----------------------------------------------------------------------------
// Chunk-callback-backed line source — used by WASM where JS pulls bytes
// from File.stream() and pushes them in via the callback. The callback
// returns true with bytes appended to `out`, or false for EOF.
//
// Also usable natively for tests: construct with a callback that reads from
// std::ifstream.
// ----------------------------------------------------------------------------
using ChunkCallback = std::function<bool(std::string& out)>;

class ChunkLineSource : public LineSource {
    ChunkCallback fetchMore;
    std::string buf;        // unconsumed bytes (including partial final line)
    size_t bufPos = 0;
    bool eof = false;

public:
    explicit ChunkLineSource(ChunkCallback cb) : fetchMore(std::move(cb)) {}

    bool hasNext() override {
        if (bufPos < buf.size()) return true;
        if (eof) return false;
        refill();
        return bufPos < buf.size();
    }

    std::string nextLine() override {
        std::string line;
        while (true) {
            if (bufPos >= buf.size()) {
                if (eof) break;
                refill();
                if (bufPos >= buf.size()) break;
            }
            char c = buf[bufPos++];
            if (c == '\n') break;
            line += c;
        }
        return line;
    }

private:
    void refill() {
        // Compact the buffer to drop already-consumed bytes.
        if (bufPos > 0) {
            buf.erase(0, bufPos);
            bufPos = 0;
        }
        std::string chunk;
        bool got = fetchMore(chunk);
        if (!got || chunk.empty()) {
            eof = true;
            return;
        }
        buf += chunk;
    }
};

// ----------------------------------------------------------------------------
// File — public type used by all parsers. Wraps a LineSource so the call
// sites (VariantParser, SampleParser, BEDParser, InputProcess) don't need to
// know which backing mode is active.
//
// Default-constructed + open(path) preserves existing behavior (mmap).
// openStream(cb) is the new entry point for WASM / test stubs.
// ----------------------------------------------------------------------------
struct File {
    std::unique_ptr<LineSource> src;
    uint64_t lineNumber = 0;

    void open(const std::string& directory) {
        auto mm = std::make_unique<MmapLineSource>();
        mm->open(directory);
        src = std::move(mm);
        lineNumber = 0;
    }

    void openStream(ChunkCallback cb) {
        src = std::make_unique<ChunkLineSource>(std::move(cb));
        lineNumber = 0;
    }

    void close() {
        src.reset();
    }

    bool hasNext() {
        return src && src->hasNext();
    }

    std::string nextLine() {
        lineNumber++;
        return src->nextLine();
    }

    int getLineNumber() {
        return static_cast<int>(lineNumber);
    }
};
