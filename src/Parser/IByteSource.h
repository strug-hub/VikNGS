#pragma once
// IByteSource is the abstraction between the VCF/sample/BED parsers and
// whatever physical byte stream they're consuming.
//
// Native CLI builds use a memory-mapped file (MmapLineSource in File.h).
// WASM builds use a chunk-callback source where JS supplies bytes lazily
// (ChunkLineSource in File.h).
//
// The parsers only ever need one line at a time and never seek backward, so
// the interface is intentionally minimal: hasNext() + nextLine(). This keeps
// both backends simple and makes the streaming WASM case trivially supportable.
#include <string>
#include <cstdint>

class LineSource {
public:
    virtual ~LineSource() = default;

    // True if there may be more bytes to read; false once EOF is reached
    // AND the internal buffer is empty.
    virtual bool hasNext() = 0;

    // Returns the next line (without the trailing '\n'). Consumes the newline.
    // If the source ends without a trailing newline, the final partial line
    // is returned and hasNext() becomes false on the next call.
    virtual std::string nextLine() = 0;
};
