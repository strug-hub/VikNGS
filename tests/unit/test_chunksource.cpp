// Tests for ChunkLineSource in src/Parser/File.h.
// Critical path: this is the line reader the WASM build will use, and it
// has to correctly handle chunk boundaries that fall mid-line.
#include "catch_amalgamated.hpp"
#include "Parser/File.h"

#include <string>
#include <vector>

static ChunkCallback feederFromChunks(std::vector<std::string> chunks) {
    auto state = std::make_shared<size_t>(0);
    auto data = std::make_shared<std::vector<std::string>>(std::move(chunks));
    return [state, data](std::string& out) -> bool {
        if (*state >= data->size()) return false;
        out = (*data)[*state];
        (*state)++;
        return true;
    };
}

TEST_CASE("ChunkLineSource handles single-chunk full input", "[chunksource]") {
    ChunkLineSource src(feederFromChunks({"line1\nline2\nline3\n"}));
    REQUIRE(src.hasNext());
    REQUIRE(src.nextLine() == "line1");
    REQUIRE(src.nextLine() == "line2");
    REQUIRE(src.nextLine() == "line3");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("ChunkLineSource handles line split across chunk boundary", "[chunksource]") {
    // "hello\nworld\n" split as ["hel", "lo\nwor", "ld\n"]
    ChunkLineSource src(feederFromChunks({"hel", "lo\nwor", "ld\n"}));
    REQUIRE(src.nextLine() == "hello");
    REQUIRE(src.nextLine() == "world");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("ChunkLineSource handles newline exactly at chunk boundary", "[chunksource]") {
    ChunkLineSource src(feederFromChunks({"abc\n", "def\n"}));
    REQUIRE(src.nextLine() == "abc");
    REQUIRE(src.nextLine() == "def");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("ChunkLineSource returns final line without trailing newline", "[chunksource]") {
    ChunkLineSource src(feederFromChunks({"a\nb\nc"}));
    REQUIRE(src.nextLine() == "a");
    REQUIRE(src.nextLine() == "b");
    REQUIRE(src.hasNext());
    REQUIRE(src.nextLine() == "c");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("ChunkLineSource handles empty lines", "[chunksource]") {
    ChunkLineSource src(feederFromChunks({"a\n\nb\n"}));
    REQUIRE(src.nextLine() == "a");
    REQUIRE(src.nextLine() == "");
    REQUIRE(src.nextLine() == "b");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("ChunkLineSource survives many small chunks", "[chunksource]") {
    // Each byte as its own chunk — worst case for buffer compaction.
    std::vector<std::string> bytes;
    std::string full = "alpha\nbeta\ngamma\n";
    for (char c : full) bytes.push_back(std::string(1, c));
    ChunkLineSource src(feederFromChunks(std::move(bytes)));
    REQUIRE(src.nextLine() == "alpha");
    REQUIRE(src.nextLine() == "beta");
    REQUIRE(src.nextLine() == "gamma");
    REQUIRE_FALSE(src.hasNext());
}

TEST_CASE("File::openStream produces same behavior as mmap for a simple input", "[chunksource]") {
    File f;
    f.openStream(feederFromChunks({"x\ny\nz\n"}));
    REQUIRE(f.hasNext());
    REQUIRE(f.nextLine() == "x");
    REQUIRE(f.nextLine() == "y");
    REQUIRE(f.nextLine() == "z");
    REQUIRE_FALSE(f.hasNext());
    REQUIRE(f.getLineNumber() == 3);
}
