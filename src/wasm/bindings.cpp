// Emscripten bindings for the VikNGS core.
//
// Phase C: minimal smoke build — exposes a hello() function so we can
// confirm the WASM toolchain produces a loadable module.
// Phase D will grow this file to expose Request / runVikNGS / result rows.
#include <emscripten/bind.h>
#include <string>

static std::string hello(const std::string& name) {
    return "hello, " + name + ", from vikngs-core wasm";
}

EMSCRIPTEN_BINDINGS(vikngs_core) {
    emscripten::function("hello", &hello);
}
