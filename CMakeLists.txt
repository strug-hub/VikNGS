cmake_minimum_required(VERSION 3.16)
project(VikNGS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -----------------------------------------------------------------------------
# Core library — Qt-free sources that the CLI, sim-cli, Catch2 unit tests, and
# WASM bindings all link against. This is the single source of truth for what
# "the core" is; everything downstream consumes it.
# -----------------------------------------------------------------------------
set(VIKNGS_SRC ${CMAKE_SOURCE_DIR}/src)

add_library(vikngs-core STATIC
    ${VIKNGS_SRC}/Global.cpp
    ${VIKNGS_SRC}/Request.cpp
    ${VIKNGS_SRC}/vikNGS.cpp
    ${VIKNGS_SRC}/cmd/Log.cpp
    ${VIKNGS_SRC}/Math/GeneticsHelper.cpp
    ${VIKNGS_SRC}/Math/RandomHelper.cpp
    ${VIKNGS_SRC}/Math/StatisticsHelper.cpp
    ${VIKNGS_SRC}/Math/VectorHelper.cpp
    ${VIKNGS_SRC}/Parser/BEDParser.cpp
    ${VIKNGS_SRC}/Parser/Filter.cpp
    ${VIKNGS_SRC}/Parser/InputProcess.cpp
    ${VIKNGS_SRC}/Parser/MemoryMapped/MemoryMapped.cpp
    ${VIKNGS_SRC}/Parser/SampleParser.cpp
    ${VIKNGS_SRC}/Parser/StringTools.cpp
    ${VIKNGS_SRC}/Parser/VariantParser.cpp
    ${VIKNGS_SRC}/Test/ScoreTestFunctions.cpp
    ${VIKNGS_SRC}/Test/Test.cpp
)

target_include_directories(vikngs-core PUBLIC ${VIKNGS_SRC})
target_compile_definitions(vikngs-core PUBLIC _GNU_SOURCE)

find_package(Threads REQUIRED)
target_link_libraries(vikngs-core PUBLIC Threads::Threads)

# -----------------------------------------------------------------------------
# Simulation library — pure-C++ but lives in the GUI tree for historical
# reasons. Separate target so the CLI doesn't need to link it.
# -----------------------------------------------------------------------------
add_library(vikngs-sim STATIC
    ${VIKNGS_SRC}/gui/src/simulation/Simulation.cpp
    ${VIKNGS_SRC}/gui/src/simulation/SimulationHelper.cpp
)
target_link_libraries(vikngs-sim PUBLIC vikngs-core)

# -----------------------------------------------------------------------------
# CLI binaries — native builds only (skip under Emscripten).
# -----------------------------------------------------------------------------
if(NOT EMSCRIPTEN)
    add_executable(vikngs ${VIKNGS_SRC}/cmd/vikNGScmd.cpp)
    target_link_libraries(vikngs PRIVATE vikngs-core)

    add_executable(sim-cli ${VIKNGS_SRC}/cmd/simcmd.cpp)
    target_link_libraries(sim-cli PRIVATE vikngs-sim)

    # Match legacy bin/ output location so existing tests find the binaries.
    set_target_properties(vikngs sim-cli PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
endif()

# -----------------------------------------------------------------------------
# WASM target — only active when building with Emscripten (emcmake cmake ...).
# Produces vikngs-core.{js,wasm} in web/public/ for the browser app to load.
# -----------------------------------------------------------------------------
if(EMSCRIPTEN)
    # Core throws std::runtime_error on bad input (caught in bindings.cpp),
    # so every TU needs -fexceptions. InputProcess.cpp is already patched to
    # use std::launch::deferred under __EMSCRIPTEN__ — no pthreads required.
    target_compile_options(vikngs-core PRIVATE -fexceptions)
    target_compile_options(vikngs-sim  PRIVATE -fexceptions)

    add_executable(vikngs-core-wasm ${VIKNGS_SRC}/wasm/bindings.cpp)
    target_compile_options(vikngs-core-wasm PRIVATE -fexceptions)
    target_link_libraries(vikngs-core-wasm PRIVATE vikngs-core)

    set_target_properties(vikngs-core-wasm PROPERTIES
        OUTPUT_NAME "vikngs-core"
        SUFFIX ".js"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/web/public
    )

    target_link_options(vikngs-core-wasm PRIVATE
        "-lembind"
        "-O3"
        "SHELL:-s MODULARIZE=1"
        "SHELL:-s EXPORT_NAME=createVikNGS"
        "SHELL:-s ENVIRONMENT=worker,web,node"
        "SHELL:-s ALLOW_MEMORY_GROWTH=1"
        "SHELL:-s INITIAL_MEMORY=256MB"
        "SHELL:-s EXPORTED_RUNTIME_METHODS=['FS']"
        "SHELL:-s FORCE_FILESYSTEM=1"
        "SHELL:-s ASSERTIONS=1"
        "-fexceptions"
    )
endif()
