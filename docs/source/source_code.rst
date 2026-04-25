.. _source_code:

Source Code
==================================

VikNGS ships in two forms:

- **Browser app** (``web/``) — a vanilla TypeScript + Vite UI driving a
  WebAssembly build of the C++ core. This replaces the legacy Qt
  desktop GUI.
- **Command-line tools** (``src/cmd/``) — ``vikngs`` (analysis) and
  ``sim-cli`` (simulation), built from the same C++ core.

Compiling the Command-Line Tools
--------------------------------

The native build uses CMake (CMake 3.16+ and a C++17 compiler).
From the repo root::

   cmake -B build
   cmake --build build -j

This produces ``bin/vikngs`` and ``bin/sim-cli``.

A legacy ``bin/Makefile`` is also available::

   cd bin
   make

.. note:: g++ and C++17 or later is required.

Building the Browser App
------------------------

The browser app needs Emscripten installed (``emcmake`` in PATH).

1. Build the WASM module::

      ./emcmake-build.sh

2. Build and serve the UI::

      cd web
      npm install
      npm run dev

   then open the printed URL.

The C++ entry points used by both builds live in
``src/vikNGS.h``/``src/vikNGS.cpp`` (analysis) and
``src/Simulation/Simulation.h``/``Simulation.cpp`` (simulation).
