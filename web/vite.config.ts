import { defineConfig } from "vite";

export default defineConfig({
    // The Emscripten-generated .cjs references itself by filename at runtime
    // (to locate the .wasm). Leave the module alone — Vite must not bundle it.
    optimizeDeps: { exclude: ["./public/vikngs-core.cjs"] },
    worker: { format: "es" },
    server: {
        fs: { allow: [".."] },
    },
});
