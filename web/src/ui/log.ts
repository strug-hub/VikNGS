export function mountLog(el: HTMLElement) {
    const append = (level: "info" | "error" | "ok", text: string) => {
        const line = document.createElement("span");
        line.className = level;
        line.textContent = text + "\n";
        el.appendChild(line);
        el.scrollTop = el.scrollHeight;
    };
    return {
        info:  (t: string) => append("info",  t),
        ok:    (t: string) => append("ok",    t),
        error: (t: string) => append("error", t),
        clear: () => { el.innerHTML = ""; },
    };
}
