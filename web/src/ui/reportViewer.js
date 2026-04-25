// Standalone viewer used by the exported HTML report. Mirrors
// renderSimResults() from simResults.ts but speaks plain JS and assumes
// uPlot is available as a global (iife build inlined alongside).
//
// Input shape:
//   data = {
//     rows: [{ stepIdx, sampleSize, testIdx, statName, genotypeSource,
//              variantIdx, pvalue }, ...],
//     steps,
//     family: "binomial" | "normal",
//     groups: [{ n, nIncrement, isCase, normalMean, normalSd, meanDepth,
//                sdDepth, errorRate, readDepth }, ...],
//     summary: { variants, processingTime, evaluationTime },
//   }
//   els = { summary, power, qq, hist, sample }

(function (root) {
    "use strict";
    var COLORS = ["#3c7ccc", "#a62a2a", "#3f7a4a", "#b68a00", "#8c46c6", "#2d8ba3"];
    var ALPHA = 0.05;

    function cohortLabel(g, family) {
        if (family === "normal") return "normal";
        return g.isCase ? "case" : "control";
    }
    function fmtDepth(g) {
        return g.meanDepth.toFixed(1) + " ± " + g.sdDepth.toFixed(1) + " (" + g.readDepth + ")";
    }
    function sampleSizeAtStep(g, step) { return g.n + g.nIncrement * step; }
    function rangeText(g, steps) {
        if (g.nIncrement === 0 || steps <= 1) return String(g.n);
        return g.n + " – " + sampleSizeAtStep(g, steps - 1);
    }
    function escapeHtml(s) {
        return String(s).replace(/[&<>"']/g, function (c) {
            return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
        });
    }

    function renderSampleTable(host, groups, family, step, steps) {
        var stepLabel = step === null
            ? (steps > 1 ? "min – max" : "")
            : "step " + (step + 1);
        var html = '<h3>Sample info' + (stepLabel ? " (" + escapeHtml(stepLabel) + ")" : "")
            + '</h3><table><thead><tr><th>n</th><th>cohort</th><th>depth</th><th>err</th></tr></thead><tbody>';
        for (var i = 0; i < groups.length; i++) {
            var g = groups[i];
            var n = step === null ? rangeText(g, steps) : sampleSizeAtStep(g, step);
            html += '<tr><td>' + escapeHtml(n) + '</td><td>' + escapeHtml(cohortLabel(g, family))
                + '</td><td>' + escapeHtml(fmtDepth(g)) + '</td><td>' + escapeHtml(g.errorRate.toFixed(3))
                + '</td></tr>';
        }
        host.innerHTML = html + '</tbody></table>';
    }

    function computePower(ps, alpha) {
        if (!ps.length) return 0;
        var k = 0;
        for (var i = 0; i < ps.length; i++) if (ps[i] < alpha) k++;
        return k / ps.length;
    }

    var activePlots = [];
    function destroyAllPlots() {
        for (var i = 0; i < activePlots.length; i++) {
            try { activePlots[i].plot.destroy(); } catch (_) { /* */ }
            try { activePlots[i].ro.disconnect(); } catch (_) { /* */ }
        }
        activePlots.length = 0;
    }

    function mountPlot(host, opts, data) {
        for (var i = activePlots.length - 1; i >= 0; i--) {
            if (activePlots[i].host === host) {
                try { activePlots[i].plot.destroy(); } catch (_) { /* */ }
                try { activePlots[i].ro.disconnect(); } catch (_) { /* */ }
                activePlots.splice(i, 1);
            }
        }
        host.innerHTML = "";
        function measure() {
            return {
                width: Math.max(50, host.clientWidth),
                height: Math.max(50, host.clientHeight),
            };
        }
        var sz = measure();
        opts.width = sz.width;
        opts.height = sz.height;
        var plot = new uPlot(opts, data, host);
        var ro = new ResizeObserver(function () { plot.setSize(measure()); });
        ro.observe(host);
        activePlots.push({ host: host, plot: plot, ro: ro });
    }

    function drawQQ(el, ps, label, step, alpha) {
        if (!ps.length) return;
        var sorted = ps.slice().sort(function (a, b) { return a - b; });
        var n = sorted.length;
        var xs = [], ys = [];
        for (var i = 0; i < n; i++) {
            xs.push(-Math.log10((i + 1) / (n + 1)));
            ys.push(-Math.log10(Math.max(sorted[i], 1e-300)));
        }
        var diagY = xs.slice();
        var alphaY = xs.map(function () { return -Math.log10(alpha); });
        mountPlot(el, {
            title: "QQ — " + label + " (step " + (step + 1) + ")",
            scales: { x: { time: false }, y: { auto: true } },
            axes: [
                { stroke: "#6f6354", label: "-log10 theoretical" },
                { stroke: "#6f6354", label: "-log10 observed" },
            ],
            series: [
                {},
                { label: "obs", stroke: "#3c7ccc", points: { show: true, size: 4 }, paths: function () { return null; } },
                { label: "y=x", stroke: "#b9ab8d", dash: [4, 4] },
                { label: "α=" + alpha, stroke: "#a62a2a", dash: [2, 4] },
            ],
        }, [xs, ys, diagY, alphaY]);
    }

    function drawHist(el, ps, label, step) {
        if (!ps.length) return;
        var BINS = 20;
        var counts = new Array(BINS).fill(0);
        for (var i = 0; i < ps.length; i++) {
            var b = Math.min(BINS - 1, Math.floor(ps[i] * BINS));
            counts[b]++;
        }
        var xs = counts.map(function (_, i) { return (i + 0.5) / BINS; });
        mountPlot(el, {
            title: "Histogram — " + label + " (step " + (step + 1) + ")",
            scales: { x: { time: false, range: [0, 1] }, y: { auto: true } },
            axes: [
                { stroke: "#6f6354", label: "p-value" },
                { stroke: "#6f6354", label: "count" },
            ],
            series: [
                {},
                { label: "count", stroke: "#3c7ccc", fill: "#a8caf0",
                  paths: uPlot.paths.bars({ size: [0.9], align: 0 }) },
            ],
        }, [xs, counts]);
    }

    function render(data, els) {
        var rows = data.rows;
        var steps = data.steps;
        var groups = data.groups;
        var family = data.family;
        var summary = data.summary;

        if (els.summary) {
            els.summary.innerHTML = "<strong>" + summary.variants + "</strong> variants &middot; "
                + "processing <strong>" + summary.processingTime.toFixed(2) + "s</strong> &middot; "
                + "evaluation <strong>" + summary.evaluationTime.toFixed(2) + "s</strong> &middot; "
                + "<strong>" + rows.length + "</strong> p-values";
        }

        renderSampleTable(els.sample, groups, family, steps > 1 ? null : 0, steps);

        var byKey = new Map();
        for (var i = 0; i < rows.length; i++) {
            var r = rows[i];
            var key = r.statName + "/" + r.genotypeSource;
            var entry = byKey.get(key);
            if (!entry) {
                entry = { pointsByStep: new Map(), sampleSizeByStep: new Map() };
                byKey.set(key, entry);
            }
            var arr = entry.pointsByStep.get(r.stepIdx);
            if (!arr) { arr = []; entry.pointsByStep.set(r.stepIdx, arr); }
            arr.push(r.pvalue);
            entry.sampleSizeByStep.set(r.stepIdx, r.sampleSize);
        }
        var keys = Array.from(byKey.keys()).sort();
        var stepIndices = [];
        for (var s = 0; s < steps; s++) stepIndices.push(s);
        var anyKey = byKey.get(keys[0]);
        var xs = stepIndices.map(function (s) {
            return anyKey ? (anyKey.sampleSizeByStep.get(s) || s) : s;
        });
        var powerData = [xs];
        for (var k = 0; k < keys.length; k++) {
            var entry2 = byKey.get(keys[k]);
            powerData.push(stepIndices.map(function (s) {
                return computePower(entry2.pointsByStep.get(s) || [], ALPHA);
            }));
        }

        var lastHover = null;
        mountPlot(els.power, {
            title: "Power (alpha=" + ALPHA + ")",
            cursor: { show: true },
            scales: { x: { time: false }, y: { auto: false, range: [0, 1] } },
            axes: [
                { stroke: "#6f6354", label: "Sample size" },
                { stroke: "#6f6354", label: "Power" },
            ],
            series: [{}].concat(keys.map(function (k, i) {
                return { label: k, stroke: COLORS[i % COLORS.length], points: { show: true, size: 6 } };
            })),
            hooks: {
                setCursor: [function (u) {
                    var idx = u.cursor.idx;
                    var next = (idx == null) ? null : idx;
                    if (next !== lastHover) {
                        lastHover = next;
                        renderSampleTable(els.sample, groups, family, next, steps);
                    }
                }],
            },
        }, powerData);

        function makePanel(host, onChange) {
            host.innerHTML = "";
            var bar = document.createElement("div");
            bar.style.cssText = "display:flex;gap:6px;margin-bottom:4px;font-size:11px;";
            var testSel = document.createElement("select");
            for (var i = 0; i < keys.length; i++) {
                var o = document.createElement("option");
                o.value = keys[i]; o.textContent = keys[i];
                testSel.appendChild(o);
            }
            var stepSel = document.createElement("select");
            for (var s = 0; s < steps; s++) {
                var o2 = document.createElement("option");
                o2.value = String(s);
                o2.textContent = steps > 1 ? "step " + (s + 1) + " (n=" + xs[s] + ")" : "n=" + xs[s];
                stepSel.appendChild(o2);
            }
            function fire() { onChange(testSel.value, parseInt(stepSel.value, 10)); }
            testSel.addEventListener("change", fire);
            stepSel.addEventListener("change", fire);
            bar.appendChild(testSel);
            bar.appendChild(stepSel);
            var plot = document.createElement("div");
            plot.className = "sim-plot";
            host.appendChild(bar);
            host.appendChild(plot);
            return plot;
        }

        var qqHost = makePanel(els.qq, function (k, st) {
            drawQQ(qqHost, byKey.get(k).pointsByStep.get(st) || [], k, st, ALPHA);
        });
        var histHost = makePanel(els.hist, function (k, st) {
            drawHist(histHost, byKey.get(k).pointsByStep.get(st) || [], k, st);
        });
        var initialKey = keys[0];
        drawQQ(qqHost, byKey.get(initialKey).pointsByStep.get(0) || [], initialKey, 0, ALPHA);
        drawHist(histHost, byKey.get(initialKey).pointsByStep.get(0) || [], initialKey, 0);
    }

    root.VikngsReport = { render: render };
})(typeof window !== "undefined" ? window : globalThis);
