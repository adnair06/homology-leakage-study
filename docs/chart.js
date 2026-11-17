/* Small canvas chart. No dependencies — the page is one folder of
   static files and stays that way. */

const C_RANDOM    = "#139A43";
const C_CLUSTERED = "#0A3200";
const C_CYAN      = "#8EF9F3";
const C_RULE      = "#c3d8c3";
const C_MUTED     = "#3d6b47";

function setupCanvas(cv, cssW, cssH) {
    const dpr = window.devicePixelRatio || 1;
    cv.width  = Math.round(cssW * dpr);
    cv.height = Math.round(cssH * dpr);
    cv.style.height = cssH + "px";
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
}

/* series: [{pts:[{band,n,v}|null], noHit, color, style, label}]
   style is "solid" | "dashed" | "dotted" */
function drawBandChart(cv, opts) {
    const cssW = cv.clientWidth || 640;
    const cssH = opts.height || 370;
    const ctx  = setupCanvas(cv, cssW, cssH);

    // Narrow screens: rotate the band labels and drop the per-band counts,
    // which otherwise overlap into an unreadable smear.
    const compact = cssW < 620;
    const L = compact ? 52 : 54;
    const R = compact ? 16 : 42;
    const T = 16;
    const B = compact ? 78 : 62;
    const plotW = cssW - L - R;
    const plotH = cssH - T - B;

    const bands = opts.bands;            // x categories, excl. "no hit"
    const yMin = opts.yMin, yMax = opts.yMax;

    // "no hit" sits in its own gutter on the left, disconnected.
    const gutter = opts.showNoHit ? 80 : 0;
    const innerW = plotW - gutter;
    const step = bands.length > 1 ? innerW / (bands.length - 1) : 0;
    const xOf = i => L + gutter + i * step;
    const xNoHit = L + gutter * 0.40;
    const yOf = v => T + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

    ctx.clearRect(0, 0, cssW, cssH);

    // shaded >90% column, clipped to the plot area
    if (opts.highlightLast) {
        const x0 = xOf(bands.length - 1) - step * 0.45;
        ctx.fillStyle = C_CYAN;
        ctx.globalAlpha = 0.32;
        ctx.fillRect(x0, T, (cssW - R) - x0, plotH);
        ctx.globalAlpha = 1;
    }

    // horizontal gridlines + y labels
    ctx.font = (compact ? 10 : 11) + "px Verdana, sans-serif";
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    const ticks = opts.ticks;
    ticks.forEach(t => {
        const y = yOf(t);
        ctx.strokeStyle = "#e8efe8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(L, Math.round(y) + 0.5);
        ctx.lineTo(cssW - R, Math.round(y) + 0.5);
        ctx.stroke();
        ctx.fillStyle = C_MUTED;
        ctx.fillText(t.toFixed(2), L - (compact ? 5 : 8), y);
    });

    // reference line
    if (opts.refLine != null) {
        const y = yOf(opts.refLine);
        ctx.save();
        ctx.setLineDash([3, 3]);
        ctx.strokeStyle = C_MUTED;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(L, Math.round(y) + 0.5);
        ctx.lineTo(cssW - R, Math.round(y) + 0.5);
        ctx.stroke();
        ctx.restore();
        ctx.fillStyle = C_MUTED;
        ctx.font = "10px Verdana, sans-serif";
        ctx.textAlign = "left";
        ctx.fillText(opts.refLabel || "", L + 4, y - 8);
    }

    // axes
    ctx.strokeStyle = "#9bb49b";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(L + 0.5, T);
    ctx.lineTo(L + 0.5, T + plotH + 0.5);
    ctx.lineTo(cssW - R, T + plotH + 0.5);
    ctx.stroke();

    // gutter divider
    if (opts.showNoHit) {
        ctx.save();
        ctx.setLineDash([2, 3]);
        ctx.strokeStyle = C_RULE;
        ctx.beginPath();
        const xd = Math.round(L + gutter * 0.82) + 0.5;
        ctx.moveTo(xd, T);
        ctx.lineTo(xd, T + plotH);
        ctx.stroke();
        ctx.restore();
    }

    // x labels. Clamp the end labels so they stay inside the canvas.
    ctx.textBaseline = "top";
    if (compact) {
        bands.forEach((b, i) => {
            ctx.save();
            ctx.translate(xOf(i), T + plotH + 7);
            ctx.rotate(-Math.PI / 4);
            ctx.textAlign = "right";
            ctx.textBaseline = "middle";
            ctx.fillStyle = "#0A3200";
            ctx.font = "10px Verdana, sans-serif";
            ctx.fillText(b, 0, 0);
            ctx.restore();
        });
        if (opts.showNoHit) {
            ctx.save();
            ctx.translate(xNoHit, T + plotH + 7);
            ctx.rotate(-Math.PI / 4);
            ctx.textAlign = "right";
            ctx.textBaseline = "middle";
            ctx.fillStyle = "#0A3200";
            ctx.font = "10px Verdana, sans-serif";
            ctx.fillText("none", 0, 0);
            ctx.restore();
        }
    } else {
        bands.forEach((b, i) => {
            const last = i === bands.length - 1;
            ctx.textAlign = last ? "right" : "center";
            const x = last ? cssW - 2 : xOf(i);
            ctx.fillStyle = "#0A3200";
            ctx.font = "10.5px Verdana, sans-serif";
            ctx.fillText(b, x, T + plotH + 8);
            const ns = opts.bandN ? opts.bandN[i] : null;
            if (ns) {
                ctx.fillStyle = C_MUTED;
                ctx.font = "9px Verdana, sans-serif";
                ctx.fillText(ns, x, T + plotH + 21);
            }
        });
        if (opts.showNoHit) {
            ctx.textAlign = "center";
            ctx.fillStyle = "#0A3200";
            ctx.font = "10.5px Verdana, sans-serif";
            ctx.fillText("no", xNoHit, T + plotH + 8);
            ctx.fillText("homolog", xNoHit, T + plotH + 20);
        }
    }

    // axis titles
    ctx.fillStyle = C_MUTED;
    ctx.font = (compact ? 10 : 11) + "px Verdana, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText(
        compact ? (opts.xTitleShort || opts.xTitle || "") : (opts.xTitle || ""),
        L + gutter + innerW / 2,
        T + plotH + (compact ? 56 : 40)
    );
    ctx.save();
    ctx.translate(compact ? 5 : 13, T + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText(opts.yTitle || "", 0, 0);
    ctx.restore();

    // series
    const hit = [];
    opts.series.forEach(s => {
        ctx.strokeStyle = s.color;
        ctx.fillStyle = s.color;
        ctx.lineWidth = 2;
        ctx.save();
        if (s.style === "dashed") ctx.setLineDash([5, 4]);
        else if (s.style === "dotted") ctx.setLineDash([1.5, 3]);

        ctx.beginPath();
        let started = false;
        s.pts.forEach((p, i) => {
            if (p == null) { started = false; return; }
            const x = xOf(i), y = yOf(p.v);
            if (!started) { ctx.moveTo(x, y); started = true; }
            else ctx.lineTo(x, y);
        });
        ctx.stroke();
        ctx.restore();

        s.pts.forEach((p, i) => {
            if (p == null) return;
            const x = xOf(i), y = yOf(p.v);
            ctx.beginPath();
            ctx.arc(x, y, 3.6, 0, Math.PI * 2);
            ctx.fill();
            if (s.style !== "solid") {           // hollow centre
                ctx.save();
                ctx.fillStyle = "#fff";
                ctx.beginPath();
                ctx.arc(x, y, s.style === "dotted" ? 2.2 : 1.7, 0, Math.PI * 2);
                ctx.fill();
                ctx.restore();
            }
            hit.push({ x, y, s, p, band: bands[i] });
        });

        if (s.noHit != null) {
            const y = yOf(s.noHit.v);
            ctx.fillStyle = s.color;
            ctx.fillRect(xNoHit - 4, y - 4, 8, 8);
            ctx.strokeStyle = "#fff";
            ctx.lineWidth = 1;
            ctx.strokeRect(xNoHit - 4, y - 4, 8, 8);
            hit.push({ x: xNoHit, y, s, p: s.noHit, band: "no homolog" });
        }
    });

    cv._hit = hit;
    return hit;
}

function attachHover(cv, readout, fmt) {
    function move(ev) {
        const r = cv.getBoundingClientRect();
        const mx = ev.clientX - r.left, my = ev.clientY - r.top;
        let best = null, bd = 400;
        (cv._hit || []).forEach(h => {
            const d = (h.x - mx) ** 2 + (h.y - my) ** 2;
            if (d < bd) { bd = d; best = h; }
        });
        readout.innerHTML = best ? fmt(best) : readout.dataset.idle;
    }
    cv.addEventListener("mousemove", move);
    cv.addEventListener("mouseleave", () => {
        readout.innerHTML = readout.dataset.idle;
    });
}
