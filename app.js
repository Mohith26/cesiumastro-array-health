/* Array Health UI. All values come from core.js output (data/bundle.js is produced by scripts/build.js). */
(function () {
  'use strict';
  const D = window.AH_DATA, AH = window.ArrayHealth;
  const arr = AH.makeArray(D.telemetry.array);
  const N = arr.N;
  let R = D.run;
  const truth = D.truth;
  const COL = { healthy: '#2b2b2b', drift: '#fdb11d', degraded: '#ff8a3d', failed: '#ef5b4c', 'rx-fault': '#ffffff' };
  const LABEL = { failed: 'Failed Tx path', degraded: 'Degraded output', drift: 'Phase drift', 'rx-fault': 'Rx-only fault', healthy: 'Healthy' };
  const MODES = [['factory', 'Factory weights'], ['rephase', 'Re-phase only'], ['compensated', 'Compensated']];
  const st = { mode: 'compensated', dir: 0, el: null, show: { failed: 1, degraded: 1, drift: 1, 'rx-fault': 1, healthy: 1 } };
  const $ = (id) => document.getElementById(id);
  const fmt = (x, d = 2) => (x === null || x === undefined || !isFinite(x) ? 'n/a' : (x < 0 ? '\u2212' : '') + Math.abs(x).toFixed(d));
  const fmtS = (x, d = 2) => (x > 0 ? '+' : '') + fmt(x, d);
  const tip = $('tip');
  function showTip(e, html) { tip.innerHTML = html; tip.style.display = 'block'; tip.style.left = e.clientX + 14 + 'px'; tip.style.top = e.clientY + 12 + 'px'; }
  function hideTip() { tip.style.display = 'none'; }
  const counts = () => { const c = { failed: 0, degraded: 0, drift: 0, 'rx-fault': 0, healthy: 0 }; R.elems.forEach((e) => c[e.s]++); return c; };

  // ---------- KPIs ----------
  function kpis() {
    const c = counts(), S = R.summary, V = R.validation;
    const cost = S.compensated.eirpMeanDb - S.factory.eirpMeanDb;
    $('kpis').innerHTML = `
      <div class="kpi"><div class="lab">Failed Tx paths</div><div class="val"><span class="bad">${c.failed}</span><small>of ${N} elements</small></div>
        <div class="sub">${c.degraded} degraded, ${c.drift} phase drift, ${c['rx-fault']} Rx-only fault</div></div>
      <div class="kpi"><div class="lab">Sidelobe mask compliance</div><div class="val"><span class="${S.factory.pass === S.factory.total ? 'good' : 'bad'}">${S.factory.pass}</span><span class="arrow">\u2192</span><span class="${S.compensated.pass === S.compensated.total ? 'good' : 'warn'}">${S.compensated.pass}</span><small>/ ${S.compensated.total} directions</small></div>
        <div class="sub">Factory weights vs. compensated, \u03b8 \u2264 ${R.cfg.thetaMaxDeg}\u00b0, ${R.cfg.maskDb} dB mask</div></div>
      <div class="kpi"><div class="lab">EIRP cost of compensation</div><div class="val">${fmt(-cost, 2)}<small>dB</small></div>
        <div class="sub">Mean over scan field, same hardware, quantized commands</div></div>
      <div class="kpi"><div class="lab">Diagnosis error (sim check)</div><div class="val">${fmt(V.rmsGainDb, 2)}<small>dB</small><span class="arrow">/</span>${fmt(V.rmsPhaseDeg, 2)}<small>° RMS</small></div>
        <div class="sub">${V.txDead.tp}/${V.txDead.tp + V.txDead.fn} dead Tx and ${V.rxDead.tp}/${V.rxDead.tp + V.rxDead.fn} dead Rx found, ${V.txDead.fp + V.rxDead.fp} false alarms</div></div>`;
  }

  // ---------- aperture ----------
  const xs = arr.els.map((e) => e.x), ys = arr.els.map((e) => e.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const W = 760, H = 520, padX = 20, padY = 40;
  const sc = Math.min((W - 2 * padX) / (maxX - minX + 0.6), (H - 2 * padY) / (maxY - minY + 0.6));
  const px = (x) => W / 2 + (x - (minX + maxX) / 2) * sc, py = (y) => H / 2 - (y - (minY + maxY) / 2) * sc;
  function hexPath(cx, cy, r) { let d = ''; for (let k = 0; k < 6; k++) { const a = Math.PI / 6 + (k * Math.PI) / 3; d += (k ? 'L' : 'M') + (cx + r * Math.cos(a)).toFixed(1) + ' ' + (cy + r * Math.sin(a)).toFixed(1); } return d + 'Z'; }
  function aperture() {
    const c = counts();
    $('ap-desc').textContent = `${N} elements on a triangular lattice in four tiles. Each element's Tx and Rx change since factory calibration is solved from ${R.pairsUsed.toLocaleString()} of ${R.pairsTotal.toLocaleString()} element-pair coupling measurements.`;
    $('ap-tag').textContent = N + '-element Ka T/R aperture';
    $('filters').innerHTML = ['failed', 'degraded', 'drift', 'rx-fault', 'healthy'].map((k) => `<div class="check ${st.show[k] ? 'on' : ''}" data-k="${k}"><span class="box"></span>${LABEL[k]}<span class="n">${c[k]}</span><span class="sw" style="background:${COL[k]};${k === 'rx-fault' ? 'background:none;border:1.5px solid #fff;' : ''}${k === 'healthy' ? 'border:1px solid #555;' : ''}"></span></div>`).join('');
    document.querySelectorAll('#filters .check').forEach((el) => el.onclick = () => { st.show[el.dataset.k] ^= 1; aperture(); });
    const r = 0.6 * arr.opts.d * sc * 0.93;
    let s = '';
    // tile separators and labels
    const midX = (px(arr.els.find((e) => e.r === 0 && e.c === arr.opts.cols / 2 - 1).x) + px(arr.els.find((e) => e.r === 0 && e.c === arr.opts.cols / 2).x)) / 2;
    const midY = (py(arr.els.find((e) => e.r === arr.opts.rows / 2 - 1).y) + py(arr.els.find((e) => e.r === arr.opts.rows / 2).y)) / 2;
    const top = py(maxY) - r - 6, bot = py(minY) + r + 6, lef = px(minX) - r - 6, rig = px(maxX) + r + 6;
    s += `<line x1="${midX}" x2="${midX}" y1="${top}" y2="${bot}" stroke="#3a3a3a" stroke-dasharray="2 3"/><line x1="${lef}" x2="${rig}" y1="${midY}" y2="${midY}" stroke="#3a3a3a" stroke-dasharray="2 3"/>`;
    const lab = (x, y, t, anchor) => `<text x="${x}" y="${y}" text-anchor="${anchor}" fill="#595959" font-family="DM Mono" font-size="11" letter-spacing="1">TILE ${t}</text>`;
    s += lab(lef, top - 8, 3, 'start') + lab(rig, top - 8, 4, 'end') + lab(lef, bot + 18, 1, 'start') + lab(rig, bot + 18, 2, 'end');
    const maxDrift = Math.max(...R.elems.filter((e) => e.s === 'drift').map((e) => Math.abs(e.p)), 1);
    arr.els.forEach((e, k) => {
      const d = R.elems[k];
      const on = st.show[d.s];
      let fill = COL[d.s], stroke = 'none', op = on ? 1 : 0.12, extra = '';
      if (d.s === 'drift') op = on ? 0.35 + 0.65 * Math.abs(d.p) / maxDrift : 0.1;
      if (d.s === 'rx-fault') { fill = '#2b2b2b'; stroke = '#fff'; }
      if (d.s === 'failed') extra = `<path d="M${px(e.x) - r * 0.42} ${py(e.y) - r * 0.42}L${px(e.x) + r * 0.42} ${py(e.y) + r * 0.42}M${px(e.x) + r * 0.42} ${py(e.y) - r * 0.42}L${px(e.x) - r * 0.42} ${py(e.y) + r * 0.42}" stroke="#000" stroke-width="1.6" opacity="${op}"/>`;
      const sel = st.el === k ? `<path d="${hexPath(px(e.x), py(e.y), r + 3.5)}" fill="none" stroke="#fff" stroke-width="1.5"/>` : '';
      s += `<g data-k="${k}" style="cursor:pointer"><path d="${hexPath(px(e.x), py(e.y), r)}" fill="${fill}" fill-opacity="${op}" stroke="${stroke}" stroke-width="1.4" stroke-opacity="${op}"/>${extra}${sel}</g>`;
    });
    $('ap').innerHTML = s;
    $('ap').querySelectorAll('g[data-k]').forEach((g) => {
      const k = +g.dataset.k, e = arr.els[k], d = R.elems[k];
      g.onclick = () => { st.el = k; aperture(); detail(); weights(); };
      g.onmousemove = (ev) => showTip(ev, `<b>${e.id}</b> &middot; ${LABEL[d.s]}<br>${d.td ? 'Tx dead (' + fmt(d.tm, 1) + ' dB vs. factory)' : 'Tx ' + fmtS(d.g, 2) + ' dB, ' + fmtS(d.p, 1) + '\u00b0'}`);
      g.onmouseleave = hideTip;
    });
    $('legend').innerHTML = ['failed', 'degraded', 'drift', 'rx-fault', 'healthy'].map((k) => `<span><i style="background:${k === 'rx-fault' ? 'transparent;border:1.5px solid #fff' : COL[k]}"></i>${LABEL[k]}</span>`).join('') + '<span style="margin-left:auto;color:#737373">Phase-drift shading scales with |\u0394\u03c6|</span>';
  }

  // ---------- commands per mode ----------
  const zeros = new Float64Array(N);
  function modeWeights(mode) {
    const taper = Float64Array.from(R.taper);
    if (mode === 'factory') return { w: [taper, zeros], est: R.elems.map(() => ({ mag: 1, ph: 0 })) };
    if (mode === 'rephase') return { w: [Float64Array.from(taper.map((v, k) => (R.elems[k].td ? 0 : v))), zeros], est: R.elems.map((e) => (e.td ? { mag: 0, ph: 0 } : { mag: 1, ph: (e.p * Math.PI) / 180 })) };
    return { w: [Float64Array.from(R.syn.wr), Float64Array.from(R.syn.wi)], est: R.elems.map((e) => (e.td ? { mag: 0, ph: 0 } : { mag: Math.pow(10, e.g / 20), ph: (e.p * Math.PI) / 180 })) };
  }
  function cmdFor(mode, d) {
    const m = modeWeights(mode);
    const c = AH.commands(arr, m.w[0], m.w[1], m.est, d.u, d.v);
    return AH.quantize(c.cr, c.ci, R.cfg.phaseBits, R.cfg.attStepDb, R.cfg.attMaxDb);
  }

  // ---------- detail ----------
  function detail() {
    if (st.el === null) return;
    const k = st.el, e = arr.els[k], d = R.elems[k], dir = R.results.compensated[st.dir];
    const qf = cmdFor('factory', dir), qc = cmdFor('compensated', dir);
    const tA = 20 * Math.log10(truth.gA[k]), tP = (Math.atan2(Math.sin(truth.gP[k]), Math.cos(truth.gP[k])) * 180) / Math.PI;
    const action = d.td ? 'Switched off; neighbors re-weighted' : d.s === 'degraded' ? 'Gain limit lowered; taper re-solved' : d.s === 'drift' ? 'Phase offset corrected' : d.s === 'rx-fault' ? 'No Tx action (Rx path only)' : 'None';
    const cmd = (q) => (q.att[k] < 0 ? 'off' : `${q.att[k]} dB / ${q.phs[k]}`);
    $('detail').innerHTML = `
      <h3>${e.id}</h3><span class="tag">Tile ${e.tile} &middot; row ${e.r + 1} &middot; column ${e.c + 1}</span>
      <div style="margin-bottom:18px"><span class="pill" style="color:${d.s === 'healthy' ? '#8a8a8a' : d.s === 'rx-fault' ? '#fff' : COL[d.s]}">${LABEL[d.s]}</span></div>
      <div class="spec">
        <div class="row"><span class="k">Tx gain change</span><span class="v">${d.td ? '<span class="bad">dead</span> (' + fmt(d.tm, 1) + ' dB)' : fmtS(d.g, 2) + ' \u00b1 ' + fmt(1.96 * d.gse, 2) + ' dB'}</span></div>
        <div class="row"><span class="k">Tx phase change</span><span class="v">${d.td ? 'n/a' : fmtS(d.p, 1) + ' \u00b1 ' + fmt(1.96 * d.pse, 1) + '\u00b0'}</span></div>
        <div class="row"><span class="k">Rx gain change</span><span class="v">${d.rd ? '<span class="bad">dead</span> (' + fmt(d.rm, 1) + ' dB)' : fmtS(d.rg, 2) + ' dB'}</span></div>
        <div class="row"><span class="k">Factory command</span><span class="v">${cmd(qf)}</span></div>
        <div class="row"><span class="k">New command</span><span class="v">${cmd(qc)}</span></div>
        <div class="row"><span class="k">Action</span><span class="v">${action}</span></div>
        <div class="row"><span class="k" style="color:#737373">Sim truth (Tx)</span><span class="v" style="color:#737373">${tA < -20 ? 'dead (' + fmt(tA, 1) + ' dB)' : fmtS(tA, 2) + ' dB, ' + fmtS(tP, 1) + '\u00b0'}</span></div>
      </div>
      <p style="font-size:13px;color:#737373;margin:12px 0 0">\u00b1 is a 95% interval. Commands are attenuator dB / 6-bit phase state. Truth is shown only because the data is simulated; it is not used by the estimator.</p>`;
  }

  // ---------- certification ----------
  function tabs() {
    $('tabs').innerHTML = MODES.map(([k, l]) => `<button class="tab ${st.mode === k ? 'on' : ''}" data-k="${k}">${l}<span class="c">${R.summary[k].pass}/${R.summary[k].total}</span></button>`).join('');
    document.querySelectorAll('#tabs .tab').forEach((b) => b.onclick = () => { st.mode = b.dataset.k; tabs(); polar(); cut(); });
  }
  function certSide() {
    const S = R.summary, last = R.refine[R.refine.length - 1];
    $('cert-desc').textContent = `Each steering direction is checked with quantized commands radiated by the simulated true hardware, including the element pattern. Re-phasing fixes drift but cannot recover the mask once elements are dead.`;
    $('cert-spec').innerHTML = [
      ['Mask', `${R.cfg.maskDb} dB rel. peak`],
      ['Main-lobe guard', `${R.cfg.mainGuard}\u00d7 first null`],
      ['Scan field', `${S.compensated.total} dirs, \u03b8 \u2264 ${R.cfg.thetaMaxDeg}\u00b0`],
      ['Element pattern', `cos^${R.cfg.elemQ} \u03b8`],
      ['Quantization', `${R.cfg.phaseBits}-bit, ${R.cfg.attStepDb} dB steps`],
      ['Design margin', `${fmt(last.marginDb, 2)} dB, ${R.refine.length} round${R.refine.length > 1 ? 's' : ''}`],
      ['Worst margin', `<span class="${S.compensated.worstDb <= 0 ? 'good' : 'bad'}">${fmtS(S.compensated.worstDb, 2)} dB</span>`],
    ].map(([k, v]) => `<div class="row"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
  }
  function polar() {
    const rs = R.results[st.mode], C = 170, Rr = 150, smax = Math.sin((R.cfg.thetaMaxDeg * Math.PI) / 180);
    let s = '';
    for (const t of R.cfg.thetas.filter((t) => t > 0)) {
      const rr = (Rr * Math.sin((t * Math.PI) / 180)) / smax;
      s += `<circle cx="${C}" cy="${C}" r="${rr}" fill="none" stroke="#262626"/><text x="${C + 4}" y="${C - rr - 4}" fill="#595959" font-family="DM Mono" font-size="10">${t}\u00b0</text>`;
    }
    for (let a = 0; a < 360; a += 30) { const ar = (a * Math.PI) / 180; s += `<line x1="${C}" y1="${C}" x2="${C + Rr * Math.cos(ar)}" y2="${C - Rr * Math.sin(ar)}" stroke="#1c1c1c"/>`; }
    
    rs.forEach((r, i) => {
      const x = C + (Rr * r.u) / smax, y = C - (Rr * r.v) / smax;
      const col = r.ok ? '#39b255' : '#ef5b4c';
      const rad = 5 + Math.min(4, Math.abs(r.w));
      s += `<g data-i="${i}" style="cursor:pointer">${i === st.dir ? `<circle cx="${x}" cy="${y}" r="${rad + 4}" fill="none" stroke="#fff" stroke-width="1.5"/>` : ''}<circle cx="${x}" cy="${y}" r="${rad}" fill="${col}" fill-opacity="${r.ok ? 0.9 : 0.85}"/></g>`;
    });
    $('polar').innerHTML = s;
    $('polar-tag').innerHTML = `<span class="${R.summary[st.mode].pass === R.summary[st.mode].total ? 'good' : 'bad'}">${R.summary[st.mode].pass}/${R.summary[st.mode].total} pass</span>`;
    $('polar').querySelectorAll('g[data-i]').forEach((g) => {
      const i = +g.dataset.i, r = rs[i];
      g.onclick = () => { st.dir = i; polar(); cut(); detail(); weights(); };
      g.onmousemove = (ev) => showTip(ev, `<b>\u03b8 ${r.t}\u00b0, \u03c6 ${r.f}\u00b0</b><br>Worst sidelobe ${fmtS(r.w, 2)} dB vs. mask`);
      g.onmouseleave = hideTip;
    });
  }
  function cut() {
    const d = R.results[st.mode][st.dir], fac = R.results.factory[st.dir];
    const q = R.cfg.elemQ;
    // radial cut through the beam and the factory pattern's worst sidelobe
    let du = fac.at[0] - d.u, dv = fac.at[1] - d.v; const nn = Math.hypot(du, dv) || 1; du /= nn; dv /= nn;
    const series = [];
    for (const mode of st.mode === 'factory' ? ['factory'] : ['factory', st.mode]) {
      const qz = cmdFor(mode, d);
      const f = AH.radiate(truth, qz.qr, qz.qi);
      series.push({ mode, c: AH.cutLine(arr, f.fr, f.fi, d.u, d.v, du, dv, q, 801) });
    }
    const ts = series[0].c.t.filter((t, k) => series[0].c.db[k] !== null);
    const tmin = Math.min(...ts), tmax = Math.max(...ts);
    const Wc = 640, Hc = 300, l = 44, r = 10, t = 14, b = 30, lo = -50;
    const X = (x) => l + ((x - tmin) / (tmax - tmin)) * (Wc - l - r), Y = (db) => t + (Math.min(0, Math.max(lo, db)) / lo) * (Hc - t - b);
    let s = '';
    for (let db = 0; db >= lo; db -= 10) s += `<line x1="${l}" x2="${Wc - r}" y1="${Y(db)}" y2="${Y(db)}" stroke="#1f1f1f"/><text x="${l - 8}" y="${Y(db) + 4}" text-anchor="end" fill="#595959" font-family="DM Mono" font-size="10">${db}</text>`;
    const step = (tmax - tmin) > 1.2 ? 0.5 : 0.25;
    for (let x = Math.ceil(tmin / step) * step; x <= tmax + 1e-9; x += step) s += `<text x="${X(x)}" y="${Hc - 10}" text-anchor="middle" fill="#595959" font-family="DM Mono" font-size="10">${x.toFixed(2)}</text>`;
    // main-lobe guard (approximate along the cut)
    const gl = Math.hypot(R.mainU * du, R.mainV * dv) > 0 ? Math.min(R.mainU / Math.max(Math.abs(du), 1e-9), R.mainV / Math.max(Math.abs(dv), 1e-9)) : R.mainU;
    s += `<rect x="${X(-gl)}" y="${t}" width="${X(gl) - X(-gl)}" height="${Hc - t - b}" fill="#ffffff" fill-opacity="0.035"/>`;
    s += `<line x1="${X(tmin)}" x2="${X(-gl)}" y1="${Y(R.cfg.maskDb)}" y2="${Y(R.cfg.maskDb)}" stroke="#fdb11d" stroke-dasharray="5 4" stroke-width="1.3"/><line x1="${X(gl)}" x2="${X(tmax)}" y1="${Y(R.cfg.maskDb)}" y2="${Y(R.cfg.maskDb)}" stroke="#fdb11d" stroke-dasharray="5 4" stroke-width="1.3"/>`;
    series.forEach((sr, j) => {
      const isTop = j === series.length - 1;
      let p = '', pen = false;
      sr.c.t.forEach((x, k) => { const v = sr.c.db[k]; if (v === null) { pen = false; return; } p += (pen ? 'L' : 'M') + X(x).toFixed(1) + ' ' + Y(v).toFixed(1); pen = true; });
      s += `<path d="${p}" fill="none" stroke="${isTop && st.mode !== 'factory' ? '#ffffff' : '#ef5b4c'}" stroke-opacity="${isTop ? 1 : 0.6}" stroke-width="${isTop ? 1.6 : 1.2}"/>`;
    });
    s += `<text x="${Wc - r}" y="${Y(R.cfg.maskDb) - 6}" text-anchor="end" fill="#fdb11d" font-family="DM Mono" font-size="10">MASK ${R.cfg.maskDb} dB</text>`;
    s += `<text x="${l + 4}" y="${Hc - 10}" fill="#3d3d3d" font-family="DM Mono" font-size="10"></text>`;
    $('cut').innerHTML = s;
    $('cut-tag').textContent = `Cut through worst factory sidelobe \u00b7 \u03b8 ${d.t}\u00b0, \u03c6 ${d.f}\u00b0`;
    $('cut-axis').textContent = `dB rel. peak`;
    const h = R.results.healthy[st.dir];
    $('cut-meta').innerHTML = `
      <div><b class="${d.ok ? 'good' : 'bad'}">${fmtS(d.w, 2)} dB</b>worst vs. mask (2D)</div>
      <div><b class="${fac.ok ? 'good' : 'bad'}">${fmtS(fac.w, 2)} dB</b>factory, same direction</div>
      <div><b>${fmt(d.pk - h.pk, 2)} dB</b>peak vs. healthy</div>
      <div style="margin-left:auto;display:flex;gap:16px;align-items:center"><span><i style="display:inline-block;width:16px;height:2px;background:#ef5b4c;vertical-align:middle;margin-right:6px;opacity:.7"></i>Factory</span>${st.mode !== 'factory' ? `<span><i style="display:inline-block;width:16px;height:2px;background:#fff;vertical-align:middle;margin-right:6px"></i>${st.mode === 'rephase' ? 'Re-phase only' : 'Compensated'}</span>` : ''}</div>`;
  }

  // ---------- weight table ----------
  function weights() {
    const d = R.results.compensated[st.dir];
    const qf = cmdFor('factory', d), qc = cmdFor('compensated', d);
    $('wt-desc').textContent = `Commands for \u03b8 ${d.t}\u00b0, \u03c6 ${d.f}\u00b0. The ${R.elems.filter((e) => e.td).length} failed elements are switched off; rows are the 12 live elements whose commands changed most. Export writes all ${N}.`;
    const rows = arr.els.map((e, k) => {
      const af = qf.att[k] < 0 ? 99 : qf.att[k], ac = qc.att[k] < 0 ? 99 : qc.att[k];
      let dp = ((qc.phs[k] - qf.phs[k] + 96) % 64) - 32;
      return { k, e, s: R.elems[k].s, af, ac, pf: qf.phs[k], pc: qc.phs[k], dp, score: R.elems[k].td ? -1 : Math.abs(ac - af) + Math.abs(dp) * 0.25 };
    }).sort((a, b) => b.score - a.score).slice(0, 12);
    if (st.el !== null && !rows.some((r) => r.k === st.el)) rows[rows.length - 1] = (() => { const k = st.el, af = qf.att[k] < 0 ? 99 : qf.att[k], ac = qc.att[k] < 0 ? 99 : qc.att[k]; return { k, e: arr.els[k], s: R.elems[k].s, af, ac, pf: qf.phs[k], pc: qc.phs[k], dp: ((qc.phs[k] - qf.phs[k] + 96) % 64) - 32 }; })();
    const att = (a) => (a === 99 ? 'off' : a + ' dB');
    $('wt').innerHTML = `<thead><tr><th>Element</th><th>Tile</th><th>Status</th><th class="num">Factory att.</th><th class="num">Factory phase</th><th class="num">New att.</th><th class="num">New phase</th><th class="num">\u0394 phase</th></tr></thead><tbody>` +
      rows.map((r) => `<tr data-k="${r.k}" class="${st.el === r.k ? 'sel' : ''}"><td style="color:#fff">${r.e.id}</td><td>${r.e.tile}</td><td><span style="color:${r.s === 'healthy' ? '#8a8a8a' : r.s === 'rx-fault' ? '#fff' : COL[r.s]}">${LABEL[r.s]}</span></td><td class="num">${att(r.af)}</td><td class="num">${r.pf}</td><td class="num">${att(r.ac)}</td><td class="num">${r.ac === 99 ? '\u2013' : r.pc}</td><td class="num">${r.ac === 99 ? '\u2013' : (r.dp > 0 ? '+' : '') + (r.dp * 5.625).toFixed(1) + '\u00b0'}</td></tr>`).join('') + '</tbody>';
    $('wt').querySelectorAll('tbody tr').forEach((tr) => tr.onclick = () => { st.el = +tr.dataset.k; aperture(); detail(); weights(); });
  }
  function csv() {
    const d = R.results.compensated[st.dir], qc = cmdFor('compensated', d);
    const lines = ['element,tile,row,col,status,attenuation_db,phase_state,phase_deg'];
    arr.els.forEach((e, k) => lines.push([e.id, e.tile, e.r + 1, e.c + 1, R.elems[k].s, qc.att[k] < 0 ? 'off' : qc.att[k], qc.att[k] < 0 ? '' : qc.phs[k], qc.att[k] < 0 ? '' : (qc.phs[k] * 5.625).toFixed(3)].join(',')));
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([lines.join('\n')], { type: 'text/csv' }));
    a.download = `array-health-weights-theta${d.t}-phi${d.f}.csv`; a.click();
  }

  // ---------- method ----------
  function method(extra) {
    const c = counts();
    $('method').innerHTML = [
      ['Telemetry', `${R.pairsTotal.toLocaleString()} element-pair coupling measurements (two neighbor rings), each divided by its factory baseline so the coupling itself cancels.`],
      ['Dead-path test', `Robust median of each element's ratios as transmitter and as receiver; below ${-15} dB flags the path dead. ${c.failed} Tx, ${R.elems.filter((e) => e.rd).length} Rx.`],
      ['Estimator', `Weighted least squares on log-amplitude and unwrapped phase, ${R.pairsUsed.toLocaleString()} pairs, Huber reweighting, 95% intervals from the residual-scaled normal matrix.`],
      ['Synthesis', `Minimal change to the ${R.cfg.taperSllDb} dB Taylor taper with adaptive sidelobe weights (${R.syn.iters} iterations). One taper covers every steering direction through a scan-aware effective mask.`],
      ['Certification', `Quantized commands are first checked on the estimated hardware model, raising the design margin until all pass (${R.refine.map((x) => x.failsOnModel).join(' \u2192 ')} failing). Results above use the simulated true hardware.`],
      ['Calibration check', `Measured 95% interval coverage: ${(100 * R.validation.coverGain).toFixed(0)}% gain, ${(100 * R.validation.coverPhase).toFixed(0)}% phase.`],
      ['Runtime', extra || `${(R.ms / 1000).toFixed(1)} s in Node for the full pipeline (scripts/build.js).`],
    ].map(([k, v]) => `<div class="row"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
  }

  function renderAll(extra) {
    kpis(); aperture(); detail(); tabs(); certSide(); polar(); cut(); weights(); method(extra);
    $('foot-run').textContent = `Seed ${R.seed} \u00b7 built ${R.generated}`;
  }

  // default selection: first failed element in a beamformer-IC block, worst compensated direction at max scan
  st.el = R.elems.findIndex((e) => e.s === 'failed');
  st.dir = R.results.factory.reduce((bi, r, i, a) => (r.w > a[bi].w ? i : bi), 0);
  renderAll();

  $('btn-csv').onclick = csv;
  $('btn-stage').onclick = (ev) => { const b = ev.currentTarget; b.firstChild.textContent = 'Staged for next pass (demo) '; };
  $('btn-rerun').onclick = () => {
    const b = $('btn-rerun'); b.textContent = 'Running\u2026'; document.body.classList.add('busy');
    setTimeout(() => {
      const t0 = performance.now();
      const scn = { seed: D.telemetry.seed, arr, telemetry: D.telemetry.pairs, truth };
      const res = AH.run(scn), val = AH.validate(scn, res.dx), ms = performance.now() - t0;
      const same = res.summary.compensated.pass === R.summary.compensated.pass && Math.abs(res.summary.compensated.eirpMeanDb - R.summary.compensated.eirpMeanDb) < 1e-6;
      R = Object.assign({}, R, {
        elems: res.dx.elems.map((e) => ({ s: e.status, td: e.txDead, rd: e.rxDead, g: e.txGainDb, gse: e.txGainSeDb, p: e.txPhaseDeg, pse: e.txPhaseSeDeg, rg: e.rxGainDb, rp: e.rxPhaseDeg, tm: e.txMedDb, rm: e.rxMedDb })),
        syn: { wr: Array.from(res.syn.wr), wi: Array.from(res.syn.wi), iters: res.syn.history.length }, refine: res.refine,
        results: Object.fromEntries(Object.entries(res.results).map(([k, v]) => [k, v.map((x) => ({ t: x.dir.theta, f: x.dir.phi, u: x.dir.u, v: x.dir.v, w: x.worstDb, pk: x.peakDb, ok: x.pass, at: x.at }))])),
        summary: res.summary, validation: val,
      });
      document.body.classList.remove('busy'); b.textContent = 'Re-run in browser';
      renderAll(`${(ms / 1000).toFixed(1)} s in this browser just now; ${same ? 'matches' : 'differs from'} the committed build.`);
    }, 30);
  };
})();
// Self-check hook: index.html#rerun recomputes everything in the browser on load.
if (location.hash === '#rerun') setTimeout(() => document.getElementById('btn-rerun').click(), 50);
