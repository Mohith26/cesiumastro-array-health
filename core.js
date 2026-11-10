/*
 * Array Health core: on-orbit element diagnosis from mutual-coupling telemetry
 * and fault-compensated beam synthesis for a planar active phased array.
 *
 * Unofficial concept demo by Mohith Gajjela. Not affiliated with CesiumAstro.
 *
 * Coordinates: element positions in wavelengths, directions in sine space (u, v).
 * Complex vectors are stored as separate Float64Array re / im pairs.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.ArrayHealth = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const TAU = 2 * Math.PI;
  const DB = (x) => 10 * Math.log10(Math.max(x, 1e-30));

  // ---------- random ----------
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(rng) {
    let u = 0;
    while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * rng());
  }

  // ---------- geometry ----------
  // Triangular lattice, rows x cols elements, four tiles (quadrants).
  function makeArray(opts) {
    const o = Object.assign({ rows: 16, cols: 18, d: 0.56 }, opts || {});
    const dy = (o.d * Math.sqrt(3)) / 2;
    const els = [];
    for (let r = 0; r < o.rows; r++) {
      for (let c = 0; c < o.cols; c++) {
        const x = (c - (o.cols - 1) / 2) * o.d + (r % 2 ? o.d / 4 : -o.d / 4);
        const y = (r - (o.rows - 1) / 2) * dy;
        const tile = 1 + (r >= o.rows / 2 ? 2 : 0) + (c >= o.cols / 2 ? 1 : 0);
        const id = 'R' + String(r + 1).padStart(2, '0') + 'C' + String(c + 1).padStart(2, '0');
        els.push({ i: els.length, r, c, x, y, tile, id });
      }
    }
    return { opts: o, dy, els, N: els.length };
  }

  // Neighbor pairs within radius (in units of lattice spacing d), both directions.
  function neighborPairs(arr, maxR) {
    const d = arr.opts.d, out = [];
    const R = (maxR || 1.8) * d;
    for (const a of arr.els) for (const b of arr.els) {
      if (a.i === b.i) continue;
      const dist = Math.hypot(a.x - b.x, a.y - b.y);
      if (dist <= R + 1e-9) out.push({ tx: a.i, rx: b.i, ring: dist < 1.2 * d ? 1 : 2 });
    }
    return out;
  }

  // ---------- Taylor taper ----------
  function taylorLine(xs, L, sllDb, nbar) {
    const A = Math.acosh(Math.pow(10, sllDb / 20)) / Math.PI;
    const sig2 = (nbar * nbar) / (A * A + (nbar - 0.5) * (nbar - 0.5));
    const F = [];
    for (let m = 1; m < nbar; m++) {
      let num = 1, den = 1;
      for (let n = 1; n < nbar; n++) {
        num *= 1 - (m * m) / (sig2 * (A * A + (n - 0.5) * (n - 0.5)));
        if (n !== m) den *= 1 - (m * m) / (n * n);
      }
      F.push((Math.pow(-1, m + 1) * num) / (2 * den));
    }
    return xs.map((x) => {
      let w = 1;
      for (let m = 1; m < nbar; m++) w += 2 * F[m - 1] * Math.cos((TAU * m * x) / L);
      return w;
    });
  }
  function taylorTaper(arr, sllDb, nbar) {
    const o = arr.opts;
    const Lx = o.cols * o.d + o.d / 2, Ly = o.rows * arr.dy;
    const wx = taylorLine(arr.els.map((e) => e.x), Lx, sllDb, nbar);
    const wy = taylorLine(arr.els.map((e) => e.y), Ly, sllDb, nbar);
    const w = wx.map((v, k) => v * wy[k]);
    const mx = Math.max(...w);
    return Float64Array.from(w.map((v) => v / mx));
  }

  // ---------- pattern evaluation ----------
  // Grid of sample points. Returns flat arrays su, sv.
  function diskGrid(radius, step) {
    const su = [], sv = [], iu = [], iv = [];
    const n = Math.floor(radius / step);
    for (let a = -n; a <= n; a++) for (let b = -n; b <= n; b++) {
      const u = a * step, v = b * step;
      if (u * u + v * v <= radius * radius + 1e-12) { su.push(u); sv.push(v); iu.push(a + n); iv.push(b + n); }
    }
    return { su: Float64Array.from(su), sv: Float64Array.from(sv), iu: Int32Array.from(iu), iv: Int32Array.from(iv), half: n, n: su.length, step };
  }

  // Fast |AF|^2 on a lattice grid from diskGrid: row-factorized,
  // AF(u,v) = sum_rows exp(j2pi y_r v) * sum_{k in row} w_k exp(j2pi x_k u).
  function patternGrid(arr, wr, wi, grid, out) {
    const rows = arr.opts.rows, N = arr.N, nu = 2 * grid.half + 1, step = grid.step;
    const ou = grid.ou || 0, ov = grid.ov || 0;
    const Sr = new Float64Array(rows * nu), Si = new Float64Array(rows * nu);
    for (let k = 0; k < N; k++) {
      if (wr[k] === 0 && wi[k] === 0) continue;
      const e = arr.els[k], base = e.r * nu;
      for (let a = 0; a < nu; a++) {
        const ph = TAU * e.x * (ou + (a - grid.half) * step), c = Math.cos(ph), s = Math.sin(ph);
        Sr[base + a] += wr[k] * c - wi[k] * s; Si[base + a] += wr[k] * s + wi[k] * c;
      }
    }
    const Yr = new Float64Array(rows * nu), Yi = new Float64Array(rows * nu), yrow = new Float64Array(rows);
    for (const e of arr.els) yrow[e.r] = e.y;
    for (let r = 0; r < rows; r++) for (let b = 0; b < nu; b++) {
      const ph = TAU * yrow[r] * (ov + (b - grid.half) * step); Yr[r * nu + b] = Math.cos(ph); Yi[r * nu + b] = Math.sin(ph);
    }
    const res = out || new Float64Array(grid.n);
    for (let m = 0; m < grid.n; m++) {
      const a = grid.iu[m], b = grid.iv[m];
      let ar = 0, ai = 0;
      for (let r = 0; r < rows; r++) {
        const xr = Sr[r * nu + a], xi = Si[r * nu + a], yr = Yr[r * nu + b], yi = Yi[r * nu + b];
        ar += xr * yr - xi * yi; ai += xr * yi + xi * yr;
      }
      res[m] = ar * ar + ai * ai;
    }
    return res;
  }

  // |AF(s)|^2 for complex weights (wr, wi), AF(s) = sum w_i exp(j 2pi (x_i su + y_i sv)).
  function patternPower(arr, wr, wi, su, sv, out) {
    const N = arr.N, M = su.length;
    const res = out || new Float64Array(M);
    const xs = new Float64Array(N), ys = new Float64Array(N);
    let nz = 0; const idx = new Int32Array(N);
    for (let k = 0; k < N; k++) if (wr[k] !== 0 || wi[k] !== 0) idx[nz++] = k;
    for (let k = 0; k < N; k++) { xs[k] = arr.els[k].x; ys[k] = arr.els[k].y; }
    for (let m = 0; m < M; m++) {
      const u = su[m], v = sv[m];
      let ar = 0, ai = 0;
      for (let q = 0; q < nz; q++) {
        const k = idx[q];
        const ph = TAU * (xs[k] * u + ys[k] * v);
        const c = Math.cos(ph), s = Math.sin(ph);
        ar += wr[k] * c - wi[k] * s;
        ai += wr[k] * s + wi[k] * c;
      }
      res[m] = ar * ar + ai * ai;
    }
    return res;
  }

  // Element power pattern cos^q(theta) in sine space.
  function elementPattern(u, v, q) {
    const c2 = 1 - u * u - v * v;
    return c2 <= 0 ? 0 : Math.pow(c2, q / 2);
  }

  // ---------- synthetic scenario ----------
  // Truth: per-element Tx change g (complex) and Rx change h since factory calibration.
  function makeScenario(seed, opts) {
    const o = Object.assign({ snrRing1Db: 40, snrRing2Db: 33, factoryGainDb: 10 }, opts || {});
    const rng = mulberry32(seed);
    const arr = makeArray();
    const N = arr.N;
    const gA = new Float64Array(N).fill(1), gP = new Float64Array(N);
    const hA = new Float64Array(N).fill(1), hP = new Float64Array(N);
    const cause = new Array(N).fill('');
    const dB = (x) => Math.pow(10, x / 20);
    for (let k = 0; k < N; k++) {
      gA[k] = dB(0.2 * gauss(rng)); gP[k] = (2.0 * Math.PI / 180) * gauss(rng);
      hA[k] = dB(0.2 * gauss(rng)); hP[k] = (2.0 * Math.PI / 180) * gauss(rng);
    }
    const at = (r, c) => arr.els.find((e) => e.r === r && e.c === c).i;
    // 1) Thermal gradient across tile 2 (hot corner near the bus): phase walk + gain droop.
    const t2 = arr.els.filter((e) => e.tile === 2);
    const cx = Math.max(...t2.map((e) => e.x)), cy = Math.min(...t2.map((e) => e.y));
    const span = Math.max(...t2.map((e) => Math.hypot(e.x - cx, e.y - cy)));
    for (const e of t2) {
      const f = 1 - Math.hypot(e.x - cx, e.y - cy) / span;
      gP[e.i] += (38 * f * Math.PI) / 180; gA[e.i] *= dB(-1.2 * f);
      if (f > 0.15) cause[e.i] = 'thermal';
    }
    // 2) Two failed 4-channel beamformer ICs (2x2 blocks): Tx and Rx dead.
    for (const [r0, c0] of [[9, 4], [5, 11]]) {
      for (const [dr, dc] of [[0, 0], [0, 1], [1, 0], [1, 1]]) {
        const k = at(r0 + dr, c0 + dc);
        gA[k] = dB(-42 + 3 * gauss(rng)); gP[k] = TAU * rng();
        hA[k] = dB(-42 + 3 * gauss(rng)); hP[k] = TAU * rng();
        cause[k] = 'bfic';
      }
    }
    // 3) Isolated dead PAs (Tx only).
    for (const [r, c] of [[2, 6], [12, 15], [7, 1]]) {
      const k = at(r, c); gA[k] = dB(-40 + 3 * gauss(rng)); gP[k] = TAU * rng(); cause[k] = 'pa';
    }
    // 4) Degraded PAs (low output).
    for (const [r, c, db] of [[3, 13, -3.4], [10, 8, -2.6], [14, 3, -4.8], [6, 16, -2.1]]) {
      const k = at(r, c); gA[k] *= dB(db); cause[k] = 'degraded';
    }
    // 5) One dead LNA (Rx only; Tx healthy).
    { const k = at(13, 10); hA[k] = dB(-41); hP[k] = TAU * rng(); cause[k] = 'lna'; }
    // 6) Phase shifter drift on isolated elements.
    for (const [r, c, deg] of [[1, 2, 24], [11, 13, -31], [4, 7, 19]]) {
      const k = at(r, c); gP[k] += (deg * Math.PI) / 180; cause[k] = 'phase';
    }

    // Telemetry: for each neighbor pair, factory baseline m0 and on-orbit m.
    const pairs = neighborPairs(arr, 1.8);
    const tel = [];
    for (const p of pairs) {
      const a = arr.els[p.tx], b = arr.els[p.rx];
      const dist = Math.hypot(a.x - b.x, a.y - b.y);
      // nominal coupling (unknown to the estimator): depends on separation + edge effects
      const edge = (a.r === 0 || a.r === 15 || a.c === 0 || a.c === 17 ? 1 : 0) + (b.r === 0 || b.r === 15 || b.c === 0 || b.c === 17 ? 1 : 0);
      const cm = dB((p.ring === 1 ? -17 : -24) + 0.6 * edge * gauss(rng));
      const cp = -TAU * dist + 0.4 + 0.1 * edge * gauss(rng) + 0.15 * Math.atan2(b.y - a.y, b.x - a.x);
      const snr = p.ring === 1 ? o.snrRing1Db : o.snrRing2Db;
      const sigma = (dB(p.ring === 1 ? -17 : -24) / dB(snr)) / Math.SQRT2; // per quadrature
      const sigma0 = sigma / dB(o.factoryGainDb);
      // factory: t0 r0 = 1 (calibrated reference state)
      const m0r = cm * Math.cos(cp) + sigma0 * gauss(rng), m0i = cm * Math.sin(cp) + sigma0 * gauss(rng);
      const amp = cm * gA[p.tx] * hA[p.rx], ph = cp + gP[p.tx] + hP[p.rx];
      const mr = amp * Math.cos(ph) + sigma * gauss(rng), mi = amp * Math.sin(ph) + sigma * gauss(rng);
      tel.push({ tx: p.tx, rx: p.rx, ring: p.ring, m0: [m0r, m0i], m: [mr, mi], sigma, sigma0 });
    }
    return {
      seed, arr, telemetry: tel,
      truth: { gA: Array.from(gA), gP: Array.from(gP), hA: Array.from(hA), hP: Array.from(hP), cause },
    };
  }

  // ---------- dense linear algebra ----------
  // Real symmetric positive definite: Cholesky in place (lower), returns L as Float64Array n*n.
  function cholesky(A, n) {
    const L = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      let s = A[j * n + j];
      for (let k = 0; k < j; k++) s -= L[j * n + k] * L[j * n + k];
      if (s <= 0) throw new Error('matrix not positive definite at ' + j);
      const d = Math.sqrt(s); L[j * n + j] = d;
      for (let i = j + 1; i < n; i++) {
        let t = A[i * n + j];
        for (let k = 0; k < j; k++) t -= L[i * n + k] * L[j * n + k];
        L[i * n + j] = t / d;
      }
    }
    return L;
  }
  function cholSolve(L, n, b) {
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k]; y[i] = s / L[i * n + i]; }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k]; x[i] = s / L[i * n + i]; }
    return x;
  }
  // Complex Hermitian PD (re, im n*n): Cholesky A = L L^H.
  function cholC(Ar, Ai, n) {
    const Lr = new Float64Array(n * n), Li = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      let s = Ar[j * n + j];
      for (let k = 0; k < j; k++) s -= Lr[j * n + k] ** 2 + Li[j * n + k] ** 2;
      if (s <= 0) throw new Error('not PD');
      const d = Math.sqrt(s); Lr[j * n + j] = d;
      for (let i = j + 1; i < n; i++) {
        let tr = Ar[i * n + j], ti = Ai[i * n + j];
        for (let k = 0; k < j; k++) {
          // L[i,k] * conj(L[j,k])
          const ar = Lr[i * n + k], ai = Li[i * n + k], br = Lr[j * n + k], bi = -Li[j * n + k];
          tr -= ar * br - ai * bi; ti -= ar * bi + ai * br;
        }
        Lr[i * n + j] = tr / d; Li[i * n + j] = ti / d;
      }
    }
    return { Lr, Li };
  }
  function cholCSolve(C, n, br, bi) {
    const { Lr, Li } = C;
    const yr = new Float64Array(n), yi = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let sr = br[i], si = bi[i];
      for (let k = 0; k < i; k++) { const lr = Lr[i * n + k], li = Li[i * n + k]; sr -= lr * yr[k] - li * yi[k]; si -= lr * yi[k] + li * yr[k]; }
      const d = Lr[i * n + i]; yr[i] = sr / d; yi[i] = si / d;
    }
    const xr = new Float64Array(n), xi = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let sr = yr[i], si = yi[i];
      for (let k = i + 1; k < n; k++) { const lr = Lr[k * n + i], li = -Li[k * n + i]; sr -= lr * xr[k] - li * xi[k]; si -= lr * xi[k] + li * xr[k]; }
      const d = Lr[i * n + i]; xr[i] = sr / d; xi[i] = si / d;
    }
    return { xr, xi };
  }

  // ---------- diagnosis ----------
  // Estimates per-element Tx and Rx complex change since factory calibration from
  // ratios z = m / m0 (coupling cancels). Dead paths are found by robust medians,
  // then log-amplitude and phase are solved by weighted least squares with
  // Huber reweighting; standard errors come from the residual-scaled inverse normal matrix.
  function diagnose(N, telemetry, opts) {
    const o = Object.assign({ deadDb: -15, z: 1.96 }, opts || {});
    const P = telemetry.length;
    const la = new Float64Array(P), ph = new Float64Array(P), vr = new Float64Array(P);
    for (let p = 0; p < P; p++) {
      const t = telemetry[p];
      const [mr, mi] = t.m, [ar, ai] = t.m0;
      const d = ar * ar + ai * ai;
      const zr = (mr * ar + mi * ai) / d, zi = (mi * ar - mr * ai) / d;
      la[p] = 0.5 * Math.log(zr * zr + zi * zi);
      ph[p] = Math.atan2(zi, zr);
      // delta-method variance of log|z| (same for phase): sigma^2/|m|^2 + sigma0^2/|m0|^2 (per quadrature)
      vr[p] = (t.sigma * t.sigma) / Math.max(mr * mr + mi * mi, 1e-30) + (t.sigma0 * t.sigma0) / d;
    }
    const median = (a) => { const s = a.slice().sort((x, y) => x - y); const n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2])) : NaN; };
    const txL = Array.from({ length: N }, () => []), rxL = Array.from({ length: N }, () => []);
    for (let p = 0; p < P; p++) { txL[telemetry[p].tx].push(la[p]); rxL[telemetry[p].rx].push(la[p]); }
    const ln20 = 20 / Math.LN10;
    const txMed = txL.map((a) => median(a) * ln20), rxMed = rxL.map((a) => median(a) * ln20);
    const txDead = txMed.map((v) => v < o.deadDb), rxDead = rxMed.map((v) => v < o.deadDb);

    // index unknowns: alpha_i (Tx) for live Tx, beta_j (Rx) for live Rx
    const txIdx = new Int32Array(N).fill(-1), rxIdx = new Int32Array(N).fill(-1);
    let n = 0;
    for (let k = 0; k < N; k++) if (!txDead[k]) txIdx[k] = n++;
    const nTx = n;
    for (let k = 0; k < N; k++) if (!rxDead[k]) rxIdx[k] = n++;
    const use = [];
    for (let p = 0; p < P; p++) { const t = telemetry[p]; if (!txDead[t.tx] && !rxDead[t.rx]) use.push(p); }

    function solve(yv, wts) {
      const A = new Float64Array(n * n), b = new Float64Array(n);
      for (let q = 0; q < use.length; q++) {
        const p = use[q], i = txIdx[telemetry[p].tx], j = rxIdx[telemetry[p].rx], w = wts[q];
        A[i * n + i] += w; A[j * n + j] += w; A[i * n + j] += w; A[j * n + i] += w;
        b[i] += w * yv[p]; b[j] += w * yv[p];
      }
      // gauge: sum of Rx unknowns = 0 (heavy constraint row)
      const G = 1e4;
      for (let a = nTx; a < n; a++) for (let c = nTx; c < n; c++) A[a * n + c] += G;
      const L = cholesky(A, n);
      return { x: cholSolve(L, n, b), L };
    }
    function irls(yv) {
      let wts = use.map((p) => 1 / vr[p]);
      let sol, res, s2 = 1;
      for (let it = 0; it < 4; it++) {
        sol = solve(yv, wts);
        res = use.map((p) => yv[p] - sol.x[txIdx[telemetry[p].tx]] - sol.x[rxIdx[telemetry[p].rx]]);
        // robust scale of standardized residuals (MAD)
        const st = res.map((r, q) => r / Math.sqrt(vr[use[q]]));
        s2 = Math.pow(1.4826 * median(st.map(Math.abs)), 2);
        const c = 2.5 * Math.sqrt(s2);
        wts = st.map((r, q) => (1 / vr[use[q]]) * (Math.abs(r) <= c ? 1 : c / Math.abs(r)));
      }
      // standard errors: s2 * diag(A^-1)
      const se = new Float64Array(n);
      for (let k = 0; k < n; k++) {
        const e = new Float64Array(n); e[k] = 1;
        se[k] = Math.sqrt(Math.max(s2, 1) * cholSolve(sol.L, n, e)[k]);
      }
      return { x: sol.x, se, s2 };
    }
    // phase: unwrap each pair around the current coarse estimate (two passes)
    const amp = irls(la);
    let phv = Float64Array.from(ph);
    let phs = irls(phv);
    for (let pass = 0; pass < 2; pass++) {
      phv = Float64Array.from(ph);
      for (const p of use) {
        const pred = phs.x[txIdx[telemetry[p].tx]] + phs.x[rxIdx[telemetry[p].rx]];
        phv[p] = ph[p] + TAU * Math.round((pred - ph[p]) / TAU);
      }
      phs = irls(phv);
    }
    // Report Tx relative to the median live Tx element (beam shape is gauge-free).
    const liveTx = []; for (let k = 0; k < N; k++) if (!txDead[k]) liveTx.push(k);
    const refA = median(liveTx.map((k) => amp.x[txIdx[k]]));
    const refP = median(liveTx.map((k) => phs.x[txIdx[k]]));
    const wrap = (a) => Math.atan2(Math.sin(a), Math.cos(a));
    const elems = [];
    for (let k = 0; k < N; k++) {
      const e = { i: k, txDead: txDead[k], rxDead: rxDead[k], txMedDb: txMed[k], rxMedDb: rxMed[k] };
      if (!txDead[k]) {
        const a = txIdx[k];
        e.txGainDb = (amp.x[a] - refA) * ln20; e.txGainSeDb = amp.se[a] * ln20;
        e.txPhaseDeg = (wrap(phs.x[a] - refP) * 180) / Math.PI; e.txPhaseSeDeg = (phs.se[a] * 180) / Math.PI;
      }
      if (!rxDead[k]) {
        const a = rxIdx[k];
        e.rxGainDb = amp.x[a] * ln20; e.rxGainSeDb = amp.se[a] * ln20;
        e.rxPhaseDeg = (wrap(phs.x[a]) * 180) / Math.PI; e.rxPhaseSeDeg = (phs.se[a] * 180) / Math.PI;
      }
      e.status = classify(e, o.z);
      elems.push(e);
    }
    return { elems, pairsUsed: use.length, pairsTotal: P, scaleAmp: Math.sqrt(amp.s2), scalePh: Math.sqrt(phs.s2) };
  }
  function classify(e, z) {
    if (e.txDead) return 'failed';
    const g = e.txGainDb, gs = z * e.txGainSeDb, p = e.txPhaseDeg, ps = z * e.txPhaseSeDeg;
    if (g + gs < -1.0) return 'degraded';
    if (Math.abs(p) - ps > 8) return 'drift';
    if (e.rxDead) return 'rx-fault';
    return 'healthy';
  }

  // ---------- mask ----------
  // Main-lobe exclusion: rectangle at the first nulls of the nominal pattern along u and v.
  function firstNull(arr, w, axis) {
    const N = arr.N, wr = w, wi = new Float64Array(N);
    let prev = Infinity;
    for (let s = 0.002; s < 0.6; s += 0.002) {
      const p = patternPower(arr, wr, wi, Float64Array.of(axis === 'u' ? s : 0), Float64Array.of(axis === 'v' ? s : 0))[0];
      if (p > prev) return s - 0.002;
      prev = p;
    }
    return 0.3;
  }
  // Effective relative mask (linear, relative to beam peak) for every relative offset s,
  // covering every steering direction d with |d| <= sin(thetaMax): the sidelobe at d+s must
  // satisfy P(s) * EP(d+s)/EP(d) <= M, so M_eff(s) = M * min_d EP(d)/EP(d+s).
  function effectiveMask(grid, maskDb, thetaMaxDeg, q, mainU, mainV, dStep) {
    const M = Math.pow(10, maskDb / 10);
    const smax = Math.sin((thetaMaxDeg * Math.PI) / 180);
    const D = diskGrid(smax, dStep || 0.05);
    const out = new Float64Array(grid.n), inMain = new Uint8Array(grid.n);
    for (let m = 0; m < grid.n; m++) {
      const u = grid.su[m], v = grid.sv[m];
      if (Math.abs(u) < mainU && Math.abs(v) < mainV) { inMain[m] = 1; out[m] = Infinity; continue; }
      let best = Infinity, any = false;
      for (let k = 0; k < D.n; k++) {
        const pu = D.su[k] + u, pv = D.sv[k] + v;
        if (pu * pu + pv * pv >= 1) continue;
        any = true;
        const r = elementPattern(D.su[k], D.sv[k], q) / elementPattern(pu, pv, q);
        if (r < best) best = r;
      }
      out[m] = any ? M * best : Infinity;
    }
    return { mask: out, inMain };
  }

  // ---------- compensation synthesis ----------
  // Minimal-perturbation sidelobe correction with adaptive (Olen-Compton style) sample weights:
  //   min ||w - w0||^2 + sum_s alpha_s |AF_s(w)|^2   s.t. sum(w) = sum(w0)
  // then project |w_i| <= limit_i. alpha grows on samples above the mask (minus margin).
  function synthesize(arr, w0, limit, grid, emask, opts) {
    const o = Object.assign({ marginDb: 0.6, iters: 60, alpha0: 0.02 }, opts || {});
    const N = arr.N;
    const live = []; for (let k = 0; k < N; k++) if (limit[k] > 0) live.push(k);
    const n = live.length;
    const Rr = new Float64Array(n * n), Ri = new Float64Array(n * n);
    for (let a = 0; a < n; a++) Rr[a * n + a] = 1;
    const alpha = new Float64Array(grid.n);
    const xs = live.map((k) => arr.els[k].x), ys = live.map((k) => arr.els[k].y);
    const target = new Float64Array(grid.n);
    const mlin = Math.pow(10, -o.marginDb / 10);
    for (let m = 0; m < grid.n; m++) target[m] = emask.mask[m] * mlin;
    const w0r = Float64Array.from(live.map((k) => w0[k])), w0i = new Float64Array(n);
    const s0 = w0r.reduce((s, v) => s + v, 0);
    let wr = Float64Array.from(w0r), wi = new Float64Array(n);
    const full = (vr, vi) => { const fr = new Float64Array(N), fi = new Float64Array(N); live.forEach((k, a) => { fr[k] = vr[a]; fi[k] = vi[a]; }); return [fr, fi]; };
    const history = [];
    let P = null, worst = 0;
    for (let it = 0; it <= o.iters; it++) {
      const [fr, fi] = full(wr, wi);
      let pk = 0; for (let a = 0; a < n; a++) pk += wr[a];
      const peak = pk * pk + wi.reduce((s, v) => s + v, 0) ** 2;
      P = patternGrid(arr, fr, fi, grid, P);
      worst = -Infinity; let nviol = 0;
      const upd = [];
      for (let m = 0; m < grid.n; m++) {
        if (emask.inMain[m] || !isFinite(target[m])) continue;
        const ratio = P[m] / peak / target[m];
        const mdb = DB(ratio);
        if (mdb > worst) worst = mdb;
        if (ratio > 1.02) {
          nviol++;
          const na = alpha[m] === 0 ? o.alpha0 : alpha[m] * Math.min(ratio, 4);
          upd.push([m, na - alpha[m]]); alpha[m] = na;
        }
      }
      history.push({ it, worstMarginDb: worst + o.marginDb, nviol });
      if (nviol === 0 || it === o.iters) break;
      // rank-1 updates R += dA * conj(b) b^T, b_a = exp(j 2pi (x u + y v)) / sqrt(peak) normalization folded into alpha
      for (const [m, da] of upd) {
        const u = grid.su[m], v = grid.sv[m];
        const scale = da / peak * s0 * s0; // keep alpha dimensionless vs. the nominal peak
        const br = new Float64Array(n), bi = new Float64Array(n);
        for (let a = 0; a < n; a++) { const ph = TAU * (xs[a] * u + ys[a] * v); br[a] = Math.cos(ph); bi[a] = Math.sin(ph); }
        for (let a = 0; a < n; a++) {
          const car = br[a], cai = -bi[a];
          for (let c = a; c < n; c++) {
            // conj(b_a) * b_c
            const re = car * br[c] - cai * bi[c], im = car * bi[c] + cai * br[c];
            Rr[a * n + c] += scale * re; Ri[a * n + c] += scale * im;
            if (c !== a) { Rr[c * n + a] += scale * re; Ri[c * n + a] -= scale * im; }
          }
        }
      }
      const C = cholC(Rr, Ri, n);
      const x0 = cholCSolve(C, n, w0r, w0i);
      const ones = cholCSolve(C, n, new Float64Array(n).fill(1), new Float64Array(n));
      // mu so that 1^T w = s0 (complex): w = R^-1 w0 + mu R^-1 1
      let a1r = 0, a1i = 0, b1r = 0, b1i = 0;
      for (let a = 0; a < n; a++) { a1r += x0.xr[a]; a1i += x0.xi[a]; b1r += ones.xr[a]; b1i += ones.xi[a]; }
      const nr = s0 - a1r, ni = -a1i, dd = b1r * b1r + b1i * b1i;
      const mur = (nr * b1r + ni * b1i) / dd, mui = (ni * b1r - nr * b1i) / dd;
      for (let a = 0; a < n; a++) {
        wr[a] = x0.xr[a] + mur * ones.xr[a] - mui * ones.xi[a];
        wi[a] = x0.xi[a] + mur * ones.xi[a] + mui * ones.xr[a];
        const k = live[a], mag = Math.hypot(wr[a], wi[a]);
        if (mag > limit[k]) { wr[a] *= limit[k] / mag; wi[a] *= limit[k] / mag; }
      }
    }
    const [fr, fi] = full(wr, wi);
    return { wr: fr, wi: fi, history, worstMarginDb: worst + o.marginDb, converged: history[history.length - 1].nviol === 0 };
  }

  // ---------- commands, quantization, certification ----------
  function quantize(cr, ci, phaseBits, attStepDb, attMaxDb) {
    const N = cr.length, qr = new Float64Array(N), qi = new Float64Array(N), att = new Float64Array(N), phs = new Int32Array(N);
    const step = TAU / (1 << phaseBits);
    for (let k = 0; k < N; k++) {
      const mag = Math.hypot(cr[k], ci[k]);
      if (mag < 1e-9) { att[k] = -1; phs[k] = 0; continue; } // element off
      let a = -20 * Math.log10(Math.min(mag, 1));
      a = Math.min(attMaxDb, Math.max(0, Math.round(a / attStepDb) * attStepDb)) + 0;
      let s = Math.round(Math.atan2(ci[k], cr[k]) / step); s = ((s % (1 << phaseBits)) + (1 << phaseBits)) % (1 << phaseBits);
      const m = Math.pow(10, -a / 20);
      qr[k] = m * Math.cos(s * step); qi[k] = m * Math.sin(s * step); att[k] = a; phs[k] = s;
    }
    return { qr, qi, att, phs };
  }
  // Commands for steering direction (u0, v0): c = (w / g_est) * exp(-j 2pi (x u0 + y v0)).
  function commands(arr, wr, wi, est, u0, v0) {
    const N = arr.N, cr = new Float64Array(N), ci = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      const g = est[k];
      if (!g || g.mag === 0 || (wr[k] === 0 && wi[k] === 0)) continue;
      // divide by estimated gain
      const dr = (wr[k] * Math.cos(-g.ph) - wi[k] * Math.sin(-g.ph)) / g.mag;
      const di = (wr[k] * Math.sin(-g.ph) + wi[k] * Math.cos(-g.ph)) / g.mag;
      const ph = -TAU * (arr.els[k].x * u0 + arr.els[k].y * v0);
      cr[k] = dr * Math.cos(ph) - di * Math.sin(ph); ci[k] = dr * Math.sin(ph) + di * Math.cos(ph);
    }
    return { cr, ci };
  }
  // Radiated field with TRUE hardware state: f = g_true * c.
  function radiate(truth, qr, qi) {
    const N = qr.length, fr = new Float64Array(N), fi = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      const a = truth.gA[k], p = truth.gP[k];
      fr[k] = a * (qr[k] * Math.cos(p) - qi[k] * Math.sin(p));
      fi[k] = a * (qr[k] * Math.sin(p) + qi[k] * Math.cos(p));
    }
    return { fr, fi };
  }
  // Certify one steering direction: worst sidelobe margin (dB above mask; <=0 passes) over
  // the visible region with the element pattern, plus peak gain (linear field units).
  function certify(arr, fr, fi, u0, v0, rel, maskDb, q, mainU, mainV) {
    // Mask is defined in beam-relative coordinates; sample offsets s on `rel` (a diskGrid of
    // radius 2), shifted by half a step so certification points interleave the synthesis grid.
    const sh = rel.step / 2;
    const g = { iu: rel.iu, iv: rel.iv, half: rel.half, step: rel.step, n: rel.n, ou: u0 + sh, ov: v0 + sh };
    const P = patternGrid(arr, fr, fi, g);
    let peak = 0, worst = -Infinity, wu = 0, wv = 0;
    const ep0 = elementPattern(u0, v0, q);
    const pc = patternPower(arr, fr, fi, Float64Array.of(u0), Float64Array.of(v0))[0] * ep0;
    peak = pc;
    for (let m = 0; m < rel.n; m++) {
      const su = rel.su[m] + sh, sv = rel.sv[m] + sh;
      if (Math.abs(su) < mainU && Math.abs(sv) < mainV) {
        const u = u0 + su, v = v0 + sv;
        if (u * u + v * v < 1) { const p = P[m] * elementPattern(u, v, q); if (p > peak) peak = p; }
      }
    }
    for (let m = 0; m < rel.n; m++) {
      const su = rel.su[m] + sh, sv = rel.sv[m] + sh;
      if (Math.abs(su) < mainU && Math.abs(sv) < mainV) continue;
      const u = u0 + su, v = v0 + sv;
      if (u * u + v * v >= 1) continue;
      const r = DB((P[m] * elementPattern(u, v, q)) / peak) - maskDb;
      if (r > worst) { worst = r; wu = u; wv = v; }
    }
    return { worstDb: worst, at: [wu, wv], peakDb: DB(peak), pass: worst <= 0 };
  }
  function scanDirections(thetas, nPhi) {
    const out = [];
    for (const t of thetas) {
      if (t === 0) { out.push({ theta: 0, phi: 0, u: 0, v: 0 }); continue; }
      for (let k = 0; k < nPhi; k++) {
        const phi = (360 * k) / nPhi, s = Math.sin((t * Math.PI) / 180);
        out.push({ theta: t, phi, u: s * Math.cos((phi * Math.PI) / 180), v: s * Math.sin((phi * Math.PI) / 180) });
      }
    }
    return out;
  }

  // ---------- full pipeline ----------
  function defaultConfig() {
    return {
      taperSllDb: 32, taperNbar: 5, maskDb: -24, thetaMaxDeg: 58.5, elemQ: 1.3,
      phaseBits: 6, attStepDb: 1.0, attMaxDb: 31, synthStep: 0.02, visStep: 0.02, marginDb: 0.6, mainGuard: 1.25,
      thetas: [0, 15, 30, 45, 58.5], nPhi: 12,
    };
  }
  function run(scenario, cfgIn, onProgress) {
    const cfg = Object.assign(defaultConfig(), cfgIn || {});
    const arr = scenario.arr, N = arr.N, truth = scenario.truth;
    const tick = onProgress || (() => {});
    tick('diagnose');
    const dx = diagnose(N, scenario.telemetry);
    const taper = taylorTaper(arr, cfg.taperSllDb, cfg.taperNbar);
    const mainU = cfg.mainGuard * firstNull(arr, taper, 'u'), mainV = cfg.mainGuard * firstNull(arr, taper, 'v');
    // estimated gains for commanding
    const estAll = dx.elems.map((e) => e.txDead ? { mag: 0, ph: 0 } : { mag: Math.pow(10, e.txGainDb / 20), ph: (e.txPhaseDeg * Math.PI) / 180 });
    const factoryEst = dx.elems.map(() => ({ mag: 1, ph: 0 }));
    // re-phase only: correct phase drift and dead-element turn-off but keep the factory taper
    const rephaseEst = dx.elems.map((e) => e.txDead ? { mag: 0, ph: 0 } : { mag: 1, ph: (e.txPhaseDeg * Math.PI) / 180 });
    tick('mask');
    const grid = diskGrid(1 + Math.sin((cfg.thetaMaxDeg * Math.PI) / 180), cfg.synthStep);
    const emask = effectiveMask(grid, cfg.maskDb, cfg.thetaMaxDeg, cfg.elemQ, mainU, mainV);
    tick('synthesize');
    const limit = new Float64Array(N), w0 = new Float64Array(N);
    for (let k = 0; k < N; k++) { limit[k] = estAll[k].mag; w0[k] = Math.min(taper[k], limit[k]); }
    const zeros = new Float64Array(N);
    const dirs = scanDirections(cfg.thetas, cfg.nPhi);
    const vis = diskGrid(2.0, cfg.visStep);
    // Synthesize, then certify the QUANTIZED commands against the ESTIMATED hardware model
    // (never the truth). If any direction fails, raise the design margin and resynthesize.
    const estHw = { gA: estAll.map((g) => g.mag), gP: estAll.map((g) => g.ph) };
    let syn, margin = cfg.marginDb;
    const refine = [];
    for (let round = 0; round < 5; round++) {
      syn = synthesize(arr, w0, limit, grid, emask, { marginDb: margin });
      let fails = 0, worst = -Infinity;
      for (const d of dirs) {
        const c = commands(arr, syn.wr, syn.wi, estAll, d.u, d.v);
        const qz = quantize(c.cr, c.ci, cfg.phaseBits, cfg.attStepDb, cfg.attMaxDb);
        const f = radiate(estHw, qz.qr, qz.qi);
        const r = certify(arr, f.fr, f.fi, d.u, d.v, vis, cfg.maskDb, cfg.elemQ, mainU, mainV);
        if (!r.pass) fails++;
        if (r.worstDb > worst) worst = r.worstDb;
      }
      refine.push({ round, marginDb: margin, failsOnModel: fails, worstOnModelDb: worst });
      if (fails === 0) break;
      margin += Math.max(0.25, worst + 0.1);
    }
    const modes = {
      factory: { w: [taper, zeros], est: factoryEst },
      rephase: { w: [Float64Array.from(taper.map((v, k) => (estAll[k].mag === 0 ? 0 : v))), zeros], est: rephaseEst },
      compensated: { w: [syn.wr, syn.wi], est: estAll },
    };
    // healthy reference: nominal taper on ideal hardware (for EIRP reference)
    const ideal = { gA: new Array(N).fill(1), gP: new Array(N).fill(0) };
    const results = {};
    tick('certify');
    for (const [name, md] of Object.entries(modes)) {
      results[name] = dirs.map((d) => {
        const c = commands(arr, md.w[0], md.w[1], md.est, d.u, d.v);
        const qz = quantize(c.cr, c.ci, cfg.phaseBits, cfg.attStepDb, cfg.attMaxDb);
        const f = radiate(truth, qz.qr, qz.qi);
        return Object.assign({ dir: d }, certify(arr, f.fr, f.fi, d.u, d.v, vis, cfg.maskDb, cfg.elemQ, mainU, mainV));
      });
    }
    const healthy = dirs.map((d) => {
      const c = commands(arr, taper, zeros, factoryEst, d.u, d.v);
      const qz = quantize(c.cr, c.ci, cfg.phaseBits, cfg.attStepDb, cfg.attMaxDb);
      const f = radiate(ideal, qz.qr, qz.qi);
      return Object.assign({ dir: d }, certify(arr, f.fr, f.fi, d.u, d.v, vis, cfg.maskDb, cfg.elemQ, mainU, mainV));
    });
    results.healthy = healthy;
    const summary = {};
    for (const name of Object.keys(results)) {
      const rs = results[name];
      const passN = rs.filter((r) => r.pass).length;
      const eirp = rs.map((r, k) => r.peakDb - healthy[k].peakDb);
      summary[name] = {
        pass: passN, total: rs.length, passFrac: passN / rs.length,
        worstDb: Math.max(...rs.map((r) => r.worstDb)),
        eirpMeanDb: eirp.reduce((s, v) => s + v, 0) / eirp.length,
        eirpMinDb: Math.min(...eirp),
      };
    }
    return { cfg, dx, taper, mainU, mainV, syn, refine, dirs, results, summary, grid, emask };
  }

  // Compare diagnosis with simulation truth (only possible because the data is synthetic).
  function validate(scenario, dx) {
    const tr = scenario.truth, N = scenario.arr.N, z = 1.96;
    const med = (a) => { const t = a.slice().sort((x, y) => x - y); return t[t.length >> 1]; };
    const live = []; for (let k = 0; k < N; k++) if (!dx.elems[k].txDead) live.push(k);
    const tA = med(live.map((k) => 20 * Math.log10(tr.gA[k]))), tP = med(live.map((k) => tr.gP[k]));
    let ca = 0, cp = 0, sa = 0, sp = 0, ma = 0, mp = 0;
    for (const k of live) {
      const e = dx.elems[k];
      const ea = e.txGainDb - (20 * Math.log10(tr.gA[k]) - tA);
      const ep = e.txPhaseDeg - ((Math.atan2(Math.sin(tr.gP[k] - tP), Math.cos(tr.gP[k] - tP)) * 180) / Math.PI);
      sa += ea * ea; sp += ep * ep; ma = Math.max(ma, Math.abs(ea)); mp = Math.max(mp, Math.abs(ep));
      if (Math.abs(ea) <= z * e.txGainSeDb) ca++;
      if (Math.abs(ep) <= z * e.txPhaseSeDeg) cp++;
    }
    const deadTrue = tr.gA.map((a) => 20 * Math.log10(a) < -20);
    let tp = 0, fp = 0, fn = 0;
    for (let k = 0; k < N; k++) { if (dx.elems[k].txDead && deadTrue[k]) tp++; else if (dx.elems[k].txDead) fp++; else if (deadTrue[k]) fn++; }
    const rxDeadTrue = tr.hA.map((a) => 20 * Math.log10(a) < -20);
    let rtp = 0, rfp = 0, rfn = 0;
    for (let k = 0; k < N; k++) { if (dx.elems[k].rxDead && rxDeadTrue[k]) rtp++; else if (dx.elems[k].rxDead) rfp++; else if (rxDeadTrue[k]) rfn++; }
    return {
      n: live.length, rmsGainDb: Math.sqrt(sa / live.length), rmsPhaseDeg: Math.sqrt(sp / live.length),
      maxGainErrDb: ma, maxPhaseErrDeg: mp, coverGain: ca / live.length, coverPhase: cp / live.length,
      txDead: { tp, fp, fn }, rxDead: { tp: rtp, fp: rfp, fn: rfn },
    };
  }

  // Pattern cut through the beam for plotting (relative dB, with element pattern).
  function cut(arr, fr, fi, u0, v0, axis, q, n) {
    const pts = n || 721, su = new Float64Array(pts), sv = new Float64Array(pts), t = [];
    for (let k = 0; k < pts; k++) {
      const s = -1 + (2 * k) / (pts - 1);
      su[k] = axis === 'u' ? s : u0; sv[k] = axis === 'v' ? s : v0; t.push(s);
    }
    const P = patternPower(arr, fr, fi, su, sv);
    let pk = 0; const vals = [];
    for (let k = 0; k < pts; k++) { const e = su[k] ** 2 + sv[k] ** 2 < 1 ? elementPattern(su[k], sv[k], q) : 0; vals.push(P[k] * e); if (P[k] * e > pk) pk = P[k] * e; }
    return { s: t, db: vals.map((v) => (v > 0 ? DB(v / pk) : null)) };
  }

  // Radial cut through the beam center (u0, v0) along unit direction (du, dv), in sine space,
  // normalized to the beam peak (max within |t| < 0.15). Points outside the visible region are null.
  function cutLine(arr, fr, fi, u0, v0, du, dv, q, pts) {
    const n = pts || 641, su = new Float64Array(n), sv = new Float64Array(n), t = [];
    for (let k = 0; k < n; k++) { const tt = -2 + (4 * k) / (n - 1); t.push(tt); su[k] = u0 + tt * du; sv[k] = v0 + tt * dv; }
    const P = patternPower(arr, fr, fi, su, sv);
    const val = [];
    let pk = 0;
    for (let k = 0; k < n; k++) {
      const vis = su[k] * su[k] + sv[k] * sv[k] < 1;
      const v = vis ? P[k] * elementPattern(su[k], sv[k], q) : null;
      val.push(v);
      if (v !== null && Math.abs(t[k]) < 0.15 && v > pk) pk = v;
    }
    return { t, db: val.map((v) => (v === null ? null : DB(v / pk))) };
  }

  return {
    cutLine, mulberry32, gauss, makeArray, patternGrid, neighborPairs, taylorTaper, diskGrid, patternPower, elementPattern,
    makeScenario, cholesky, cholSolve, cholC, cholCSolve, diagnose, classify, firstNull, effectiveMask,
    synthesize, validate, quantize, commands, radiate, certify, scanDirections, defaultConfig, run, cut, DB,
  };
});
