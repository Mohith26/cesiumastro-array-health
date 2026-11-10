const test = require('node:test');
const assert = require('node:assert');
const AH = require('../core.js');

test('real Cholesky solves an SPD system', () => {
  const rng = AH.mulberry32(1), n = 12, B = Array.from({ length: n * n }, () => rng() - 0.5);
  const A = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { let s = i === j ? 1 : 0; for (let k = 0; k < n; k++) s += B[i * n + k] * B[j * n + k]; A[i * n + j] = s; }
  const x = Float64Array.from({ length: n }, () => rng()), b = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) b[i] += A[i * n + j] * x[j];
  const y = AH.cholSolve(AH.cholesky(A, n), n, b);
  for (let i = 0; i < n; i++) assert.ok(Math.abs(y[i] - x[i]) < 1e-9);
});

test('complex Cholesky solves a Hermitian PD system', () => {
  const rng = AH.mulberry32(2), n = 10;
  const Br = Array.from({ length: n * n }, () => rng() - 0.5), Bi = Array.from({ length: n * n }, () => rng() - 0.5);
  const Ar = new Float64Array(n * n), Ai = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let sr = i === j ? 1 : 0, si = 0;
    for (let k = 0; k < n; k++) { // B[i,k] * conj(B[j,k])
      sr += Br[i * n + k] * Br[j * n + k] + Bi[i * n + k] * Bi[j * n + k];
      si += Bi[i * n + k] * Br[j * n + k] - Br[i * n + k] * Bi[j * n + k];
    }
    Ar[i * n + j] = sr; Ai[i * n + j] = si;
  }
  const xr = Float64Array.from({ length: n }, () => rng()), xi = Float64Array.from({ length: n }, () => rng());
  const br = new Float64Array(n), bi = new Float64Array(n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { br[i] += Ar[i * n + j] * xr[j] - Ai[i * n + j] * xi[j]; bi[i] += Ar[i * n + j] * xi[j] + Ai[i * n + j] * xr[j]; }
  const s = AH.cholCSolve(AH.cholC(Ar, Ai, n), n, br, bi);
  for (let i = 0; i < n; i++) { assert.ok(Math.abs(s.xr[i] - xr[i]) < 1e-9); assert.ok(Math.abs(s.xi[i] - xi[i]) < 1e-9); }
});

test('fast row-factorized grid pattern matches direct summation', () => {
  const arr = AH.makeArray(), rng = AH.mulberry32(3);
  const wr = Float64Array.from({ length: arr.N }, () => rng()), wi = Float64Array.from({ length: arr.N }, () => rng() - 0.5);
  const g = AH.diskGrid(0.5, 0.05); g.ou = 0.013; g.ov = -0.021;
  const fast = AH.patternGrid(arr, wr, wi, g);
  const su = Float64Array.from(g.su, (u) => u + 0.013), sv = Float64Array.from(g.sv, (v) => v - 0.021);
  const slow = AH.patternPower(arr, wr, wi, su, sv);
  for (let m = 0; m < g.n; m++) assert.ok(Math.abs(fast[m] - slow[m]) < 1e-6 * Math.max(1, slow[m]));
});

test('uniform steered array peaks at N^2 in the steering direction', () => {
  const arr = AH.makeArray(), N = arr.N;
  const est = Array.from({ length: N }, () => ({ mag: 1, ph: 0 }));
  const c = AH.commands(arr, new Float64Array(N).fill(1), new Float64Array(N), est, 0.31, -0.22);
  const p = AH.patternPower(arr, c.cr, c.ci, Float64Array.of(0.31), Float64Array.of(-0.22))[0];
  assert.ok(Math.abs(p - N * N) < 1e-6);
});

test('quantizer uses 6-bit phase states and whole-dB attenuator steps', () => {
  const q = AH.quantize(Float64Array.of(0.5, 0, 1), Float64Array.of(0.5, 0, 0), 6, 1, 31);
  assert.ok(q.phs[0] >= 0 && q.phs[0] < 64);
  assert.strictEqual(q.phs[0], 8); // 45 deg / 5.625 deg
  assert.strictEqual(q.att[0], 3); // |0.707| -> 3.01 dB -> 3 dB
  assert.strictEqual(q.att[1], -1); // off element
  assert.strictEqual(q.att[2], 0);
});

test('diagnosis finds every dead path and recovers drift to well under the injected size', () => {
  for (const seed of [20261105, 7]) {
    const sc = AH.makeScenario(seed);
    const dx = AH.diagnose(sc.arr.N, sc.telemetry);
    const v = AH.validate(sc, dx);
    assert.deepStrictEqual(v.txDead, { tp: 11, fp: 0, fn: 0 });
    assert.deepStrictEqual(v.rxDead, { tp: 9, fp: 0, fn: 0 });
    assert.ok(v.rmsGainDb < 0.1, 'rms gain ' + v.rmsGainDb);
    assert.ok(v.rmsPhaseDeg < 0.5, 'rms phase ' + v.rmsPhaseDeg);
    assert.ok(v.coverGain > 0.8 && v.coverPhase > 0.85, 'coverage ' + v.coverGain + ' ' + v.coverPhase);
    // the Rx-only fault is not confused with a Tx fault
    const lna = sc.truth.cause.indexOf('lna');
    assert.ok(!dx.elems[lna].txDead && dx.elems[lna].rxDead);
    // degraded PAs are flagged degraded
    sc.truth.cause.forEach((c, k) => { if (c === 'degraded') assert.strictEqual(dx.elems[k].status, 'degraded'); });
  }
});

test('effective mask is tighter than the raw mask wherever scanning can bring a sidelobe toward broadside', () => {
  const g = AH.diskGrid(1.85, 0.1);
  const em = AH.effectiveMask(g, -24, 58.5, 1.3, 0.2, 0.25);
  const M = Math.pow(10, -24 / 10);
  let tightest = Infinity;
  for (let m = 0; m < g.n; m++) if (!em.inMain[m] && isFinite(em.mask[m]) && Math.hypot(g.su[m], g.sv[m]) <= 1.2) { assert.ok(em.mask[m] <= M * 1.0000001); tightest = Math.min(tightest, em.mask[m]); }
  // worst case ratio cos^q(58.5deg) ~ -3.7 dB for q = 1.3
  assert.ok(10 * Math.log10(tightest / M) < -3);
});

test('steering is a pure shift in sine space (no element pattern)', () => {
  const arr = AH.makeArray(), N = arr.N, rng = AH.mulberry32(5);
  const wr = Float64Array.from({ length: N }, () => rng()), wi = new Float64Array(N);
  const est = Array.from({ length: N }, () => ({ mag: 1, ph: 0 }));
  const c = AH.commands(arr, wr, wi, est, 0.4, 0.1);
  const a = AH.patternPower(arr, c.cr, c.ci, Float64Array.of(0.4 + 0.07), Float64Array.of(0.1 - 0.12))[0];
  const b = AH.patternPower(arr, wr, wi, Float64Array.of(0.07), Float64Array.of(-0.12))[0];
  assert.ok(Math.abs(a - b) < 1e-8 * b);
});

test('full pipeline: factory beams break the mask, compensated beams pass every scan direction', () => {
  const sc = AH.makeScenario(20261105);
  const res = AH.run(sc);
  assert.strictEqual(res.summary.healthy.pass, res.summary.healthy.total);
  assert.ok(res.summary.factory.pass < res.summary.factory.total / 2);
  assert.strictEqual(res.summary.compensated.pass, res.summary.compensated.total);
  assert.ok(res.summary.compensated.eirpMeanDb - res.summary.factory.eirpMeanDb > -1.0, 'compensation EIRP cost under 1 dB');
});

test('compensation never reads the simulation truth', () => {
  const sc = AH.makeScenario(20261105);
  const a = AH.run(sc, { thetas: [0, 45], nPhi: 4 });
  const scrambled = Object.assign({}, sc, { truth: { gA: sc.truth.gA.map(() => 0.5), gP: sc.truth.gP.map(() => 1), hA: sc.truth.hA, hP: sc.truth.hP, cause: [] } });
  const b = AH.run(scrambled, { thetas: [0, 45], nPhi: 4 });
  for (let k = 0; k < sc.arr.N; k++) { assert.strictEqual(a.syn.wr[k], b.syn.wr[k]); assert.strictEqual(a.syn.wi[k], b.syn.wi[k]); }
});
