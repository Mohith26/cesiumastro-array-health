// Generates the synthetic on-orbit telemetry (inputs) and the simulation truth (kept separate;
// the estimator never reads truth.json). Usage: node scripts/make_data.js [seed]
const fs = require('fs'), path = require('path');
const AH = require('../core.js');
const seed = Number(process.argv[2] || 20261105);
const sc = AH.makeScenario(seed);
const r6 = (x) => Math.round(x * 1e9) / 1e9;
const tel = sc.telemetry.map((t) => ({ tx: t.tx, rx: t.rx, ring: t.ring, m0: t.m0.map(r6), m: t.m.map(r6), sigma: r6(t.sigma), sigma0: r6(t.sigma0) }));
const out = path.join(__dirname, '..', 'data');
fs.writeFileSync(path.join(out, 'telemetry.json'), JSON.stringify({ seed, array: sc.arr.opts, note: 'Synthetic element-pair coupling telemetry: factory baseline m0 and on-orbit m (complex, linear volts). Not CesiumAstro data.', pairs: tel }));
fs.writeFileSync(path.join(out, 'truth.json'), JSON.stringify({ seed, note: 'Simulation truth: per-element Tx (g) and Rx (h) change since factory calibration. Used only to score the estimator and to radiate the true hardware.', gA: sc.truth.gA.map(r6), gP: sc.truth.gP.map(r6), hA: sc.truth.hA.map(r6), hP: sc.truth.hP.map(r6), cause: sc.truth.cause }));
console.log('wrote', tel.length, 'pairs, seed', seed);
