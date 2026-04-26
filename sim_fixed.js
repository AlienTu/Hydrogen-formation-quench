const els = {
  canvas: document.getElementById('simCanvas'),
  flavourCanvas: document.getElementById('flavourCanvas'),
  playPauseBtn: document.getElementById('playPauseBtn'),
  resetBtn: document.getElementById('resetBtn'),
  defaultPresetBtn: document.getElementById('defaultPresetBtn'),
  nearThresholdPresetBtn: document.getElementById('nearThresholdPresetBtn'),
  bcPeriodicBtn: document.getElementById('bcPeriodicBtn'),
  bcDirichletBtn: document.getElementById('bcDirichletBtn'),
  bcNeumannBtn: document.getElementById('bcNeumannBtn'),
  g: document.getElementById('g'), gNumber: document.getElementById('gNumber'),
  kappa1: document.getElementById('kappa1'), kappa1Number: document.getElementById('kappa1Number'),
  kappa2: document.getElementById('kappa2'), kappa2Number: document.getElementById('kappa2Number'),
  L: document.getElementById('L'), LNumber: document.getElementById('LNumber'),
  N: document.getElementById('N'), NNumber: document.getElementById('NNumber'),
  dt: document.getElementById('dt'), dtNumber: document.getElementById('dtNumber'),
  topCharge: document.getElementById('topCharge'), topChargeNumber: document.getElementById('topChargeNumber'),
  relaxSteps: document.getElementById('relaxSteps'), relaxStepsNumber: document.getElementById('relaxStepsNumber'),
  substeps: document.getElementById('substeps'), substepsNumber: document.getElementById('substepsNumber'),
  viewMode: document.getElementById('viewMode'),
  bcFormula: document.getElementById('bcFormula'),
  initialFormula: document.getElementById('initialFormula'),
  bcReadout: document.getElementById('bcReadout'),
  topChargeReadout: document.getElementById('topChargeReadout'),
  timeReadout: document.getElementById('timeReadout'),
  energyReadout: document.getElementById('energyReadout'),
  dxReadout: document.getElementById('dxReadout'),
  phiCReadout: document.getElementById('phiCReadout'),
  phiSReadout: document.getElementById('phiSReadout'),
  qSReadout: document.getElementById('qSReadout'),
  consistencyReadout: document.getElementById('consistencyReadout'),
  messageReadout: document.getElementById('messageReadout')
};

const ctx = els.canvas.getContext('2d');
const fctx = els.flavourCanvas.getContext('2d');
const SQRT2 = Math.sqrt(2);
const SQRT_PI = Math.sqrt(Math.PI);
const ALPHA = Math.sqrt(2 * Math.PI);
const DELTA_CS = ALPHA;
const DELTA_12 = SQRT_PI;
const BETA = 2 * SQRT_PI;

let sim = null;
let running = false;
let bcType = 'periodic';

const pairs = [
  ['g', 'gNumber'], ['kappa1', 'kappa1Number'], ['kappa2', 'kappa2Number'],
  ['L', 'LNumber'], ['N', 'NNumber'], ['dt', 'dtNumber'],
  ['topCharge', 'topChargeNumber'], ['relaxSteps', 'relaxStepsNumber'],
  ['substeps', 'substepsNumber']
];

function params() {
  return { g: Number(els.gNumber.value), kappa1: Number(els.kappa1Number.value), kappa2: Number(els.kappa2Number.value), L: Number(els.LNumber.value), N: Math.round(Number(els.NNumber.value)), dt: Number(els.dtNumber.value), z: Math.round(Number(els.topChargeNumber.value)), relaxSteps: Math.round(Number(els.relaxStepsNumber.value)), substeps: Math.round(Number(els.substepsNumber.value)), viewMode: els.viewMode.value, bcType };
}

function syncPair(from, to) { els[to].value = els[from].value; }

function setPreset() {
  els.g.value = els.gNumber.value = '1.0';
  els.kappa1.value = els.kappa1Number.value = '5.0';
  els.kappa2.value = els.kappa2Number.value = '0.5';
  els.L.value = els.LNumber.value = '200';
  els.N.value = els.NNumber.value = '2000';
  els.dt.value = els.dtNumber.value = '0.012';
  els.topCharge.value = els.topChargeNumber.value = '0';
  els.relaxSteps.value = els.relaxStepsNumber.value = '160';
  els.substeps.value = els.substepsNumber.value = '2';
  els.viewMode.value = 'both';
}

function setNearThresholdPreset() {
  els.g.value = els.gNumber.value = '1.58';
  els.kappa1.value = els.kappa1Number.value = '1.0';
  els.kappa2.value = els.kappa2Number.value = '1.0';
  els.L.value = els.LNumber.value = '240';
  els.N.value = els.NNumber.value = '2800';
  els.dt.value = els.dtNumber.value = '0.008';
  els.topCharge.value = els.topChargeNumber.value = '1';
  els.relaxSteps.value = els.relaxStepsNumber.value = '120';
  els.substeps.value = els.substepsNumber.value = '2';
  els.viewMode.value = 'both';
}

function setBoundaryCondition(nextBC) {
  bcType = nextBC;
  els.bcPeriodicBtn.classList.toggle('active', nextBC === 'periodic');
  els.bcDirichletBtn.classList.toggle('active', nextBC === 'dirichlet');
  els.bcNeumannBtn.classList.toggle('active', nextBC === 'neumann');
  reinitialize();
}

function renderFormulaBlocks() {
  let bcLatex;
  if (bcType === 'periodic') bcLatex = '\\[\\phi_a(-L/2,t)=\\phi_a(L/2,t),\\qquad \\partial_x\\phi_a(-L/2,t)=\\partial_x\\phi_a(L/2,t),\\quad a=c,s.\\]';
  else if (bcType === 'dirichlet') bcLatex = '\\[\\phi_c(\\pm L/2,t)=0,\\qquad \\phi_s(-L/2,t)=0,\\qquad \\phi_s(L/2,t)=z\\sqrt{2\\pi}.\\]';
  else bcLatex = '\\[\\partial_x\\phi_c(\\pm L/2,t)=0,\\qquad \\partial_x\\phi_s(\\pm L/2,t)=0.\\]';
  els.bcFormula.innerHTML = bcLatex;
  els.initialFormula.innerHTML = '\\[\\dot\\phi_c(x,0)=\\dot\\phi_s(x,0)=\\dot\\phi_1(x,0)=\\dot\\phi_2(x,0)=0.\\]' + '\\[z=0:\\quad \\phi_2(x,0)=0,\\qquad \\phi_1(-L/4,0)=\\phi_1(L/4,0)=\\sqrt\\pi/2.\\]' + '\\[z\\neq0:\\quad \\phi_c(x,0)=0,\\qquad \\phi_s(0,0)=\\frac z2\\sqrt{2\\pi}.\\]';
  if (window.MathJax && typeof window.MathJax.typesetPromise === 'function') window.MathJax.typesetPromise([els.bcFormula, els.initialFormula]).catch(() => {});
}

function reinitialize() { running = false; els.playPauseBtn.textContent = 'Play'; initialize(); }

function nearestIndex(x, target) { let best = 0, err = Math.abs(x[0] - target); for (let i = 1; i < x.length; i++) { const e = Math.abs(x[i] - target); if (e < err) { err = e; best = i; } } return best; }
function kinkStep(x, x0, width, amplitude) { return amplitude * (2 / Math.PI) * Math.atan(Math.exp((x - x0) / width)); }
function applyFixed(phi, fixed) { for (let i = 0; i < phi.length; i++) if (fixed[i] !== null) phi[i] = fixed[i]; }

function relaxPhi1Pinned(phi1, fixed, p) {
  const n = phi1.length, dx = p.L / p.N, dtau = 0.14 * dx * dx, old = new Float64Array(n);
  applyFixed(phi1, fixed);
  for (let step = 0; step < p.relaxSteps; step++) {
    old.set(phi1);
    for (let i = 1; i < n - 1; i++) {
      if (fixed[i] !== null) continue;
      const lap = (old[i + 1] - 2 * old[i] + old[i - 1]) / (dx * dx);
      const force = lap - 0.5 * p.g * p.g * old[i] - SQRT_PI * p.kappa1 * Math.sin(BETA * old[i]);
      phi1[i] += dtau * force;
    }
    applyFixed(phi1, fixed);
  }
}

function relaxPhiSPinned(phiS, fixed, p) {
  const n = phiS.length, dx = p.L / p.N, dtau = 0.14 * dx * dx, old = new Float64Array(n);
  applyFixed(phiS, fixed);
  for (let step = 0; step < p.relaxSteps; step++) {
    old.set(phiS);
    for (let i = 1; i < n - 1; i++) {
      if (fixed[i] !== null) continue;
      const lap = (old[i + 1] - 2 * old[i] + old[i - 1]) / (dx * dx);
      const force = lap - 0.5 * ALPHA * (p.kappa1 + p.kappa2) * Math.sin(ALPHA * old[i]);
      phiS[i] += dtau * force;
    }
    applyFixed(phiS, fixed);
  }
}

function buildInitialProfile(p) {
  const n = p.N, dx = p.L / p.N, x = new Float64Array(n), phi1 = new Float64Array(n), phi2 = new Float64Array(n), phiC = new Float64Array(n), phiS = new Float64Array(n);
  for (let i = 0; i < n; i++) x[i] = -p.L / 2 + i * dx;
  const width1 = 1.8 / Math.sqrt(Math.max(1e-8, 0.5 * p.g * p.g + 2 * Math.PI * Math.max(p.kappa1, 1e-4)));
  const widthS = 1.8 / Math.sqrt(Math.max(1e-8, Math.PI * Math.max(p.kappa1 + p.kappa2, 1e-4)));

  if (p.bcType === 'periodic' || p.z === 0) {
    for (let i = 0; i < n; i++) { phi1[i] = kinkStep(x[i], -p.L / 4, width1, SQRT_PI) - kinkStep(x[i], p.L / 4, width1, SQRT_PI); phi2[i] = 0; }
    if (p.bcType !== 'periodic') {
      const fixed1 = new Array(n).fill(null);
      fixed1[0] = 0; fixed1[n - 1] = 0; fixed1[nearestIndex(x, -p.L / 4)] = SQRT_PI / 2; fixed1[nearestIndex(x, p.L / 4)] = SQRT_PI / 2;
      relaxPhi1Pinned(phi1, fixed1, p);
    }
    for (let i = 0; i < n; i++) { phiC[i] = (phi1[i] + phi2[i]) / SQRT2; phiS[i] = (phi1[i] - phi2[i]) / SQRT2; }
  } else {
    const z = p.z;
    for (let i = 0; i < n; i++) { phiC[i] = 0; phiS[i] = z * kinkStep(x[i], 0, widthS, DELTA_CS); }
    if (p.bcType === 'dirichlet') {
      const fixedS = new Array(n).fill(null);
      fixedS[0] = 0; fixedS[n - 1] = z * DELTA_CS; fixedS[nearestIndex(x, 0)] = 0.5 * z * DELTA_CS;
      relaxPhiSPinned(phiS, fixedS, p);
    }
    if (p.bcType === 'neumann') { phiS[0] = phiS[1]; phiS[n - 1] = phiS[n - 2]; }
    for (let i = 0; i < n; i++) { phi1[i] = (phiC[i] + phiS[i]) / SQRT2; phi2[i] = (phiC[i] - phiS[i]) / SQRT2; }
  }

  if (p.bcType === 'dirichlet') {
    phiC[0] = 0; phiC[n - 1] = 0; phiS[0] = 0; phiS[n - 1] = p.z * DELTA_CS;
    for (const i of [0, n - 1]) { phi1[i] = (phiC[i] + phiS[i]) / SQRT2; phi2[i] = (phiC[i] - phiS[i]) / SQRT2; }
  }
  if (p.bcType === 'neumann') {
    phiC[0] = phiC[1]; phiC[n - 1] = phiC[n - 2]; phiS[0] = phiS[1]; phiS[n - 1] = phiS[n - 2]; phi1[0] = phi1[1]; phi1[n - 1] = phi1[n - 2]; phi2[0] = phi2[1]; phi2[n - 1] = phi2[n - 2];
  }
  return { x, phiC, phiS, phi1, phi2 };
}

function initialize() {
  const p = params(), init = buildInitialProfile(p), zero = new Float64Array(p.N);
  sim = { ...p, dx: p.L / p.N, x: init.x, phiC: init.phiC, phiS: init.phiS, piC: zero.slice(), piS: zero.slice(), phi1: init.phi1, phi2: init.phi2, pi1: zero.slice(), pi2: zero.slice(), t: 0 };
  enforceBoundaryCS(); enforceBoundary12(); renderFormulaBlocks(); updateReadout('initialized'); draw(); drawFlavourCheck();
}

function laplacian(field) {
  const n = field.length, out = new Float64Array(n), inv = 1 / (sim.dx * sim.dx);
  for (let i = 0; i < n; i++) {
    let im = i - 1, ip = i + 1;
    if (sim.bcType === 'periodic') { if (i === 0) im = n - 1; if (i === n - 1) ip = 0; }
    else if (sim.bcType === 'neumann') { if (i === 0) im = 1; if (i === n - 1) ip = n - 2; }
    else { if (i === 0) im = 0; if (i === n - 1) ip = n - 1; }
    out[i] = (field[ip] - 2 * field[i] + field[im]) * inv;
  }
  return out;
}

function enforceBoundaryCS() {
  const n = sim.phiC.length;
  if (sim.bcType === 'periodic') return;
  if (sim.bcType === 'dirichlet') { sim.phiC[0] = 0; sim.phiC[n - 1] = 0; sim.phiS[0] = 0; sim.phiS[n - 1] = sim.z * DELTA_CS; sim.piC[0] = 0; sim.piC[n - 1] = 0; sim.piS[0] = 0; sim.piS[n - 1] = 0; }
  else { sim.phiC[0] = sim.phiC[1]; sim.phiC[n - 1] = sim.phiC[n - 2]; sim.phiS[0] = sim.phiS[1]; sim.phiS[n - 1] = sim.phiS[n - 2]; sim.piC[0] = sim.piC[1]; sim.piC[n - 1] = sim.piC[n - 2]; sim.piS[0] = sim.piS[1]; sim.piS[n - 1] = sim.piS[n - 2]; }
}

function enforceBoundary12() {
  const n = sim.phi1.length;
  if (sim.bcType === 'periodic') return;
  if (sim.bcType === 'dirichlet') { sim.phi1[0] = 0; sim.phi2[0] = 0; sim.phi1[n - 1] = sim.z * SQRT_PI; sim.phi2[n - 1] = -sim.z * SQRT_PI; sim.pi1[0] = 0; sim.pi1[n - 1] = 0; sim.pi2[0] = 0; sim.pi2[n - 1] = 0; }
  else { sim.phi1[0] = sim.phi1[1]; sim.phi1[n - 1] = sim.phi1[n - 2]; sim.phi2[0] = sim.phi2[1]; sim.phi2[n - 1] = sim.phi2[n - 2]; sim.pi1[0] = sim.pi1[1]; sim.pi1[n - 1] = sim.pi1[n - 2]; sim.pi2[0] = sim.pi2[1]; sim.pi2[n - 1] = sim.pi2[n - 2]; }
}

function accelerationsCS(phiC, phiS) {
  const n = phiC.length, lapC = laplacian(phiC), lapS = laplacian(phiS), aC = new Float64Array(n), aS = new Float64Array(n);
  for (let i = 0; i < n; i++) { const arg1 = ALPHA * (phiC[i] + phiS[i]), arg2 = ALPHA * (phiC[i] - phiS[i]); aC[i] = lapC[i] - sim.g * sim.g * phiC[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1) - 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2); aS[i] = lapS[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1) + 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2); }
  return { aC, aS };
}

function accelerations12(phi1, phi2) {
  const n = phi1.length, lap1 = laplacian(phi1), lap2 = laplacian(phi2), a1 = new Float64Array(n), a2 = new Float64Array(n);
  for (let i = 0; i < n; i++) { const chargeTerm = 0.5 * sim.g * sim.g * (phi1[i] + phi2[i]); a1[i] = lap1[i] - chargeTerm - SQRT_PI * sim.kappa1 * Math.sin(BETA * phi1[i]); a2[i] = lap2[i] - chargeTerm - SQRT_PI * sim.kappa2 * Math.sin(BETA * phi2[i]); }
  return { a1, a2 };
}

function stepLeapfrog() {
  const n = sim.phiC.length, dt = sim.dt, cs = accelerationsCS(sim.phiC, sim.phiS), fl = accelerations12(sim.phi1, sim.phi2);
  const pcH = new Float64Array(n), psH = new Float64Array(n), p1H = new Float64Array(n), p2H = new Float64Array(n), cNew = new Float64Array(n), sNew = new Float64Array(n), f1New = new Float64Array(n), f2New = new Float64Array(n);
  for (let i = 0; i < n; i++) { pcH[i] = sim.piC[i] + 0.5 * dt * cs.aC[i]; psH[i] = sim.piS[i] + 0.5 * dt * cs.aS[i]; p1H[i] = sim.pi1[i] + 0.5 * dt * fl.a1[i]; p2H[i] = sim.pi2[i] + 0.5 * dt * fl.a2[i]; cNew[i] = sim.phiC[i] + dt * pcH[i]; sNew[i] = sim.phiS[i] + dt * psH[i]; f1New[i] = sim.phi1[i] + dt * p1H[i]; f2New[i] = sim.phi2[i] + dt * p2H[i]; }
  sim.phiC = cNew; sim.phiS = sNew; sim.phi1 = f1New; sim.phi2 = f2New; enforceBoundaryCS(); enforceBoundary12();
  const cs2 = accelerationsCS(sim.phiC, sim.phiS), fl2 = accelerations12(sim.phi1, sim.phi2);
  for (let i = 0; i < n; i++) { sim.piC[i] = pcH[i] + 0.5 * dt * cs2.aC[i]; sim.piS[i] = psH[i] + 0.5 * dt * cs2.aS[i]; sim.pi1[i] = p1H[i] + 0.5 * dt * fl2.a1[i]; sim.pi2[i] = p2H[i] + 0.5 * dt * fl2.a2[i]; }
  enforceBoundaryCS(); enforceBoundary12(); sim.t += dt;
}

function maxAbs(arr) { let m = 0; for (let i = 0; i < arr.length; i++) m = Math.max(m, Math.abs(arr[i])); return m; }
function minValue(arr) { let m = Infinity; for (let i = 0; i < arr.length; i++) m = Math.min(m, arr[i]); return m; }
function maxValue(arr) { let m = -Infinity; for (let i = 0; i < arr.length; i++) m = Math.max(m, arr[i]); return m; }
function convertedFields() { const n = sim.phi1.length, phiCFrom12 = new Float64Array(n), phiSFrom12 = new Float64Array(n), diffC = new Float64Array(n), diffS = new Float64Array(n); for (let i = 0; i < n; i++) { phiCFrom12[i] = (sim.phi1[i] + sim.phi2[i]) / SQRT2; phiSFrom12[i] = (sim.phi1[i] - sim.phi2[i]) / SQRT2; diffC[i] = phiCFrom12[i] - sim.phiC[i]; diffS[i] = phiSFrom12[i] - sim.phiS[i]; } return { phiCFrom12, phiSFrom12, diffC, diffS }; }
function consistencyError() { const d = convertedFields(); return Math.max(maxAbs(d.diffC), maxAbs(d.diffS)); }

function energyDensity() {
  const n = sim.phiC.length, ed = new Float64Array(n);
  for (let i = 0; i < n; i++) { let ip = i + 1, im = i - 1; if (sim.bcType === 'periodic') { if (i === 0) im = n - 1; if (i === n - 1) ip = 0; } else { if (i === 0) im = 0; if (i === n - 1) ip = n - 1; } const gradC = (sim.phiC[ip] - sim.phiC[im]) / (2 * sim.dx), gradS = (sim.phiS[ip] - sim.phiS[im]) / (2 * sim.dx), arg1 = ALPHA * (sim.phiC[i] + sim.phiS[i]), arg2 = ALPHA * (sim.phiC[i] - sim.phiS[i]); const potential = 0.5 * sim.g * sim.g * sim.phiC[i] * sim.phiC[i] - 0.5 * sim.kappa1 * Math.cos(arg1) - 0.5 * sim.kappa2 * Math.cos(arg2); ed[i] = 0.5 * sim.piC[i] * sim.piC[i] + 0.5 * sim.piS[i] * sim.piS[i] + 0.5 * gradC * gradC + 0.5 * gradS * gradS + potential; }
  return ed;
}
function totalEnergy(ed) { let sum = 0; for (let i = 0; i < ed.length; i++) sum += ed[i]; return sum * sim.dx; }
function updateReadout(message) { const ed = energyDensity(); els.bcReadout.textContent = sim.bcType; els.topChargeReadout.textContent = sim.z.toString(); els.timeReadout.textContent = sim.t.toFixed(2); els.energyReadout.textContent = totalEnergy(ed).toFixed(4); els.dxReadout.textContent = sim.dx.toFixed(4); els.phiCReadout.textContent = (maxAbs(sim.phiC) / DELTA_CS).toFixed(3); els.phiSReadout.textContent = (maxAbs(sim.phiS) / DELTA_CS).toFixed(3); els.qSReadout.textContent = ((sim.phiS[sim.phiS.length - 1] - sim.phiS[0]) / DELTA_CS).toFixed(3); els.consistencyReadout.textContent = consistencyError().toExponential(3); els.messageReadout.textContent = message; }
function formatAxisValue(v) { const av = Math.abs(v); if (av >= 100) return v.toFixed(0); if (av >= 10) return v.toFixed(1); if (av >= 1e-2) return v.toFixed(2); return v.toExponential(1); }
function boundsForArrays(arrays, sectorStep = null) { let yMinRaw = Infinity, yMaxRaw = -Infinity; arrays.forEach(arr => { yMinRaw = Math.min(yMinRaw, minValue(arr)); yMaxRaw = Math.max(yMaxRaw, maxValue(arr)); }); if (Math.abs(yMaxRaw - yMinRaw) < 1e-12) { yMinRaw -= 1; yMaxRaw += 1; } if (sectorStep === null) { const pad = 0.18 * (yMaxRaw - yMinRaw); return { ymin: yMinRaw - pad, ymax: yMaxRaw + pad, sectorStep: null, nLow: 1, nHigh: 0 }; } let nLow = Math.ceil(yMinRaw / sectorStep), nHigh = Math.floor(yMaxRaw / sectorStep); if (nLow > nHigh) { const nearest = Math.round(0.5 * (yMinRaw + yMaxRaw) / sectorStep); nLow = nearest; nHigh = nearest; } return { ymin: (nLow - 1) * sectorStep, ymax: (nHigh + 1) * sectorStep, sectorStep, nLow, nHigh }; }

function drawAxesOn(c, x0, y0, w, h, bounds, title) { c.strokeStyle = '#cbd5e1'; c.lineWidth = 1; c.strokeRect(x0, y0, w, h); c.fillStyle = '#0f172a'; c.font = '20px sans-serif'; c.fillText(title, x0, y0 - 12); c.fillStyle = '#475569'; c.font = '15px sans-serif'; c.fillText(formatAxisValue(bounds.ymax), x0 - 44, y0 + 14); c.fillText(formatAxisValue(bounds.ymin), x0 - 44, y0 + h); c.fillText((-sim.L / 2).toFixed(0), x0 - 4, y0 + h + 24); c.fillText((sim.L / 2).toFixed(0), x0 + w - 34, y0 + h + 24); if (bounds.sectorStep !== null) { for (let n = bounds.nLow; n <= bounds.nHigh; n++) { const yy = y0 + h - ((n * bounds.sectorStep - bounds.ymin) / (bounds.ymax - bounds.ymin)) * h; c.setLineDash([7, 7]); c.strokeStyle = '#94a3b8'; c.beginPath(); c.moveTo(x0, yy); c.lineTo(x0 + w, yy); c.stroke(); c.setLineDash([]); } } }
function drawCurveOn(c, arr, x0, y0, w, h, ymin, ymax, color, lineWidth = 2.4, dash = []) { c.strokeStyle = color; c.lineWidth = lineWidth; c.setLineDash(dash); c.beginPath(); for (let i = 0; i < arr.length; i++) { const xx = x0 + (i / (arr.length - 1)) * w, yy = y0 + h - ((arr[i] - ymin) / (ymax - ymin)) * h; if (i === 0) c.moveTo(xx, yy); else c.lineTo(xx, yy); } c.stroke(); c.setLineDash([]); }
function drawLegendOn(c, entries, x, y) { c.font = '14px sans-serif'; entries.forEach((entry, idx) => { const yy = y + idx * 20; c.strokeStyle = entry.color; c.lineWidth = 3; c.setLineDash(entry.dash || []); c.beginPath(); c.moveTo(x, yy); c.lineTo(x + 20, yy); c.stroke(); c.setLineDash([]); c.fillStyle = '#0f172a'; c.fillText(entry.label, x + 28, yy + 5); }); }
function drawSinglePanel(c, arr, x0, y0, w, h, title, color, sectorStep) { const bounds = boundsForArrays([arr], sectorStep); drawAxesOn(c, x0, y0, w, h, bounds, title); drawCurveOn(c, arr, x0, y0, w, h, bounds.ymin, bounds.ymax, color); }
function drawComparePanel(c, arrA, arrB, x0, y0, w, h, title, colorA, colorB, sectorStep, labelA, labelB) { const bounds = boundsForArrays([arrA, arrB], sectorStep); drawAxesOn(c, x0, y0, w, h, bounds, title); drawCurveOn(c, arrA, x0, y0, w, h, bounds.ymin, bounds.ymax, colorA); drawCurveOn(c, arrB, x0, y0, w, h, bounds.ymin, bounds.ymax, colorB, 2.2, [10, 5]); drawLegendOn(c, [{ color: colorA, label: labelA }, { color: colorB, label: labelB, dash: [10, 5] }], x0 + w - 172, y0 + 22); }
function drawResidualPanel(c, arr, x0, y0, w, h, title, color) { const bounds = boundsForArrays([arr], null); drawAxesOn(c, x0, y0, w, h, bounds, title); drawCurveOn(c, arr, x0, y0, w, h, bounds.ymin, bounds.ymax, color); }
function draw() { ctx.clearRect(0, 0, els.canvas.width, els.canvas.height); ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, els.canvas.width, els.canvas.height); const padX = 64, padY = 52, gapX = 38, gapY = 52, plotW = (els.canvas.width - 2 * padX - gapX) / 2, plotH = (els.canvas.height - 2 * padY - gapY) / 2; if (sim.viewMode === 'both' || sim.viewMode === 'charge_only') { drawSinglePanel(ctx, sim.phiC, padX, padY, plotW, plotH, 'phi_c(x,t)', '#2563eb', DELTA_CS); drawSinglePanel(ctx, sim.phiS, padX + plotW + gapX, padY, plotW, plotH, 'phi_s(x,t)', '#dc2626', DELTA_CS); } if (sim.viewMode === 'both' || sim.viewMode === 'flavour_only') { const y0 = padY + plotH + gapY; drawSinglePanel(ctx, sim.phi1, padX, y0, plotW, plotH, 'phi_1(x,t)', '#7c3aed', DELTA_12); drawSinglePanel(ctx, sim.phi2, padX + plotW + gapX, y0, plotW, plotH, 'phi_2(x,t)', '#059669', DELTA_12); } }
function drawFlavourCheck() { fctx.clearRect(0, 0, els.flavourCanvas.width, els.flavourCanvas.height); fctx.fillStyle = '#ffffff'; fctx.fillRect(0, 0, els.flavourCanvas.width, els.flavourCanvas.height); const d = convertedFields(), padX = 64, padY = 52, gapX = 38, gapY = 54, plotW = (els.flavourCanvas.width - 2 * padX - gapX) / 2, plotH = (els.flavourCanvas.height - 2 * padY - gapY) / 2; drawComparePanel(fctx, sim.phiC, d.phiCFrom12, padX, padY, plotW, plotH, 'phi_c: direct vs reconstructed', '#2563eb', '#f97316', DELTA_CS, 'direct', 'reconstructed'); drawComparePanel(fctx, sim.phiS, d.phiSFrom12, padX + plotW + gapX, padY, plotW, plotH, 'phi_s: direct vs reconstructed', '#dc2626', '#0ea5e9', DELTA_CS, 'direct', 'reconstructed'); const y0 = padY + plotH + gapY; drawResidualPanel(fctx, d.diffC, padX, y0, plotW, plotH, 'Delta phi_c', '#111827'); drawResidualPanel(fctx, d.diffS, padX + plotW + gapX, y0, plotW, plotH, 'Delta phi_s', '#ef4444'); }
function tick() { if (running) { for (let k = 0; k < sim.substeps; k++) stepLeapfrog(); updateReadout('running'); draw(); drawFlavourCheck(); } requestAnimationFrame(tick); }
function bindPair(rangeId, numberId) { els[rangeId].addEventListener('input', () => { syncPair(rangeId, numberId); reinitialize(); }); els[numberId].addEventListener('input', () => { syncPair(numberId, rangeId); reinitialize(); }); }

els.playPauseBtn.addEventListener('click', () => { running = !running; els.playPauseBtn.textContent = running ? 'Pause' : 'Play'; updateReadout(running ? 'running' : 'paused'); });
els.resetBtn.addEventListener('click', reinitialize);
els.defaultPresetBtn.addEventListener('click', () => { setPreset(); reinitialize(); });
els.nearThresholdPresetBtn.addEventListener('click', () => { setNearThresholdPreset(); reinitialize(); });
els.bcPeriodicBtn.addEventListener('click', () => setBoundaryCondition('periodic'));
els.bcDirichletBtn.addEventListener('click', () => setBoundaryCondition('dirichlet'));
els.bcNeumannBtn.addEventListener('click', () => setBoundaryCondition('neumann'));
pairs.forEach(([rangeId, numberId]) => bindPair(rangeId, numberId));
els.viewMode.addEventListener('input', reinitialize);

setPreset();
initialize();
tick();
