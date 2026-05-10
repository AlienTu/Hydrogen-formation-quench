// Patch the quench page to use the Schwinger convention
// V_charge = g^2 (phi1+phi2)^2 / 2 = g^2 phi_c^2.
// Therefore the EOM contain +2 g^2 phi_c in charge basis and +g^2(phi1+phi2) in flavour basis.
// This file also patches the z=0 initial kink-antikink separation so it is user-adjustable.

(function(){
  if (typeof window === 'undefined') return;

  function numericValue(id, fallback) {
    const el = document.getElementById(id);
    if (!el) return fallback;
    const v = Number(el.value);
    return Number.isFinite(v) ? v : fallback;
  }

  function kinkDistanceDefault(g, kappa2) {
    const gg = Math.max(Math.abs(Number(g)), 1e-8);
    const k2 = Math.max(Number(kappa2) || 0, 0);
    return (2 / Math.PI) * Math.sqrt(gg * gg + 2 * Math.PI * k2) / (gg * gg);
  }

  function kinkDistanceDefaultFromInputs() {
    return kinkDistanceDefault(numericValue('gNumber', 1.0), numericValue('kappa2Number', 0.5));
  }

  function kinkDistanceMaxFromInputs() {
    return Math.max(1.0, 0.98 * Math.max(numericValue('LNumber', 200), 1.0));
  }

  function clampKinkDistance(D, L) {
    const maxD = Math.max(0, 0.98 * Math.max(Number(L) || 0, 0));
    if (!Number.isFinite(D)) return Math.min(kinkDistanceDefaultFromInputs(), maxD || Infinity);
    return Math.min(Math.max(D, 0), maxD || D);
  }

  function updateKinkDistanceRangeLimit() {
    const range = document.getElementById('kinkDistance');
    const number = document.getElementById('kinkDistanceNumber');
    if (!range || !number) return;
    const maxD = kinkDistanceMaxFromInputs();
    range.max = maxD.toFixed(3);
    number.max = maxD.toFixed(3);
  }

  function setKinkDistanceValue(D) {
    const range = document.getElementById('kinkDistance');
    const number = document.getElementById('kinkDistanceNumber');
    if (!range || !number) return;
    updateKinkDistanceRangeLimit();
    const maxD = Number(range.max) || kinkDistanceMaxFromInputs();
    const val = Math.min(Math.max(Number(D), 0), maxD);
    const text = val.toFixed(3);
    range.value = text;
    number.value = text;
  }

  function readKinkDistance(p) {
    const el = document.getElementById('kinkDistanceNumber');
    const raw = el ? Number(el.value) : kinkDistanceDefault(p.g, p.kappa2);
    return clampKinkDistance(raw, p.L);
  }

  function installKinkDistanceControls() {
    const topCharge = document.getElementById('topCharge');
    const controls = document.querySelector('.controls-horizontal');
    if (!controls || !topCharge) return;

    let range = document.getElementById('kinkDistance');
    let number = document.getElementById('kinkDistanceNumber');
    if (!range || !number) {
      const holder = document.createElement('div');
      holder.className = 'control triple';
      holder.innerHTML = '<label for="kinkDistance">D <small>z=0 kink-antikink separation</small></label>' +
        '<input id="kinkDistance" type="range" min="0" max="196" step="0.05" value="1.296" />' +
        '<input id="kinkDistanceNumber" type="number" min="0" max="196" step="0.05" value="1.296" />';
      const anchor = topCharge.closest('.control.triple');
      if (anchor && anchor.parentNode) anchor.parentNode.insertBefore(holder, anchor.nextSibling);
      else controls.appendChild(holder);
      range = document.getElementById('kinkDistance');
      number = document.getElementById('kinkDistanceNumber');
    }

    let readout = document.getElementById('kinkDistanceReadout');
    if (!readout) {
      const statusCard = document.querySelector('.status-card.status');
      if (statusCard) {
        const row = document.createElement('div');
        row.innerHTML = '<strong>kink separation D:</strong> <span id="kinkDistanceReadout">0.000</span>';
        const topChargeReadout = document.getElementById('topChargeReadout');
        const topChargeRow = topChargeReadout ? topChargeReadout.closest('div') : null;
        if (topChargeRow && topChargeRow.parentNode) topChargeRow.parentNode.insertBefore(row, topChargeRow.nextSibling);
        else statusCard.appendChild(row);
        readout = document.getElementById('kinkDistanceReadout');
      }
    }

    try {
      if (typeof els !== 'undefined') {
        els.kinkDistance = range;
        els.kinkDistanceNumber = number;
        els.kinkDistanceReadout = readout;
      }
    } catch (_) {}

    let manualD = false;
    const applyManualD = (fromRange) => {
      manualD = true;
      if (fromRange) number.value = range.value;
      else range.value = number.value;
      setKinkDistanceValue(Number(number.value));
      if (typeof reinitialize === 'function') reinitialize();
    };
    range.addEventListener('input', () => applyManualD(true));
    number.addEventListener('input', () => applyManualD(false));

    const maybeRefreshDefaultD = () => {
      updateKinkDistanceRangeLimit();
      if (!manualD) setKinkDistanceValue(kinkDistanceDefaultFromInputs());
      else setKinkDistanceValue(Number(number.value));
      if (typeof reinitialize === 'function') reinitialize();
    };
    ['g', 'gNumber', 'kappa2', 'kappa2Number', 'L', 'LNumber'].forEach((id) => {
      const el = document.getElementById(id);
      if (el) el.addEventListener('input', maybeRefreshDefaultD);
    });

    const defaultBtn = document.getElementById('defaultPresetBtn');
    const nearBtn = document.getElementById('nearThresholdPresetBtn');
    if (defaultBtn) defaultBtn.addEventListener('click', () => { manualD = false; setTimeout(maybeRefreshDefaultD, 0); });
    if (nearBtn) nearBtn.addEventListener('click', () => { manualD = false; setTimeout(maybeRefreshDefaultD, 0); });

    manualD = false;
    setKinkDistanceValue(kinkDistanceDefaultFromInputs());
    if (typeof reinitialize === 'function') reinitialize();
  }

  if (document.readyState === 'loading') window.addEventListener('DOMContentLoaded', installKinkDistanceControls);
  else installKinkDistanceControls();

  window.addEventListener('DOMContentLoaded', () => {
    const existing = document.getElementById('staticSpectrumTopLink');
    if (existing) {
      existing.href = 'spectrum/';
      existing.target = '_top';
      existing.textContent = 'Open static spectrum solver';
    }

    const target = document.querySelector('.top-buttons') || document.querySelector('.hero-card');
    if (target && !document.getElementById('spectrumLinkBtn')) {
      const a = document.createElement('a');
      a.id = 'spectrumLinkBtn';
      a.href = 'spectrum/';
      a.target = '_top';
      a.textContent = 'Static spectrum solver';
      a.style.display = 'inline-block';
      a.style.margin = '8px 0';
      a.style.padding = '8px 10px';
      a.style.borderRadius = '8px';
      a.style.background = '#111827';
      a.style.color = 'white';
      a.style.textDecoration = 'none';
      target.parentNode.insertBefore(a, target);
    }
  });

  try {
    const baseParams = typeof params === 'function' ? params : null;
    if (baseParams) {
      params = function() {
        const p = baseParams();
        p.kinkDistance = readKinkDistance(p);
        return p;
      };
    }

    const baseSetPreset = typeof setPreset === 'function' ? setPreset : null;
    if (baseSetPreset) {
      setPreset = function() {
        baseSetPreset();
        setKinkDistanceValue(kinkDistanceDefaultFromInputs());
      };
    }

    const baseSetNearThresholdPreset = typeof setNearThresholdPreset === 'function' ? setNearThresholdPreset : null;
    if (baseSetNearThresholdPreset) {
      setNearThresholdPreset = function() {
        baseSetNearThresholdPreset();
        setKinkDistanceValue(kinkDistanceDefaultFromInputs());
      };
    }

    renderFormulaBlocks = function() {
      let bcLatex;
      if (bcType === 'periodic') bcLatex = '\\[\\phi_a(-L/2,t)=\\phi_a(L/2,t),\\qquad \\partial_x\\phi_a(-L/2,t)=\\partial_x\\phi_a(L/2,t),\\quad a=c,s.\\]';
      else if (bcType === 'dirichlet') bcLatex = '\\[\\phi_c(\\pm L/2,t)=0,\\qquad \\phi_s(-L/2,t)=0,\\qquad \\phi_s(L/2,t)=z\\sqrt{2\\pi}.\\]';
      else bcLatex = '\\[\\partial_x\\phi_c(\\pm L/2,t)=0,\\qquad \\partial_x\\phi_s(\\pm L/2,t)=0.\\]';
      els.bcFormula.innerHTML = bcLatex;
      els.initialFormula.innerHTML = '\\[\\dot\\phi_c(x,0)=\\dot\\phi_s(x,0)=\\dot\\phi_1(x,0)=\\dot\\phi_2(x,0)=0.\\]' +
        '\\[z=0:\\quad x_K=-\\frac{D}{2},\\qquad x_{\\bar K}=+\\frac{D}{2},\\qquad \\phi_2(x,0)=0,\\qquad \\phi_1(\\pm D/2,0)=\\sqrt\\pi/2.\\]' +
        '\\[D_0=\\frac{2}{\\pi}\\,\\frac{\\sqrt{g^2+2\\pi\\kappa_2}}{g^2}.\\]' +
        '\\[z\\neq0:\\quad \\phi_c(x,0)=0,\\qquad \\phi_s(0,0)=\\frac z2\\sqrt{2\\pi}.\\]';
      if (window.MathJax && typeof window.MathJax.typesetPromise === 'function') window.MathJax.typesetPromise([els.bcFormula, els.initialFormula]).catch(() => {});
    };

    buildInitialProfile = function(p) {
      const n = p.N, dx = p.L / p.N, x = new Float64Array(n), phi1 = new Float64Array(n), phi2 = new Float64Array(n), phiC = new Float64Array(n), phiS = new Float64Array(n);
      for (let i = 0; i < n; i++) x[i] = -p.L / 2 + i * dx;
      const width1 = 1.8 / Math.sqrt(Math.max(1e-8, 0.5 * p.g * p.g + 2 * Math.PI * Math.max(p.kappa1, 1e-4)));
      const widthS = 1.8 / Math.sqrt(Math.max(1e-8, Math.PI * Math.max(p.kappa1 + p.kappa2, 1e-4)));

      if (p.bcType === 'periodic' || p.z === 0) {
        const D = clampKinkDistance(p.kinkDistance, p.L);
        const xK = -0.5 * D;
        const xAK = 0.5 * D;
        for (let i = 0; i < n; i++) {
          phi1[i] = kinkStep(x[i], xK, width1, SQRT_PI) - kinkStep(x[i], xAK, width1, SQRT_PI);
          phi2[i] = 0;
        }
        if (p.bcType !== 'periodic') {
          const fixed1 = new Array(n).fill(null);
          fixed1[0] = 0; fixed1[n - 1] = 0;
          fixed1[nearestIndex(x, xK)] = SQRT_PI / 2;
          fixed1[nearestIndex(x, xAK)] = SQRT_PI / 2;
          relaxPhi1Pinned(phi1, fixed1, p);
        }
        for (let i = 0; i < n; i++) {
          phiC[i] = (phi1[i] + phi2[i]) / SQRT2;
          phiS[i] = (phi1[i] - phi2[i]) / SQRT2;
        }
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
    };

    const baseUpdateReadout = typeof updateReadout === 'function' ? updateReadout : null;
    if (baseUpdateReadout) {
      updateReadout = function(message) {
        baseUpdateReadout(message);
        const readout = document.getElementById('kinkDistanceReadout');
        if (readout && sim && Number.isFinite(sim.kinkDistance)) readout.textContent = sim.kinkDistance.toFixed(3);
      };
    }

    accelerationsCS = function(phiC, phiS) {
      const n = phiC.length;
      const lapC = laplacian(phiC);
      const lapS = laplacian(phiS);
      const aC = new Float64Array(n);
      const aS = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const arg1 = ALPHA * (phiC[i] + phiS[i]);
        const arg2 = ALPHA * (phiC[i] - phiS[i]);
        aC[i] = lapC[i] - 2 * sim.g * sim.g * phiC[i]
          - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1)
          - 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2);
        aS[i] = lapS[i]
          - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1)
          + 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2);
      }
      return { aC, aS };
    };

    accelerations12 = function(phi1, phi2) {
      const n = phi1.length;
      const lap1 = laplacian(phi1);
      const lap2 = laplacian(phi2);
      const a1 = new Float64Array(n);
      const a2 = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        const chargeTerm = sim.g * sim.g * (phi1[i] + phi2[i]);
        a1[i] = lap1[i] - chargeTerm - SQRT_PI * sim.kappa1 * Math.sin(BETA * phi1[i]);
        a2[i] = lap2[i] - chargeTerm - SQRT_PI * sim.kappa2 * Math.sin(BETA * phi2[i]);
      }
      return { a1, a2 };
    };

    energyDensity = function() {
      const n = sim.phiC.length;
      const ed = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        let ip = i + 1, im = i - 1;
        if (sim.bcType === 'periodic') {
          if (i === 0) im = n - 1;
          if (i === n - 1) ip = 0;
        } else {
          if (i === 0) im = 0;
          if (i === n - 1) ip = n - 1;
        }
        const gradC = (sim.phiC[ip] - sim.phiC[im]) / (2 * sim.dx);
        const gradS = (sim.phiS[ip] - sim.phiS[im]) / (2 * sim.dx);
        const arg1 = ALPHA * (sim.phiC[i] + sim.phiS[i]);
        const arg2 = ALPHA * (sim.phiC[i] - sim.phiS[i]);
        const potential = sim.g * sim.g * sim.phiC[i] * sim.phiC[i]
          - 0.5 * sim.kappa1 * Math.cos(arg1)
          - 0.5 * sim.kappa2 * Math.cos(arg2);
        ed[i] = 0.5 * sim.piC[i] * sim.piC[i] + 0.5 * sim.piS[i] * sim.piS[i]
          + 0.5 * gradC * gradC + 0.5 * gradS * gradS + potential;
      }
      return ed;
    };
  } catch (e) {
    console.warn('Schwinger g / kink separation patch failed:', e);
  }
})();
