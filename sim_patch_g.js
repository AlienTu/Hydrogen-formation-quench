// Patch the quench page to use the Schwinger convention
// V_charge = g^2 (phi1+phi2)^2 / 2 = g^2 phi_c^2.
// Therefore the EOM contain +2 g^2 phi_c in charge basis and +g^2(phi1+phi2) in flavour basis.

(function(){
  if (typeof window === 'undefined') return;

  window.addEventListener('DOMContentLoaded', () => {
    const target = document.querySelector('.top-buttons') || document.querySelector('.hero-card');
    if (target && !document.getElementById('spectrumLinkBtn')) {
      const a = document.createElement('a');
      a.id = 'spectrumLinkBtn';
      a.href = 'spectrum.html';
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
    console.warn('Schwinger g patch failed:', e);
  }
})();
