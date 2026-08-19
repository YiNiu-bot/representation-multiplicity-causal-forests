"use strict";

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function choose(n, k) {
  if (k < 0 || n < k) return 0;
  k = Math.min(k, n - k);
  let out = 1;
  for (let j = 1; j <= k; j += 1) out *= (n - k + j) / j;
  return out;
}

function matVec(v, K) {
  const out = Array(v.length).fill(0);
  for (let i = 0; i < v.length; i += 1) {
    for (let j = 0; j < v.length; j += 1) out[j] += v[i] * K[i][j];
  }
  return out;
}

function evolve(v, K, steps) {
  let out = v.slice();
  for (let j = 0; j < steps; j += 1) out = matVec(out, K);
  return out;
}

function hitWithin(K, success, state, steps) {
  const n = K.length;
  let value = Array(n).fill(0).map((_, i) => (success.has(i) ? 1 : 0));
  for (let h = 0; h < steps; h += 1) {
    const next = value.slice();
    for (let i = 0; i < n; i += 1) {
      if (!success.has(i)) {
        next[i] = K[i].reduce((sum, p, j) => sum + p * value[j], 0);
      }
    }
    value = next;
  }
  return value[state];
}

// Null actions make the class-level kernel exactly stochastic.
{
  const M = 10;
  const q = 3;
  const denom = choose(M, q);
  const p1 = (choose(M, q) - choose(M - 2, q)) / denom;
  const p2 = (choose(M - 2, q) - choose(M - 2 - 3, q)) / denom;
  const pNull = choose(5, q) / denom;
  assert(Math.abs(p1 + p2 + pNull - 1) < 1e-12, "Null kernel does not sum to one");
}

// The integrated leaf bound remains valid although the true indicator is not
// pointwise constant in a child split at t != 1/2.
for (const t of [0.49, 0.53]) {
  const n = 200000;
  const rows = [];
  const leafTotals = new Map();
  for (let i = 0; i < n; i += 1) {
    const x = (i + 0.5) / n;
    const b = x >= 0.5 ? 1 : 0;
    const bhat = x >= t ? 1 : 0;
    let leaf;
    if (x < t) leaf = x < 0.2 ? 0 : 1;
    else leaf = x < 0.7 ? 2 : 3;
    const z = leafTotals.get(leaf) || { n: 0, b: 0 };
    z.n += 1;
    z.b += b;
    leafTotals.set(leaf, z);
    rows.push({ b, bhat, leaf });
  }
  let lhs = 0;
  let mismatch = 0;
  for (const row of rows) {
    const z = leafTotals.get(row.leaf);
    lhs += Math.abs(z.b / z.n - row.b) / n;
    mismatch += Math.abs(row.bhat - row.b) / n;
  }
  assert(lhs <= 2 * mismatch + 2 / n, "Integrated leaf-contamination inequality failed");
}

// Adaptive resolution is an integrated statement. A first split at t leaves
// exactly |t - 1/2| units of weighted minority mass, and arbitrary descendant
// refinements cannot increase that mass or the integrated residual CART gain.
{
  const intervalMasses = (lo, hi) => ({
    zero: Math.max(0, Math.min(hi, 0.5) - lo),
    one: Math.max(0, hi - Math.max(lo, 0.5))
  });
  const minorityMass = intervals => intervals.reduce((sum, [lo, hi]) => {
    const z = intervalMasses(lo, hi);
    return sum + Math.min(z.zero, z.one);
  }, 0);
  const cartGain = (lo, hi, t, a) => {
    const left = intervalMasses(lo, t);
    const right = intervalMasses(t, hi);
    const nL = left.zero + left.one;
    const nR = right.zero + right.one;
    if (nL === 0 || nR === 0) return 0;
    const muL = left.one / nL;
    const muR = right.one / nR;
    return a * a * (nL / (hi - lo)) * (nR / (hi - lo))
      * (muL - muR) ** 2;
  };

  for (const t0 of [0.47, 0.49, 0.53, 0.58]) {
    const firstChildren = [[0, t0], [t0, 1]];
    assert(
      Math.abs(minorityMass(firstChildren) - Math.abs(t0 - 0.5)) < 1e-12,
      "First-split weighted impurity identity failed"
    );

    const cuts = [...new Set([0, t0, 0.12, 0.31, 0.495, 0.505, 0.71, 0.9, 1])]
      .sort((a, b) => a - b);
    const refined = cuts.slice(0, -1).map((lo, i) => [lo, cuts[i + 1]]);
    assert(
      minorityMass(refined) <= Math.abs(t0 - 0.5) + 1e-12,
      "Refinement increased weighted minority mass"
    );

    const coefficient = 0.35;
    let integratedMaxGain = 0;
    for (const [lo, hi] of refined) {
      let best = 0;
      for (let g = 1; g < 100; g += 1) {
        const t = lo + (hi - lo) * g / 100;
        best = Math.max(best, cartGain(lo, hi, t, coefficient));
      }
      integratedMaxGain += (hi - lo) * best;
    }
    assert(
      integratedMaxGain <= coefficient ** 2 * minorityMass(refined) + 1e-12,
      "Integrated residual-gain bound failed"
    );
  }
}

// A favorable route does not imply recovery when a competing transition goes
// to an unresolved trap. The corrected global hitting probability detects it.
{
  const G = 10000;
  const p = 1 / Math.sqrt(G);
  // States: initial, route predecessor, success, unresolved trap.
  const K = [
    [0, p, 0, 1 - p],
    [0, 1 - p, p, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
  ];
  const success = new Set([2]);
  const globalHit = Math.min(
    hitWithin(K, success, 0, 2),
    hitWithin(K, success, 1, 2),
    hitWithin(K, success, 3, 2)
  );
  assert(globalHit === 0, "Unresolved trap should force the global hitting bound to zero");
  const eventual = evolve([1, 0, 0, 0], K, 100000)[2];
  assert(eventual <= p + 1e-8, "Favorable-route counterexample unexpectedly recovered");
}

// Under closed progress, a uniform G^{-1/2} progress hazard gives the stated
// G/H^2 phase for finitely many required progress moves.
{
  const G = 10000;
  const p = 1 / Math.sqrt(G);
  const K = [
    [1, 0, 0],
    [p, 1 - p, 0],
    [0, p, 1 - p]
  ];
  const HRecovery = Math.floor(G ** 0.75);
  const recovery = evolve([0, 0, 1], K, HRecovery)[0];
  assert(recovery > 0.999, "Closed-progress recovery phase failed");

  const HOmission = Math.floor(G ** 0.25);
  const resolved = evolve([0, 1], [[1, 0], [p, 1 - p]], HOmission)[0];
  assert(resolved < 0.11, "Omission phase failed");
}

console.log("PASS: null kernel, adaptive resolution, leaf contamination, global hitting, and closed-progress phase checks");
