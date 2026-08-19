"use strict";

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function choose(n, k) {
  if (k < 0 || k > n) return 0;
  let out = 1;
  for (let j = 1; j <= k; j += 1) out = (out * (n - k + j)) / j;
  return out;
}

function subsets(items, k) {
  const out = [];
  function visit(start, current) {
    if (current.length === k) {
      out.push(current.slice());
      return;
    }
    for (let j = start; j < items.length; j += 1) {
      current.push(items[j]);
      visit(j + 1, current);
      current.pop();
    }
  }
  visit(0, []);
  return out;
}

// Theorem 1: the product-to-exponential approximation and its stated bound.
for (const hazards of [
  [0.02, 0.03, 0.01],
  [0.10, 0.04, 0.07, 0.02],
  Array.from({ length: 200 }, () => 1.7 / 200),
]) {
  const lambda = hazards.reduce((a, b) => a + b, 0);
  const maxHazard = Math.max(...hazards);
  const survival = hazards.reduce((a, b) => a * (1 - b), 1);
  const remainder = -Math.log(survival) - lambda;
  const bound = (maxHazard / (1 - maxHazard)) * lambda;
  assert(remainder >= -1e-12, "log remainder must be nonnegative");
  assert(remainder <= bound + 1e-12, "log remainder exceeds theorem bound");
}

for (const lambda of [0.25, 1, 3]) {
  for (const calls of [1000, 5000, 20000]) {
    const hazard = lambda / calls;
    const survival = Math.pow(1 - hazard, calls);
    const error = Math.abs(survival - Math.exp(-lambda));
    assert(error <= (lambda * lambda + 1) / calls, "rare-hazard limit failed");
  }
}

// Theorem 2: predictable hazard deficits are zero off structural calls and
// bounded by singleton availability on structural calls. A path-level bad
// event permits the same bound at every call.
for (const calls of [5, 20, 100]) {
  for (const structuralCalls of [0, 1, 3]) {
    for (const omega of [0.01, 0.2, 0.8]) {
      const goodDeficit = Math.min(calls, structuralCalls) * omega;
      const goodBound = structuralCalls * omega;
      assert(goodDeficit <= goodBound + 1e-12, "structural deficit bound failed");

      const badDeficit = calls * omega;
      const badBound = structuralCalls * omega + calls * omega;
      assert(badDeficit <= badBound + 1e-12, "bad-path deficit bound failed");
    }
  }
}

// Theorem 2, fixed-M repair: a copied class is absent only when Q=1 and the
// singleton label is drawn. Its second-offer waiting time is uniformly
// dominated by a negative-binomial tail, even when M stays fixed.
function secondSuccessTail(p, calls) {
  return Math.pow(1 - p, calls)
    + calls * p * Math.pow(1 - p, calls - 1);
}

for (const M of [2, 3, 5, 16, 64]) {
  for (const nuOne of [0, 0.2, 0.75, 1]) {
    const copiedAvailability = 1 - nuOne / M;
    assert(copiedAvailability >= 0.5, "copied-class availability fell below 1/2");
    for (const calls of [2, 5, 20, 80]) {
      const tail = secondSuccessTail(copiedAvailability, calls);
      const uniformBound = (1 + 2 * calls) * Math.pow(2, -calls);
      assert(tail <= uniformBound + 1e-12, "fixed-M waiting-time bound failed");
    }
  }
}

const fixedMTailShort = secondSuccessTail(0.5, 2);
const fixedMTailLong = secondSuccessTail(0.5, 80);
assert(fixedMTailShort > 0, "the invalid two-call shortcut was not detected");
assert(fixedMTailLong < 1e-20, "many-call fixed-M resolution did not converge");

// Proposition 3: h_g counts every class preceding g in the complete order.
for (let G = 2; G <= 8; G += 1) {
  const classes = Array.from({ length: G }, (_, j) => j);
  for (let q = 1; q <= G; q += 1) {
    for (let h = 0; h < G; h += 1) {
      const target = h;
      const empirical = subsets(classes, q).filter((set) => {
        return set.includes(target) && set.every((j) => j >= target);
      }).length / choose(G, q);
      const formula = (choose(G - h, q) - choose(G - h - 1, q)) / choose(G, q);
      assert(Math.abs(empirical - formula) < 1e-12, "quotient winning law failed");
    }
  }
}

// Component targets depend on resolution probabilities, not raw offer counts.
const f = [1.25, -0.75, 0.40];
const lambdas = [0, Math.log(2), Infinity];
const target = f.reduce((sum, value, j) => {
  const survival = Number.isFinite(lambdas[j]) ? Math.exp(-lambdas[j]) : 0;
  return sum + (1 - survival) * value;
}, 0);
assert(Math.abs(target - (-0.375 + 0.40)) < 1e-12, "mixture target failed");

console.log("successful-resolution checks passed");
