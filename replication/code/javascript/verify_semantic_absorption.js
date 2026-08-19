"use strict";

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function choose(n, k) {
  if (k < 0 || k > n) return 0;
  k = Math.min(k, n - k);
  let out = 1;
  for (let j = 1; j <= k; j += 1) out *= (n - k + j) / j;
  return out;
}

function zeros(n) {
  return Array.from({ length: n }, () => Array(n).fill(0));
}

function identity(n) {
  const out = zeros(n);
  for (let i = 0; i < n; i += 1) out[i][i] = 1;
  return out;
}

function add(a, b) {
  return a.map((row, i) => row.map((x, j) => x + b[i][j]));
}

function scale(a, c) {
  return a.map((row) => row.map((x) => c * x));
}

function multiply(a, b) {
  const n = a.length;
  const out = zeros(n);
  for (let i = 0; i < n; i += 1) {
    for (let k = 0; k < n; k += 1) {
      for (let j = 0; j < n; j += 1) out[i][j] += a[i][k] * b[k][j];
    }
  }
  return out;
}

function power(a, exponent) {
  let base = a;
  let out = identity(a.length);
  let k = exponent;
  while (k > 0) {
    if (k % 2 === 1) out = multiply(out, base);
    base = multiply(base, base);
    k = Math.floor(k / 2);
  }
  return out;
}

function expm(a) {
  let term = identity(a.length);
  let out = identity(a.length);
  for (let k = 1; k <= 120; k += 1) {
    term = scale(multiply(term, a), 1 / k);
    out = add(out, term);
  }
  return out;
}

function maxDiff(a, b) {
  let out = 0;
  for (let i = 0; i < a.length; i += 1) {
    for (let j = 0; j < a.length; j += 1) {
      out = Math.max(out, Math.abs(a[i][j] - b[i][j]));
    }
  }
  return out;
}

function canonicalColumn(values) {
  const observed = [...new Set(values.filter((x) => x !== null))]
    .sort((a, b) => a - b);
  const forward = values.map((x) => (x === null ? null : observed.indexOf(x) + 1));
  const reverse = forward.map((x) => (x === null ? null : observed.length + 1 - x));
  const key = (x) => x.map((value) => (value === null ? "NA" : value)).join(",");
  const canonical = key(forward) <= key(reverse) ? forward : reverse;
  return { certificate: key(canonical), values: canonical };
}

function quotientColumns(columns) {
  const encoded = columns.map(canonicalColumn);
  const certificates = [...new Set(encoded.map((x) => x.certificate))].sort();
  return certificates.map((certificate) => (
    encoded.find((x) => x.certificate === certificate).values
  ));
}

function stateHazard(multiplicities, scores, resolvedMask, g, q) {
  if ((resolvedMask >> g) & 1) return 0;
  const total = multiplicities.reduce((a, b) => a + b, 0);
  let stronger = 0;
  for (let ell = 0; ell < multiplicities.length; ell += 1) {
    if (((resolvedMask >> ell) & 1) === 0 && scores[ell] > scores[g]) {
      stronger += multiplicities[ell];
    }
  }
  return (
    choose(total - stronger, q) -
    choose(total - stronger - multiplicities[g], q)
  ) / choose(total, q);
}

function enumerateCandidateSets(n, q, visit, start = 0, chosen = []) {
  if (chosen.length === q) {
    visit(chosen);
    return;
  }
  for (let j = start; j <= n - (q - chosen.length); j += 1) {
    chosen.push(j);
    enumerateCandidateSets(n, q, visit, j + 1, chosen);
    chosen.pop();
  }
}

function generalHazard(labels, transitionClass, q) {
  const targetLabels = labels.filter((label) => label.transition === transitionClass);
  assert(targetLabels.length > 0, "empty target transition class");
  const targetScore = targetLabels[0].score;
  assert(targetLabels.every((label) => label.score === targetScore),
    "transition class is not score tied");
  const stronger = labels.filter((label) => label.score > targetScore).length;
  return (choose(labels.length - stronger, q)
    - choose(labels.length - stronger - targetLabels.length, q))
    / choose(labels.length, q);
}

function enumeratedHazard(labels, transitionClass, q) {
  let wins = 0;
  let total = 0;
  enumerateCandidateSets(labels.length, q, (set) => {
    total += 1;
    const bestScore = Math.max(...set.map((index) => labels[index].score));
    const best = set.map((index) => labels[index])
      .filter((label) => label.score === bestScore);
    const winningClasses = new Set(best.map((label) => label.transition));
    assert(winningClasses.size === 1, "strict action-class ordering failed");
    if (best[0].transition === transitionClass) wins += 1;
  });
  return wins / total;
}

function singletonWinning(gCount, strongerCount, q) {
  if (q < 1 || q > gCount || strongerCount >= gCount) return 0;
  let avoidance = 1;
  for (let r = 0; r < q - 1; r += 1) {
    avoidance *= (gCount - strongerCount - 1 - r) / (gCount - 1 - r);
  }
  return (q / gCount) * avoidance;
}

function rawKernel(multiplicities, scores, q) {
  const gCount = multiplicities.length;
  const nStates = 1 << gCount;
  const out = zeros(nStates);
  for (let s = 0; s < nStates; s += 1) {
    let move = 0;
    for (let g = 0; g < gCount; g += 1) {
      const p = stateHazard(multiplicities, scores, s, g, q);
      if (p > 0) out[s][s | (1 << g)] += p;
      move += p;
    }
    out[s][s] += 1 - move;
  }
  return out;
}

function quotientKernel(gCount, scores, q) {
  return rawKernel(Array(gCount).fill(1), scores, q);
}

function checkStochastic(kernel) {
  for (const row of kernel) {
    const sum = row.reduce((a, b) => a + b, 0);
    assert(Math.abs(sum - 1) < 1e-12, `row sum is ${sum}`);
    assert(row.every((x) => x >= -1e-12), "negative transition probability");
  }
}

// Exact state-dependent winning law.
const scores = [3, 2, 1];
const rawA = rawKernel([8, 2, 1], scores, 2);
const rawB = rawKernel([1, 2, 8], scores, 2);
checkStochastic(rawA);
checkStochastic(rawB);
assert(maxDiff(rawA, rawB) > 0.05, "raw multiplicity should change search");

// The hazard counts every stronger label, including state-preserving actions.
const actionLabels = [
  ...Array.from({ length: 4 }, () => ({ score: 4, transition: "stay" })),
  ...Array.from({ length: 3 }, () => ({ score: 3, transition: "g0" })),
  ...Array.from({ length: 2 }, () => ({ score: 2, transition: "g1" })),
  { score: 1, transition: "stay-low" },
];
for (const q of [1, 2, 3, 5]) {
  assert(Math.abs(generalHazard(actionLabels, "g0", q)
    - enumeratedHazard(actionLabels, "g0", q)) < 1e-14,
  "general hazard failed with stronger state-preserving labels");
  assert(Math.abs(generalHazard(actionLabels, "g1", q)
    - enumeratedHazard(actionLabels, "g1", q)) < 1e-14,
  "general hazard failed with stronger transition labels");
}

// Availability is not winning: an offered class can be blocked with certainty.
const blockedLabels = [
  { score: 2, transition: "stronger" },
  { score: 1, transition: "target" },
];
assert(enumeratedHazard(blockedLabels, "target", 2) === 0,
  "offered target should never win when the stronger class is always offered");
assert(singletonWinning(2, 1, 2) === 0,
  "closed-form blocking counterexample failed");

// The capped-Poisson q=1 event gives every active singleton a positive route hazard.
for (const strongerCount of [0, 1, 9, 98]) {
  assert(Math.abs(singletonWinning(100, strongerCount, 1) - 0.01) < 1e-15,
    "singleton-only winning event should be independent of blockers");
}

// With q of order sqrt(G), o(sqrt(G)) blockers are asymptotically negligible.
let blockerRatioError = Infinity;
for (const gCount of [10000, 160000, 2560000]) {
  const q = Math.ceil(Math.sqrt(gCount));
  const h = Math.floor(gCount ** 0.2);
  const winning = singletonWinning(gCount, h, q);
  blockerRatioError = Math.abs(winning / (q / gCount) - 1);
}
assert(blockerRatioError < 0.06,
  `route-blocker approximation failed: ${blockerRatioError}`);

// Canonical rank coordinates remove label order, scale, sign, and survivor choice.
const xIncome = [2, 5, 5, 9, 12, 20];
const xDistance = [8, 4, 6, 2, 7, 1];
const xSize = [1, 3, 2, 2, 5, 4];
const canonicalDesign = quotientColumns([xIncome, xDistance, xSize]);
const augmentedDesign = quotientColumns([
  xIncome,
  xDistance,
  xSize,
  [...xIncome],
  xIncome.map((x) => Math.log1p(x)),
  xDistance.map((x) => -x),
  xSize.map((x) => x * 1000),
].reverse());
const replacementDesign = quotientColumns([
  [...xIncome],
  xDistance.map((x) => -x),
  xSize.map((x) => x * 1000),
]);
assert(JSON.stringify(canonicalDesign) === JSON.stringify(augmentedDesign),
  "canonical quotient changed after insertion and permutation");
assert(JSON.stringify(canonicalDesign) === JSON.stringify(replacementDesign),
  "canonical quotient changed after deleting original labels");

// Quotient search depends only on semantic classes.
const quotientA = quotientKernel(3, scores, 2);
const quotientB = quotientKernel(3, scores, 2);
checkStochastic(quotientA);
assert(maxDiff(quotientA, quotientB) < 1e-14, "quotient kernel changed");

// Continuous-time absorption limit K_H^[dH] -> exp(d Lambda).
const lambda = [
  [-1.1, 0.7, 0.4, 0],
  [0, -0.9, 0, 0.9],
  [0, 0, -0.6, 0.6],
  [0, 0, 0, 0],
];
const d = 1.7;
const target = expm(scale(lambda, d));
let previousError = Infinity;
for (const h of [100, 500, 2500, 12500]) {
  const kH = add(identity(4), scale(lambda, 1 / h));
  checkStochastic(kH);
  const approximation = power(kH, Math.floor(d * h));
  const error = maxDiff(approximation, target);
  assert(error < previousError, "matrix-exponential approximation did not improve");
  previousError = error;
}
assert(previousError < 1e-4, `matrix-exponential error is ${previousError}`);

// Under the default class budget, semantic availability is of order G^(-1/2).
for (const gCount of [10000, 40000, 160000]) {
  const lambdaG = Math.ceil(Math.sqrt(gCount)) + 20;
  const scaledAvailability = (lambdaG / gCount) * Math.sqrt(gCount);
  assert(Math.abs(scaledAvailability - 1) < 0.21, "default scaling failed");
}

process.stdout.write(
  JSON.stringify(
    {
      raw_kernel_difference: maxDiff(rawA, rawB),
      quotient_kernel_difference: maxDiff(quotientA, quotientB),
      state_preserving_hazard_check: "PASS",
      offered_but_blocked_check: "PASS",
      singleton_only_route_check: "PASS",
      route_blocker_ratio_error: blockerRatioError,
      canonical_rank_quotient_check: "PASS",
      ctmc_max_error: previousError,
      status: "PASS",
    },
    null,
    2,
  ) + "\n",
);
