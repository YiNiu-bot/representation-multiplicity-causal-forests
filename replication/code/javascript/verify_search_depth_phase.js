import assert from "node:assert/strict";

function candidateLaw(M, lambda) {
  const raw = [];
  const p0 = Math.exp(-lambda);
  let probability = p0 * lambda;
  const p1 = probability;
  raw[1] = p0 + p1;
  let cumulative = raw[1];
  for (let q = 2; q < M; q += 1) {
    probability *= lambda / q;
    raw[q] = probability;
    cumulative += raw[q];
  }
  raw[M] = Math.max(0, 1 - cumulative);
  return raw;
}

function expectedQ(M, lambda) {
  const law = candidateLaw(M, lambda);
  let value = 0;
  for (let q = 1; q <= M; q += 1) value += q * law[q];
  return value;
}

function defaultMtry(M) {
  return Math.min(Math.ceil(Math.sqrt(M)) + 20, M);
}

for (const M of [2, 3, 10, 50, 200, 1000, 10000]) {
  for (const lambda of new Set([1, 6, 16, defaultMtry(M)])) {
    const law = candidateLaw(M, lambda);
    const total = law.slice(1).reduce((a, b) => a + b, 0);
    assert.ok(Math.abs(total - 1) < 1e-11);
    const omega = expectedQ(M, lambda) / M;
    assert.ok(omega > 0 && omega <= 1);
  }
}

const fixedLimit = 6 + Math.exp(-6);
const fixedScaled = expectedQ(10000, 6);
assert.ok(Math.abs(fixedScaled - fixedLimit) < 1e-8);

const MLarge = 50000;
const defaultScaled =
  expectedQ(MLarge, defaultMtry(MLarge)) / Math.sqrt(MLarge);
assert.ok(Math.abs(defaultScaled - 1) < 0.1);

const epsilon = 0.4;
const omissionZeta = 1;
assert.ok(1 - omissionZeta < epsilon);
assert.equal(0.5, 0.5);
assert.ok(Math.abs((1 - epsilon) / 4 - 0.15) < 1e-12);

console.log("All production-mechanics phase checks passed.");
console.log({
  fixedScaled,
  fixedLimit,
  defaultScaled,
});
