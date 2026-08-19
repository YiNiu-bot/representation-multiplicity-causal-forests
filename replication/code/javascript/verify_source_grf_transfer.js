"use strict";

// Dependency-free checks for the exact identities used in Theorem 1.

let state = 20260808 >>> 0;
function uniform() {
  state = (1664525 * state + 1013904223) >>> 0;
  return state / 4294967296;
}

function assertClose(actual, expected, tolerance, label) {
  const error = Math.abs(actual - expected);
  if (!Number.isFinite(actual) || error > tolerance) {
    throw new Error(`${label}: error=${error}, actual=${actual}, expected=${expected}`);
  }
  return error;
}

let maxResponseError = 0;
let maxCriterionError = 0;
let maxPredictionError = 0;
let maxInformationError = 0;

for (let trial = 0; trial < 10000; trial += 1) {
  const n = 8 + Math.floor(24 * uniform());
  const z = Array.from({ length: n }, (_, i) =>
    i === 0 ? -0.5 : i === 1 ? 0.5 : (uniform() < 0.5 ? -0.5 : 0.5));
  const y = Array.from({ length: n }, () => 4 * uniform() - 2);
  const mean = (values) => values.reduce((a, b) => a + b, 0) / values.length;
  const delta = mean(z);
  const eta = mean(y);
  const zy = mean(z.map((value, i) => value * y[i]));
  const variance = 0.25 - delta * delta;
  const theta = (zy - delta * eta) / variance;
  const response = z.map((value, i) =>
    (value - delta) * (y[i] - eta - theta * (value - delta)));

  for (let i = 0; i < n; i += 1) {
    const gamma = 4 * z[i] * y[i];
    const remainder = -4 * z[i] * eta - 4 * delta * y[i]
      + 4 * delta * eta + 8 * theta * z[i] * delta
      - 4 * theta * delta * delta;
    maxResponseError = Math.max(
      maxResponseError,
      assertClose(4 * response[i], gamma - theta + remainder, 2e-12,
        "source relabeling identity")
    );
  }

  const nLeft = 2 + Math.floor((n - 3) * uniform());
  const sumLeft = response.slice(0, nLeft).reduce((a, b) => a + b, 0);
  const sumRight = response.slice(nLeft).reduce((a, b) => a + b, 0);
  assertClose(sumLeft + sumRight, 0, 2e-11, "parent response centering");
  const sourceCriterion = (16 / n) *
    (sumLeft * sumLeft / nLeft + sumRight * sumRight / (n - nLeft));
  const leftMean = 4 * sumLeft / nLeft;
  const rightMean = 4 * sumRight / (n - nLeft);
  const cartCriterion = (nLeft / n) * ((n - nLeft) / n) *
    (leftMean - rightMean) ** 2;
  maxCriterionError = Math.max(
    maxCriterionError,
    assertClose(sourceCriterion, cartCriterion, 3e-11,
      "source/CART criterion identity")
  );

  const weightsRaw = Array.from({ length: n }, () => 0.1 + uniform());
  const weightSum = weightsRaw.reduce((a, b) => a + b, 0);
  const alpha = weightsRaw.map((value) => value / weightSum);
  const barZ = alpha.reduce((sum, value, i) => sum + value * z[i], 0);
  const barY = alpha.reduce((sum, value, i) => sum + value * y[i], 0);
  const aMoment = alpha.reduce((sum, value, i) => sum + value * z[i] * y[i], 0);
  const sourcePrediction = (aMoment - barZ * barY) / (0.25 - barZ * barZ);
  const linearPrediction = 4 * aMoment;
  const predictionRemainder = 4 * barZ * (4 * aMoment * barZ - barY)
    / (1 - 4 * barZ * barZ);
  maxPredictionError = Math.max(
    maxPredictionError,
    assertClose(sourcePrediction - linearPrediction, predictionRemainder, 2e-12,
      "source prediction identity")
  );

  const informationDirect = z.reduce((sum, value) => sum + value * value, 0)
    - z.reduce((sum, value) => sum + value, 0) ** 2 / n;
  const informationClosed = n * (0.25 - delta * delta);
  maxInformationError = Math.max(
    maxInformationError,
    assertClose(informationDirect, informationClosed, 2e-12,
      "binary treatment-information identity")
  );
}

process.stdout.write(JSON.stringify({
  trials: 10000,
  maxResponseError,
  maxCriterionError,
  maxPredictionError,
  maxInformationError,
  status: "PASS"
}, null, 2) + "\n");
