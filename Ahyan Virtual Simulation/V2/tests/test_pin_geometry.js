const assert = require('assert');
const { ellipseExtent } = require('../web/two-cell-attachment/pin_geometry.js');
let seed = 20261014;
const random = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296);
for (let i = 0; i < 500; i++) {
  const x = i % 7 ? 4 * random() - 2 : 0;
  const y = i % 11 ? 4 * random() - 2 : 0;
  const b = .1 + 2 * random(), a = b + (i % 13 ? 3 * random() : 0);
  const exact = ellipseExtent(x, y, a, b);
  let sampled = 0;
  for (let j = 0; j < 8192; j++) {
    const t = j * 2 * Math.PI / 8192;
    sampled = Math.max(sampled, Math.hypot(x + a * Math.cos(t), y + b * Math.sin(t)));
  }
  assert.ok(exact.radius >= sampled - 1e-12, 'maximum must contain sampled ellipse');
  assert.ok(exact.radius - sampled < 2e-6, 'analytical maximum agrees with dense independent samples');
  assert.ok(Math.abs(ellipseExtent(-x, -y, a, b).radius - exact.radius) < 1e-12);
}
assert.throws(() => ellipseExtent(0, 0, 1, 2), RangeError);
console.log('500 exact ellipse containment cases passed');
