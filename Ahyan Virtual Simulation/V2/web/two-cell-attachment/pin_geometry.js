(function (root) {
  "use strict";

  // Farthest point of an axis-aligned ellipse translated by (x,y).
  // In the reflected first quadrant the stationary equation is monotone:
  // a*x*tan(t) + (a*a-b*b)*sin(t) = b*y, for a >= b.
  function ellipseExtent(x, y, a, b) {
    if (![x, y, a, b].every(Number.isFinite) || a < b || b < 0) {
      throw new RangeError("Ellipse axes must be finite and a >= b >= 0");
    }
    const qx = Math.abs(x), qy = Math.abs(y);
    if (Math.abs(a - b) < 1e-14) {
      const length = Math.hypot(x, y);
      return { radius: length + a, x: length ? x * (1 + a / length) : a,
        y: length ? y * (1 + a / length) : 0 };
    }
    let lo = 0, hi = Math.PI / 2;
    if (qx === 0 && b * qy >= a * a - b * b) lo = hi;
    else for (let i = 0; i < 48; i += 1) {
      const t = (lo + hi) / 2;
      if (a * qx * Math.tan(t) + (a * a - b * b) * Math.sin(t) < b * qy) lo = t;
      else hi = t;
    }
    const t = (lo + hi) / 2;
    const px = (x < 0 ? -1 : 1) * (qx + a * Math.cos(t));
    const py = (y < 0 ? -1 : 1) * (qy + b * Math.sin(t));
    return { radius: Math.hypot(px, py), x: px, y: py };
  }

  const api = { ellipseExtent };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.RADPinGeometry = api;
})(typeof window === "undefined" ? globalThis : window);
