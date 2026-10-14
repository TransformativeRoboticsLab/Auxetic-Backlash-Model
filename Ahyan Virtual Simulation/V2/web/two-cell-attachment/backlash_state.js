(function (root) {
  "use strict";

  // Rate-independent play approximation, not a rigid-body force solver.
  // Previously accepted passive angles stay put while the clearance permits it.
  function advance(previous, sourceRow, sourceCol, command, gap) {
    if (!Number.isFinite(command) || !Number.isFinite(gap) || gap < 0 || !previous[sourceRow] || !Number.isFinite(previous[sourceRow][sourceCol])) {
      throw new RangeError("Invalid backlash state or command");
    }
    return previous.map((row, r) => row.map((angle, c) => {
      const radius = gap * (Math.abs(r - sourceRow) + Math.abs(c - sourceCol));
      return Math.max(command - radius, Math.min(command + radius, angle));
    }));
  }
  const api = { advance };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.RADBacklashState = api;
})(typeof window === "undefined" ? globalThis : window);
