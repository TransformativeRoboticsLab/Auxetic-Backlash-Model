const assert = require("assert");
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const root = path.resolve(__dirname, "..");
const source = fs.readFileSync(path.join(root, "web", "two-cell-attachment", "two_cell_attachment.js"), "utf8");
const start = source.indexOf("function playOperator");
const end = source.indexOf("function resetContactMemory");
assert.ok(start >= 0 && end > start, "could not locate hysteresis function block");

const context = { Math };
vm.createContext(context);
vm.runInContext(`${source.slice(start, end)}
this.playOperator = playOperator;
this.playTrace = playTrace;
this.backlashHysteresisDiagnostic = backlashHysteresisDiagnostic;
this.contactSide = contactSide;`, context);

const fieldStart = source.indexOf("function neighborCells");
const fieldEnd = source.indexOf("function contactDriveField");
assert.ok(fieldStart >= 0 && fieldEnd > fieldStart, "could not locate play-field propagation block");
vm.runInContext(`${source.slice(fieldStart, fieldEnd)}
this.propagatePlayField = propagatePlayField;`, context);
const plain = (value) => JSON.parse(JSON.stringify(value));

assert.strictEqual(context.playOperator(0, 0, 2), 0, "follower should float inside centered clearance");
assert.strictEqual(context.playOperator(5, 0, 2), 3, "positive motion should transmit after clearance");
assert.strictEqual(context.playOperator(-5, 3, 2), -3, "reversal should cross the opposite clearance wall");

const forward = context.playTrace([0, 1, 2, 3, 4, 5], 2, 0);
assert.deepStrictEqual(
  forward.map((entry) => entry.follower),
  [0, 0, 0, 1, 2, 3],
  "forward trace should delay follower until clearance is taken up"
);

const reverse = context.playTrace([5, 4, 3, 2, 1, 0, -1], 2, 3);
assert.deepStrictEqual(
  reverse.map((entry) => entry.follower),
  [3, 3, 3, 3, 3, 2, 1],
  "reverse trace should retain prior contact until opposite clearance is taken up"
);

const diagnostic = context.backlashHysteresisDiagnostic({ backlashDeg: 2, amplitudeDeg: 8, steps: 32 });
assert.ok(diagnostic.hysteretic, "diagnostic should detect a nonzero hysteresis loop");
assert.ok(diagnostic.loopGap > 1, "forward and reverse paths should differ at the same driver angle");
assert.notStrictEqual(
  Math.sign(diagnostic.downMid.follower),
  Math.sign(diagnostic.upMid.follower),
  "opposite traversal directions should land on opposite sides of the backlash loop"
);

const emptyString = Array.from({ length: 1 }, () => Array.from({ length: 15 }, () => 0));
const forwardField = context.propagatePlayField({
  rows: 1,
  cols: 15,
  sourceRow: 0,
  sourceCol: 0,
  sourceDelta: 5,
  previousDeltaField: emptyString,
  backlash: 2,
});
assert.deepStrictEqual(
  plain(forwardField.deltaField[0].slice(0, 6)),
  [5, 3, 1, 0, 0, 0],
  "1x15 field should transmit only after each local clearance is taken up"
);
assert.deepStrictEqual(
  plain(Object.keys(forwardField.connectionField).sort()),
  ["0,0|0,1", "0,1|0,2"],
  "field propagation should store backlash memory on physical neighbor connections"
);
assert.deepStrictEqual(
  plain(forwardField.connectionSideField),
  { "0,0|0,1": 1, "0,1|0,2": 1 },
  "field propagation should expose physical neighbor connection contact direction"
);

const reverseToNeutral = context.propagatePlayField({
  rows: 1,
  cols: 15,
  sourceRow: 0,
  sourceCol: 0,
  sourceDelta: 0,
  previousDeltaField: forwardField.deltaField,
  previousConnectionField: forwardField.connectionField,
  backlash: 2,
});
assert.deepStrictEqual(
  plain(reverseToNeutral.deltaField[0].slice(0, 6)),
  [0, 2, 1, 0, 0, 0],
  "returning the source to neutral should update through the opposite clearance, not freeze the old field"
);
assert.notDeepStrictEqual(
  plain(reverseToNeutral.deltaField[0].slice(0, 6)),
  plain(context.propagatePlayField({
    rows: 1,
    cols: 15,
    sourceRow: 0,
    sourceCol: 0,
    sourceDelta: 0,
    previousDeltaField: emptyString,
    backlash: 2,
  }).deltaField[0].slice(0, 6)),
  "the same source angle should produce a different lattice state when reached from a different history"
);

function zeroField(rows, cols) {
  return Array.from({ length: rows }, () => Array.from({ length: cols }, () => 0));
}

function commandSweep(min, max, step) {
  const up = [];
  for (let value = min; value <= max + 1e-9; value += step) up.push(Number(value.toFixed(6)));
  const down = [];
  for (let value = max - step; value >= min - 1e-9; value -= step) down.push(Number(value.toFixed(6)));
  return up.concat(down);
}

function maxFieldJump(a, b) {
  let jump = 0;
  for (let row = 0; row < a.length; row += 1) {
    for (let col = 0; col < a[row].length; col += 1) {
      jump = Math.max(jump, Math.abs((a[row][col] || 0) - (b[row][col] || 0)));
    }
  }
  return jump;
}

function simulateString({ rows, cols, sourceRow, sourceCol, backlash, commands }) {
  let deltaField = zeroField(rows, cols);
  let connectionField = {};
  const frames = [];
  commands.forEach((sourceDelta) => {
    const result = context.propagatePlayField({
      rows,
      cols,
      sourceRow,
      sourceCol,
      sourceDelta,
      previousDeltaField: deltaField,
      previousConnectionField: connectionField,
      backlash,
    });
    frames.push(result);
    deltaField = result.deltaField;
    connectionField = result.connectionField;
  });
  return frames;
}

const backlashSweep = Array.from({ length: 21 }, (_, index) => Number((index * 0.25).toFixed(2)));
const commands = commandSweep(-6, 6, 0.5);
[
  { rows: 1, cols: 15, sourceRow: 0, sourceCol: 0, name: "1x15" },
  { rows: 15, cols: 1, sourceRow: 0, sourceCol: 0, name: "15x1" },
].forEach((shape) => {
  backlashSweep.forEach((backlash) => {
    const frames = simulateString({ ...shape, backlash, commands });
    for (let index = 1; index < frames.length; index += 1) {
      const jump = maxFieldJump(frames[index - 1].deltaField, frames[index].deltaField);
      assert.ok(
        jump <= 0.500001,
        `${shape.name} backlash ${backlash} has discontinuous field jump ${jump} at frame ${index}`
      );
    }
    if (backlash > 0) {
      const freshNeutral = context.propagatePlayField({
        rows: shape.rows,
        cols: shape.cols,
        sourceRow: shape.sourceRow,
        sourceCol: shape.sourceCol,
        sourceDelta: 0,
        previousDeltaField: zeroField(shape.rows, shape.cols),
        backlash,
      });
      const returnedNeutral = frames[frames.length - 1];
      assert.notDeepStrictEqual(
        plain(returnedNeutral.deltaField),
        plain(freshNeutral.deltaField),
        `${shape.name} backlash ${backlash} should retain path history at the same neutral command`
      );
    }
  });
});

console.log("two-cell hysteresis validation passed");
