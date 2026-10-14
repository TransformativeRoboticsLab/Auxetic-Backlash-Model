const assert = require('assert');
const { advance } = require('../web/two-cell-attachment/backlash_state.js');
const rest = [Array(10).fill(40)];
let field = advance(rest, 0, 0, 41, 2);
assert.strictEqual(field[0][1], 40, 'neighbor waits inside clearance');
field = advance(field, 0, 0, 60, 2);
assert.strictEqual(field[0][1], 58);
const reversed = advance(field, 0, 0, 59, 2);
assert.strictEqual(reversed[0][1], 58, 'reversal must not immediately drag neighbor');
assert.strictEqual(advance(rest, 0, 0, 59, 2)[0][1], 57, 'same command has a different loading path');
assert.deepStrictEqual(advance(field, 0, 6, field[0][6], 2), field, 'selecting another cell does not move it');
assert.ok(advance(field, 0, 3, 51, 0)[0].every(a => a === 51));
assert.ok(rest[0].every(a => a === 40), 'trial must not mutate accepted state');
let sheet = Array.from({length:5}, () => Array(5).fill(40));
const rotate = a => a[0].map((_,c) => a.map(row=>row[c]).reverse());
for (const [r,c,q] of [[1,2,55],[4,0,32],[3,3,57],[0,1,42]]) {
  const next = advance(sheet,r,c,q,2);
  assert.deepStrictEqual(rotate(next),advance(rotate(sheet),c,4-r,q,2));
  for(let i=0;i<5;i++) for(let j=0;j<5;j++) {
    if(i) assert.ok(Math.abs(next[i][j]-next[i-1][j])<=2);
    if(j) assert.ok(Math.abs(next[i][j]-next[i][j-1])<=2);
  }
  sheet=next;
}
console.log('Backlash reversal, sequential selection, zero clearance, and quarter-turn tests passed');
