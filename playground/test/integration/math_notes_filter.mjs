#!/usr/bin/env node
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import { fileURLToPath } from 'node:url';
import { readFileSync } from 'node:fs';

const quarto = process.argv[2] ?? 'quarto';
const filter = fileURLToPath(new URL('./math_notes_filter.lua', import.meta.url));
const math = String.raw`$$L=\frac{\cssId{imag}{\Im(Z)}}{\omega}$$`;
const note = '::: {.lcm-math-note target="imag"}\n\n**Imaginary part**\n\nAn explanation with $\\omega$.\n\n:::';
const render = source => spawnSync(quarto, ['pandoc', '-f', 'markdown', '-t', 'revealjs',
  '--slide-level=2', '--mathjax', '--lua-filter=' + filter], { input: source, encoding: 'utf8' });
const source = '## Equation\n\n' + math + '\n\n' + note;
const result = render(source);
assert.equal(result.status, 0, result.stderr);
assert.match(result.stdout, /data-lcm-math-target="imag"/);
assert.match(result.stdout, /hidden=""/);
assert.match(result.stdout, /\\cssId\{imag\}/);
assert.match(result.stdout, /<strong>Imaginary part<\/strong>/);

const invalid = [
  [note, 'missing'],
  [source + '\n\n' + note, 'duplicate explanation'],
  [source + '\n\n' + math, 'duplicate equation anchor'],
  [source.replace('target="imag"', 'target="missing"'), 'missing'],
  [source.replace('\n\n' + note, '\n\n## Another slide\n\n' + note), 'same slide'],
  [source.replace('target="imag"', 'target="imag" trigger="hover"'), 'click-toggle'],
  [source.replace('**Imaginary part**\n\nAn explanation with $\\omega$.', ''), 'empty explanation'],
  [source.replace('**Imaginary part**', '<input>'), 'raw markup'],
  [source.replace('**Imaginary part**', '### Heading'), 'headings'],
  [source.replace('**Imaginary part**', '$\\cssId{nested}{x}$'), 'inside math notes'],
  ['## Equation {#imag}\n\n' + math + '\n\n' + note, 'duplicate equation anchor'],
  [source.replace('\\cssId{imag}', '\\cssId{invalid id}'), 'literal ID'],
  ['## Equation\n\n$$x % \\cssId{imag}{x}\n$$\n\n' + note, 'missing'],
];
for (const [input, message] of invalid) {
  const failure = render(input);
  assert.notEqual(failure.status, 0, `accepted invalid authoring: ${input}`);
  assert.ok(failure.stderr.includes(message), failure.stderr);
}
const lists = render('## Lists\n\n::: {.incremental}\n\n- Parent\n  - Child\n- Next\n\n:::\n\n' + source);
assert.equal(lists.status, 0, lists.stderr);
assert.equal((lists.stdout.match(/class="fragment"/g) ?? []).length, 3);
const guide = readFileSync(new URL('../../presentations/math-notes.qmd', import.meta.url), 'utf8');
const example = guide.match(/```\{\.markdown shortcodes=false\}\n([\s\S]*?)```/)[1];
const copied = render(example);
assert.equal(copied.status, 0, copied.stderr);
assert.match(copied.stdout, /data-lcm-math-target="impedance-imag"/);
assert.match(copied.stdout, /data-lcm-layout="full-canvas"/);
console.log(`Math note authoring contract passed (${invalid.length} invalid cases, valid markup, native nested fragments)`);
