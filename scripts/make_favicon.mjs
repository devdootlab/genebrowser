// Rasterise favicon.svg's avocado to a real 32x32 PNG, for browsers that will not take an SVG
// favicon. Written because the alternative was pasting a base64 blob nobody had rendered -- an
// icon you have not looked at is indistinguishable from a corrupt one until it ships.
//
// No canvas, no image library: the shape is drawn from its own geometry into an RGBA buffer and
// encoded with zlib, which is in Node core. 4x4 supersampling for edges that survive at tab size.
//
//   node scripts/make_favicon.mjs        -> favicon.png (+ prints the data URI for the <link>)
import fs from 'node:fs';
import zlib from 'node:zlib';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const N = 32, SS = 4;                       // output size, supersample factor

// --- the same geometry as favicon.svg, in its 64-unit space ------------------------------------
// half-width of the avocado body at vertical position y, narrow at the top and fat at the bottom
function halfWidth(y) {
  const u = (y - 4) / 56;                   // 0 at the tip, 1 at the base
  if (u < 0 || u > 1) return 0;
  return 17 * Math.sqrt(Math.max(0, 1 - (2 * u - 1) ** 2)) * (0.5 + 0.6 * u);
}
const SKIN = [0x3f, 0x6b, 0x34], FLESH = [0xcf, 0xe0, 0x8d], INNER = [0xe6, 0xef, 0xb8],
      PIT = [0x8a, 0x5a, 0x30], RIM = [0x6d, 0x44, 0x23], LIT = [0xff, 0xff, 0xff],
      DOT = [0x2b, 0x1a, 0x0d];

function sample(x, y) {                     // -> [r,g,b,a] at a point in the 64-unit space
  const inPit = (x - 32) ** 2 + (y - 37) ** 2;
  if (inPit <= 3.2 ** 2) {
    const l = (x - 27.5) ** 2 + (y - 32.5) ** 2;
    if (l <= 1.6 ** 2) return [...DOT, 255];
  }
  if ((x - 27.5) ** 2 + (y - 32.5) ** 2 <= 4.2 ** 2 && inPit <= 12.5 ** 2) {
    const l = (x - 27.5) ** 2 + (y - 32.5) ** 2;
    return l <= 1.6 ** 2 ? [...DOT, 255] : [...LIT, 217];
  }
  if (inPit <= 12.5 ** 2) return inPit >= 11.5 ** 2 ? [...RIM, 255] : [...PIT, 255];
  const hw = halfWidth(y), d = Math.abs(x - 32);
  if (hw <= 0 || d > hw) return [0, 0, 0, 0];
  if (d <= hw * 0.62) return [...INNER, 255];
  if (d <= hw * 0.80) return [...FLESH, 255];
  return [...SKIN, 255];
}

// --- rasterise with supersampling ---------------------------------------------------------------
const px = Buffer.alloc(N * N * 4);
for (let py = 0; py < N; py++) {
  for (let pxx = 0; pxx < N; pxx++) {
    let r = 0, g = 0, b = 0, a = 0;
    for (let sy = 0; sy < SS; sy++) for (let sx = 0; sx < SS; sx++) {
      const gx = (pxx + (sx + 0.5) / SS) * (64 / N);
      const gy = (py + (sy + 0.5) / SS) * (64 / N);
      const c = sample(gx, gy);
      r += c[0] * c[3]; g += c[1] * c[3]; b += c[2] * c[3]; a += c[3];
    }
    const n = SS * SS, o = (py * N + pxx) * 4;
    px[o]     = a ? Math.round(r / a) : 0;    // un-premultiply
    px[o + 1] = a ? Math.round(g / a) : 0;
    px[o + 2] = a ? Math.round(b / a) : 0;
    px[o + 3] = Math.round(a / n);
  }
}

// --- PNG container -------------------------------------------------------------------------------
const CRC = (() => { const t = []; for (let n = 0; n < 256; n++) { let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
const crc32 = buf => { let c = 0xFFFFFFFF;
  for (const b of buf) c = CRC[(c ^ b) & 0xFF] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; };
function chunk(type, data) {
  const len = Buffer.alloc(4); len.writeUInt32BE(data.length);
  const td = Buffer.concat([Buffer.from(type, 'ascii'), data]);
  const c = Buffer.alloc(4); c.writeUInt32BE(crc32(td));
  return Buffer.concat([len, td, c]);
}
const raw = Buffer.alloc(N * (N * 4 + 1));
for (let y = 0; y < N; y++) {                       // filter byte 0 (None) per scanline
  raw[y * (N * 4 + 1)] = 0;
  px.copy(raw, y * (N * 4 + 1) + 1, y * N * 4, (y + 1) * N * 4);
}
const ihdr = Buffer.alloc(13);
ihdr.writeUInt32BE(N, 0); ihdr.writeUInt32BE(N, 4);
ihdr[8] = 8; ihdr[9] = 6; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;   // 8-bit RGBA
const png = Buffer.concat([
  Buffer.from([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A]),
  chunk('IHDR', ihdr),
  chunk('IDAT', zlib.deflateSync(raw, { level: 9 })),
  chunk('IEND', Buffer.alloc(0))
]);
fs.writeFileSync(R + 'favicon.png', png);
console.log(`wrote favicon.png  ${N}x${N}  ${png.length} bytes`);
console.log(`data:image/png;base64,${png.toString('base64')}`);
