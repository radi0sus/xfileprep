// laue-3d.js
// 3D companion to the layer viewer ("Laue" tab): the whole measured
// reciprocal lattice as a point cloud inside the resolution sphere
// (radius d*max = 2·(sinθ/λ)max), drawn in the TRUE reciprocal metric.
// The layer currently shown by the layer viewer is drawn as a plane
// (disc = plane ∩ sphere) and its points are highlighted, so paging
// through the layers moves the plane through the sphere. The principal
// zones 0kl / h0l / hk0 can be shown as planes through the origin.
//
//  • Point size, colours and intensity reference are the ones of the 2D
//    plots (laueRadius, getPlotColors, Iref of the layer data), so the
//    current layer looks like the 2D layer picture seen in perspective.
//  • Red rings = gaps, exactly as in the layer viewer, but for the whole
//    sphere (blind region, beamstop shadow and missing scans show up as
//    shapes). Computed lazily, only while "mark gaps" is on.
//  • Plain WebGL 1, no library: index.html is opened from file://, where
//    ES-module builds of 3D libraries cannot be loaded.
//
// Called from layer-viewer.js: initLaue3d(host, ctx) when the layer viewer
// is (re)built, laue3dSync(state) after every layer redraw.
// Needs: laueRadius, laueFillSymmetryEnabled, symmetryFilledReflections
// (laue-viz.js), canonicalHKL (laue-merge.js), isForbiddenByHeldConditions
// (absences.js), sinThetaOverLambda (d-spacing.js), reflectionInfoHtml /
// showReflectionTip (reflection-marks.js), getPlotColors /
// onPlotThemeChange (plot-theme.js), LAYER_AXES (layer-viewer.js).

const L3D_AXIS_COLORS = { h: '#c0504d', k: '#4f9a3a', l: '#3b6fc4' }; // a*/0kl, b*/h0l, c*/hk0
const L3D_PLANE_ZONES = [ // principal zones: fixed index = 0
  { axis: 'h', name: '0kl' },
  { axis: 'k', name: 'h0l' },
  { axis: 'l', name: 'hk0' },
];

const _l3d = {
  canvas: null, overlay: null, gl: null, prog: null, host: null, els: null, ro: null,
  ctx: null, state: null,
  basis: null, R: 1,               // reciprocal basis (Cartesian, 1/Å), sphere radius
  pointBuf: null, nPoints: 0, builtFor: null, builtMode: null, builtFill: null, builtMarks: null, cpu: null, counts: null,
  gapsFor: null, gaps: null,
  lineBuf: null,
  rot: null, zoom: 1, others: 0.2, marksOthers: 0.9, size: 0.5, // opacity of the other layers' points / marks
  show: { sphere: true, h: false, k: false, l: false },
  frame: 0, failed: false,
};

// ---------- geometry ----------

// Cartesian reciprocal basis: a* along x, b* in the xy plane. Also the
// direct axes in the same frame (a = (b*×c*)/V* etc.), which are the
// normals of the layers h = n, k = n, l = n.
function l3dBasis(cell) {
  const rad = Math.PI / 180;
  const [ca, cb, cg] = [cell.alpha, cell.beta, cell.gamma].map(d => Math.cos(d * rad));
  const [sa, sb, sg] = [cell.alpha, cell.beta, cell.gamma].map(d => Math.sin(d * rad));
  const V = cell.a * cell.b * cell.c * Math.sqrt(Math.max(1e-12, 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg));
  const as = cell.b * cell.c * sa / V, bs = cell.a * cell.c * sb / V, cs = cell.a * cell.b * sg / V;
  const cas = (cb * cg - ca) / (sb * sg), cbs = (ca * cg - cb) / (sa * sg), cgs = (ca * cb - cg) / (sa * sb);
  const sgs = Math.sqrt(Math.max(1e-12, 1 - cgs * cgs));
  const cy = (cas - cbs * cgs) / sgs;
  const A = [as, 0, 0];
  const B = [bs * cgs, bs * sgs, 0];
  const C = [cs * cbs, cs * cy, cs * Math.sqrt(Math.max(0, 1 - cbs * cbs - cy * cy))];
  const Vs = l3dDot(A, l3dCross(B, C));
  const direct = { h: l3dScale(l3dCross(B, C), 1 / Vs), k: l3dScale(l3dCross(C, A), 1 / Vs), l: l3dScale(l3dCross(A, B), 1 / Vs) };
  return { A, B, C, star: { h: A, k: B, l: C }, direct };
}

const l3dDot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const l3dCross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const l3dScale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const l3dLen = a => Math.sqrt(l3dDot(a, a));
const l3dNorm = a => l3dScale(a, 1 / (l3dLen(a) || 1));
const l3dSub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

function l3dPos(h, k, l) {
  const { A, B, C } = _l3d.basis;
  return [h * A[0] + k * B[0] + l * C[0], h * A[1] + k * B[1] + l * C[1], h * A[2] + k * B[2] + l * C[2]];
}

// Rotation = 3 rows (screen x, screen y, towards the viewer) in world coords.
function l3dMul(A, B) {
  return A.map(row => [0, 1, 2].map(j => row[0] * B[0][j] + row[1] * B[1][j] + row[2] * B[2][j]));
}
function l3dOrtho(R) {
  const x = l3dNorm(R[0]);
  const y = l3dNorm(l3dSub(R[1], l3dScale(x, l3dDot(R[1], x))));
  return [x, y, l3dCross(x, y)];
}

// View straight onto a layer family, oriented like the 2D layer picture:
// its u axis to the right, v upwards.
function l3dViewOntoLayer(axis) {
  const cfg = LAYER_AXES[axis];
  const keys = ['h', 'k', 'l'];
  const U = _l3d.basis.star[keys[cfg.uIdx]], V = _l3d.basis.star[keys[cfg.vIdx]];
  const x = l3dNorm(U), z = l3dNorm(l3dCross(U, V));
  return [x, l3dCross(z, x), z];
}

// Oblique default: c* up, looked at from slightly above and to the side.
function l3dViewOblique() {
  const up = l3dNorm(_l3d.basis.C);
  const right = l3dNorm(l3dSub(_l3d.basis.B, l3dScale(up, l3dDot(_l3d.basis.B, up))));
  const base = [right, up, l3dCross(right, up)];
  const ry = -0.65, rx = 0.38;
  const Ry = [[Math.cos(ry), 0, Math.sin(ry)], [0, 1, 0], [-Math.sin(ry), 0, Math.cos(ry)]];
  const Rx = [[1, 0, 0], [0, Math.cos(rx), -Math.sin(rx)], [0, Math.sin(rx), Math.cos(rx)]];
  return l3dMul(Rx, l3dMul(Ry, base));
}

// ---------- data ----------

// One vertex per unique hkl (mean over exact repeats, like the 2D plots),
// plus the gaps and the outlier / forbidden marks (reflection-marks.js).
// Layout per vertex: pos(3) hkl(3) radius(1) kind(1);
// kind 0 measured, 1 filled in by symmetry, 2 gap, 3 outlier ring, 4 forbidden cross.
function l3dBuildPoints() {
  const { ctx, state } = _l3d;
  const data = state.data;
  const modeEl = document.getElementById('laueScaleSelect');
  const mode = modeEl ? modeEl.value : 'area';
  const fill = typeof laueFillSymmetryEnabled === 'function' && laueFillSymmetryEnabled();
  const marks = state.marks || null;
  if (_l3d.builtFor === data && _l3d.builtMode === mode && _l3d.builtFill === fill && _l3d.gapsFor === (state.showGaps ? data : null) && _l3d.builtMarks === marks) return;

  let refl = ctx.reflections;
  if (fill) refl = refl.concat(symmetryFilledReflections(ctx.reflections, ctx.group));
  const map = new Map();
  for (const r of refl) {
    const id = r.h + ',' + r.k + ',' + r.l;
    let e = map.get(id);
    if (!e) { e = { h: r.h, k: r.k, l: r.l, sum: 0, n: 0, generated: !!r.generated }; map.set(id, e); }
    e.sum += r.f2; e.n++;
  }
  const b = _l3d.basis;
  const rmax = 0.5 * Math.min(l3dLen(b.A), l3dLen(b.B), l3dLen(b.C));
  const gaps = state.showGaps ? l3dGaps(data) : [];
  const nMarks = marks ? marks.outliers.size + marks.forbidden.size : 0;
  const n = map.size + gaps.length + nMarks;
  const arr = new Float32Array(n * 8);
  let i = 0;
  const put = (h, k, l, rad, kind) => {
    const p = l3dPos(h, k, l);
    arr[i++] = p[0]; arr[i++] = p[1]; arr[i++] = p[2];
    arr[i++] = h; arr[i++] = k; arr[i++] = l;
    arr[i++] = rad; arr[i++] = kind;
  };
  for (const g of gaps) put(g[0], g[1], g[2], rmax * 0.6, 2);
  // laueRadius floors at 0.8 (meant as pixels), so ask it for a fraction of
  // a large rmax and scale that to world units; I <= 0 gives 0 (drawn at the
  // minimum pixel size by the shader).
  const frac = I => (I > 0 ? laueRadius(I, data.Iref, 1000, mode) / 1000 : 0);
  for (const e of map.values()) put(e.h, e.k, e.l, rmax * frac(e.sum / e.n), e.generated ? 1 : 0);
  // marks last, so they are drawn on top of their points
  if (marks) {
    const radAt = key => { const e = map.get(key); return e ? rmax * frac(e.sum / e.n) : 0; };
    for (const key of marks.outliers.keys()) { const [h, k, l] = key.split(',').map(Number); put(h, k, l, radAt(key), 3); }
    for (const key of marks.forbidden.keys()) { const [h, k, l] = key.split(',').map(Number); put(h, k, l, radAt(key), 4); }
  }
  _l3d.cpu = arr; // for hover / click picking

  const gl = _l3d.gl;
  if (!_l3d.pointBuf) _l3d.pointBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, _l3d.pointBuf);
  gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW);
  _l3d.nPoints = n;
  _l3d.builtFor = data; _l3d.builtMode = mode; _l3d.builtFill = fill; _l3d.builtMarks = marks;
  _l3d.gapsFor = state.showGaps ? data : null;
  l3dStatus(map.size, gaps.length);
}

// Totals for the whole sphere; marks follow their checkboxes.
function l3dStatus(nPoints, nGaps) {
  const st = _l3d.state, mk = st.mark || {}, marks = st.marks;
  if (nPoints !== undefined) _l3d.counts = { nPoints, nGaps };
  const c = _l3d.counts || { nPoints: 0, nGaps: 0 };
  let html = `${c.nPoints} lattice points`;
  if (st.showGaps) html += ` · ${c.nGaps} gaps`;
  if (marks && mk.outliers) html += ` · <span class="reflMarkOutlier">${marks.outliers.size} outliers</span>`;
  if (marks && mk.forbidden) html += ` · <span class="reflMarkForbidden">${marks.forbidden.size} forbidden but observed</span>`;
  _l3d.els.status.innerHTML = `<small>${html} (whole sphere)</small>`;
}

// Same definition as computeLayerInfo (layer-viewer.js), for all layers at
// once: allowed by the observed lattice, inside the resolution sphere, not
// systematically forbidden, and no member of its symmetry class measured.
function l3dGaps(data) {
  if (_l3d.gaps && _l3d.gaps.data === data) return _l3d.gaps.list;
  const { ctx } = _l3d;
  const ax = data.axes.l; // u = h, v = k, n = l
  const lim = data.sMax + 1e-9;
  const list = [];
  for (let h = ax.uMin - 1; h <= ax.uMax + 1; h++) {
    for (let k = ax.vMin - 1; k <= ax.vMax + 1; k++) {
      for (let l = ax.nMin - 1; l <= ax.nMax + 1; l++) {
        if (h === 0 && k === 0 && l === 0) continue;
        if (sinThetaOverLambda(h, k, l, data.recip) > lim) continue;
        if (ctx.latticeTest && !ctx.latticeTest({ h, k, l })) continue;
        if (ctx.scan && isForbiddenByHeldConditions(ctx.scan, h, k, l)) continue;
        const c = canonicalHKL(h, k, l, ctx.group);
        if (!data.classes.has(c[0] + ',' + c[1] + ',' + c[2])) list.push([h, k, l]);
      }
    }
  }
  _l3d.gaps = { data, list };
  return list;
}

// ---------- WebGL ----------

const L3D_VS_POINTS = `
attribute vec3 aPos; attribute vec3 aHkl; attribute float aRad; attribute float aKind;
uniform mat4 uMvp; uniform float uPx; uniform float uSize; uniform float uMinPx;
uniform vec3 uAxis; uniform float uLayer; uniform float uPass; uniform vec3 uMarks; // gaps, outliers, forbidden
varying float vKind;
void main() {
  float on = abs(dot(aHkl, uAxis) - uLayer) < 0.5 ? 1.0 : 0.0;
  float markOff = aKind > 3.5 ? 1.0 - uMarks.z : aKind > 2.5 ? 1.0 - uMarks.y : aKind > 1.5 ? 1.0 - uMarks.x : 0.0;
  bool hide = (uPass < 0.5 && on < 0.5) || (uPass > 0.5 && on > 0.5) || markOff > 0.5;
  gl_Position = hide ? vec4(2.0, 2.0, 2.0, 1.0) : uMvp * vec4(aPos, 1.0);
  float d = 2.0 * aRad * uPx * uSize;
  float u = uMinPx / 1.3; // device pixels per CSS pixel
  if (aKind > 3.5) gl_PointSize = max(9.0 * u, d * 1.15 + 3.0 * u);
  else if (aKind > 2.5) gl_PointSize = max(9.0 * u, d + 6.0 * u);
  else if (aKind > 1.5) gl_PointSize = max(uMinPx * 3.0, d);
  else gl_PointSize = max(uMinPx, d);
  vKind = aKind;
}`;
const L3D_FS_POINTS = `
precision mediump float;
uniform vec4 uColM; uniform vec4 uColG; uniform vec4 uColGap; uniform vec4 uColOut; uniform vec4 uColForb;
uniform float uAlpha; uniform float uGapAlpha;
varying float vKind;
void main() {
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(p, p);
  if (r2 > 1.0) discard;
  vec4 c;
  float a = uAlpha;
  if (vKind > 3.5) { // cross
    if (abs(p.x - p.y) > 0.42 && abs(p.x + p.y) > 0.42) discard;
    c = uColForb; a = uGapAlpha;
  }
  else if (vKind > 2.5) { if (r2 < 0.6) discard; c = uColOut; a = uGapAlpha; } // ring around the point
  else if (vKind > 1.5) { if (r2 < 0.5) discard; c = uColGap; a = uGapAlpha; }
  else { c = vKind > 0.5 ? uColG : uColM; c.rgb *= 0.82 + 0.18 * (1.0 - r2); }
  if (a <= 0.0) discard;
  gl_FragColor = vec4(c.rgb, a);
}`;
const L3D_VS_FLAT = `
attribute vec3 aPos; uniform mat4 uMvp;
void main() { gl_Position = uMvp * vec4(aPos, 1.0); }`;
const L3D_FS_FLAT = `
precision mediump float; uniform vec4 uCol;
void main() { gl_FragColor = uCol; }`;

function l3dProgram(gl, vs, fs) {
  const sh = (type, src) => {
    const s = gl.createShader(type);
    gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    return s;
  };
  const p = gl.createProgram();
  gl.attachShader(p, sh(gl.VERTEX_SHADER, vs)); gl.attachShader(p, sh(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  const loc = {};
  const nA = gl.getProgramParameter(p, gl.ACTIVE_ATTRIBUTES), nU = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < nA; i++) { const n = gl.getActiveAttrib(p, i).name; loc[n] = gl.getAttribLocation(p, n); }
  for (let i = 0; i < nU; i++) { const n = gl.getActiveUniform(p, i).name; loc[n] = gl.getUniformLocation(p, n); }
  return { p, loc };
}

function l3dInitGL() {
  if (_l3d.gl || _l3d.failed) return !!_l3d.gl;
  const opts = { antialias: true, alpha: false, premultipliedAlpha: false, preserveDrawingBuffer: false };
  const gl = _l3d.canvas.getContext('webgl', opts) || _l3d.canvas.getContext('experimental-webgl', opts);
  if (!gl) { _l3d.failed = true; return false; }
  try {
    _l3d.prog = { pts: l3dProgram(gl, L3D_VS_POINTS, L3D_FS_POINTS), flat: l3dProgram(gl, L3D_VS_FLAT, L3D_FS_FLAT) };
  } catch (err) { _l3d.failed = true; console.error('laue-3d:', err); return false; }
  _l3d.gl = gl;
  _l3d.lineBuf = gl.createBuffer();
  _l3d.canvas.addEventListener('webglcontextlost', e => { e.preventDefault(); _l3d.gl = null; _l3d.pointBuf = null; _l3d.builtFor = null; });
  _l3d.canvas.addEventListener('webglcontextrestored', () => { l3dInitGL(); l3dRequest(); });
  return true;
}

// Parses '#rgb', '#rrggbb' and 'rgba(r,g,b,a)' into [r,g,b,a] 0..1.
function l3dColor(str, alpha) {
  let c = [0, 0, 0, 1];
  const m = /rgba?\(([^)]+)\)/.exec(str);
  if (m) { const v = m[1].split(',').map(Number); c = [v[0] / 255, v[1] / 255, v[2] / 255, v.length > 3 ? v[3] : 1]; }
  else if (str[0] === '#') {
    const hex = str.length === 4 ? str.slice(1).split('').map(x => x + x).join('') : str.slice(1);
    c = [0, 2, 4].map(i => parseInt(hex.slice(i, i + 2), 16) / 255).concat(1);
  }
  if (alpha !== undefined) c[3] = alpha;
  return c;
}

function l3dMvp() {
  const R = _l3d.rot, s = _l3d.zoom / (1.25 * _l3d.R), zr = 1.4 * _l3d.R;
  return {
    s,
    m: new Float32Array([
      s * R[0][0], s * R[1][0], -R[2][0] / zr, 0,
      s * R[0][1], s * R[1][1], -R[2][1] / zr, 0,
      s * R[0][2], s * R[1][2], -R[2][2] / zr, 0,
      0, 0, 0, 1]),
  };
}

// Circle of the plane {r : r·n̂ = d} with the sphere, as a polygon.
function l3dDiscPoints(normal, d, segs = 72) {
  const R = _l3d.R;
  if (Math.abs(d) >= R) return null;
  const rho = Math.sqrt(R * R - d * d);
  const n = l3dNorm(normal);
  const t = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  const e1 = l3dNorm(l3dCross(n, t)), e2 = l3dCross(n, e1);
  const pts = [];
  for (let i = 0; i <= segs; i++) {
    const a = 2 * Math.PI * i / segs, ca = Math.cos(a) * rho, sa = Math.sin(a) * rho;
    pts.push([n[0] * d + e1[0] * ca + e2[0] * sa, n[1] * d + e1[1] * ca + e2[1] * sa, n[2] * d + e1[2] * ca + e2[2] * sa]);
  }
  return { pts, center: l3dScale(n, d) };
}

function l3dDrawFlat(mode, verts, col, mvp) {
  const gl = _l3d.gl, { p, loc } = _l3d.prog.flat;
  gl.useProgram(p);
  gl.uniformMatrix4fv(loc.uMvp, false, mvp);
  gl.uniform4fv(loc.uCol, col);
  gl.bindBuffer(gl.ARRAY_BUFFER, _l3d.lineBuf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(verts.flat()), gl.DYNAMIC_DRAW);
  gl.enableVertexAttribArray(loc.aPos);
  gl.vertexAttribPointer(loc.aPos, 3, gl.FLOAT, false, 0, 0);
  gl.drawArrays(mode, 0, verts.length);
  gl.disableVertexAttribArray(loc.aPos);
}

function l3dDrawPoints(pass, alpha, gapAlpha, size, mvp, colors) {
  const gl = _l3d.gl, { p, loc } = _l3d.prog.pts, st = _l3d.state;
  gl.useProgram(p);
  gl.uniformMatrix4fv(loc.uMvp, false, mvp.m);
  gl.uniform1f(loc.uPx, mvp.s * _l3d.canvas.width / 2);
  gl.uniform1f(loc.uSize, size);
  gl.uniform1f(loc.uMinPx, 1.3 * (window.devicePixelRatio || 1));
  gl.uniform3fv(loc.uAxis, [st.axis === 'h' ? 1 : 0, st.axis === 'k' ? 1 : 0, st.axis === 'l' ? 1 : 0]);
  gl.uniform1f(loc.uLayer, st.n);
  gl.uniform1f(loc.uPass, pass);
  const mk = st.mark || {}, on = st.marks ? 1 : 0;
  gl.uniform3fv(loc.uMarks, [st.showGaps ? 1 : 0, on && mk.outliers ? 1 : 0, on && mk.forbidden ? 1 : 0]);
  gl.uniform4fv(loc.uColM, l3dColor(colors.point));
  gl.uniform4fv(loc.uColG, l3dColor(colors.symFilled));
  gl.uniform4fv(loc.uColGap, l3dColor(colors.gap));
  gl.uniform4fv(loc.uColOut, l3dColor(colors.outlier));
  gl.uniform4fv(loc.uColForb, l3dColor(colors.forbidden));
  gl.uniform1f(loc.uAlpha, alpha);
  gl.uniform1f(loc.uGapAlpha, gapAlpha);
  gl.bindBuffer(gl.ARRAY_BUFFER, _l3d.pointBuf);
  const F = 4, stride = 8 * F;
  const attr = (name, size, off) => { gl.enableVertexAttribArray(loc[name]); gl.vertexAttribPointer(loc[name], size, gl.FLOAT, false, stride, off * F); };
  attr('aPos', 3, 0); attr('aHkl', 3, 3); attr('aRad', 1, 6); attr('aKind', 1, 7);
  gl.drawArrays(gl.POINTS, 0, _l3d.nPoints);
  for (const n of ['aPos', 'aHkl', 'aRad', 'aKind']) gl.disableVertexAttribArray(loc[n]);
}

function l3dRender() {
  _l3d.frame = 0;
  const cv = _l3d.canvas;
  if (!_l3d.state || !cv || !cv.isConnected) return;
  const cssW = cv.parentElement.clientWidth;
  if (!cssW) return; // tab not visible yet; the ResizeObserver calls again
  if (!l3dInitGL()) { l3dShowFallback(); return; }
  const dpr = window.devicePixelRatio || 1;
  const px = Math.round(cssW * dpr);
  if (cv.width !== px || cv.height !== px) { cv.width = px; cv.height = px; }
  cv.style.width = cssW + 'px'; cv.style.height = cssW + 'px';
  l3dBuildPoints();

  const gl = _l3d.gl, colors = getPlotColors(), mvp = l3dMvp(), st = _l3d.state;
  const bg = l3dColor(colors.bg);
  gl.viewport(0, 0, px, px);
  gl.clearColor(bg[0], bg[1], bg[2], 1);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL);
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

  // 1) points of the current layer, opaque and at the size of the 2D picture
  gl.depthMask(true);
  l3dDrawPoints(0, 1, 1, 1, mvp, colors);

  // 2) sphere wireframe and axes (lines)
  const R = _l3d.R;
  if (_l3d.show.sphere) {
    const gridCol = l3dColor(colors.axis, 0.28);
    const up = l3dNorm(_l3d.basis.C);
    for (let i = -2; i <= 2; i++) { // parallels around c*
      const disc = l3dDiscPoints(up, R * i / 3);
      if (disc) l3dDrawFlat(gl.LINE_STRIP, disc.pts, gridCol, mvp.m);
    }
    const ref = l3dNorm(_l3d.basis.A), side = l3dCross(up, ref);
    for (let i = 0; i < 6; i++) { // meridians
      const a = Math.PI * i / 6;
      l3dDrawFlat(gl.LINE_STRIP, l3dDiscPoints(l3dAdd(l3dScale(ref, Math.cos(a)), l3dScale(side, Math.sin(a))), 0).pts, gridCol, mvp.m);
    }
  }
  for (const key of ['h', 'k', 'l']) {
    const dir = l3dNorm(_l3d.basis.star[key]);
    l3dDrawFlat(gl.LINES, [[0, 0, 0], l3dScale(dir, R * 1.12)], l3dColor(L3D_AXIS_COLORS[key], 0.95), mvp.m);
  }

  // 3) transparent: planes, then all other points
  gl.depthMask(false);
  const planes = [];
  for (const z of L3D_PLANE_ZONES) if (_l3d.show[z.axis] && !(z.axis === st.axis && st.n === 0)) planes.push({ axis: z.axis, n: 0, fill: 0.10, rim: 0.55 });
  planes.push({ axis: st.axis, n: st.n, fill: 0.20, rim: 0.95 });
  for (const pl of planes) {
    const dvec = _l3d.basis.direct[pl.axis];
    const disc = l3dDiscPoints(dvec, pl.n / l3dLen(dvec));
    if (!disc) continue;
    const col = L3D_AXIS_COLORS[pl.axis];
    l3dDrawFlat(gl.TRIANGLE_FAN, [disc.center].concat(disc.pts), l3dColor(col, pl.fill), mvp.m);
    l3dDrawFlat(gl.LINE_STRIP, disc.pts, l3dColor(col, pl.rim), mvp.m);
  }
  // marks of the other layers have their own opacity: with the points hidden
  // and the marks full, a blind region or a cloud of outliers shows up as a
  // shape of its own
  const anyMarks = (st.showGaps || (st.marks && st.mark && (st.mark.outliers || st.mark.forbidden))) && _l3d.marksOthers > 0;
  if (_l3d.others > 0 || anyMarks) l3dDrawPoints(1, _l3d.others, anyMarks ? _l3d.marksOthers : 0, _l3d.size, mvp, colors);
  gl.depthMask(true);

  l3dDrawOverlay(cssW, colors, mvp.s);
}

const l3dAdd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];

// Axis labels and caption on a 2D canvas on top of the WebGL one.
function l3dDrawOverlay(cssW, colors, s) {
  const ov = _l3d.overlay, dpr = window.devicePixelRatio || 1;
  const px = Math.round(cssW * dpr);
  if (ov.width !== px || ov.height !== px) { ov.width = px; ov.height = px; }
  ov.style.width = cssW + 'px'; ov.style.height = cssW + 'px';
  const g = ov.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, cssW, cssW);
  const Rm = _l3d.rot;
  const proj = p => [(1 + s * l3dDot(Rm[0], p)) * cssW / 2, (1 - s * l3dDot(Rm[1], p)) * cssW / 2];
  g.font = 'bold 12px sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle';
  const names = { h: 'a*', k: 'b*', l: 'c*' };
  for (const key of ['h', 'k', 'l']) {
    const [x, y] = proj(l3dScale(l3dNorm(_l3d.basis.star[key]), _l3d.R * 1.2));
    g.fillStyle = L3D_AXIS_COLORS[key];
    g.fillText(names[key], x, y);
  }
  const st = _l3d.state, cfg = LAYER_AXES[st.axis];
  g.font = '13px sans-serif'; g.textAlign = 'right'; g.textBaseline = 'alphabetic';
  g.fillStyle = colors.text;
  g.fillText(`${cfg.fixed} = ${st.n}`, cssW - 8, 20);
  g.font = '10px sans-serif'; g.textAlign = 'left';
  g.fillText(`sphere: d*max = ${_l3d.R.toFixed(3)} Å⁻¹ (d = ${(1 / _l3d.R).toFixed(2)} Å)`, 8, cssW - 8);
}

function l3dShowFallback() {
  if (!_l3d.host) return;
  _l3d.host.innerHTML = '<p><small>3D view not available: this browser offers no WebGL. The 2D layer viewer is unaffected.</small></p>';
}

function l3dRequest() {
  if (!_l3d.frame) _l3d.frame = requestAnimationFrame(l3dRender);
}

// ---------- interaction ----------

function l3dRotate(ax, ay) {
  const cy = Math.cos(ax), sy = Math.sin(ax), cx = Math.cos(ay), sx = Math.sin(ay);
  const Ry = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]];
  const Rx = [[1, 0, 0], [0, cx, -sx], [0, sx, cx]];
  _l3d.rot = l3dOrtho(l3dMul(Rx, l3dMul(Ry, _l3d.rot)));
}

// Nearest visible point / gap under the cursor (CSS pixels of the overlay).
// Points of the current layer win over the others, then the one in front.
function l3dPick(clientX, clientY) {
  const arr = _l3d.cpu, st = _l3d.state;
  if (!arr || !st || !_l3d.rot) return null;
  const b = _l3d.overlay.getBoundingClientRect();
  if (!b.width) return null;
  const cssW = b.width, mx = clientX - b.left, my = clientY - b.top;
  const { s } = l3dMvp(), R = _l3d.rot, half = cssW / 2, pxW = s * half;
  const mk = st.mark || {};
  const showKind = [
    _l3d.others > 0, _l3d.others > 0, !!st.showGaps,
    !!(st.marks && mk.outliers), !!(st.marks && mk.forbidden),
  ];
  const axIdx = { h: 3, k: 4, l: 5 }[st.axis];
  let best = null;
  for (let i = 0; i < arr.length; i += 8) {
    const kind = arr[i + 7], on = Math.abs(arr[i + axIdx] - st.n) < 0.5;
    const visible = kind < 1.5 ? (on || showKind[kind]) : showKind[kind] && (on || _l3d.marksOthers > 0);
    if (!visible) continue;
    const x = arr[i], y = arr[i + 1], z = arr[i + 2];
    const sx = half + pxW * (R[0][0] * x + R[0][1] * y + R[0][2] * z);
    const sy = half - pxW * (R[1][0] * x + R[1][1] * y + R[1][2] * z);
    const d = Math.hypot(sx - mx, sy - my);
    const rad = arr[i + 6] * pxW * (on ? 1 : _l3d.size);
    if (d > Math.max(rad + 2, typeof touchSlop === 'function' ? touchSlop(6) : 6)) continue;
    const depth = R[2][0] * x + R[2][1] * y + R[2][2] * z;
    const score = (on ? 1e6 : 0) - d + 0.5 * depth / _l3d.R; // nearest on screen; depth only breaks ties
    if (!best || score > best.score) best = { score, h: arr[i + 3], k: arr[i + 4], l: arr[i + 5], kind };
  }
  if (!best) return null;
  return { h: best.h, k: best.k, l: best.l, kind: best.kind === 2 ? 'gap' : best.kind === 1 ? 'generated' : 'measured' };
}

function l3dTipAt(e, pin) {
  const it = l3dPick(e.clientX, e.clientY);
  if (it) showReflectionTip(reflectionInfoHtml(_l3d.ctx, it, _l3d.state.marks || null), e.clientX, e.clientY, pin);
  else if (!pin) hideReflectionTip(false);
  _l3d.overlay.style.cursor = it ? 'pointer' : '';
}

function l3dAttachPointer(el) {
  const ptrs = new Map();
  let pinch = null, start = null, hoverFrame = 0, lastMove = null, lastTap = null;
  el.addEventListener('pointerdown', e => {
    if (ptrs.size === 0) start = { x: e.clientX, y: e.clientY, moved: false };
    hideReflectionTip(false);
    el.setPointerCapture(e.pointerId);
    ptrs.set(e.pointerId, [e.clientX, e.clientY]);
    if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; pinch = { d: Math.hypot(a[0] - b[0], a[1] - b[1]), zoom: _l3d.zoom }; }
  });
  el.addEventListener('pointermove', e => {
    const prev = ptrs.get(e.pointerId);
    if (!prev) { // hovering (mouse, no button): info of the point under the cursor
      if (e.pointerType !== 'mouse') return;
      lastMove = e;
      if (!hoverFrame) hoverFrame = requestAnimationFrame(() => { hoverFrame = 0; if (lastMove) l3dTipAt(lastMove, false); });
      return;
    }
    if (start && Math.hypot(e.clientX - start.x, e.clientY - start.y) > 4) start.moved = true;
    ptrs.set(e.pointerId, [e.clientX, e.clientY]);
    if (ptrs.size === 2 && pinch) {
      const [a, b] = [...ptrs.values()];
      _l3d.zoom = Math.max(0.4, Math.min(8, pinch.zoom * Math.hypot(a[0] - b[0], a[1] - b[1]) / (pinch.d || 1)));
    } else if (ptrs.size === 1) {
      l3dRotate((e.clientX - prev[0]) * 0.01, (e.clientY - prev[1]) * 0.01);
    }
    l3dRequest();
  });
  const up = e => { ptrs.delete(e.pointerId); if (ptrs.size < 2) pinch = null; };
  el.addEventListener('pointerup', e => {
    const click = ptrs.size === 1 && start && !start.moved;
    up(e);
    if (!click) return;
    // double tap with a finger = reset (dblclick does not reach a page on iPad)
    if (e.pointerType !== 'mouse') {
      const now = Date.now();
      if (lastTap && now - lastTap.t < 350 && Math.hypot(e.clientX - lastTap.x, e.clientY - lastTap.y) < 30) {
        lastTap = null; _l3d.rot = l3dViewOblique(); _l3d.zoom = 1; l3dRequest(); return;
      }
      lastTap = { t: now, x: e.clientX, y: e.clientY };
    }
    l3dTipAt(e, true); // a click (or tap) without dragging pins the info
  });
  el.addEventListener('pointercancel', up);
  el.addEventListener('pointerleave', () => { lastMove = null; hideReflectionTip(false); });
  el.addEventListener('wheel', e => {
    e.preventDefault();
    _l3d.zoom = Math.max(0.4, Math.min(8, _l3d.zoom * Math.exp(-e.deltaY * 0.0015)));
    l3dRequest();
  }, { passive: false });
  el.addEventListener('dblclick', () => { _l3d.rot = l3dViewOblique(); _l3d.zoom = 1; l3dRequest(); });
}

// ---------- public ----------

function initLaue3d(host, ctx) {
  _l3d.host = host;
  closeViewModal(); // a new data set while enlarged: back to the page first
  _l3d.ctx = ctx;
  _l3d.state = null;
  _l3d.builtFor = null; _l3d.builtMarks = null; _l3d.gaps = null; _l3d.gapsFor = null; _l3d.cpu = null; _l3d.counts = null;
  _l3d.basis = l3dBasis(ctx.cell);
  if (!_l3d.canvas) { // one canvas (= one WebGL context) for the lifetime of the page
    _l3d.canvas = document.createElement('canvas');
    _l3d.canvas.className = 'laue3dGL';
    _l3d.overlay = document.createElement('canvas');
    _l3d.overlay.className = 'laue3dOverlay';
    l3dAttachPointer(_l3d.overlay);
  }
  if (_l3d.failed) { l3dShowFallback(); return; }
  const chk = (id, key, label, title) =>
    `<label${title ? ` title="${title}"` : ''}><input type="checkbox" id="${id}"${_l3d.show[key] ? ' checked' : ''}> ${label}</label>`;
  host.innerHTML = `
    <div class="laue3dMain">
      <div class="viewHead"><span class="ctlTitle">3D · whole sphere</span>
        <button type="button" class="viewExpand" title="enlarge the 3D view (with its controls)" aria-label="Enlarge 3D view">⤢</button></div>
      <div class="laue3dWrap"></div>
      <div class="laue3dStatus"></div>
    </div>
    <div class="laue3dSide">
      <div class="laue3dGroup">
        <span class="laue3dGroupTitle">View</span>
        <button type="button" id="l3dViewLayer" title="look straight onto the current layer, oriented like the 2D picture">⊥ layer</button>
        <button type="button" id="l3dViewOblique">oblique</button>
      </div>
      <div class="laue3dGroup">
        <span class="laue3dGroupTitle">Other layers</span>
        <label class="l3dRow" title="reflections of the other layers"><span>points</span>
          <select id="l3dOthers">
            <option value="0.2">faint</option>
            <option value="0.75">full</option>
            <option value="0">hidden</option>
          </select></label>
        <label class="l3dRow" title="gaps, outliers and forbidden-but-observed marks of the other layers (switched on under Mark); the current layer always shows its marks"><span>marks</span>
          <select id="l3dMarksOthers">
            <option value="0.3">faint</option>
            <option value="0.9">full</option>
            <option value="0">hidden</option>
          </select></label>
        <label class="l3dRow" title="spot size of the other layers; the current layer is drawn at the size of the 2D picture"><span>size</span>
          <input type="range" id="l3dSize" min="0.15" max="1" step="0.05" value="${_l3d.size}"></label>
      </div>
      <div class="laue3dGroup">
        <span class="laue3dGroupTitle">Show</span>
        ${chk('l3dSphere', 'sphere', 'sphere')}
        ${L3D_PLANE_ZONES.map(z => chk('l3dPlane_' + z.axis, z.axis, `<span style="color:${L3D_AXIS_COLORS[z.axis]}">■</span> ${z.name}`, `zone plane ${z.name} through the origin`)).join('')}
      </div>
      <p class="laue3dHint"><small>Drag to rotate · wheel or pinch to zoom · double-click (double-tap) to reset.</small></p>
    </div>`;
  const wrap = host.querySelector('.laue3dWrap');
  wrap.appendChild(_l3d.canvas);
  wrap.appendChild(_l3d.overlay);
  const q = id => host.querySelector('#' + id);
  _l3d.els = { status: host.querySelector('.laue3dStatus') };
  q('l3dOthers').value = String(_l3d.others);
  q('l3dOthers').addEventListener('change', e => { _l3d.others = parseFloat(e.target.value); l3dRequest(); });
  q('l3dMarksOthers').value = String(_l3d.marksOthers);
  q('l3dMarksOthers').addEventListener('change', e => { _l3d.marksOthers = parseFloat(e.target.value); l3dRequest(); });
  q('l3dSize').addEventListener('input', e => { _l3d.size = parseFloat(e.target.value); l3dRequest(); });
  q('l3dSphere').addEventListener('change', e => { _l3d.show.sphere = e.target.checked; l3dRequest(); });
  for (const z of L3D_PLANE_ZONES) q('l3dPlane_' + z.axis).addEventListener('change', e => { _l3d.show[z.axis] = e.target.checked; l3dRequest(); });
  q('l3dViewLayer').addEventListener('click', () => { if (_l3d.state) { _l3d.rot = l3dViewOntoLayer(_l3d.state.axis); l3dRequest(); } });
  q('l3dViewOblique').addEventListener('click', () => { _l3d.rot = l3dViewOblique(); l3dRequest(); });
  host.querySelector('.viewExpand').addEventListener('click', () => (viewModalIsOpen() ? closeViewModal() : l3dOpenModal()));
  _l3d.rot = l3dViewOblique();
  _l3d.zoom = 1;
  if (_l3d.ro) _l3d.ro.disconnect();
  if (window.ResizeObserver) { _l3d.ro = new ResizeObserver(l3dRequest); _l3d.ro.observe(wrap); }
}

// ---------- enlarged view (view-modal.js) ----------
// The 3D picture and its control column, plus the layer viewer's controls.

function l3dOpenModal() {
  if (!_l3d.host) return;
  const main = _l3d.host.querySelector('.laue3dMain'), side = _l3d.host.querySelector('.laue3dSide');
  const wrap = main.querySelector('.laue3dWrap'), btn = main.querySelector('.viewExpand');
  openViewModal({
    title: 'Layer viewer · 3D',
    top: [document.querySelector('#layerViewerContainer .layerControls')],
    body: [main, side],
    onResize: () => { wrap.style.width = viewModalSquareSize(170) + 'px'; l3dRequest(); },
    onClose: () => {
      wrap.style.width = '';
      btn.textContent = '⤢'; btn.title = 'enlarge the 3D view (with its controls)';
      l3dRequest();
    },
  });
  btn.textContent = '⤡'; btn.title = 'back to the page';
}

// state: { axis, n, data, showGaps } of the layer viewer.
function laue3dSync(state) {
  if (!_l3d.host || !_l3d.ctx) return;
  _l3d.state = state;
  _l3d.R = 2 * state.data.sMax || 1;
  if (_l3d.els && _l3d.counts) l3dStatus();
  l3dRequest();
}

window.addEventListener('resize', l3dRequest);
onPlotThemeChange(l3dRequest);
