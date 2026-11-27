// touch.js
// Tablet / touch support. No switch: two automatic levels.
//
//  • Layout: <html class="touch"> while the PRIMARY pointer is a finger
//    (matchMedia '(pointer: coarse)'), e.g. an iPad without trackpad. The CSS
//    then enlarges controls and the tab bar; info buttons (ⓘ) are added to
//    labels that carry an explanation in their title.
//  • Behaviour: decided per input. touchInput() is true when the last pointer
//    that went down was a finger or a pen. Then a tap on an element with a
//    title (or on a symbol of the space-group diagram) shows that text in the
//    info bar at the bottom of the screen instead of a hover tooltip, and the
//    reflection info of the plots goes there too (reflection-marks.js).
//    With a mouse or trackpad nothing changes, also on a touch laptop.
//
//   touchInput()            last input was touch / pen
//   touchSlop(mouseRadius)  hit radius in CSS px for picking a plot point
//   showInfoBar(html, kind) hideInfoBar(kind?)   the bottom info bar
//   attachSwipeStep(el, step)   horizontal swipe on a plot = previous/next layer

const TOUCH = {
  last: 'mouse',
  coarse: window.matchMedia ? window.matchMedia('(pointer: coarse)') : null,
  bar: null, barKind: null,
};

function touchInput() { return TOUCH.last === 'touch' || TOUCH.last === 'pen'; }
function touchSlop(mouseRadius) { return touchInput() ? 22 : mouseRadius; }

// Capture phase, so every other handler of the same tap already sees it.
document.addEventListener('pointerdown', e => { TOUCH.last = e.pointerType || 'mouse'; }, true);
// Browsers without pointer events on some old iPads: touchstart counts too.
document.addEventListener('touchstart', () => { TOUCH.last = 'touch'; }, { capture: true, passive: true });

function _touchApplyClass() {
  const on = !!(TOUCH.coarse && TOUCH.coarse.matches);
  document.documentElement.classList.toggle('touch', on);
  if (on) _touchScanSoon();
}
if (TOUCH.coarse) {
  if (TOUCH.coarse.addEventListener) TOUCH.coarse.addEventListener('change', _touchApplyClass);
  else if (TOUCH.coarse.addListener) TOUCH.coarse.addListener(_touchApplyClass);
}

// ---------- info bar ----------
function showInfoBar(html, kind = 'tip') {
  if (!TOUCH.bar) {
    const bar = document.createElement('div');
    bar.className = 'infoBar';
    bar.setAttribute('role', 'status');
    bar.setAttribute('aria-live', 'polite');
    bar.innerHTML = '<div class="infoBarText"></div><button type="button" class="infoBarClose" aria-label="Close">×</button>';
    bar.querySelector('.infoBarClose').addEventListener('click', e => { e.stopPropagation(); hideInfoBar(); });
    document.body.appendChild(bar);
    TOUCH.bar = bar;
  }
  TOUCH.bar.querySelector('.infoBarText').innerHTML = html;
  TOUCH.bar.classList.add('open');
  TOUCH.barKind = kind;
  document.body.classList.add('infoBarOpen');
}

function hideInfoBar(kind) {
  if (!TOUCH.bar || (kind && TOUCH.barKind !== kind)) return;
  TOUCH.bar.classList.remove('open');
  TOUCH.barKind = null;
  document.body.classList.remove('infoBarOpen');
  _touchClearSgdSelection();
}

const _touchEsc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

// ---------- space-group diagram: tap = select + explain ----------
function _touchClearSgdSelection() {
  document.querySelectorAll('.sgdSel').forEach(el => el.classList.remove('sgdSel'));
}

function _touchSelectSgd(g) {
  _touchClearSgdSelection();
  g.classList.add('sgdSel');
  const t = g.querySelector(':scope > title');
  const text = t ? t.textContent : '';
  // a general position names its lines of the coordinate list: "(3) …"
  if (g.classList.contains('sgdPosGroup')) {
    const box = g.closest('.sgdBox');
    const items = box ? box.querySelectorAll('.sgdList > li:not(.sgdCentring)') : [];
    for (const m of text.matchAll(/\((\d+)\)/g)) {
      const li = items[Number(m[1]) - 1];
      if (li) li.classList.add('sgdSel');
    }
  }
  showInfoBar(_touchEsc(text), 'sgd');
}

// ---------- generic: tap on an element with a title ----------
const _TOUCH_INTERACTIVE = 'button, select, input, textarea, option, a[href], summary, iframe, canvas, [role="button"], .zoomable';

function _touchTipTarget(target) {
  const el = target.closest('[title]');
  if (!el) return null;
  if (el.matches(_TOUCH_INTERACTIVE)) return null;           // a control: its tap does something else
  if (target.closest(_TOUCH_INTERACTIVE) && el.contains(target.closest(_TOUCH_INTERACTIVE))) return null;
  if (el.tagName === 'LABEL' && el.querySelector('input, select')) return null; // has its own ⓘ button
  return el;
}

document.addEventListener('click', e => {
  const t = e.target;
  if (!(t instanceof Element)) return;
  if (TOUCH.bar && TOUCH.bar.contains(t)) return;
  // ⓘ next to a label: show the label's explanation, don't toggle the switch
  const info = t.closest('.tipBtn');
  if (info) {
    e.preventDefault(); e.stopPropagation();
    const owner = info.closest('[title], [data-tip]');
    const text = owner ? (owner.getAttribute('title') || owner.dataset.tip) : '';
    if (text) showInfoBar(_touchEsc(text));
    return;
  }
  if (!touchInput()) return;
  const sgd = t.closest('.sgDiagram .sgdEl, .sgDiagram .sgdPosGroup');
  if (sgd) { _touchSelectSgd(sgd); return; }
  const el = _touchTipTarget(t);
  if (el) { showInfoBar(_touchEsc(el.getAttribute('title'))); return; }
  // a tap anywhere else closes a general tip (plot info stays until the next pick or ×)
  if (TOUCH.barKind === 'tip' || TOUCH.barKind === 'sgd') hideInfoBar();
}, true);

document.addEventListener('keydown', e => { if (e.key === 'Escape') hideInfoBar(); });

// ---------- layout level: ⓘ buttons and scrollable tables ----------
let _touchScanQueued = false;
function _touchScanSoon() {
  if (_touchScanQueued) return;
  _touchScanQueued = true;
  requestAnimationFrame(() => { _touchScanQueued = false; _touchScan(); });
}

function _touchScan() {
  if (!document.documentElement.classList.contains('touch')) return;
  // labels with a control and an explanation get an ⓘ at the end
  for (const lab of document.querySelectorAll('label[title]')) {
    if (lab.querySelector(':scope > .tipBtn') || !lab.querySelector('input, select')) continue;
    const b = document.createElement('button');
    b.type = 'button'; b.className = 'tipBtn'; b.textContent = 'i';
    b.setAttribute('aria-label', 'Explanation');
    lab.appendChild(b);
  }
  // wide result tables scroll sideways inside their own box instead of widening the page
  for (const tb of document.querySelectorAll('.tabPanel table')) {
    if (tb.closest('#cellSection, .tableScroll, .twinSecBody, .sgdBox')) continue;
    const wrap = document.createElement('div');
    wrap.className = 'tableScroll';
    tb.parentNode.insertBefore(wrap, tb);
    wrap.appendChild(tb);
  }
}

new MutationObserver(muts => {
  if (!document.documentElement.classList.contains('touch')) return;
  for (const m of muts) {
    if (m.addedNodes.length || m.attributeName === 'title') { _touchScanSoon(); return; }
  }
}).observe(document.documentElement, { childList: true, subtree: true, attributes: true, attributeFilter: ['title'] });

// ---------- swipe on a layer plot ----------
// One finger sideways = previous / next layer (one step per ~36 px, so a long
// swipe scrubs through several). Vertical moves still scroll the page
// (touch-action: pan-y). Returns a function that tells whether the last
// gesture was a swipe, so the click that may follow does not pick a point.
function attachSwipeStep(el, step) {
  el.style.touchAction = 'pan-y pinch-zoom';
  let s = null, lastSwipe = 0;
  el.addEventListener('pointerdown', e => {
    if (e.pointerType === 'mouse') return;
    s = { id: e.pointerId, x: e.clientX, y: e.clientY, done: 0, horiz: null };
  });
  el.addEventListener('pointermove', e => {
    if (!s || e.pointerId !== s.id) return;
    const dx = e.clientX - s.x, dy = e.clientY - s.y;
    if (s.horiz === null && Math.hypot(dx, dy) > 10) s.horiz = Math.abs(dx) > Math.abs(dy);
    if (!s.horiz) return;
    const want = Math.trunc(-dx / 36); // finger to the left = next layer, like turning a page
    while (s.done < want) { step(1); s.done++; }
    while (s.done > want) { step(-1); s.done--; }
    if (s.done) lastSwipe = Date.now();
  });
  const end = e => { if (s && e.pointerId === s.id) { if (s.horiz) lastSwipe = Date.now(); s = null; } };
  el.addEventListener('pointerup', end);
  el.addEventListener('pointercancel', end);
  return () => Date.now() - lastSwipe < 400;
}

_touchApplyClass();
