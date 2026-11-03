// plot-zoom.js
// Click-to-enlarge for the small canvas plots: a click on a plot card opens
// it in a centered overlay over a dimmed, blurred page. The plot is
// RE-RENDERED at a larger zoom (not a stretched bitmap), so lines, markers
// and labels stay sharp. Closes on Esc, the × button, or a click outside
// the card / on the enlarged plot itself.
//
// Usage (see absences-viz.js / laue-viz.js):
//   attachPlotZoom(wrapEl, captionHtml, baseCssSize, (canvas, zoom) => drawPlot(canvas, ..., zoom));
// `draw` must render at `zoom` times its normal CSS size.

let _plotZoomOverlay = null;
let _plotZoomRedraw = null;
let _plotZoomKeyHandler = null;
let _plotZoomPrevOverflow = '';

function closePlotZoom() {
  if (!_plotZoomOverlay) return;
  _plotZoomOverlay.remove();
  _plotZoomOverlay = null;
  _plotZoomRedraw = null;
  document.body.style.overflow = _plotZoomPrevOverflow;
  if (_plotZoomKeyHandler) document.removeEventListener('keydown', _plotZoomKeyHandler);
  _plotZoomKeyHandler = null;
}

// Largest zoom that still fits the viewport (never below 1, capped so a
// huge screen doesn't create a gigantic bitmap).
function plotZoomFactor(baseSize) {
  const fit = Math.min(window.innerWidth - 64, window.innerHeight - 120);
  return Math.max(1, Math.min(3.5, fit / baseSize));
}

function openPlotZoom(titleHtml, baseSize, draw) {
  closePlotZoom();

  const overlay = document.createElement('div');
  overlay.className = 'plotZoomOverlay';
  overlay.setAttribute('role', 'dialog');
  overlay.setAttribute('aria-modal', 'true');

  const card = document.createElement('div');
  card.className = 'plotZoomCard';

  const head = document.createElement('div');
  head.className = 'plotZoomHead';
  const title = document.createElement('span');
  title.innerHTML = titleHtml;
  const closeBtn = document.createElement('button');
  closeBtn.type = 'button';
  closeBtn.className = 'plotZoomClose';
  closeBtn.setAttribute('aria-label', 'Close');
  closeBtn.textContent = '×';
  head.appendChild(title);
  head.appendChild(closeBtn);

  const canvas = document.createElement('canvas');
  canvas.className = 'plotZoomCanvas';
  card.appendChild(head);
  card.appendChild(canvas);
  overlay.appendChild(card);

  const redraw = () => draw(canvas, plotZoomFactor(baseSize));
  redraw();

  overlay.addEventListener('click', (e) => { if (e.target === overlay) closePlotZoom(); });
  canvas.addEventListener('click', closePlotZoom);
  closeBtn.addEventListener('click', closePlotZoom);
  _plotZoomKeyHandler = (e) => { if (e.key === 'Escape') closePlotZoom(); };
  document.addEventListener('keydown', _plotZoomKeyHandler);

  _plotZoomPrevOverflow = document.body.style.overflow;
  document.body.style.overflow = 'hidden'; // no scrolling of the page behind the overlay
  document.body.appendChild(overlay);
  _plotZoomOverlay = overlay;
  _plotZoomRedraw = redraw;
  closeBtn.focus();
}

function attachPlotZoom(wrapEl, titleHtml, baseSize, draw) {
  wrapEl.classList.add('zoomable');
  wrapEl.tabIndex = 0;
  wrapEl.setAttribute('role', 'button');
  wrapEl.title = 'Click to enlarge';
  wrapEl.addEventListener('click', () => openPlotZoom(titleHtml, baseSize, draw));
  wrapEl.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); openPlotZoom(titleHtml, baseSize, draw); }
  });
}

// keep an open overlay correct when the window is resized or the OS theme flips
window.addEventListener('resize', () => { if (_plotZoomRedraw) _plotZoomRedraw(); });
onPlotThemeChange(() => { if (_plotZoomRedraw) _plotZoomRedraw(); });
