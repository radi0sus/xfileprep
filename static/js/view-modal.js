// view-modal.js
// Enlarged view for the layer viewer's 2D layer and 3D sphere: the given
// page elements are MOVED into a centered card over the dimmed, blurred page
// (same look as plot-zoom.js) and moved back on close. Moving instead of
// copying keeps canvases, the WebGL context and all event listeners.
//
//   openViewModal({ title, top: [el, …], body: [el, …], onResize, onClose })
//   closeViewModal()   viewModalIsOpen()
// `top` goes above (e.g. the layer controls), `body` side by side below.
// Closes on ×, Esc or a click outside the card.

const _viewModal = { overlay: null, moved: [], key: null, prevOverflow: '', onClose: null, onResize: null };

function viewModalIsOpen() { return !!_viewModal.overlay; }

function openViewModal({ title, top = [], body = [], onResize = null, onClose = null }) {
  closeViewModal();
  const overlay = document.createElement('div');
  overlay.className = 'plotZoomOverlay';
  overlay.setAttribute('role', 'dialog');
  overlay.setAttribute('aria-modal', 'true');
  overlay.innerHTML = `
    <div class="plotZoomCard viewModalCard">
      <div class="plotZoomHead"><span><strong>${title}</strong></span>
        <button type="button" class="plotZoomClose" aria-label="Close">×</button></div>
      <div class="viewModalTop"></div>
      <div class="viewModalBody"></div>
    </div>`;
  const move = (el, target) => {
    if (!el || !el.parentNode) return;
    const mark = document.createComment('view-modal');
    el.parentNode.insertBefore(mark, el);
    _viewModal.moved.push({ el, mark });
    target.appendChild(el);
  };
  for (const el of top) move(el, overlay.querySelector('.viewModalTop'));
  for (const el of body) move(el, overlay.querySelector('.viewModalBody'));

  overlay.addEventListener('click', e => { if (e.target === overlay) closeViewModal(); });
  overlay.querySelector('.plotZoomClose').addEventListener('click', closeViewModal);
  _viewModal.key = e => { if (e.key === 'Escape') closeViewModal(); };
  document.addEventListener('keydown', _viewModal.key);
  _viewModal.prevOverflow = document.body.style.overflow;
  document.body.style.overflow = 'hidden';
  document.body.appendChild(overlay);
  _viewModal.overlay = overlay;
  _viewModal.onClose = onClose;
  _viewModal.onResize = onResize;
  if (onResize) onResize();
  overlay.querySelector('.plotZoomClose').focus();
}

function closeViewModal() {
  if (!_viewModal.overlay) return;
  for (const { el, mark } of _viewModal.moved.reverse()) {
    if (mark.parentNode) { mark.parentNode.insertBefore(el, mark); mark.remove(); }
  }
  _viewModal.moved = [];
  _viewModal.overlay.remove();
  _viewModal.overlay = null;
  document.removeEventListener('keydown', _viewModal.key);
  document.body.style.overflow = _viewModal.prevOverflow;
  const cb = _viewModal.onClose;
  _viewModal.onClose = null; _viewModal.onResize = null;
  if (cb) cb();
}

// Edge length (CSS px) of a square picture that fits the viewport next to a
// control column of `sideWidth` (stacked below it on narrow screens).
function viewModalSquareSize(sideWidth) {
  const narrow = window.innerWidth < 760;
  const fit = narrow ? Math.min(window.innerWidth - 48, window.innerHeight - 330)
                     : Math.min(window.innerWidth - 120 - sideWidth, window.innerHeight - 290);
  return Math.max(260, Math.floor(fit));
}

window.addEventListener('resize', () => { if (_viewModal.onResize) _viewModal.onResize(); });
