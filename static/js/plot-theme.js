// plot-theme.js
// Canvas plots (absences-viz.js, e-statistics.js) draw with plain 2D
// context calls, so they can't pick up CSS the way the rest of the
// page does. This gives them a small palette that mirrors the
// light/dark variables in style.css, plus a hook so a plot re-renders
// itself if the OS-level color scheme flips while the page is open
// (no reload needed).

function plotIsDarkMode() {
  return !!(window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches);
}

// Mirrors the --bg/--text/etc. custom properties in style.css. Kept as
// a plain object (not read live from getComputedStyle) so the exact
// same values used for the CSS dark theme are guaranteed available
// here even before the stylesheet has painted anything.
function getPlotColors() {
  return plotIsDarkMode() ? {
    bg: '#1c1f26',
    present: '#4fce84',   // significant reflection measured (dark)
    absent: '#4a4f59',    // measured but weak (dark)
    axis: '#9aa1ac',
    grid: '#333844',
    text: '#c7cbd3',
    fitLine: '#ff8a80',
    point: 'rgba(157, 181, 255, 0.8)',
  } : {
    bg: '#ffffff',
    present: '#0a8a3c',
    absent: '#cfcfcf',
    axis: '#888',
    grid: '#eee',
    text: '#555',
    fitLine: '#c33',
    point: 'rgba(20,90,180,0.75)',
  };
}

// Registers a redraw callback that fires whenever the OS color scheme
// changes. Returns nothing; callers just pass a zero-arg function that
// re-runs their last render with whatever data they last used.
function onPlotThemeChange(redraw) {
  if (!window.matchMedia) return;
  const mq = window.matchMedia('(prefers-color-scheme: dark)');
  const handler = () => redraw();
  if (mq.addEventListener) mq.addEventListener('change', handler);
  else if (mq.addListener) mq.addListener(handler); // older Safari
}
