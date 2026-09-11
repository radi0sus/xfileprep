(function () {
  const tabBar = document.getElementById('tabBar');
  if (!tabBar) return;

  const buttons = Array.from(tabBar.querySelectorAll('.tabBtn'));
  const panels = Array.from(document.querySelectorAll('.tabPanel'));

  function switchTab(tabId) {
    buttons.forEach(btn => btn.classList.toggle('active', btn.dataset.tab === tabId));
    panels.forEach(panel => panel.classList.toggle('active', panel.id === tabId));
  }

  buttons.forEach(btn => btn.addEventListener('click', () => switchTab(btn.dataset.tab)));

  // Exposed so app.js can switch tabs itself (e.g. jump to the space-group
  // tab once the auto-run pipeline lands a candidate table there).
  window.switchTab = switchTab;
})();
