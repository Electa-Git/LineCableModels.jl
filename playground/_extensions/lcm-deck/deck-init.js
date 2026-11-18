// Normalize PDF entry points before Reveal reads its URL. LCM owns pagination.
(() => {
  const url = new URL(location.href);
  const print = url.searchParams.has('print-pdf') ||
    url.searchParams.get('view') === 'print' || url.searchParams.has('lcm-print');
  if (print) {
    url.searchParams.delete('print-pdf');
    url.searchParams.delete('view');
    url.searchParams.set('lcm-print', '');
    history.replaceState(history.state, '', url);
    document.documentElement.classList.add('lcm-print-view', 'lcm-static-live', 'print-pdf');
  }
})();
