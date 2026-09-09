/* Page-only wiring; theme semantics remain in the shared theme initializer. */
(() => {
  const selector = document.querySelector("[data-lcm-theme-selector]");
  const theme = globalThis.LineCableModelsTheme;
  if (!selector || !theme) return;
  selector.value = theme.preference();
  selector.addEventListener("change", event => theme.select(event.currentTarget.value));
})();
