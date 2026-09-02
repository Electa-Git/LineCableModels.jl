(() => {
  const initialize = (shell) => {
    if (shell.dataset.lcCollapsibleSidebarReady === "true") return;
    shell.dataset.lcCollapsibleSidebarReady = "true";

    const toggle = shell.querySelector("[data-sidebar-toggle]");
    const stateLabel = shell.querySelector("[data-sidebar-state-label]");
    const pageTitle = shell.querySelector("[data-demo-page-title]");
    const navigationItems = shell.querySelectorAll(
      ".lc-cs-nav-item:not(:disabled)"
    );

    const synchronize = () => {
      const collapsed = shell.dataset.sidebarState === "collapsed";
      toggle.setAttribute("aria-expanded", String(!collapsed));
      toggle.setAttribute(
        "aria-label",
        collapsed ? "Expand navigation" : "Collapse navigation"
      );
      toggle.dataset.tooltip = collapsed
        ? "Expand navigation"
        : "Collapse navigation";
      stateLabel.textContent = collapsed
        ? "Collapsed rail · 4.25rem reserved"
        : "Expanded navigation · 18rem reserved";
    };

    toggle.addEventListener("click", () => {
      shell.dataset.sidebarState =
        shell.dataset.sidebarState === "collapsed" ? "expanded" : "collapsed";
      synchronize();
    });

    navigationItems.forEach((item) => {
      item.addEventListener("click", () => {
        navigationItems.forEach((candidate) => {
          candidate.classList.remove("is-active");
          candidate.removeAttribute("aria-current");
        });
        item.classList.add("is-active");
        item.setAttribute("aria-current", "page");
        pageTitle.textContent = item.dataset.demoTitle;
      });
    });

    synchronize();
  };

  const initializeAll = () => {
    document
      .querySelectorAll("[data-lc-collapsible-sidebar]")
      .forEach(initialize);
  };

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initializeAll, { once: true });
  } else {
    initializeAll();
  }
})();
