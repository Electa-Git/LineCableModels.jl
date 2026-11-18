(() => {
    const accepted = new Set([
        'lcm:slide-enter',
        'lcm:slide-leave',
        'lcm:viewport-settling',
        'lcm:viewport-settled',
        'lcm:print-mode',
        'lcm:host-ready',
        'lcm:theme'
    ]);
    const counts = Object.create(null);
    const root = document.documentElement;

    window.addEventListener('message', event => {
        if (event.source !== window.parent || event.origin !== window.location.origin) return;
        const message = event.data;
        if (!message || message.namespace !== 'lcm-deck' || !accepted.has(message.type)) return;
        root.classList.add('lcm-deck-embedded');
        if (message.type === 'lcm:host-ready') {
            window.parent.postMessage({
                namespace: 'lcm-deck', type: 'lcm:child-ready'
            }, window.location.origin);
            return;
        }
        if (message.type === 'lcm:theme') {
            const detail = message.detail || {};
            if (!['dark', 'light', 'system'].includes(detail.preference) ||
                !['dark', 'light'].includes(detail.resolved)) return;
            root.dataset.lcmTheme = detail.preference;
            root.dataset.lcmResolvedTheme = detail.resolved;
            window.dispatchEvent(new CustomEvent('lcm:theme-changed', { detail }));
            return;
        }
        counts[message.type] = (counts[message.type] || 0) + 1;
        const eventName = message.type.replace(/^lcm:/, '').split('-')
            .map(part => part.charAt(0).toUpperCase() + part.slice(1)).join('');
        root.dataset.lcmDeckEvent = message.type;
        root.dataset.lcmDeckEventCount = String(counts[message.type]);
        root.dataset[`lcmDeck${eventName}Count`] = String(counts[message.type]);
        if (message.type === 'lcm:slide-enter') root.classList.add('lcm-deck-visible');
        if (message.type === 'lcm:slide-leave') root.classList.remove('lcm-deck-visible');
        if (message.type === 'lcm:viewport-settling') root.classList.add('lcm-deck-settling');
        if (message.type === 'lcm:viewport-settled') root.classList.remove('lcm-deck-settling');
        window.dispatchEvent(new CustomEvent(message.type, { detail: message.detail || {} }));
    });

    if (window.parent !== window) {
        window.parent.postMessage({
            namespace: 'lcm-deck',
            type: 'lcm:child-ready'
        }, window.location.origin);
    }
})();
