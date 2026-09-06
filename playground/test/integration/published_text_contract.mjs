// Shared assertion for Quarto documents and decks; never injected into widgets.
export async function assertPublishedText(devtools, selector) {
  const evaluate = async expression => {
    const result = await devtools.command("Runtime.evaluate", {
      expression, returnByValue: true, awaitPromise: true,
    });
    if (result.exceptionDetails) throw new Error(result.exceptionDetails.text);
    return result.result.value;
  };
  const assert = (condition, message) => {
    if (!condition) throw new Error(message);
  };
  const sample = await evaluate(`(() => {
    const heading = document.querySelector(${JSON.stringify(selector)});
    const walker = document.createTreeWalker(heading, NodeFilter.SHOW_TEXT);
    const range = document.createRange();
    range.selectNodeContents(walker.nextNode());
    const box = range.getClientRects()[0];
    return {
      caret: getComputedStyle(heading).caretColor,
      x1: box.left + 2, x2: Math.min(box.right - 2, box.left + 180),
      y: box.top + box.height / 2,
      loaded: [...document.styleSheets].some(sheet => sheet.href?.endsWith('/published-text.css'))
    };
  })()`);
  assert(sample.loaded, "published surface omitted the shared text policy");
  assert(sample.caret === "rgba(0, 0, 0, 0)", "published heading paints a browsing caret");
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mousePressed", x: sample.x1, y: sample.y, button: "left", buttons: 1, clickCount: 1,
  });
  for (let step = 1; step <= 5; step++) {
    await devtools.command("Input.dispatchMouseEvent", {
      type: "mouseMoved", x: sample.x1 + (sample.x2 - sample.x1) * step / 5,
      y: sample.y, button: "left", buttons: 1,
    });
  }
  await devtools.command("Input.dispatchMouseEvent", {
    type: "mouseReleased", x: sample.x2, y: sample.y, button: "left", buttons: 0, clickCount: 1,
  });
  assert(await evaluate("getSelection().toString().trim().length > 0"),
    "mouse dragging no longer selects published text");

  // Ordinary page inputs and the three valid contenteditable modes retain their
  // caret. Read-only descendants must not inherit an editor's visible caret.
  try {
    await evaluate(`(() => {
      const fixture = document.createElement('div');
      fixture.id = 'published-text-regression';
      fixture.style.cssText = 'position:fixed;top:10px;left:10px;z-index:10000';
      fixture.innerHTML = '<p>Read-only prose <code>code</code></p>' +
        '<button>Action</button><a href="#">Navigation</a>' +
        '<input aria-label="Test text"><textarea></textarea>' +
        '<input readonly value="Read only"><input disabled>' +
        '<div contenteditable="true"><span>Editor</span>' +
          '<span contenteditable="false">Read only island</span></div>' +
        '<div contenteditable="">Editor</div><div contenteditable="plaintext-only">Editor</div>';
      document.body.append(fixture);
      fixture.querySelector('input').focus();
    })()`);
    const controls = await evaluate(`(() => {
      const fixture = document.querySelector('#published-text-regression');
      const carets = selector => [...fixture.querySelectorAll(selector)]
        .map(node => getComputedStyle(node).caretColor);
      return {
        reading: carets('p, code, button, a, input[readonly], input[disabled], [contenteditable="false"]'),
        editing: carets('input:not([readonly]):not([disabled]), textarea, [contenteditable]:not([contenteditable="false"])')
      };
    })()`);
    assert(controls.reading.every(value => value === "rgba(0, 0, 0, 0)"),
      "non-input published content paints a caret");
    assert(controls.editing.every(value => value !== "rgba(0, 0, 0, 0)"),
      "published text-entry control lost its caret");
    await devtools.command("Input.insertText", { text: "entry" });
    assert(await evaluate("document.querySelector('#published-text-regression input').value === 'entry'"),
      "published text input stopped accepting input");
  } finally {
    await evaluate("document.querySelector('#published-text-regression')?.remove(); getSelection().removeAllRanges()");
  }
}
