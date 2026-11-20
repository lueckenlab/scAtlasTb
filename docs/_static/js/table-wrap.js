// Add line-break opportunities (<wbr>) so that long entries wrap at natural positions
// instead of overflowing the page or widening table columns:
// - in table cells: after "/", "_", ",", ".", "(", "=", ":", "~" and "-" (e.g. bio_conservation)
// - in inline code anywhere in the text (long file paths, slot names)
// Code blocks (<pre>) are not touched.
const BREAK_AFTER = /(?<=[\/_,.(=:~-])(?=\S)/;

function addBreaks(element) {
  const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
  const nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  nodes.forEach((node) => {
    if (node.parentElement.closest("pre")) return;
    const parts = node.textContent.split(BREAK_AFTER);
    if (parts.length < 2) return;
    const fragment = document.createDocumentFragment();
    parts.forEach((part, i) => {
      fragment.append(part);
      if (i < parts.length - 1) fragment.append(document.createElement("wbr"));
    });
    node.replaceWith(fragment);
  });
}

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("article table td, article table th").forEach(addBreaks);
  document.querySelectorAll("article code.literal").forEach((code) => {
    if (!code.closest("table")) addBreaks(code);
  });
});
