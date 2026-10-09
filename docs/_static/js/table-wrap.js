// Add line-break opportunities after "/", "_", ",", ".", "(", "=" and ":" in table cells,
// so that long paths, slot names and identifiers (e.g. bio_conservation) wrap at natural
// positions instead of widening their column or breaking mid-word.
document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("article table td, article table th").forEach((cell) => {
    const walker = document.createTreeWalker(cell, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach((node) => {
      const parts = node.textContent.split(/(?<=[\/_,.(=:])(?=\S)/);
      if (parts.length < 2) return;
      const fragment = document.createDocumentFragment();
      parts.forEach((part, i) => {
        fragment.append(part);
        if (i < parts.length - 1) fragment.append(document.createElement("wbr"));
      });
      node.replaceWith(fragment);
    });
  });
});
