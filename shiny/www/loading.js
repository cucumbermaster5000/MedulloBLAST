// Do not show readiness while Shiny is rendering visible outputs.
let medulloRendering = false;
let medulloRenderCheck;
function finishMedulloRendering() {
  if (medulloRendering) return;
  const pending = Array.from(document.querySelectorAll('.shiny-plot-output img'))
    .some(img => img.getClientRects().length && !img.complete);
  if (pending) {
    medulloRenderCheck = setTimeout(finishMedulloRendering, 100);
    return;
  }
  requestAnimationFrame(function () {
    if (!medulloRendering) document.documentElement.classList.remove('analysis-rendering-busy');
  });
}
$(document).on('shiny:busy', function () {
  medulloRendering = true;
  clearTimeout(medulloRenderCheck);
  document.documentElement.classList.add('analysis-rendering-busy');
}).on('shiny:idle', function () {
  medulloRendering = false;
  finishMedulloRendering();
});
