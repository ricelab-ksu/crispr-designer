/* ==========================================================================
 * www/pjd_assets/pjd_map.js — interaction layer for the plasmid map
 * ==========================================================================
 * Vanilla JS, no dependencies. Adds to the server-rendered SVG:
 *   - hover tooltip            (data-pjd-* attributes written by R)
 *   - wheel / button zoom+pan  (svg viewBox)
 *   - click -> Shiny input     ("pjd_map_selected_feature")
 * Everything degrades gracefully: with JS disabled the SVG still renders and
 * its native <title> elements show tooltips.
 * ========================================================================== */
(function () {
  'use strict';
  if (window.__pjdMapInit) { return; }
  window.__pjdMapInit = true;

  function featureFrom(target) {
    return target && target.closest ? target.closest('[data-pjd-map="feature"]') : null;
  }
  function svgOf(target) {
    return target && target.closest ? target.closest('svg.pjd-map-svg') : null;
  }
  function boxOf(svg) {
    return svg && svg.closest ? svg.closest('.pjd-interactive-map') : null;
  }

  function tooltip(box) {
    return box ? box.querySelector('.pjd-map-tooltip') : null;
  }

  function showTooltip(el, ev) {
    var svg = svgOf(el), box = boxOf(svg), tip = tooltip(box);
    if (!tip) { return; }
    tip.innerHTML =
      '<strong>' + (el.getAttribute('data-label') || '') + '</strong>' +
      '<span>' + (el.getAttribute('data-type') || '') + '</span>' +
      '<span>' + (el.getAttribute('data-start') || '?') + '..' +
      (el.getAttribute('data-end') || '?') + ' (' +
      (el.getAttribute('data-strand') || '+') + ')</span>' +
      '<span>' + (el.getAttribute('data-length') || '?') + ' bp</span>';
    tip.hidden = false;
    moveTooltip(box, tip, ev);
  }

  function moveTooltip(box, tip, ev) {
    if (!tip || tip.hidden || !box) { return; }
    var r = box.getBoundingClientRect();
    var x = ev.clientX - r.left + 14;
    var y = ev.clientY - r.top + 14;
    var maxX = r.width - tip.offsetWidth - 8;
    var maxY = r.height - tip.offsetHeight - 8;
    tip.style.left = Math.max(8, Math.min(x, maxX)) + 'px';
    tip.style.top = Math.max(8, Math.min(y, maxY)) + 'px';
  }

  document.addEventListener('mouseover', function (ev) {
    var f = featureFrom(ev.target);
    if (f) { showTooltip(f, ev); }
  });

  document.addEventListener('mousemove', function (ev) {
    var f = featureFrom(ev.target);
    if (!f) { return; }
    var tip = tooltip(boxOf(svgOf(f)));
    moveTooltip(boxOf(svgOf(f)), tip, ev);
  });

  document.addEventListener('mouseout', function (ev) {
    var f = featureFrom(ev.target);
    if (!f) { return; }
    var tip = tooltip(boxOf(svgOf(f)));
    if (tip) { tip.hidden = true; }
  });

  document.addEventListener('click', function (ev) {
    var f = featureFrom(ev.target);
    if (f && window.Shiny && window.Shiny.setInputValue) {
      window.Shiny.setInputValue('pjd_map_selected_feature', {
        label: f.getAttribute('data-label'),
        type: f.getAttribute('data-type'),
        start: parseInt(f.getAttribute('data-start'), 10),
        end: parseInt(f.getAttribute('data-end'), 10),
        strand: f.getAttribute('data-strand'),
        length: parseInt(f.getAttribute('data-length'), 10)
      }, { priority: 'event' });
    }
  });

  function setViewBox(svg, vb) {
    svg.setAttribute('viewBox', vb.join(' '));
  }
  function currentViewBox(svg) {
    var raw = svg.getAttribute('data-pjd-original-viewbox') || svg.getAttribute('viewBox');
    return raw.split(/[\s,]+/).map(Number);
  }
  function zoom(svg, factor, cx, cy) {
    if (!svg.getAttribute('data-pjd-original-viewbox')) {
      svg.setAttribute('data-pjd-original-viewbox', svg.getAttribute('viewBox'));
    }
    var vb = currentViewBox(svg);
    var nw = vb[2] * factor, nh = vb[3] * factor;
    var px = (cx == null) ? vb[2] / 2 : cx;
    var py = (cy == null) ? vb[3] / 2 : cy;
    vb[0] = vb[0] + (vb[2] - nw) * (px / vb[2]);
    vb[1] = vb[1] + (vb[3] - nh) * (py / vb[3]);
    vb[2] = nw; vb[3] = nh;
    setViewBox(svg, vb);
  }

  document.addEventListener('click', function (ev) {
    var btn = ev.target.closest ? ev.target.closest('[data-pjd-action]') : null;
    if (!btn) { return; }
    var svg = boxOf(btn) && boxOf(btn).querySelector('svg.pjd-map-svg');
    if (!svg) { return; }
    var action = btn.getAttribute('data-pjd-action');
    if (action === 'zoom-in') { zoom(svg, 0.8); }
    else if (action === 'zoom-out') { zoom(svg, 1.25); }
    else if (action === 'reset') {
      var orig = svg.getAttribute('data-pjd-original-viewbox');
      if (orig) { setViewBox(svg, orig.split(/[\s,]+/).map(Number)); }
    }
  });

  document.addEventListener('wheel', function (ev) {
    var svg = svgOf(ev.target);
    if (!svg) { return; }
    ev.preventDefault();
    var rect = svg.getBoundingClientRect();
    var cx = (ev.clientX - rect.left) / rect.width * 100;
    var cy = (ev.clientY - rect.top) / rect.height * 100;
    zoom(svg, ev.deltaY > 0 ? 1.12 : 0.89, cx, cy);
  }, { passive: false });
})();
