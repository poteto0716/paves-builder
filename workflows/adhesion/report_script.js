(function () {
  const D = JSON.parse(document.getElementById('data').textContent);
  const $ = id => document.getElementById(id);
  $('eyebrow').textContent = D.meta.eyebrow;
  $('title').textContent = D.meta.title;
  $('lede').innerHTML = D.meta.lede;
  $('hint').textContent = D.meta.hint;
  $('metrics').innerHTML = D.metrics.map(m =>
    `<div class="metric"><div class="k">${m.k}</div><div class="v">${m.v}<small>${m.u || ''}</small></div><div class="d">${m.d || ''}</div></div>`).join('');
  $('notes').innerHTML = D.notes.map(n => `<div><h4>${n.h}</h4><p>${n.p}</p></div>`).join('');

  // Displayed atoms: D.shown[k] is the original index of display atom k.
  const shown = D.shown, NSI = D.n_si, elem = D.elements;
  const pulled = new Set(D.pulled), fixed = new Set(D.fixed);
  const firstPolymer = shown.findIndex(i => i >= NSI);
  const COLORS = { Si: '#D9B24C', Osi: '#D8665A', C: '#8C98A3', O: '#E0544A', H: '#E9EEF1',
                   pulled: '#4FA3F0', fixed: '#8A6420', contact: '#39C27A' };
  const RADIUS = { Si: 0.9, O: 0.7, C: 0.72, H: 0.42 };
  const ELEM = { S: 'Si', O: 'O', C: 'C', H: 'H' };
  const legend = [['Si', COLORS.Si], ['O (silica)', COLORS.Osi], ['C', COLORS.C], ['O (PMMA)', COLORS.O]];
  if (D.has_h) legend.push(['H', COLORS.H]);
  legend.push(['within 3 Å of silica', COLORS.contact]);
  if (D.pulled.length) legend.push(['pulled group', COLORS.pulled]);
  legend.push(['fixed atoms', COLORS.fixed]);
  $('legend').innerHTML = legend.map(([n, c]) => `<span><i class="swatch" style="background:${c}"></i>${n}</span>`).join('');
  if (!D.pulled.length) $('t-pulled').closest('label').hidden = true;

  const decode = b64 => {
    const s = atob(b64), u = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
    const q = new Int16Array(u.buffer), out = new Float32Array(q.length);
    for (let i = 0; i < q.length; i++) out[i] = q[i] * D.quantum;
    return out;
  };
  D.stages.forEach(st => st.frames.forEach(f => {
    f.xyz = decode(f.x); delete f.x;
    f.contact = f.c ? new Set(f.c) : new Set(); delete f.c;
  }));

  const opts = () => ({ silica: $('t-silica').checked, pulled: $('t-pulled').checked, fixed: $('t-fixed').checked,
                        box: $('t-box').checked, contact: $('t-contact').checked });
  let stageIdx = 0, frameIdx = 0, timer = null, view = null, lastStage = -1;
  const viewer = $3Dmol.createViewer($('viewer'), { backgroundColor: '#0B1114', antialias: true });

  function colorOf(i, st, fr, o) {
    if (o.pulled && st.pulled && pulled.has(i)) return COLORS.pulled;
    if (o.contact && fr.contact.has(i)) return COLORS.contact;
    if (o.fixed && i < NSI && fixed.has(i)) return COLORS.fixed;
    const e = ELEM[elem[i]];
    if (e === 'O' && i < NSI) return COLORS.Osi;
    return COLORS[e];
  }

  function draw() {
    const st = D.stages[stageIdx], fr = st.frames[frameIdx], o = opts();
    const k0 = st.subset === 'polymer' ? firstPolymer : 0;
    const idx = [], lines = [];
    for (let k = k0; k < shown.length; k++) {
      const i = shown[k], j = k - k0;
      if (!o.silica && i < NSI) continue;
      idx.push(i);
      lines.push(`${ELEM[elem[i]]} ${fr.xyz[3 * j].toFixed(2)} ${fr.xyz[3 * j + 1].toFixed(2)} ${fr.xyz[3 * j + 2].toFixed(2)}`);
    }
    if (lastStage === stageIdx) view = viewer.getView();
    viewer.removeAllModels(); viewer.removeAllShapes();
    const m = viewer.addModel(`${lines.length}\n${st.id}\n${lines.join('\n')}`, 'xyz');
    for (const e of ['Si', 'O', 'C', 'H'])
      m.setStyle({ elem: e }, { sphere: { radius: RADIUS[e], colorfunc: a => colorOf(idx[a.index], st, fr, o) } });
    if (o.box) viewer.addBox({ corner: { x: 0, y: 0, z: 0 }, dimensions: { w: st.box[0], h: st.box[1], d: st.box[2] }, wireframe: true, color: '#5E7380' });
    if (lastStage !== stageIdx) {
      viewer.zoomTo(); viewer.rotate(-90, { x: 1, y: 0, z: 0 }); viewer.rotate(-25, { x: 0, y: 0, z: 1 }); viewer.zoom(0.95);
      lastStage = stageIdx;
    } else if (view) viewer.setView(view);
    viewer.render();
    $('ov-title').textContent = st.label;
    $('ov-sub').textContent = st.sub;
    $('ov-frame').innerHTML = `${lines.length.toLocaleString()} atoms shown<br>${st.box.map(v => v.toFixed(1)).join(' × ')} Å`;
    $('frame').max = st.frames.length - 1; $('frame').value = frameIdx;
    $('frame-label').textContent = `frame ${frameIdx + 1} / ${st.frames.length} · ${fr.t}`;
    $('play').disabled = st.frames.length < 2;
    document.querySelectorAll('.steps button').forEach((b, k) => b.setAttribute('aria-current', k === stageIdx ? 'true' : 'false'));
  }

  $('steps').innerHTML = D.stages.map((st, k) =>
    `<li><button type="button" data-k="${k}"><span class="n">${st.num}</span><span class="t">${st.label}</span><span class="f">${st.frames.length > 1 ? st.frames.length + ' frames' : ''}</span><span class="s">${st.sub}</span></button></li>`).join('');
  $('steps').addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; stop(); stageIdx = +b.dataset.k; frameIdx = 0; draw(); });
  $('frame').addEventListener('input', e => { frameIdx = +e.target.value; draw(); });
  ['t-silica', 't-pulled', 't-fixed', 't-box', 't-contact'].forEach(id => $(id).addEventListener('change', draw));
  function stop() { if (timer) { clearInterval(timer); timer = null; $('play').textContent = 'Play'; } }
  $('play').addEventListener('click', () => {
    if (timer) return stop();
    $('play').textContent = 'Pause';
    timer = setInterval(() => { frameIdx = (frameIdx + 1) % D.stages[stageIdx].frames.length; draw(); }, 900);
  });
  stageIdx = D.start_stage; frameIdx = D.stages[stageIdx].frames.length - 1;
  draw();

  // --- charts ------------------------------------------------------------------
  function nice(lo, hi, n) {
    const span = hi - lo || Math.abs(hi) || 1, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map(k => k * mag).find(s => s >= raw);
    const a = Math.floor(lo / step) * step, b = Math.ceil(hi / step) * step, t = [];
    for (let v = a; v <= b + step / 2; v += step) t.push(+v.toFixed(10));
    return { a, b, t, step };
  }
  function chart(el, cfg) {
    const W = 320, H = 190, L = 46, R = cfg.series.length > 1 ? 62 : 14, T = 8, B = 30;
    const xs = cfg.series[0].pts.map(p => p[0]);
    const all = cfg.series.flatMap(s => s.pts.map(p => p[1])).concat(cfg.ref != null ? [cfg.ref] : []);
    const X = nice(Math.min(...xs), Math.max(...xs), 4), Y = nice(Math.min(...all), Math.max(...all), 4);
    const sx = v => L + (v - X.a) / (X.b - X.a) * (W - L - R), sy = v => T + (1 - (v - Y.a) / (Y.b - Y.a)) * (H - T - B);
    const fmt = (v, s) => Math.abs(s) >= 1 ? v.toFixed(0) : v.toFixed(Math.min(4, Math.ceil(-Math.log10(s))));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${cfg.title}"><g class="grid">`;
    Y.t.forEach(v => svg += `<line x1="${L}" x2="${W - R}" y1="${sy(v)}" y2="${sy(v)}"/>`);
    svg += `</g><g class="axis"><line x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}"/></g>`;
    Y.t.forEach(v => svg += `<text x="${L - 6}" y="${sy(v) + 3.5}" text-anchor="end">${fmt(v, Y.step)}</text>`);
    X.t.forEach(v => svg += `<text x="${sx(v)}" y="${H - B + 15}" text-anchor="middle">${fmt(v, X.step)}</text>`);
    svg += `<text x="${(L + W - R) / 2}" y="${H - 2}" text-anchor="middle">${cfg.x}</text>`;
    if (cfg.ref != null) svg += `<line x1="${L}" x2="${W - R}" y1="${sy(cfg.ref)}" y2="${sy(cfg.ref)}" stroke="var(--ink)" stroke-dasharray="5 4" stroke-width="1"/><text x="${W - R}" y="${sy(cfg.ref) - 5}" text-anchor="end" style="fill:var(--ink)">${cfg.refLabel}</text>`;
    cfg.series.forEach(s => {
      const d = s.pts.map((p, k) => `${k ? 'L' : 'M'}${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`).join('');
      svg += `<path d="${d}" fill="none" stroke="var(${s.color})" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>`;
      if (s.pts.length <= 12) s.pts.forEach(p => svg += `<circle cx="${sx(p[0])}" cy="${sy(p[1])}" r="4" fill="var(${s.color})" stroke="var(--panel)" stroke-width="2"/>`);
      if (cfg.series.length > 1) { const p = s.pts[s.pts.length - 1]; svg += `<text x="${sx(p[0]) + 8}" y="${sy(p[1]) + (s.dy || 0) + 3.5}" style="fill:var(--ink)">${s.short}</text>`; }
    });
    svg += `<line class="xh" x1="0" x2="0" y1="${T}" y2="${H - B}" stroke="var(--faint)" stroke-dasharray="3 3" visibility="hidden"/>`;
    cfg.series.forEach((s, k) => svg += `<circle class="hd${k}" r="4.5" fill="var(${s.color})" stroke="var(--panel)" stroke-width="2" visibility="hidden"/>`);
    svg += `<rect x="${L}" y="${T}" width="${W - L - R}" height="${H - T - B}" fill="transparent"/></svg>`;
    const keys = cfg.series.length > 1 ? `<div class="keys">${cfg.series.map(s => `<span><i style="background:var(${s.color})"></i>${s.name}</span>`).join('')}</div>` : '';
    el.innerHTML = `<h3>${cfg.title}</h3><p class="sub">${cfg.sub}</p>${keys}${svg}<div class="tip" hidden></div>`;
    const sv = el.querySelector('svg'), tip = el.querySelector('.tip'), xh = sv.querySelector('.xh');
    sv.addEventListener('mousemove', ev => {
      const r = sv.getBoundingClientRect(), x = (ev.clientX - r.left) / r.width * W;
      let k = 0, best = Infinity; xs.forEach((v, j) => { const dd = Math.abs(sx(v) - x); if (dd < best) { best = dd; k = j; } });
      xh.setAttribute('x1', sx(xs[k])); xh.setAttribute('x2', sx(xs[k])); xh.setAttribute('visibility', 'visible');
      cfg.series.forEach((s, j) => { const c = sv.querySelector('.hd' + j); c.setAttribute('cx', sx(s.pts[k][0])); c.setAttribute('cy', sy(s.pts[k][1])); c.setAttribute('visibility', 'visible'); });
      tip.innerHTML = `${cfg.x.split(' (')[0]} ${xs[k].toFixed(cfg.xd ?? 2)}<br>` + cfg.series.map(s => `${s.short}: ${s.pts[k][1].toFixed(cfg.yd ?? 3)} ${cfg.unit}`).join('<br>');
      tip.hidden = false; tip.style.left = (sx(xs[k]) / W * r.width + sv.offsetLeft) + 'px'; tip.style.top = (sv.offsetTop + sy(cfg.series[0].pts[k][1]) / H * r.height) + 'px';
    });
    sv.addEventListener('mouseleave', () => { tip.hidden = true; xh.setAttribute('visibility', 'hidden'); sv.querySelectorAll('[class^=hd]').forEach(c => c.setAttribute('visibility', 'hidden')); });
  }
  D.charts.forEach(cfg => { const el = document.createElement('div'); el.className = 'chart'; $('charts').appendChild(el); chart(el, cfg); });
})();
