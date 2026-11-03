/* OmicLearningHub · Ask
 * 1. Search: instant keyword (BM25) search over assets/ask/index.json, built from the site's pages.
 * 2. Optional answer: a small open model (Qwen2.5 via WebLLM) runs on the visitor's own GPU (WebGPU).
 *    It downloads only when someone clicks "Generate answer"; then the browser caches it.
 *    Questions never leave the browser. The model is told to answer only from the matched sections.
 */
(function () {
  "use strict";
  var SCRIPT = document.currentScript;
  var BASE = SCRIPT ? new URL(".", SCRIPT.src).href : "assets/ask/";
  var WEBLLM_URL = "https://cdn.jsdelivr.net/npm/@mlc-ai/web-llm@0.2.86/+esm";
  var MODELS = [
    { key: "1.5b", label: "Qwen2.5 1.5B · better answers (~1 GB download)", f16: "Qwen2.5-1.5B-Instruct-q4f16_1-MLC", f32: "Qwen2.5-1.5B-Instruct-q4f32_1-MLC" },
    { key: "0.5b", label: "Qwen2.5 0.5B · faster, simpler (~0.4 GB)", f16: "Qwen2.5-0.5B-Instruct-q4f16_1-MLC", f32: "Qwen2.5-0.5B-Instruct-q4f32_1-MLC" }
  ];
  var EXAMPLES = [
    "Why use TMM normalisation?",
    "How does Harmony integration work?",
    "Why is the mouse the replicate, not the cell?",
    "What does GSEA rank genes by?"
  ];
  var TOP_K = 6, CONTEXT_K = 4, CONTEXT_CHARS = 900;
  var STOP = new Set(("a an and are as at be but by can do does for from how i if in into is it its me my of on or " +
    "so than that the their them then there these this to was we what when where which who why will with you your " +
    "should would could use using used").split(" "));

  var docs = null, idf = null, avgdl = 0, loading = null;
  var engine = null, engineKey = null, busy = false;
  var lastHits = [], lastQuestion = "";

  // ---------- helpers ----------
  function esc(s) { return String(s).replace(/[&<>"']/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]; }); }
  function tokens(s) {
    return (s.toLowerCase().match(/[a-z0-9][a-z0-9+\-_.]*[a-z0-9+]|[a-z0-9]/g) || [])
      .map(function (t) { return t.replace(/^(.{4,}?)(es|s)$/, "$1"); })      // light plural stemming
      .filter(function (t) { return !STOP.has(t); });
  }
  function store(k, v) { try { if (v === undefined) return localStorage.getItem(k); localStorage.setItem(k, v); } catch (e) { return null; } }
  function el(tag, attrs, html) {
    var n = document.createElement(tag);
    if (attrs) for (var k in attrs) n.setAttribute(k, attrs[k]);
    if (html != null) n.innerHTML = html;
    return n;
  }

  // ---------- index + BM25 ----------
  function loadIndex() {
    if (loading) return loading;
    loading = fetch(BASE + "index.json").then(function (r) {
      if (!r.ok) throw new Error("index " + r.status);
      return r.json();
    }).then(function (list) {
      var df = Object.create(null), total = 0;
      docs = list.map(function (d) {
        var toks = tokens(d.title + " " + d.title + " " + d.text);       // title counts twice
        var tf = Object.create(null);
        toks.forEach(function (t) { tf[t] = (tf[t] || 0) + 1; });
        Object.keys(tf).forEach(function (t) { df[t] = (df[t] || 0) + 1; });
        total += toks.length;
        return { d: d, tf: tf, len: toks.length, tt: new Set(tokens(d.title)) };
      });
      avgdl = total / Math.max(docs.length, 1);
      idf = Object.create(null);
      var N = docs.length;
      Object.keys(df).forEach(function (t) { idf[t] = Math.log(1 + (N - df[t] + 0.5) / (df[t] + 0.5)); });
      return docs;
    });
    return loading;
  }
  function search(q, k) {
    var qt = Array.from(new Set(tokens(q))), k1 = 1.4, b = 0.75;
    if (!qt.length) return [];
    var here = location.pathname.split("/").pop() || "index.html";
    return docs.map(function (x) {
      var s = 0, hit = 0, th = 0;
      qt.forEach(function (t) {
        var f = x.tf[t]; if (!f) return;
        hit++; if (x.tt.has(t)) th++;
        s += idf[t] * (f * (k1 + 1)) / (f + k1 * (1 - b + b * x.len / avgdl));
      });
      if (s > 0) {
        s *= 1 + 0.15 * (hit - 1);                              // reward matching several words
        s *= 1 + 0.6 * th / qt.length;                          // and words that appear in the heading
        if (x.d.url.split("#")[0] === here) s *= 1.1;           // slight preference for the current page
      }
      return { s: s, d: x.d };
    }).filter(function (r) { return r.s > 0; })
      .sort(function (a, b) { return b.s - a.s; })
      .slice(0, k).map(function (r) { return r.d; });
  }
  function snippet(text, q) {
    var qt = tokens(q), lower = text.toLowerCase(), pos = -1;
    for (var i = 0; i < qt.length && pos < 0; i++) pos = lower.indexOf(qt[i]);
    var start = Math.max(0, pos - 60), s = text.slice(start, start + 220).replace(/\s+/g, " ");
    s = (start > 0 ? "…" : "") + esc(s) + (start + 220 < text.length ? "…" : "");
    qt.filter(function (t) { return t.length > 2; }).forEach(function (t) {
      s = s.replace(new RegExp("(" + t.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + "[a-z]*)", "gi"), "<mark>$1</mark>");
    });
    return s;
  }

  // ---------- UI ----------
  var fab = el("button", { "class": "ask-fab", type: "button", "aria-haspopup": "dialog", "aria-expanded": "false", "aria-controls": "ask-panel" },
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M21 12a8 8 0 0 1-11.6 7.1L4 20l1-4.6A8 8 0 1 1 21 12z"/><path d="M9.5 9.5a2.5 2.5 0 0 1 4.8 1c0 1.7-2.3 2-2.3 3.5"/><path d="M12 17h.01"/></svg><span>Ask</span>');
  var panel = el("div", { "class": "ask-panel", id: "ask-panel", role: "dialog", "aria-modal": "false", "aria-labelledby": "ask-title", hidden: "" });
  panel.innerHTML =
    '<div class="ask-head"><h2 id="ask-title">Ask the hub</h2><button class="ask-x" type="button" aria-label="Close">×</button></div>' +
    '<div class="ask-body">' +
    '<form class="ask-form" autocomplete="off"><input type="search" name="q" placeholder="e.g. why use TMM normalisation?" aria-label="Your question" maxlength="300">' +
    '<button class="ask-btn" type="submit">Search</button></form>' +
    '<p class="ask-hint">Searches all four pages. Then, optionally, a small AI model running in your browser writes an answer from the sections found.</p>' +
    '<div class="ask-chips"></div>' +
    '<div class="ask-out" aria-live="polite"></div>' +
    "</div>";
  var form = panel.querySelector("form"), input = form.q, out = panel.querySelector(".ask-out"), chips = panel.querySelector(".ask-chips");
  EXAMPLES.forEach(function (ex) {
    var c = el("button", { "class": "ask-chip", type: "button" }, esc(ex));
    c.onclick = function () { input.value = ex; run(ex); };
    chips.appendChild(c);
  });

  function open() {
    panel.hidden = false; fab.setAttribute("aria-expanded", "true");
    loadIndex().catch(function () {});
    setTimeout(function () { input.focus(); }, 0);
  }
  function close() { panel.hidden = true; fab.setAttribute("aria-expanded", "false"); fab.focus(); }
  fab.onclick = function () { panel.hidden ? open() : close(); };
  panel.querySelector(".ask-x").onclick = close;
  document.addEventListener("keydown", function (e) { if (e.key === "Escape" && !panel.hidden) close(); });
  form.onsubmit = function (e) { e.preventDefault(); run(input.value); };

  function run(q) {
    q = (q || "").trim();
    if (!q) return;
    lastQuestion = q;
    out.innerHTML = '<p class="ask-note">Searching…</p>';
    loadIndex().then(function () {
      lastHits = search(q, TOP_K);
      render(q);
    }).catch(function () {
      out.innerHTML = '<p class="ask-note ask-err">Could not load the search index. Check your connection and try again.</p>';
    });
  }

  function render(q) {
    chips.hidden = true;
    if (!lastHits.length) {
      out.innerHTML = '<div class="ask-sec"><p class="ask-note">Nothing on the site matches those words. Try other terms (for example a method name such as <i>voom</i>, <i>Harmony</i> or <i>Banksy</i>).</p></div>';
      return;
    }
    var html = aiBox() + '<div class="ask-sec"><h3>Relevant sections</h3><ol class="ask-res">';
    lastHits.forEach(function (d, i) {
      html += '<li><span class="n">[' + (i + 1) + ']</span><a href="' + esc(d.url) + '">' + esc(d.title) + '</a> <span class="pg">· ' + esc(d.page) + "</span>" +
        "<p>" + snippet(d.text, q) + "</p></li>";
    });
    html += "</ol></div>";
    out.innerHTML = html;
    out.querySelectorAll(".ask-res a").forEach(function (a) { a.addEventListener("click", function () { if (window.innerWidth < 700) close(); }); });
    wireAi();
  }

  // ---------- on-device model ----------
  function aiBox() {
    if (!("gpu" in navigator)) {
      return '<div class="ask-ai"><p class="ask-note" style="margin:0">AI answers need a browser with WebGPU (recent Chrome or Edge on a laptop or desktop; Safari 26+). ' +
        "This browser doesn't have it, so use the sections below.</p></div>";
    }
    var saved = store("ask-model") || "1.5b", opts = MODELS.map(function (m) {
      return '<option value="' + m.key + '"' + (m.key === saved ? " selected" : "") + ">" + esc(m.label) + "</option>";
    }).join("");
    return '<div class="ask-ai"><div class="ask-ai-row">' +
      '<button class="ask-btn ask-gen" type="button">Generate answer</button>' +
      '<select class="ask-model" aria-label="Model">' + opts + "</select></div>" +
      '<p class="ask-note">Runs on your device; your question is not sent anywhere. The first use downloads the model once, then it is cached. ' +
      "Small models make mistakes: check the numbered sections below.</p>" +
      '<div class="ask-prog" hidden><div class="ask-bar"><i></i></div><p class="ask-note ask-prog-t"></p></div>' +
      '<div class="ask-answer"></div></div>';
  }
  function wireAi() {
    var gen = out.querySelector(".ask-gen");
    if (!gen) return;
    var sel = out.querySelector(".ask-model");
    sel.onchange = function () { store("ask-model", sel.value); };
    gen.onclick = function () {
      if (busy && engine) { try { engine.interruptGenerate(); } catch (e) {} return; }
      answer(sel.value, gen);
    };
  }

  function pickModel(key, hasF16) {
    var m = MODELS.filter(function (x) { return x.key === key; })[0] || MODELS[0];
    return hasF16 ? m.f16 : m.f32;
  }
  function getEngine(key, progress) {
    return navigator.gpu.requestAdapter().then(function (adapter) {
      if (!adapter) throw new Error("No WebGPU adapter: your GPU or browser settings block WebGPU.");
      var id = pickModel(key, adapter.features && adapter.features.has("shader-f16"));
      if (engine && engineKey === id) return engine;
      return import(WEBLLM_URL).then(function (webllm) {
        var prev = engine ? engine.unload() : Promise.resolve();
        return prev.then(function () {
          engine = null;
          return webllm.CreateMLCEngine(id, { initProgressCallback: progress });
        }).then(function (eng) { engine = eng; engineKey = id; return eng; });
      });
    });
  }

  function buildMessages(q, hits) {
    var ctx = hits.slice(0, CONTEXT_K).map(function (d, i) {
      var t = d.text.length > CONTEXT_CHARS ? d.text.slice(0, CONTEXT_CHARS) + " …" : d.text;
      return "[" + (i + 1) + "] " + d.page + " — " + d.title + "\n" + t;
    }).join("\n\n");
    return [
      { role: "system", content:
        "You are the study assistant for OmicLearningHub, a bioinformatics tutorial site (bulk RNA-seq, single-cell RNA-seq, Xenium spatial). " +
        "Answer the question using ONLY the numbered excerpts. Cite the excerpts you used like [1] or [2]. " +
        "If the excerpts do not contain the answer, say that the site does not cover it and do not guess. " +
        "Be concise: at most about 150 words, plain language, and short code only if the excerpts contain it." },
      { role: "user", content: "Excerpts:\n\n" + ctx + "\n\nQuestion: " + q }
    ];
  }

  function md(text, hits) {
    var parts = text.split(/```[a-zA-Z]*\n?/), html = "";
    parts.forEach(function (p, i) {
      if (i % 2) { html += "<pre><code>" + esc(p.replace(/\n$/, "")) + "</code></pre>"; return; }
      p.split(/\n{2,}/).forEach(function (para) {
        if (!para.trim()) return;
        var h = esc(para.trim())
          .replace(/`([^`]+)`/g, "<code>$1</code>")
          .replace(/\*\*([^*]+)\*\*/g, "<b>$1</b>")
          .replace(/\n/g, "<br>")
          .replace(/\[(\d+)\]/g, function (m, n) {
            var d = hits[+n - 1];
            return d ? '<a class="cite" href="' + esc(d.url) + '" title="' + esc(d.title) + '">[' + n + "]</a>" : m;
          });
        html += "<p>" + h + "</p>";
      });
    });
    return html;
  }

  function answer(key, gen) {
    var box = out.querySelector(".ask-answer"), prog = out.querySelector(".ask-prog"),
      bar = prog.querySelector("i"), pt = prog.querySelector(".ask-prog-t"),
      sel = out.querySelector(".ask-model"), q = lastQuestion, hits = lastHits.slice();
    busy = true; gen.textContent = "Stop"; sel.disabled = true; box.innerHTML = "";
    prog.hidden = false; pt.textContent = "Loading the model…"; bar.style.width = "2%";

    getEngine(key, function (r) {
      bar.style.width = Math.max(2, Math.round((r.progress || 0) * 100)) + "%";
      pt.textContent = r.text || "Loading the model…";
    }).then(function (eng) {
      prog.hidden = true;
      box.innerHTML = '<p class="ask-note">Thinking…</p>';
      return eng.chat.completions.create({
        messages: buildMessages(q, hits), stream: true, temperature: 0.2, top_p: 0.9, max_tokens: 400
      }).then(function (stream) {
        var text = "", it = stream[Symbol.asyncIterator]();
        function step() {
          return it.next().then(function (r) {
            if (r.done) return text;
            var delta = r.value.choices && r.value.choices[0] && r.value.choices[0].delta;
            if (delta && delta.content) { text += delta.content; box.innerHTML = md(text, hits); }
            return step();
          });
        }
        return step();
      });
    }).then(function (text) {
      if (!text || !text.trim()) box.innerHTML = '<p class="ask-note">No answer was produced. Use the sections below.</p>';
      else box.insertAdjacentHTML("beforeend", '<p class="ask-note">AI-generated from the sections below. Check them before relying on it.</p>');
    }).catch(function (err) {
      prog.hidden = true;
      var msg = String((err && err.message) || err);
      if (/memory|OOM|allocate|device lost|buffer/i.test(msg)) msg += " — try the smaller 0.5B model, or close other tabs.";
      box.innerHTML = '<p class="ask-note ask-err">The AI model could not run here: ' + esc(msg) + " The sections below still answer most questions.</p>";
    }).then(function () {
      busy = false; gen.textContent = "Generate answer"; sel.disabled = false;
    });
  }

  function mount() { document.body.appendChild(fab); document.body.appendChild(panel); }
  if (document.body) mount(); else document.addEventListener("DOMContentLoaded", mount);

  // test hook (used by the repository's checks; harmless in production)
  window.__omicAsk = { loadIndex: loadIndex, search: function (q, k) { return search(q, k || TOP_K); }, buildMessages: buildMessages, md: md, run: run, open: open };
})();
