#!/usr/bin/env python3
"""Build assets/ask/index.json: the site's pages split into small, linkable
sections that the Ask widget searches. Standard library only, so it also runs
in the GitHub Pages workflow before every deploy.

Usage: python3 tools/build_ask_index.py
"""
import json
import re
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAGES = [("index.html", "Home"), ("rnaseq.html", "Module 1 · Bulk RNA-seq"),
         ("scrnaseq.html", "Module 2 · Single-cell RNA-seq"), ("xenium.html", "Module 3 · Spatial (Xenium)")]
SKIP = {"script", "style", "nav", "header", "footer", "svg", "button", "canvas", "figcaption", "label"}
HEADINGS = {"h1", "h2", "h3", "h4"}
VOID = {"br", "img", "input", "meta", "link", "hr", "source", "wbr"}
MAX_CHARS = 1600


class Splitter(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []            # (tag, had_id)
        self.skip = 0
        self.anchor = "top"
        self.heading = ""
        self.in_heading = False
        self.head_buf = []
        self.buf = []
        self.chunks = []

    def flush(self):
        text = re.sub(r"[ \t]+", " ", "".join(self.buf))
        text = re.sub(r"\s*\n\s*", "\n", text).strip()
        if len(text) > 40:
            self.chunks.append({"anchor": self.anchor, "heading": self.heading, "text": text})
        self.buf = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in VOID:
            if tag == "br":
                self.buf.append("\n")
            return
        self.stack.append(tag)
        if tag in SKIP or self.skip:
            self.skip += 1
            return
        if a.get("id") and tag in ("section", "article", "div", "details"):
            self.flush()
            self.anchor = a["id"]
        if tag in HEADINGS:
            self.flush()
            self.in_heading = True
            self.head_buf = []
        if tag in ("p", "li", "tr", "pre", "div", "dd", "dt", "summary"):
            self.buf.append("\n")
        if tag in ("td", "th"):
            self.buf.append(" | ")

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        while self.stack:
            t = self.stack.pop()
            if self.skip:
                self.skip -= 1
            if t == tag:
                break
        if tag == "b" and not self.skip:
            self.buf.append(" ")                     # callout labels: "Decide" + text
        if tag in HEADINGS and self.in_heading:
            self.in_heading = False
            self.heading = re.sub(r"\s+", " ", "".join(self.head_buf)).strip()

    def handle_data(self, data):
        if self.skip:
            return
        if self.in_heading:
            self.head_buf.append(data)
        else:
            self.buf.append(data)


def split_long(text, limit=MAX_CHARS):
    if len(text) <= limit:
        return [text]
    parts, cur = [], ""
    for para in text.split("\n"):
        if len(cur) + len(para) + 1 > limit and cur:
            parts.append(cur.strip())
            cur = ""
        cur += para + "\n"
    if cur.strip():
        parts.append(cur.strip())
    return parts


def main():
    out = []
    for fname, page_title in PAGES:
        html = (ROOT / fname).read_text(encoding="utf-8")
        body = html[html.find("<main"):] if "<main" in html else html
        sp = Splitter()
        sp.feed(body)
        sp.flush()
        for c in sp.chunks:
            for i, part in enumerate(split_long(c["text"])):
                out.append({
                    "id": len(out),
                    "page": page_title,
                    "title": c["heading"] or page_title,
                    "url": f"{fname}#{c['anchor']}" if c["anchor"] != "top" else fname,
                    "text": part,
                })
    dest = ROOT / "assets" / "ask" / "index.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"{len(out)} sections from {len(PAGES)} pages -> {dest.relative_to(ROOT)} "
          f"({dest.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
