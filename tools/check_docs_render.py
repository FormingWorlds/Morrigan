"""Compare the built documentation against its sources.

Markdown fails quietly: a continuation paragraph that is not indented ends the
list it belongs to, and the items after it are printed as literal "6." in a
paragraph. Maths fails quietly too, leaving raw TeX in the text. Both look fine
in the source and wrong on the page, so this compares the two.

This reads the built HTML; it does not execute the page. It therefore catches a
maths renderer that is configured wrongly, but not one whose script fails to
load at run time, and it cannot see anything a stylesheet or a script decides.
Read a page in a browser as well before trusting a clean run.

Usage: python tools/check_docs_render.py [site_dir]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

DOCS = Path('docs')

# A list marker that reached the page as text rather than as a list item.
LEAKED_MARKER = re.compile(r'(?<![\w.>])(\d{1,2})\.\s+(?=[A-Z*_`])')

# Raw TeX that never became maths.
RAW_TEX = re.compile(
    r'\\(?:frac|mathrm|sqrt|left|right|alpha|beta|Delta|sum)\b|\$\$?[^$\n]{2,}\$\$?'
)


def strip_tags(html: str) -> str:
    html = re.sub(r'<(script|style|pre|code)\b.*?</\1>', ' ', html, flags=re.S)
    # Maths reaches the page either already typeset, or as the \( \) that
    # arithmatex writes for the renderer to pick up in the browser. Both count as
    # maths; only TeX that sits outside them is a failure.
    html = re.sub(r'<mjx-container.*?</mjx-container>', ' MATH ', html, flags=re.S)
    html = re.sub(r'<(span|div) class="arithmatex">.*?</\1>', ' MATH ', html, flags=re.S)
    html = re.sub(r'<[^>]+>', ' ', html)
    return re.sub(r'\s+', ' ', html)


def source_list_counts(md: Path) -> list[int]:
    """Length of every top-level ordered list in a markdown source."""
    counts, run = [], 0
    for line in md.read_text().splitlines():
        if re.match(r'^\d{1,2}\.\s', line):
            run += 1
        elif line.strip() and not line.startswith(('    ', '\t')):
            if run:
                counts.append(run)
            run = 0
    if run:
        counts.append(run)
    return counts


def rendered_list_counts(html: str) -> list[int]:
    body = html[html.index('<article') :] if '<article' in html else html
    body = re.sub(r'<(script|style)\b.*?</\1>', '', body, flags=re.S)
    # the footnote definitions are an ordered list of their own, not content
    body = re.sub(r'<div class="footnote">.*?</div>', '', body, flags=re.S)
    return [
        len(re.findall(r'<li\b', m.group(1)))
        for m in re.finditer(r'<ol\b[^>]*>(.*?)</ol>', body, re.S)
    ]


# A figure pulled from another repository at a branch ref stops resolving when
# that branch is deleted. Reported so a pin taken while a pull request is open is
# not left behind once it merges.
PINNED_ASSET = re.compile(
    r'(?:cdn\.jsdelivr\.net/gh|raw\.githubusercontent\.com)/'
    r'FormingWorlds/(\w+)[@/](?!main\b)([\w./-]+?)/docs/assets/([\w.-]+)'
)


def check_math_renderer(html: str) -> str | None:
    """The maths renderer has to be able to style what it produces.

    KaTeX hides its MathML-and-TeX layer with its own stylesheet. Load its script
    without that stylesheet and the layer becomes visible text, printing the TeX
    source next to every equation.
    """
    if 'class="arithmatex"' not in html:
        return None  # no maths on this page, so no renderer to get wrong
    if 'katex' in html and 'katex.min.css' not in html:
        return (
            'KaTeX is loaded without katex.min.css, so its own output is left '
            'unstyled and renders twice'
        )
    if not re.search(r'mathjax|katex', html, re.I):
        return 'the page has maths but loads no renderer'
    return None


def main(site='site'):
    site = Path(site)
    problems = 0
    for page in sorted(site.rglob('*.html')):
        if page.name == '404.html':
            continue
        html = page.read_text()
        body = html[html.index('<article') :] if '<article' in html else html
        # a heading may legitimately begin with a number, so headings are not
        # searched for leaked list markers
        text = strip_tags(re.sub(r'<h[1-6]\b.*?</h[1-6]>', ' ', body, flags=re.S))

        for m in LEAKED_MARKER.finditer(text):
            # a numbered marker in running prose is a list that broke
            ctx = text[max(0, m.start() - 60) : m.end() + 50]
            print(
                f'{page.relative_to(site)}: list marker "{m.group(1)}." in prose\n    ...{ctx}...'
            )
            problems += 1

        for m in RAW_TEX.finditer(text):
            ctx = text[max(0, m.start() - 50) : m.end() + 30]
            print(f'{page.relative_to(site)}: raw TeX on the page\n    ...{ctx}...')
            problems += 1

        problem = check_math_renderer(html)
        if problem:
            print(f'{page.relative_to(site)}: {problem}')
            problems += 1

        md = DOCS / page.relative_to(site).with_suffix('.md')
        if md.exists():
            src, out = source_list_counts(md), rendered_list_counts(html)
            if src and sum(src) != sum(out):
                print(
                    f'{page.relative_to(site)}: ordered-list items {sum(src)} in source, '
                    f'{sum(out)} on the page (source {src}, page {out})'
                )
                problems += 1

    pins = {}
    for page in sorted(site.rglob('*.html')):
        for repo, ref, asset in PINNED_ASSET.findall(page.read_text()):
            pins.setdefault((repo, ref), set()).add(asset)
    for (repo, ref), assets in pins.items():
        print(
            f'note: {repo} assets pinned to "{ref}" rather than main: '
            f'{", ".join(sorted(assets))}'
        )

    print(f'\n{problems} problem(s)')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
