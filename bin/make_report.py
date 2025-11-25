#!/usr/bin/env python3
import os
import glob
import html
import argparse

def build_report(sample: str, jaffal_dir: str, out_html: str):
    final_dir = os.path.join(jaffal_dir, 'final')
    files = []
    if os.path.isdir(final_dir):
        for ext in ('*.txt','*.tsv','*.csv'):
            files.extend(glob.glob(os.path.join(final_dir, ext)))

    rows = []
    for fp in files:
        name = os.path.basename(fp)
        try:
            with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
                lines = []
                for i, line in enumerate(f):
                    if i >= 10:
                        break
                    lines.append(line.rstrip('\n'))
            preview = html.escape('\n'.join(lines))
        except Exception as e:
            preview = html.escape(str(e))
        rows.append((name, preview))

    parts = []
    parts.append(f"<html><head><meta charset='utf-8'><title>{html.escape(sample)} JAFFAL Report</title></head><body>")
    parts.append(f"<h1>Sample: {html.escape(sample)}</h1>")
    parts.append(f"<p>JAFFAL output directory: {html.escape(jaffal_dir)}</p>")
    if rows:
        parts.append("<h2>Final Outputs</h2>")
        parts.append("<ul>")
        for name, _ in rows:
            rel = os.path.join('jaffal', 'final', name)
            parts.append(f"<li><a href='{html.escape(rel)}'>{html.escape(name)}</a></li>")
        parts.append("</ul>")
        parts.append("<h2>Previews (first 10 lines)</h2>")
        for name, preview in rows:
            parts.append(f"<h3>{html.escape(name)}</h3>")
            parts.append("<pre style='background:#f7f7f7;border:1px solid #ddd;padding:10px;'>" + preview + "</pre>")
    else:
        parts.append("<p>No files found under JAFFAL final/ directory.</p>")
    parts.append("</body></html>")

    with open(out_html, 'w', encoding='utf-8') as out:
        out.write('\n'.join(parts))


def main():
    p = argparse.ArgumentParser(description='Generate JAFFAL HTML report')
    p.add_argument('--sample', required=True)
    p.add_argument('--jaffal_dir', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    build_report(args.sample, args.jaffal_dir, args.out)

if __name__ == '__main__':
    main()

