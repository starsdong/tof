#!/usr/bin/env python3
"""Tile PNGs into a contact sheet. usage: contact_sheet.py out.png img1.png img2.png ... [--cols N] [--w W]"""
import sys
from PIL import Image, ImageDraw

args = sys.argv[1:]
cols, w = 4, 400
if "--cols" in args:
    i = args.index("--cols"); cols = int(args[i + 1]); del args[i:i + 2]
if "--w" in args:
    i = args.index("--w"); w = int(args[i + 1]); del args[i:i + 2]
out, files = args[0], args[1:]
ims = []
for f in files:
    im = Image.open(f).convert("RGB")
    h = int(im.height * w / im.width)
    ims.append((f, im.resize((w, h))))
hmax = max(im.height for _, im in ims) + 16
rows = (len(ims) + cols - 1) // cols
sheet = Image.new("RGB", (cols * w, rows * hmax), "white")
d = ImageDraw.Draw(sheet)
for k, (f, im) in enumerate(ims):
    x, y = (k % cols) * w, (k // cols) * hmax
    sheet.paste(im, (x, y + 16))
    d.text((x + 4, y + 2), f.split("/")[-1], fill="black")
sheet.save(out)
print(out, len(ims))
