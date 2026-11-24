
"""Render the submission workflow figure as a vector PDF."""
from pathlib import Path

from reportlab.lib.colors import HexColor
from reportlab.pdfgen.canvas import Canvas

OUTPUT = Path(__file__).resolve().parents[1] / "fig1_pipeline.pdf"
canvas = Canvas(str(OUTPUT), pagesize=(462, 350))
canvas.setTitle("TEAnnotate workflow")
INK = HexColor("#152C2D")
MUTED = HexColor("#305351")
BORDER = HexColor("#ACBFBA")
FILL = HexColor("#F2F7F4")
SIDE = HexColor("#FAF3E4")
ACCENT = HexColor("#15584D")


def text(x, y, lines, size=9.1, heading=False):
    canvas.setFillColor(INK if heading else MUTED)
    canvas.setFont("Helvetica-Bold" if heading else "Helvetica", size)
    for i, line in enumerate(lines):
        canvas.drawString(x, y - i * 12, line)


def box(x, y, width, height, fill=FILL):
    canvas.setFillColor(fill)
    canvas.setStrokeColor(BORDER)
    canvas.setLineWidth(0.65)
    canvas.roundRect(x, y, width, height, 5, fill=1, stroke=1)


def arrow(x1, y1, x2, y2):
    canvas.setStrokeColor(ACCENT)
    canvas.setLineWidth(0.9)
    canvas.line(x1, y1, x2, y2)
    if x1 == x2:
        canvas.line(x2 - 3, y2 + 4, x2, y2)
        canvas.line(x2 + 3, y2 + 4, x2, y2)
    else:
        canvas.line(x2 - 4, y2 - 3, x2, y2)
        canvas.line(x2 - 4, y2 + 3, x2, y2)


box(2, 305, 336, 43)
text(12, 333, ["1   Validate inputs"], 11, True)
text(12, 318, ["Genome assembly, lineage and optional user library"])
arrow(170, 305, 170, 293)

box(2, 228, 336, 65)
text(12, 278, ["2   Discover transposable elements"], 11, True)
text(12, 263, ["Animals and fungi     Earl Grey + HiTE",
               "Plants                         EDTA + HELIANO + HiTE",
               "Outputs                      Pooled consensi and tool provenance"])
box(353, 234, 107, 55, SIDE)
text(361, 274, ["Starfish report"], 9.1, True)
text(361, 260, ["Fungal captains", "from the assembly"], 8.8)
arrow(338, 261, 353, 261)
arrow(170, 228, 170, 216)

box(2, 124, 336, 92)
text(12, 201, ["3   Construct the library"], 11, True)
text(12, 185, ["Class-aware deduplication"], 9.1, True)
text(130, 185, ["CD-HIT-EST or MMseqs2-based search"])
text(12, 166, ["Evidence gate"], 9.1, True)
text(130, 166, ["TE domains, tool agreement or copy number"])
text(12, 147, ["Optional TEtrimmer"], 9.1, True)
text(130, 147, ["Restore supported inputs and refine the library"])
arrow(170, 124, 170, 112)

box(2, 2, 336, 110)
text(12, 97, ["4   Classify and annotate"], 11, True)
text(12, 80, ["Domain-based relabelling"], 9.1, True)
text(12, 62, ["Two-step vote for remaining unknowns"], 9.1, True)
text(12, 50, ["Order then superfamily, restricted to DNA transposons and Helitrons by default"], 8.6)
text(12, 31, ["RepeatMasker"], 9.1, True)
text(12, 19, ["Masked assembly, repeat coordinates, libraries and summary plots"])
box(353, 12, 107, 55, SIDE)
text(361, 53, ["RIP report"], 9.1, True)
text(361, 39, ["Fungal RIP-like", "signatures after masking"], 8.8)
arrow(338, 28, 353, 28)

canvas.showPage()
canvas.save()

