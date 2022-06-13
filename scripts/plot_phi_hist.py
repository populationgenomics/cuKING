#!/usr/bin/env python3

import hail as hl
from bokeh.io import show
from bokeh.models import Span

hl.init()

mt = hl.read_matrix_table('king_1kg.mt')
# Filter to avoid huge outlier bin.
mt = mt.filter_entries(mt.phi >= 0.)

p = hl.plot.histogram(mt.phi, bins=100)
p.legend.visible = False

# Mark relatedness degree ranges (https://www.kingrelatedness.com/manual.shtml)
vlines = []
for x in [0.0884, 0.0442, 0.177, 0.354]:
    vlines.append(Span(location=x, dimension='height', line_color='red', line_width=1))
p.renderers.extend(vlines)

show(p)

