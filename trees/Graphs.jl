#Julia has built-in support for "symmetric views" of graphs

using CairoMakie
using NetworkLayout
using GraphMakie, Graphs

g = watts_strogatz(1000, 5, 0.03; seed=5)
layout = Spectral(dim=2)
f, ax, p = graphplot(g, layout=layout, node_size=0.0, edge_width=1.0)