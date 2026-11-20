# Geometric states

`examples/graph_states/build.py` compacts a CFA into a ToCUMG (CDBG) and a coloured graph tensor whose node features include mean longitude and latitude.

```bash
python examples/graph_states/build.py --cfa CFA_DIR --out OUT_DIR --node-namespace NODE_NS --edge-namespaces NS1,NS2 [--draw-edge-namespace NS] [--epochs 80] [--seed 0]
```

It writes `cdbg/`, `cgt/`, `states.pdf`, and `gcn/` under `OUT_DIR`.

The PDF has three pages, in order: CFA, ToCUMG, CGT. The script always draws longitude and latitude; it has no layout option. Fruchterman–Reingold coordinates (`Layout x`, `Layout y`) are available only from Python, as `layout="fr"` in `metametro.viz.states.plot_geometric_states`. A graph with a `coverage` column uses that column and one shared YlGnBu scale. Otherwise edges are coloured by `--draw-edge-namespace` (Set1) and nodes by degree.

The GCN reads only longitude and latitude. It predicts one node class (softmax) and an independent 0/1 probability for every edge colour in `--edge-namespaces` (sigmoid, threshold 0.5). Metrics in `gcn/metrics.json` are computed on the seed-0 split.
