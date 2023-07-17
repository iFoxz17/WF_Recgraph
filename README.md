# WF_RecGraph
WF_Recgraph contains an optimizazion of [RecGraph](https://github.com/AlgoLab/RecGraph) for **variation graph** to **sequence alignment**, implementing an a **multithreading** approach based on the ***wavefront*** **algorithm**.

## Optimizated complexity
- $T(n, m, \overline{d}, p) = O(\min{(n, m)} \cdot \overline{d} \cdot p)$
- $M(n, m, \overline{d}, p) = O(\max{(n, m)} \cdot \overline{d} \cdot p)$

## Future improvements
- Extend **weighted edit distance** to **gap pairwise alignment**;
- Reduce memory complexity of ***wf_vec*** implementation to $O(\overline{d}^2 \cdot p)$;
- Implements ***wavefront heuristics***.
