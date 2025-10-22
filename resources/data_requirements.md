Identity Time Series Data Requirements
=====================================

Required columns
----------------
- `ID`: Person identifier (integer or string)
- `week`: Context identifier (e.g., 1–3 for three weekly bursts)
- `sday`: Study day within context (1–7 if daily data)
- `GrL`: Leader identity measure (numeric)
- `GrF`: Follower identity measure (numeric)

Events (choose one option)
-------------------------
1) Provide component items and the script will compute an aggregate:
   - `e2`, `e3`, `e6` (numeric)

2) Provide a precomputed event strength composite:
   - `event_agg` (numeric)

Example CSV (header + first two rows)
-------------------------------------
```
ID,week,sday,GrL,GrF,e2,e3,e6
101,1,1,5.2,4.8,0.1,0.2,0.0
101,1,2,5.0,4.7,0.4,0.1,0.3
```

Notes
-----
- `sday` should be evenly spaced within each `week` for stable derivative estimates.
- If your design uses a single context, you may still name it `week = 1`.
- If you already have `event_agg`, you may omit `e2,e3,e6`.
- Scales: ensure `GrL`/`GrF` are comparable across people/contexts (e.g., same Likert scale).

