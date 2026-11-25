# Synthetic instances for IbexOpt (Section 6.3 of the paper, Table 4)

| Folder | Instances | Generator |
|---|---|---|
| `synq/` | 78: n ∈ {3,5,8,10}, ρ ∈ {0.2,0.5,0.9}, diagonal pos/posneg (+2 extra configurations), 3 seeds | `generators/generar_sinteticas.py` |
| `synt/` | 24: n ∈ {3,5,8,10}, ρ ∈ {0.5,0.9}, 3 seeds | `generators/generar_sinteticas_funciones.py --familia synt` |
| `synp/` | 18: n ∈ {3,5,8}, 6 seeds | `generators/gen_synp.py <out_dir>` |

Each folder of the first two has a `manifest.csv` with the parameters of every instance. The COCONUT set of Table 4
is the optimization benchmark distributed with Ibex (see `ibex/README.md`).
