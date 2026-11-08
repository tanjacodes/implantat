import json
from run_optim import evaluate
combos = [("4-hole, bar 2.0", dict(holes_p=[-5.5, 5.5], wb=2.0)),
          ("4-hole, bar 2.4, t 0.8", dict(holes_p=[-5.5, 5.5], wb=2.4, tp=0.8)),
          ("4-hole, bar 2.4", dict(holes_p=[-5.5, 5.5], wb=2.4)),
          ("8-hole outer struts, bar 2.2, t 0.8", dict(struts="outer", wb=2.2, tp=0.8))]
out = [evaluate(n, p) for n, p in combos]
json.dump(out, open("results/optim_combos.json", "w"), indent=1)
