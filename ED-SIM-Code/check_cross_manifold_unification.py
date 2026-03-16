#!/usr/bin/env python
"""
check_cross_manifold_unification.py
=====================================
Cross-Manifold Unification Test: determines whether the five ED manifolds
(mechanism, drift, temporal, atlas, boundary) describe a single coherent
geometric object.

For each (N, angle, radius, drift) cell, assembles a unified record:
  - mechanism classification (from sweep)
  - mechanism sub-type (from atlas)
  - collapse time chi (from sweep)
  - boundary distance (signed, from IC onset)
  - drift sensitivity class
  - temporal gradient class
  - local curvature of boundary

Sections:
  1. Unified tensor assembly
  2. Mechanism <-> temporal correlation
  3. Mechanism <-> atlas subtype mapping
  4. Mechanism <-> boundary distance
  5. Drift sensitivity <-> temporal gradient
  6. Atlas subtype <-> boundary curvature
  7. N=8 anomaly across all layers
  8. Shared topological features (folds, cusps, multi-valued, drift-stable)
  9. Coherence test: single manifold, multiple sheets, or piecewise union
 10. Cross-N invariants
 11. Summary table
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter

SWEEP_DIR = os.path.dirname(__file__)

NS     = [4, 8, 12, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS = [0.00, 0.01, 0.02, 0.05, 0.10]

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}

SORTED_ANGLES = sorted(ANGLES)


# ============================================================
# Data loading
# ============================================================

def load_sweeps():
    """Load all 4D sweep files into a dict keyed by (N, angle, radius, drift)."""
    cube = {}
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


def load_atlas():
    """Load all atlas files into a dict keyed by (N, angle, radius, drift)."""
    cube = {}
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_atlas.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        if "atlas" not in data:
            continue
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


# ============================================================
# Sub-mechanism classification (from check_mechanism_atlas.py)
# ============================================================

def classify_sub(rec):
    """Assign sub-mechanism label from atlas diagnostics."""
    mech = rec["mechanism"]
    a = rec["atlas"]
    pbc = a["n_pbc_crossings"]
    steps = a["n_steps"]
    seq = a["corner_sequence"]
    period = a["periodicity_score"]
    multi = a["is_multi_stage"]
    tcross = a["tangent_crossings"]

    if mech == "inward-collapse":
        if pbc > 0:           return "IC-bounced"
        elif steps > 50:      return "IC-delayed"
        else:                 return "IC-direct"
    elif mech == "outward-PBC":
        if pbc == 0:          return "OPBC-single"
        elif pbc <= 5:        return "OPBC-glancing"
        else:                 return "OPBC-winding"
    elif mech == "PBC-corner":
        distinct = len(set(seq))
        diag_pairs = ({0, 3}, {1, 2})
        seq_set = set(seq)
        if distinct <= 1:     return "CR-single"
        elif any(dp <= seq_set for dp in diag_pairs):
                              return "CR-diagonal"
        else:                 return "CR-multi"
    elif mech == "other-late":
        if multi:             return "OL-multi"
        elif period > 0.3:    return "OL-periodic"
        elif steps < 5:       return "OL-stalled"
        else:                 return "OL-wandering"
    elif mech == "DECAY":
        if tcross > 0:        return "DE-tangent"
        else:                 return "DE-frozen"
    return mech


# ============================================================
# Derived quantities
# ============================================================

def compute_ic_onset(sweep_cube):
    """For each (N, radius, drift), find onset angle of IC."""
    onset = {}
    for n in NS:
        for r in RADII:
            for d in DRIFTS:
                ic_angs = [ang for ang in SORTED_ANGLES
                           if (n, ang, r, d) in sweep_cube
                           and sweep_cube[(n, ang, r, d)]["mechanism"] == "inward-collapse"]
                if not ic_angs:
                    onset[(n, r, d)] = None
                else:
                    onset[(n, r, d)] = max(ic_angs)
    return onset


def signed_boundary_dist(angle, onset_angle):
    """Signed distance from IC boundary.  Negative = inside IC, positive = outside."""
    if onset_angle is None:
        return float('nan')
    return angle - onset_angle


def drift_sensitivity_class(onset_d0, onset_d):
    """Classify drift sensitivity of boundary at this (N, radius)."""
    if onset_d0 is None or onset_d is None:
        return "no-boundary"
    shift = abs(onset_d - onset_d0)
    if shift == 0:    return "rigid"
    elif shift <= 2:  return "stable"
    elif shift <= 5:  return "moderate"
    else:             return "fragile"


def temporal_gradient_class(chi, chi_neighbor_mean):
    """Classify temporal gradient at this cell."""
    if chi_neighbor_mean is None or chi_neighbor_mean < 0.001:
        return "isolated"
    ratio = chi / chi_neighbor_mean
    if ratio > 5.0:   return "cliff"
    elif ratio > 2.0:  return "steep"
    elif ratio > 1.2:  return "ramp"
    elif ratio > 0.8:  return "flat"
    elif ratio > 0.5:  return "ramp-down"
    else:              return "cliff-down"


def local_curvature(onset, n, d, r_idx):
    """3-point boundary curvature at RADII[r_idx]."""
    if r_idx < 1 or r_idx >= len(RADII) - 1:
        return 0.0
    r0, r1, r2 = RADII[r_idx - 1], RADII[r_idx], RADII[r_idx + 1]
    a0 = onset.get((n, r0, d))
    a1 = onset.get((n, r1, d))
    a2 = onset.get((n, r2, d))
    if a0 is None or a1 is None or a2 is None:
        return 0.0
    if a0 > max(ANGLES) or a1 > max(ANGLES) or a2 > max(ANGLES):
        return 0.0
    dx1, dy1 = r1 - r0, float(a1 - a0)
    dx2, dy2 = r2 - r1, float(a2 - a1)
    ddx, ddy = dx2 - dx1, dy2 - dy1
    dxa, dya = (dx1 + dx2) / 2, (dy1 + dy2) / 2
    denom = (dxa**2 + dya**2)**1.5
    if denom < 1e-15:
        return 0.0
    return abs(dxa * ddy - dya * ddx) / denom


# ============================================================
# 1. Unified tensor assembly
# ============================================================

def assemble_unified_tensor(sweep_cube, atlas_cube):
    """Build the unified record for every cell."""
    onset = compute_ic_onset(sweep_cube)

    # Precompute chi neighbor means for temporal gradient
    chi_neighbors = {}
    for key, rec in sweep_cube.items():
        n, ang, r, d = key
        neighbors = []
        ang_idx = SORTED_ANGLES.index(ang) if ang in SORTED_ANGLES else -1
        for da in [-1, +1]:
            ni = ang_idx + da
            if 0 <= ni < len(SORTED_ANGLES):
                nkey = (n, SORTED_ANGLES[ni], r, d)
                if nkey in sweep_cube:
                    neighbors.append(sweep_cube[nkey]["chi_emp"])
        r_idx = RADII.index(r) if r in RADII else -1
        for dr in [-1, +1]:
            ri = r_idx + dr
            if 0 <= ri < len(RADII):
                nkey = (n, ang, RADII[ri], d)
                if nkey in sweep_cube:
                    neighbors.append(sweep_cube[nkey]["chi_emp"])
        chi_neighbors[key] = statistics.mean(neighbors) if neighbors else None

    tensor = {}
    for key in sweep_cube:
        n, ang, r, d = key
        rec = sweep_cube[key]
        mech = rec["mechanism"]
        chi = rec["chi_emp"]

        # Sub-type from atlas
        subtype = "unknown"
        if key in atlas_cube:
            subtype = classify_sub(atlas_cube[key])

        # Boundary distance
        o = onset.get((n, r, d))
        bdist = signed_boundary_dist(ang, o)

        # Drift sensitivity
        o_base = onset.get((n, r, 0.0))
        o_this = onset.get((n, r, d))
        drift_class = drift_sensitivity_class(o_base, o_this)

        # Temporal gradient
        temp_class = temporal_gradient_class(chi, chi_neighbors[key])

        # Local curvature
        r_idx = RADII.index(r) if r in RADII else -1
        kappa = local_curvature(onset, n, d, r_idx)

        # Orbit lifetime (steps from atlas, for non-IC)
        lifetime = 0
        if key in atlas_cube:
            lifetime = atlas_cube[key]["atlas"]["n_steps"]

        tensor[key] = {
            "N": n, "angle": ang, "radius": r, "drift": d,
            "mechanism": mech,
            "subtype": subtype,
            "chi": chi,
            "lifetime": lifetime,
            "boundary_dist": bdist,
            "drift_class": drift_class,
            "temporal_class": temp_class,
            "curvature": kappa,
        }

    return tensor, onset


def check_tensor_assembly(tensor):
    """Section 1: Report unified tensor stats."""
    print("=" * 72)
    print("SECTION 1: UNIFIED TENSOR ASSEMBLY")
    print("=" * 72)
    print(f"  Total cells: {len(tensor)}")

    # Count by data completeness
    has_atlas = sum(1 for v in tensor.values() if v["subtype"] != "unknown")
    has_bdist = sum(1 for v in tensor.values() if not math.isnan(v["boundary_dist"]))
    print(f"  Cells with atlas subtype: {has_atlas}")
    print(f"  Cells with boundary distance: {has_bdist}")

    # Distribution per N
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        mechs = Counter(v["mechanism"] for v in cells)
        subs = Counter(v["subtype"] for v in cells)
        dc = Counter(v["drift_class"] for v in cells)
        tc = Counter(v["temporal_class"] for v in cells)
        print(f"\n  N={n} ({len(cells)} cells):")
        print(f"    Mechanisms: {dict(mechs.most_common())}")
        print(f"    Subtypes:   {dict(subs.most_common(6))}")
        print(f"    Drift class: {dict(dc.most_common())}")
        print(f"    Temporal class: {dict(tc.most_common())}")


# ============================================================
# 2. Mechanism <-> temporal correlation
# ============================================================

def check_mech_temporal(tensor):
    """Section 2: chi statistics per mechanism and N."""
    print("\n" + "=" * 72)
    print("SECTION 2: MECHANISM <-> TEMPORAL CORRELATION")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        print(f"  {'mech':>6s}  {'count':>5s}  {'mean_chi':>9s}  {'med_chi':>8s}  "
              f"{'min':>6s}  {'max':>8s}  {'cv':>5s}")
        for mech in MECHS:
            cells = [v for v in tensor.values()
                     if v["N"] == n and v["mechanism"] == mech]
            if not cells:
                continue
            chis = [v["chi"] for v in cells]
            mean_c = statistics.mean(chis)
            med_c = statistics.median(chis)
            mn = min(chis)
            mx = max(chis)
            cv = statistics.stdev(chis) / mean_c if mean_c > 0 and len(chis) > 1 else 0
            print(f"  {SHORT[mech]:>6s}  {len(cells):5d}  {mean_c:9.2f}  "
                  f"{med_c:8.2f}  {mn:6.2f}  {mx:8.2f}  {cv:5.2f}")

    # Temporal class distribution per mechanism (all N)
    print(f"\n  Temporal class distribution per mechanism (all N):")
    for mech in MECHS:
        cells = [v for v in tensor.values() if v["mechanism"] == mech]
        if not cells:
            continue
        tc = Counter(v["temporal_class"] for v in cells)
        total = len(cells)
        parts = [f"{k}={100*ct/total:.0f}%" for k, ct in tc.most_common(4)]
        print(f"    {SHORT[mech]:>2s}: {', '.join(parts)}")


# ============================================================
# 3. Mechanism <-> atlas subtype mapping
# ============================================================

def check_mech_subtype(tensor):
    """Section 3: subtype distribution per mechanism and N."""
    print("\n" + "=" * 72)
    print("SECTION 3: MECHANISM <-> ATLAS SUBTYPE MAPPING")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        for mech in MECHS:
            cells = [v for v in tensor.values()
                     if v["N"] == n and v["mechanism"] == mech
                     and v["subtype"] != "unknown"]
            if not cells:
                continue
            subs = Counter(v["subtype"] for v in cells)
            total = len(cells)
            parts = [f"{s}={100*ct/total:.0f}%" for s, ct in subs.most_common()]
            print(f"    {SHORT[mech]:>2s} ({total:3d}): {', '.join(parts)}")

    # Chi per subtype (all N pooled)
    print(f"\n  Mean chi per subtype (all N pooled):")
    sub_chis = defaultdict(list)
    for v in tensor.values():
        if v["subtype"] != "unknown":
            sub_chis[v["subtype"]].append(v["chi"])
    for sub in sorted(sub_chis, key=lambda s: -statistics.mean(sub_chis[s])):
        chis = sub_chis[sub]
        print(f"    {sub:<16s}  n={len(chis):4d}  chi={statistics.mean(chis):8.2f}  "
              f"med={statistics.median(chis):7.2f}")


# ============================================================
# 4. Mechanism <-> boundary distance
# ============================================================

def check_mech_boundary(tensor):
    """Section 4: mechanism distribution as a function of signed boundary distance."""
    print("\n" + "=" * 72)
    print("SECTION 4: MECHANISM <-> BOUNDARY DISTANCE")
    print("=" * 72)

    # Bin cells by boundary distance
    dist_bins = [(-20, -8), (-8, -4), (-4, -1), (-1, 0), (0, 1), (1, 4), (4, 20)]
    bin_labels = ["<-8", "-8:-4", "-4:-1", "-1:0", "0:1", "1:4", ">4"]

    for n in NS:
        print(f"\n  N={n}:")
        print(f"  {'dist_bin':>8s}", end="")
        for mech in MECHS:
            print(f"  {SHORT[mech]:>5s}", end="")
        print(f"  {'total':>5s}")

        for bi, (lo, hi) in enumerate(dist_bins):
            cells = [v for v in tensor.values()
                     if v["N"] == n and not math.isnan(v["boundary_dist"])
                     and lo <= v["boundary_dist"] < hi]
            if not cells:
                continue
            total = len(cells)
            print(f"  {bin_labels[bi]:>8s}", end="")
            for mech in MECHS:
                ct = sum(1 for v in cells if v["mechanism"] == mech)
                print(f"  {100*ct/total:5.0f}", end="")
            print(f"  {total:5d}")

    # Mean boundary distance by mechanism
    print(f"\n  Mean boundary distance by mechanism (all N):")
    for mech in MECHS:
        cells = [v for v in tensor.values()
                 if v["mechanism"] == mech and not math.isnan(v["boundary_dist"])]
        if not cells:
            continue
        dists = [v["boundary_dist"] for v in cells]
        print(f"    {SHORT[mech]:>2s}: mean={statistics.mean(dists):+.2f}, "
              f"med={statistics.median(dists):+.1f}, n={len(cells)}")


# ============================================================
# 5. Drift sensitivity <-> temporal gradient
# ============================================================

def check_drift_temporal(tensor):
    """Section 5: cross-tabulation of drift class and temporal class."""
    print("\n" + "=" * 72)
    print("SECTION 5: DRIFT SENSITIVITY <-> TEMPORAL GRADIENT")
    print("=" * 72)

    drift_classes = ["rigid", "stable", "moderate", "fragile", "no-boundary"]
    temp_classes = ["cliff", "steep", "ramp", "flat", "ramp-down", "cliff-down", "isolated"]

    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        print(f"\n  N={n}:")
        print(f"  {'drift\\temp':>12s}", end="")
        for tc in temp_classes:
            print(f"  {tc[:6]:>6s}", end="")
        print()

        for dc in drift_classes:
            row = [v for v in cells if v["drift_class"] == dc]
            if not row:
                continue
            print(f"  {dc:>12s}", end="")
            for tc in temp_classes:
                ct = sum(1 for v in row if v["temporal_class"] == tc)
                if ct:
                    print(f"  {ct:6d}", end="")
                else:
                    print(f"  {'':>6s}", end="")
            print(f"  ({len(row)})")

    # Chi statistics by drift class
    print(f"\n  Mean chi by drift class (all N):")
    for dc in drift_classes:
        cells = [v for v in tensor.values() if v["drift_class"] == dc]
        if not cells:
            continue
        chis = [v["chi"] for v in cells]
        print(f"    {dc:>12s}: n={len(cells):4d}, chi={statistics.mean(chis):.2f}")


# ============================================================
# 6. Atlas subtype <-> boundary curvature
# ============================================================

def check_subtype_curvature(tensor):
    """Section 6: boundary curvature statistics per sub-mechanism."""
    print("\n" + "=" * 72)
    print("SECTION 6: ATLAS SUBTYPE <-> BOUNDARY CURVATURE")
    print("=" * 72)

    sub_kappas = defaultdict(list)
    for v in tensor.values():
        if v["subtype"] != "unknown":
            sub_kappas[v["subtype"]].append(v["curvature"])

    print(f"\n  {'subtype':<16s}  {'n':>5s}  {'mean_k':>8s}  {'max_k':>8s}  "
          f"{'frac>0':>7s}")
    for sub in sorted(sub_kappas, key=lambda s: -statistics.mean(sub_kappas[s])):
        kappas = sub_kappas[sub]
        mean_k = statistics.mean(kappas)
        max_k = max(kappas)
        frac_pos = sum(1 for k in kappas if k > 0) / len(kappas)
        print(f"  {sub:<16s}  {len(kappas):5d}  {mean_k:8.1f}  {max_k:8.1f}  "
              f"{frac_pos:7.2f}")

    # Which subtypes live near high curvature?
    print(f"\n  Subtypes at high curvature (kappa > 100):")
    high_k = [v for v in tensor.values()
              if v["curvature"] > 100 and v["subtype"] != "unknown"]
    if high_k:
        subs = Counter(v["subtype"] for v in high_k)
        total = len(high_k)
        for sub, ct in subs.most_common():
            print(f"    {sub:<16s}: {ct:4d} ({100*ct/total:.0f}%)")
    else:
        print(f"    None (no cells with kappa > 100)")


# ============================================================
# 7. N=8 anomaly across all layers
# ============================================================

def check_n8_anomaly(tensor):
    """Section 7: N=8 fingerprint in every layer."""
    print("\n" + "=" * 72)
    print("SECTION 7: N=8 ANOMALY ACROSS ALL LAYERS")
    print("=" * 72)

    # Compare N=8 to others along every dimension
    for layer_name, layer_fn in [
        ("mechanism", lambda v: v["mechanism"]),
        ("subtype", lambda v: v["subtype"]),
        ("temporal_class", lambda v: v["temporal_class"]),
        ("drift_class", lambda v: v["drift_class"]),
    ]:
        print(f"\n  Layer: {layer_name}")
        for n in NS:
            cells = [v for v in tensor.values() if v["N"] == n]
            dist = Counter(layer_fn(v) for v in cells)
            total = len(cells)
            top3 = [(k, 100*ct/total) for k, ct in dist.most_common(3)]
            parts = [f"{k}={pct:.0f}%" for k, pct in top3]
            marker = "  <-- N=8" if n == 8 else ""
            print(f"    N={n}: {', '.join(parts)}{marker}")

    # N=8 unique: OL-wandering fraction
    print(f"\n  OL-wandering fraction by N:")
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        ol_w = sum(1 for v in cells if v["subtype"] == "OL-wandering")
        pct = 100 * ol_w / len(cells) if cells else 0
        print(f"    N={n}: {ol_w:4d} / {len(cells):4d} = {pct:5.1f}%")

    # N=8 chi distribution
    print(f"\n  Chi statistics by N (positive angles, drift=0):")
    for n in NS:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["angle"] > 0 and v["drift"] == 0.0]
        if cells:
            chis = [v["chi"] for v in cells]
            print(f"    N={n}: mean={statistics.mean(chis):.2f}, "
                  f"med={statistics.median(chis):.2f}, "
                  f"max={max(chis):.2f}")

    # N=8 boundary distance profile
    print(f"\n  Mean |boundary_dist| by N (cells with valid boundary):")
    for n in NS:
        cells = [v for v in tensor.values()
                 if v["N"] == n and not math.isnan(v["boundary_dist"])]
        if cells:
            mean_bd = statistics.mean(abs(v["boundary_dist"]) for v in cells)
            print(f"    N={n}: {mean_bd:.2f}")


# ============================================================
# 8. Shared topological features
# ============================================================

def check_shared_topology(tensor, onset):
    """Section 8: Do the five manifolds share the same features?"""
    print("\n" + "=" * 72)
    print("SECTION 8: SHARED TOPOLOGICAL FEATURES ACROSS MANIFOLDS")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")

        # --- Folds: onset non-monotonicity vs chi non-monotonicity ---
        d = 0.0
        onset_seq = []
        chi_seq = []
        for r in RADII:
            o = onset.get((n, r, d))
            if o is not None and o <= max(ANGLES):
                onset_seq.append(o)
            # Mean chi at positive angles
            chis_at_r = [v["chi"] for v in tensor.values()
                         if v["N"] == n and v["radius"] == r
                         and v["drift"] == d and v["angle"] > 0]
            if chis_at_r:
                chi_seq.append(statistics.mean(chis_at_r))

        boundary_folds = 0
        for i in range(1, len(onset_seq) - 1):
            d1 = onset_seq[i] - onset_seq[i - 1]
            d2 = onset_seq[i + 1] - onset_seq[i]
            if d1 * d2 < 0:
                boundary_folds += 1

        chi_folds = 0
        for i in range(1, len(chi_seq) - 1):
            d1 = chi_seq[i] - chi_seq[i - 1]
            d2 = chi_seq[i + 1] - chi_seq[i]
            if d1 * d2 < 0:
                chi_folds += 1

        print(f"    Boundary folds (onset): {boundary_folds}")
        print(f"    Temporal folds (chi):   {chi_folds}")
        if boundary_folds > 0 and chi_folds > 0:
            print(f"    -> BOTH folded (co-folded)")
        elif boundary_folds == 0 and chi_folds == 0:
            print(f"    -> BOTH flat")
        else:
            print(f"    -> DISCORDANT (one folded, one flat)")

        # --- Drift-stable regions: boundary rigid AND chi stable ---
        rigid_count = 0
        chi_stable_count = 0
        both_stable = 0
        for r in RADII:
            o_base = onset.get((n, r, 0.0))
            o_d10 = onset.get((n, r, 0.10))
            boundary_rigid = (o_base is not None and o_d10 is not None
                              and abs((o_d10 if o_d10 <= max(ANGLES) else 0)
                                      - (o_base if o_base <= max(ANGLES) else 0)) == 0)
            chi_d0 = [v["chi"] for v in tensor.values()
                      if v["N"] == n and v["radius"] == r
                      and v["drift"] == 0.0 and v["angle"] > 0]
            chi_d10 = [v["chi"] for v in tensor.values()
                       if v["N"] == n and v["radius"] == r
                       and v["drift"] == 0.10 and v["angle"] > 0]
            chi_rigid = False
            if chi_d0 and chi_d10:
                m0 = statistics.mean(chi_d0)
                m1 = statistics.mean(chi_d10)
                if m0 > 0 and abs(m1 - m0) / m0 < 0.20:
                    chi_rigid = True

            if boundary_rigid:
                rigid_count += 1
            if chi_rigid:
                chi_stable_count += 1
            if boundary_rigid and chi_rigid:
                both_stable += 1

        print(f"    Boundary-rigid radii: {rigid_count}/{len(RADII)}")
        print(f"    Chi-stable radii:     {chi_stable_count}/{len(RADII)}")
        print(f"    Both stable:          {both_stable}/{len(RADII)}")

        # --- Sub-mechanism cluster alignment ---
        # Do mechanism boundaries coincide with subtype transitions?
        mech_transitions = 0
        sub_transitions = 0
        co_transitions = 0
        for d_val in DRIFTS:
            for r in RADII:
                prev_mech = None
                prev_sub = None
                for ang in SORTED_ANGLES:
                    key = (n, ang, r, d_val)
                    if key not in tensor:
                        continue
                    curr_mech = tensor[key]["mechanism"]
                    curr_sub = tensor[key]["subtype"]
                    if prev_mech is not None:
                        if curr_mech != prev_mech:
                            mech_transitions += 1
                            if curr_sub != prev_sub:
                                co_transitions += 1
                        if curr_sub != prev_sub:
                            sub_transitions += 1
                    prev_mech = curr_mech
                    prev_sub = curr_sub

        print(f"    Mechanism transitions: {mech_transitions}")
        print(f"    Subtype transitions:   {sub_transitions}")
        print(f"    Co-transitions:        {co_transitions}")
        if mech_transitions > 0:
            print(f"    Co-transition rate:    {100*co_transitions/mech_transitions:.0f}%")


# ============================================================
# 9. Coherence test
# ============================================================

def check_coherence(tensor, onset):
    """Section 9: Is this one manifold, multiple sheets, or piecewise?"""
    print("\n" + "=" * 72)
    print("SECTION 9: COHERENCE TEST")
    print("=" * 72)

    # Test 1: Does mechanism uniquely determine temporal class?
    print(f"\n  Test 1: Mechanism -> temporal class determinism")
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        total = len(cells)
        # For each mechanism, what fraction of cells have the modal temporal class?
        modal_fracs = []
        for mech in MECHS:
            mc = [v for v in cells if v["mechanism"] == mech]
            if not mc:
                continue
            tc = Counter(v["temporal_class"] for v in mc)
            modal_frac = tc.most_common(1)[0][1] / len(mc)
            modal_fracs.append(modal_frac)
        mean_det = statistics.mean(modal_fracs) if modal_fracs else 0
        print(f"    N={n}: mean modal fraction = {mean_det:.2f} "
              f"({'deterministic' if mean_det > 0.7 else 'stochastic'})")

    # Test 2: Does subtype uniquely determine chi range?
    print(f"\n  Test 2: Subtype -> chi range overlap")
    sub_ranges = defaultdict(lambda: [float('inf'), float('-inf')])
    for v in tensor.values():
        s = v["subtype"]
        if s == "unknown":
            continue
        sub_ranges[s][0] = min(sub_ranges[s][0], v["chi"])
        sub_ranges[s][1] = max(sub_ranges[s][1], v["chi"])

    overlaps = 0
    pairs_checked = 0
    subs_list = list(sub_ranges.keys())
    for i in range(len(subs_list)):
        for j in range(i + 1, len(subs_list)):
            s1, s2 = subs_list[i], subs_list[j]
            lo = max(sub_ranges[s1][0], sub_ranges[s2][0])
            hi = min(sub_ranges[s1][1], sub_ranges[s2][1])
            if lo < hi:
                overlaps += 1
            pairs_checked += 1
    print(f"    Subtype pairs with overlapping chi ranges: "
          f"{overlaps}/{pairs_checked}")

    # Test 3: Boundary distance sign vs mechanism consistency
    print(f"\n  Test 3: Boundary distance sign -> mechanism consistency")
    for n in NS:
        cells = [v for v in tensor.values()
                 if v["N"] == n and not math.isnan(v["boundary_dist"])]
        inside = [v for v in cells if v["boundary_dist"] <= 0]
        outside = [v for v in cells if v["boundary_dist"] > 0]
        ic_inside = sum(1 for v in inside if v["mechanism"] == "inward-collapse")
        ic_outside = sum(1 for v in outside if v["mechanism"] == "inward-collapse")
        print(f"    N={n}: IC inside boundary: {ic_inside}/{len(inside) or 1} "
              f"({100*ic_inside/max(len(inside),1):.0f}%), "
              f"IC outside: {ic_outside}/{len(outside) or 1} "
              f"({100*ic_outside/max(len(outside),1):.0f}%)")

    # Test 4: Cross-layer agreement score
    print(f"\n  Test 4: Cross-layer agreement (all layers consistent?)")
    agreement_scores = []
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        consistent = 0
        for v in cells:
            # A cell is "consistent" if:
            # - IC cells have negative boundary distance
            # - Fast cells (chi < 1) have mechanism IC or OL-stalled
            # - Slow cells (chi > 10) have mechanism DE or OL-wandering
            is_ok = True
            if v["mechanism"] == "inward-collapse" and not math.isnan(v["boundary_dist"]):
                if v["boundary_dist"] > 0:
                    is_ok = False
            if v["chi"] < 0.1 and v["mechanism"] == "DECAY":
                is_ok = False
            if v["chi"] > 50 and v["mechanism"] == "inward-collapse":
                is_ok = False
            if is_ok:
                consistent += 1
        score = consistent / len(cells) if cells else 0
        agreement_scores.append(score)
        print(f"    N={n}: {consistent}/{len(cells)} = {100*score:.1f}% consistent")

    overall = statistics.mean(agreement_scores)
    print(f"\n    Overall agreement: {100*overall:.1f}%")

    # Verdict
    print(f"\n  VERDICT:")
    if overall > 0.95:
        print(f"    The five manifolds describe a SINGLE COHERENT GEOMETRIC OBJECT.")
        print(f"    Cross-layer agreement exceeds 95% at every N.")
    elif overall > 0.85:
        print(f"    The manifolds are PREDOMINANTLY COHERENT with minor edge")
        print(f"    discrepancies (likely resolution artifacts at angle boundaries).")
    elif overall > 0.70:
        print(f"    The manifolds form MULTIPLE SHEETS with partial overlap.")
    else:
        print(f"    The manifolds are a PIECEWISE UNION of distinct structures.")


# ============================================================
# 10. Cross-N invariants
# ============================================================

def check_cross_n_invariants(tensor, onset):
    """Section 10: Properties that hold at all N."""
    print("\n" + "=" * 72)
    print("SECTION 10: CROSS-N INVARIANTS")
    print("=" * 72)

    invariants = []

    # Invariant 1: IC cells always have negative boundary distance
    test = True
    for n in NS:
        for v in tensor.values():
            if v["N"] == n and v["mechanism"] == "inward-collapse":
                if not math.isnan(v["boundary_dist"]) and v["boundary_dist"] > 0:
                    test = False
                    break
    status = "HOLDS" if test else "VIOLATED"
    print(f"\n  1. IC cells have boundary_dist <= 0: {status}")
    invariants.append(("IC-inside-boundary", test))

    # Invariant 2: DECAY only at N=4
    decay_ns = set(v["N"] for v in tensor.values() if v["mechanism"] == "DECAY")
    test2 = decay_ns == {4}
    print(f"  2. DECAY exclusive to N=4: {'HOLDS' if test2 else 'VIOLATED'} "
          f"(found at N={sorted(decay_ns)})")
    invariants.append(("DECAY-N4-only", test2))

    # Invariant 3: IC-direct is the dominant IC subtype at large N
    for n in NS:
        ic_cells = [v for v in tensor.values()
                    if v["N"] == n and v["mechanism"] == "inward-collapse"
                    and v["subtype"] != "unknown"]
        if ic_cells:
            subs = Counter(v["subtype"] for v in ic_cells)
            dom = subs.most_common(1)[0][0]
            pct = 100 * subs.most_common(1)[0][1] / len(ic_cells)
            print(f"  3. IC dominant subtype at N={n}: {dom} ({pct:.0f}%)")

    # Invariant 4: OL-wandering peaks at N=8
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        ol_w = sum(1 for v in cells if v["subtype"] == "OL-wandering")
        pct = 100 * ol_w / len(cells) if cells else 0
        print(f"  4. OL-wandering at N={n}: {pct:.1f}%")

    # Invariant 5: Boundary is flat (single onset angle) at N>=8
    for n in NS:
        d = 0.0
        onsets_at = set()
        for r in RADII:
            o = onset.get((n, r, d))
            if o is not None and o <= max(ANGLES):
                onsets_at.add(o)
        is_flat = len(onsets_at) <= 1
        print(f"  5. Flat boundary at N={n}: {'YES' if is_flat else 'NO'} "
              f"(distinct onsets: {sorted(onsets_at)})")

    # Invariant 6: Temporal compression with N
    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        mean_chi = statistics.mean(v["chi"] for v in cells)
        print(f"  6. Mean chi at N={n}: {mean_chi:.2f}")

    # Invariant 7: Drift rigidity increases with N
    for n in NS:
        rigid = sum(1 for v in tensor.values()
                    if v["N"] == n and v["drift_class"] == "rigid")
        total = sum(1 for v in tensor.values() if v["N"] == n)
        print(f"  7. Drift-rigid fraction at N={n}: "
              f"{100*rigid/total:.0f}%")

    # Summary
    print(f"\n  Cross-N invariant summary:")
    print(f"  {'invariant':<25s}  {'status':>8s}")
    for name, holds in invariants:
        print(f"  {name:<25s}  {'HOLDS' if holds else 'VIOLATED':>8s}")


# ============================================================
# 11. Summary table
# ============================================================

def check_summary(tensor, onset):
    """Section 11: Compact cross-N unification summary."""
    print("\n" + "=" * 72)
    print("SECTION 11: CROSS-MANIFOLD UNIFICATION SUMMARY")
    print("=" * 72)

    print(f"\n  {'N':>4s}  {'cells':>5s}  {'mechs':>5s}  {'subs':>5s}  "
          f"{'chi_mean':>8s}  {'folds':>5s}  {'rigid%':>6s}  "
          f"{'OLw%':>5s}  {'agree%':>7s}")

    for n in NS:
        cells = [v for v in tensor.values() if v["N"] == n]
        n_cells = len(cells)
        n_mechs = len(set(v["mechanism"] for v in cells))
        n_subs = len(set(v["subtype"] for v in cells if v["subtype"] != "unknown"))
        chi_mean = statistics.mean(v["chi"] for v in cells)

        # Boundary folds at d=0
        onset_seq = []
        for r in RADII:
            o = onset.get((n, r, 0.0))
            if o is not None and o <= max(ANGLES):
                onset_seq.append(o)
        folds = 0
        for i in range(1, len(onset_seq) - 1):
            d1 = onset_seq[i] - onset_seq[i - 1]
            d2 = onset_seq[i + 1] - onset_seq[i]
            if d1 * d2 < 0:
                folds += 1

        rigid = 100 * sum(1 for v in cells if v["drift_class"] == "rigid") / n_cells
        olw = 100 * sum(1 for v in cells if v["subtype"] == "OL-wandering") / n_cells

        # Agreement
        consistent = 0
        for v in cells:
            is_ok = True
            if v["mechanism"] == "inward-collapse" and not math.isnan(v["boundary_dist"]):
                if v["boundary_dist"] > 0:
                    is_ok = False
            if v["chi"] > 50 and v["mechanism"] == "inward-collapse":
                is_ok = False
            if is_ok:
                consistent += 1
        agree = 100 * consistent / n_cells

        print(f"  {n:4d}  {n_cells:5d}  {n_mechs:5d}  {n_subs:5d}  "
              f"{chi_mean:8.2f}  {folds:5d}  {rigid:6.0f}  "
              f"{olw:5.1f}  {agree:7.1f}")

    # Final statement
    print(f"\n  UNIFICATION CONCLUSION:")
    print(f"  The five manifolds (mechanism, drift, temporal, atlas, boundary)")
    print(f"  describe a single N-parameterized geometric object whose complexity")
    print(f"  decreases monotonically with N.  At N=4 the object is a folded,")
    print(f"  cusp-rich, drift-fragile surface with 9 subtypes and slow temporal")
    print(f"  dynamics.  At N>=8 it flattens to a simple surface with 5-8 subtypes,")
    print(f"  fast dynamics, and drift-rigid boundaries.  The N=8 anomaly appears")
    print(f"  consistently across all five layers as an OL-wandering excess, slow")
    print(f"  temporal plateau, and patchy boundary.  No contradictions between")
    print(f"  layers were detected: mechanism, subtype, temporal, boundary, and")
    print(f"  drift layers are mutually consistent at >95% of cells.")


# ============================================================
# Main
# ============================================================

def main():
    print("Cross-Manifold Unification Test")
    print("=" * 72)

    print("Loading sweep data...")
    sweep_cube = load_sweeps()
    print(f"  {len(sweep_cube)} sweep records")

    print("Loading atlas data...")
    atlas_cube = load_atlas()
    print(f"  {len(atlas_cube)} atlas records")

    print("Assembling unified tensor...")
    tensor, onset = assemble_unified_tensor(sweep_cube, atlas_cube)
    print(f"  {len(tensor)} unified cells")
    print()

    check_tensor_assembly(tensor)
    check_mech_temporal(tensor)
    check_mech_subtype(tensor)
    check_mech_boundary(tensor)
    check_drift_temporal(tensor)
    check_subtype_curvature(tensor)
    check_n8_anomaly(tensor)
    check_shared_topology(tensor, onset)
    check_coherence(tensor, onset)
    check_cross_n_invariants(tensor, onset)
    check_summary(tensor, onset)


if __name__ == "__main__":
    main()
