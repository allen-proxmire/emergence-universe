#!/usr/bin/env python
"""
check_mechanism_atlas.py
=========================
Mechanism Atlas analysis: resolves the internal structure of OPBC,
PBC-corner, and other-late into sub-mechanism families using trajectory
diagnostics from atlas sweep files.

Sub-mechanism classification rules:
  IC sub-types:
    IC-direct     : n_steps <= 50, no PBC crossings
    IC-delayed    : n_steps > 50, no PBC crossings
    IC-bounced    : any PBC crossings before collapse

  OPBC sub-types:
    OPBC-single   : 0 PBC crossings (edge collapse without wrapping)
    OPBC-glancing : 1-5 PBC crossings (single pass through boundary)
    OPBC-winding  : 6+ PBC crossings (multi-wrap orbit)

  PBC-corner sub-types:
    CR-single     : only 1 corner visited
    CR-multi      : 2+ distinct corners visited
    CR-diagonal   : corner_sequence contains [0,3] or [1,2] (diagonal pair)

  Other-late sub-types:
    OL-stalled    : n_steps < 5 (immediate catch-all)
    OL-wandering  : n_steps >= 5, no periodicity
    OL-periodic   : periodicity_score > 0.3
    OL-multi      : is_multi_stage = True

  DECAY sub-types:
    DE-tangent    : tangent_crossings > 0 (oscillating around tangent)
    DE-frozen     : tangent_crossings = 0 (static orbit)
"""

import json
import os
import glob
import re
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


# ============================================================
# Data loading
# ============================================================

def load_all_atlas():
    """Load all atlas files. Returns list of dicts."""
    records = []
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_atlas.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        if "atlas" not in data:
            continue
        records.append(data)
    return records


# ============================================================
# Sub-mechanism classification
# ============================================================

def classify_sub(rec):
    """Assign a sub-mechanism label based on atlas diagnostics."""
    mech = rec["mechanism"]
    a = rec["atlas"]
    pbc = a["n_pbc_crossings"]
    steps = a["n_steps"]
    corners = a["corner_approaches"]
    seq = a["corner_sequence"]
    period = a["periodicity_score"]
    multi = a["is_multi_stage"]
    tcross = a["tangent_crossings"]

    if mech == "inward-collapse":
        if pbc > 0:
            return "IC-bounced"
        elif steps > 50:
            return "IC-delayed"
        else:
            return "IC-direct"

    elif mech == "outward-PBC":
        if pbc == 0:
            return "OPBC-single"
        elif pbc <= 5:
            return "OPBC-glancing"
        else:
            return "OPBC-winding"

    elif mech == "PBC-corner":
        distinct_corners = sum(1 for c in corners if c > 0)
        if distinct_corners <= 1:
            return "CR-single"
        else:
            # Check for diagonal pair
            has_diag = False
            for i in range(len(seq) - 1):
                pair = {seq[i], seq[i+1]}
                if pair == {0, 3} or pair == {1, 2}:
                    has_diag = True
                    break
            if has_diag:
                return "CR-diagonal"
            return "CR-multi"

    elif mech == "other-late":
        if multi:
            return "OL-multi"
        elif period > 0.3:
            return "OL-periodic"
        elif steps < 5:
            return "OL-stalled"
        else:
            return "OL-wandering"

    elif mech == "DECAY":
        if tcross > 0:
            return "DE-tangent"
        else:
            return "DE-frozen"

    return mech


# ============================================================
# Section 1: Overall taxonomy
# ============================================================

def print_taxonomy(records):
    print("=" * 90)
    print("  MECHANISM ATLAS: SUB-MECHANISM TAXONOMY")
    print("=" * 90)

    # Classify all records
    for rec in records:
        rec["sub_mech"] = classify_sub(rec)

    # Count by parent mechanism and sub-mechanism
    by_parent = defaultdict(lambda: Counter())
    for rec in records:
        parent = SHORT[rec["mechanism"]]
        by_parent[parent][rec["sub_mech"]] += 1

    total = len(records)
    print(f"\n  Total records: {total}")

    for parent in ["IC", "OP", "DE", "CR", "OL"]:
        subs = by_parent[parent]
        if not subs:
            continue
        parent_total = sum(subs.values())
        print(f"\n  {parent} ({parent_total} cells, "
              f"{100*parent_total/total:.1f}%):")
        for sub, count in sorted(subs.items(), key=lambda x: -x[1]):
            pct = 100 * count / parent_total
            bar = "#" * max(1, int(pct / 3))
            print(f"    {sub:>15}: {count:5d}  ({pct:5.1f}%)  {bar}")


# ============================================================
# Section 2: Sub-mechanism profiles
# ============================================================

def print_sub_profiles(records):
    print("\n" + "=" * 90)
    print("  SUB-MECHANISM PROFILES: diagnostic fingerprints")
    print("=" * 90)

    by_sub = defaultdict(list)
    for rec in records:
        by_sub[rec["sub_mech"]].append(rec)

    print(f"\n  {'sub_mech':>15}  {'n':>4}  {'chi':>8}  {'steps':>7}  "
          f"{'pbc':>5}  {'edge':>5}  {'R_rng':>7}  {'tcross':>6}  {'period':>7}")
    print(f"  " + "-" * 72)

    for sub in sorted(by_sub):
        recs = by_sub[sub]
        chis = [r["chi_emp"] for r in recs]
        steps = [r["atlas"]["n_steps"] for r in recs]
        pbcs = [r["atlas"]["n_pbc_crossings"] for r in recs]
        edges = [r["atlas"]["n_edge_approaches"] for r in recs]
        r_ranges = [r["atlas"]["R_range"] for r in recs]
        tcrosses = [r["atlas"]["tangent_crossings"] for r in recs]
        periods = [r["atlas"]["periodicity_score"] for r in recs]

        mn = lambda vs: statistics.mean(vs) if vs else 0
        print(f"  {sub:>15}  {len(recs):4d}  {mn(chis):8.2f}  "
              f"{mn(steps):7.0f}  {mn(pbcs):5.1f}  {mn(edges):5.0f}  "
              f"{mn(r_ranges):7.4f}  {mn(tcrosses):6.1f}  {mn(periods):7.3f}")


# ============================================================
# Section 3: Sub-mechanism distribution by N
# ============================================================

def print_sub_by_n(records):
    print("\n" + "=" * 90)
    print("  SUB-MECHANISM DISTRIBUTION BY N")
    print("=" * 90)

    for n in NS:
        n_recs = [r for r in records if r["N"] == n]
        if not n_recs:
            continue

        subs = Counter(r["sub_mech"] for r in n_recs)
        total = len(n_recs)
        print(f"\n  N = {n} ({total} cells):")

        for sub, count in sorted(subs.items(), key=lambda x: -x[1]):
            pct = 100 * count / total
            bar = "#" * max(1, int(pct / 2))
            print(f"    {sub:>15}: {count:4d}  ({pct:5.1f}%)  {bar}")


# ============================================================
# Section 4: OPBC sub-structure
# ============================================================

def print_opbc_structure(records):
    print("\n" + "=" * 90)
    print("  OPBC INTERNAL STRUCTURE")
    print("=" * 90)

    opbc = [r for r in records if r["mechanism"] == "outward-PBC"]
    if not opbc:
        print("  No OPBC cells")
        return

    print(f"\n  Total OPBC: {len(opbc)}")

    # PBC crossing distribution
    crossing_counts = [r["atlas"]["n_pbc_crossings"] for r in opbc]
    print(f"\n  PBC crossing distribution:")
    bins = [(0, 0, "zero"), (1, 2, "1-2"), (3, 5, "3-5"),
            (6, 20, "6-20"), (21, 100, "21-100"), (101, 99999, ">100")]
    for lo, hi, label in bins:
        n = sum(1 for c in crossing_counts if lo <= c <= hi)
        if n > 0:
            pct = 100 * n / len(opbc)
            print(f"    {label:>8}: {n:4d}  ({pct:5.1f}%)")

    # Winding number distribution
    print(f"\n  Winding number magnitude:")
    wind_mags = [max(abs(r["atlas"]["winding_number"][0]),
                     abs(r["atlas"]["winding_number"][1]))
                 for r in opbc]
    for thr, label in [(0.01, "< 0.01"), (0.1, "0.01-0.1"),
                       (0.5, "0.1-0.5"), (1.0, "0.5-1.0"),
                       (999, "> 1.0")]:
        prev_thr = 0 if label.startswith("<") else float(label.split("-")[0]) if "-" in label else float(label.split("> ")[1]) if "> " in label else 0
        n = sum(1 for w in wind_mags
                if (w < thr if label.startswith("<") else
                    w >= float(label.split("-")[0] if "-" in label else label.replace("> ", ""))
                    and w < thr))
        if n > 0:
            print(f"    {label:>12}: {n:4d}  ({100*n/len(opbc):5.1f}%)")

    # By N
    print(f"\n  OPBC sub-types by N:")
    print(f"  {'N':>3}  {'single':>8}  {'glancing':>10}  {'winding':>9}  {'total':>6}")
    print(f"  " + "-" * 42)
    for n in NS:
        n_opbc = [r for r in opbc if r["N"] == n]
        subs = Counter(r["sub_mech"] for r in n_opbc)
        print(f"  {n:3d}  {subs.get('OPBC-single', 0):8d}  "
              f"{subs.get('OPBC-glancing', 0):10d}  "
              f"{subs.get('OPBC-winding', 0):9d}  {len(n_opbc):6d}")


# ============================================================
# Section 5: PBC-corner sub-structure
# ============================================================

def print_corner_structure(records):
    print("\n" + "=" * 90)
    print("  PBC-CORNER INTERNAL STRUCTURE")
    print("=" * 90)

    cr = [r for r in records if r["mechanism"] == "PBC-corner"]
    if not cr:
        print("  No PBC-corner cells")
        return

    print(f"\n  Total PBC-corner: {len(cr)}")

    # Corner count distribution
    print(f"\n  Distinct corners visited:")
    for nc in range(5):
        n = sum(1 for r in cr
                if sum(1 for c in r["atlas"]["corner_approaches"] if c > 0) == nc)
        if n > 0:
            pct = 100 * n / len(cr)
            print(f"    {nc} corners: {n:4d}  ({pct:5.1f}%)")

    # Most common corner sequences
    print(f"\n  Most common corner sequences (first 3 entries):")
    seq_counter = Counter()
    for r in cr:
        seq = tuple(r["atlas"]["corner_sequence"][:3])
        if seq:
            seq_counter[seq] += 1
    for seq, count in seq_counter.most_common(10):
        print(f"    {list(seq)}: {count}")

    # By N
    print(f"\n  PBC-corner sub-types by N:")
    print(f"  {'N':>3}  {'single':>8}  {'multi':>7}  {'diagonal':>10}  {'total':>6}")
    print(f"  " + "-" * 38)
    for n in NS:
        n_cr = [r for r in cr if r["N"] == n]
        subs = Counter(r["sub_mech"] for r in n_cr)
        print(f"  {n:3d}  {subs.get('CR-single', 0):8d}  "
              f"{subs.get('CR-multi', 0):7d}  "
              f"{subs.get('CR-diagonal', 0):10d}  {len(n_cr):6d}")


# ============================================================
# Section 6: Other-late sub-structure
# ============================================================

def print_ol_structure(records):
    print("\n" + "=" * 90)
    print("  OTHER-LATE INTERNAL STRUCTURE")
    print("=" * 90)

    ol = [r for r in records if r["mechanism"] == "other-late"]
    if not ol:
        print("  No other-late cells")
        return

    print(f"\n  Total other-late: {len(ol)}")

    # Step distribution
    steps = [r["atlas"]["n_steps"] for r in ol]
    print(f"\n  Step count distribution:")
    for lo, hi, label in [(0, 1, "0-1"), (2, 5, "2-5"), (6, 50, "6-50"),
                           (51, 500, "51-500"), (501, 5000, "501-5000"),
                           (5001, 99999, ">5000")]:
        n = sum(1 for s in steps if lo <= s <= hi)
        if n > 0:
            pct = 100 * n / len(ol)
            print(f"    {label:>10}: {n:4d}  ({pct:5.1f}%)")

    # Multi-stage vs single
    multi = sum(1 for r in ol if r["atlas"]["is_multi_stage"])
    print(f"\n  Multi-stage: {multi}/{len(ol)} ({100*multi/len(ol):.1f}%)")

    # Stage transition patterns for multi-stage OL
    if multi > 0:
        patterns = Counter()
        for r in ol:
            if r["atlas"]["is_multi_stage"]:
                patterns[" -> ".join(r["atlas"]["stage_sequence"])] += 1
        print(f"  Multi-stage transition patterns:")
        for pattern, count in patterns.most_common(10):
            print(f"    {pattern}: {count}")

    # By N
    print(f"\n  OL sub-types by N:")
    print(f"  {'N':>3}  {'stalled':>9}  {'wandering':>11}  "
          f"{'periodic':>10}  {'multi':>7}  {'total':>6}")
    print(f"  " + "-" * 52)
    for n in NS:
        n_ol = [r for r in ol if r["N"] == n]
        subs = Counter(r["sub_mech"] for r in n_ol)
        print(f"  {n:3d}  {subs.get('OL-stalled', 0):9d}  "
              f"{subs.get('OL-wandering', 0):11d}  "
              f"{subs.get('OL-periodic', 0):10d}  "
              f"{subs.get('OL-multi', 0):7d}  {len(n_ol):6d}")


# ============================================================
# Section 7: N=8 anomaly — sub-mechanism fingerprint
# ============================================================

def print_n8_fingerprint(records):
    print("\n" + "=" * 90)
    print("  N=8 ANOMALY: SUB-MECHANISM FINGERPRINT")
    print("=" * 90)

    n8 = [r for r in records if r["N"] == 8]
    others = [r for r in records if r["N"] != 8]

    if not n8:
        return

    # Compare sub-mechanism distributions
    n8_subs = Counter(r["sub_mech"] for r in n8)
    other_subs = Counter(r["sub_mech"] for r in others)
    all_subs = sorted(set(list(n8_subs.keys()) + list(other_subs.keys())))

    print(f"\n  Sub-mechanism fractions: N=8 vs rest")
    print(f"  {'sub_mech':>15}  {'N=8':>8}  {'rest':>8}  {'delta':>8}")
    print(f"  " + "-" * 42)

    n8_tot = len(n8)
    ot_tot = len(others)
    unique_to_n8 = []
    for sub in all_subs:
        n8_pct = 100 * n8_subs.get(sub, 0) / n8_tot
        ot_pct = 100 * other_subs.get(sub, 0) / ot_tot
        delta = n8_pct - ot_pct
        flag = " ***" if abs(delta) > 5 else ""
        print(f"  {sub:>15}  {n8_pct:7.1f}%  {ot_pct:7.1f}%  "
              f"{delta:+7.1f}%{flag}")
        if n8_pct > 0 and ot_pct == 0:
            unique_to_n8.append(sub)

    if unique_to_n8:
        print(f"\n  Sub-mechanisms UNIQUE to N=8: {unique_to_n8}")

    # N=8 OL deep-dive
    n8_ol = [r for r in n8 if r["mechanism"] == "other-late"]
    if n8_ol:
        print(f"\n  N=8 other-late deep-dive ({len(n8_ol)} cells):")
        chis = [r["chi_emp"] for r in n8_ol]
        steps = [r["atlas"]["n_steps"] for r in n8_ol]
        pbcs = [r["atlas"]["n_pbc_crossings"] for r in n8_ol]
        edges = [r["atlas"]["n_edge_approaches"] for r in n8_ol]
        r_ranges = [r["atlas"]["R_range"] for r in n8_ol]
        tcrosses = [r["atlas"]["tangent_crossings"] for r in n8_ol]

        print(f"    chi:   mean={statistics.mean(chis):.2f}  "
              f"med={statistics.median(chis):.2f}  "
              f"range=[{min(chis):.2f}, {max(chis):.2f}]")
        print(f"    steps: mean={statistics.mean(steps):.0f}  "
              f"range=[{min(steps)}, {max(steps)}]")
        print(f"    pbc_crossings: mean={statistics.mean(pbcs):.1f}  "
              f"range=[{min(pbcs)}, {max(pbcs)}]")
        print(f"    R_range: mean={statistics.mean(r_ranges):.4f}")
        print(f"    tangent_crossings: mean={statistics.mean(tcrosses):.1f}")

        # OL sub-types at N=8
        ol_subs = Counter(r["sub_mech"] for r in n8_ol)
        print(f"    Sub-types: {dict(ol_subs)}")


# ============================================================
# Section 8: Large-N taxonomy reduction
# ============================================================

def print_large_n_reduction(records):
    print("\n" + "=" * 90)
    print("  LARGE-N TAXONOMY REDUCTION")
    print("=" * 90)

    print(f"\n  Distinct sub-mechanisms by N:")
    print(f"  {'N':>3}  {'n_sub':>5}  {'sub-mechanisms'}")
    print(f"  " + "-" * 60)

    for n in NS:
        n_recs = [r for r in records if r["N"] == n]
        subs = sorted(set(r["sub_mech"] for r in n_recs))
        # Filter to >1% prevalence
        sub_counts = Counter(r["sub_mech"] for r in n_recs)
        significant = [s for s in subs
                       if sub_counts[s] / len(n_recs) > 0.01]
        print(f"  {n:3d}  {len(significant):5d}  {', '.join(significant)}")

    # Convergence test: do the top-3 sub-mechanisms stabilize?
    print(f"\n  Top-3 sub-mechanisms by N (>10% each):")
    for n in NS:
        n_recs = [r for r in records if r["N"] == n]
        sub_counts = Counter(r["sub_mech"] for r in n_recs)
        total = len(n_recs)
        top3 = [(s, c, 100*c/total)
                for s, c in sub_counts.most_common(3)]
        parts = [f"{s}({pct:.0f}%)" for s, _, pct in top3]
        print(f"  N={n:2d}: {', '.join(parts)}")


# ============================================================
# Section 9: Cross-N invariants
# ============================================================

def print_cross_n_invariants(records):
    print("\n" + "=" * 90)
    print("  CROSS-N INVARIANTS IN SUB-MECHANISM STRUCTURE")
    print("=" * 90)

    # Check which sub-mechanisms appear at ALL N values
    all_n_subs = {}
    for n in NS:
        n_recs = [r for r in records if r["N"] == n]
        all_n_subs[n] = set(r["sub_mech"] for r in n_recs)

    universal = set.intersection(*all_n_subs.values())
    print(f"\n  Universal sub-mechanisms (present at all N): {sorted(universal)}")

    n_specific = {}
    for n in NS:
        specific = all_n_subs[n] - set.union(*[all_n_subs[m]
                                                for m in NS if m != n])
        if specific:
            n_specific[n] = specific
    if n_specific:
        print(f"\n  N-specific sub-mechanisms:")
        for n, subs in sorted(n_specific.items()):
            print(f"    N={n}: {sorted(subs)}")
    else:
        print(f"\n  No N-specific sub-mechanisms found")

    # Diagnostic invariants
    print(f"\n  Mean PBC crossings per mechanism, by N:")
    print(f"  {'N':>3}  {'IC':>6}  {'OP':>6}  {'CR':>6}  {'OL':>6}  {'DE':>6}")
    print(f"  " + "-" * 36)
    for n in NS:
        row = f"  {n:3d}"
        for mech in MECHS:
            vals = [r["atlas"]["n_pbc_crossings"]
                    for r in records
                    if r["N"] == n and r["mechanism"] == mech]
            if vals:
                row += f"  {statistics.mean(vals):5.1f}"
            else:
                row += "    ---"
        print(row)

    print(f"\n  Mean R_range per mechanism, by N:")
    print(f"  {'N':>3}  {'IC':>8}  {'OP':>8}  {'CR':>8}  {'OL':>8}  {'DE':>8}")
    print(f"  " + "-" * 46)
    for n in NS:
        row = f"  {n:3d}"
        for mech in MECHS:
            vals = [r["atlas"]["R_range"]
                    for r in records
                    if r["N"] == n and r["mechanism"] == mech]
            if vals:
                row += f"  {statistics.mean(vals):7.4f}"
            else:
                row += "      ---"
        print(row)


# ============================================================
# Section 10: Drift effect on sub-mechanism taxonomy
# ============================================================

def print_drift_effect(records):
    print("\n" + "=" * 90)
    print("  DRIFT EFFECT ON SUB-MECHANISM TAXONOMY")
    print("=" * 90)

    print(f"\n  Distinct sub-mechanism count by (N, drift):")
    print(f"  {'N':>3}", end="")
    for d in DRIFTS:
        print(f"  d={d:.2f}", end="")
    print()
    print(f"  " + "-" * (4 + 8 * len(DRIFTS)))

    for n in NS:
        row = f"  {n:3d}"
        for d in DRIFTS:
            recs = [r for r in records if r["N"] == n
                    and abs(r["drift"] - d) < 0.001]
            subs = set(r["sub_mech"] for r in recs)
            row += f"  {len(subs):6d}"
        print(row)

    # IC sub-type evolution with drift
    print(f"\n  IC sub-type fractions vs drift:")
    for n in NS:
        print(f"\n  N={n}:")
        for d in DRIFTS:
            ic_recs = [r for r in records if r["N"] == n
                       and abs(r["drift"] - d) < 0.001
                       and r["mechanism"] == "inward-collapse"]
            if not ic_recs:
                continue
            subs = Counter(r["sub_mech"] for r in ic_recs)
            total = len(ic_recs)
            parts = [f"{s}:{100*c/total:.0f}%"
                    for s, c in sorted(subs.items())]
            print(f"    d={d:.2f} (n={total:3d}): {', '.join(parts)}")


# ============================================================
# Section 11: Summary table
# ============================================================

def print_summary(records):
    print("\n" + "=" * 90)
    print("  MECHANISM ATLAS SUMMARY")
    print("=" * 90)

    by_parent = defaultdict(set)
    for rec in records:
        by_parent[rec["mechanism"]].add(rec["sub_mech"])

    print(f"\n  {'Parent':>18}  {'n_sub':>5}  sub-mechanisms")
    print(f"  " + "-" * 70)
    total_subs = 0
    for mech in MECHS:
        subs = sorted(by_parent.get(mech, set()))
        total_subs += len(subs)
        print(f"  {SHORT[mech]:>18}  {len(subs):5d}  {', '.join(subs)}")
    print(f"\n  Total sub-mechanisms: {total_subs}")

    # Taxonomy by N
    print(f"\n  Taxonomy size by N:")
    for n in NS:
        n_recs = [r for r in records if r["N"] == n]
        n_subs = len(set(r["sub_mech"] for r in n_recs))
        n_parents = len(set(r["mechanism"] for r in n_recs))
        print(f"    N={n:2d}: {n_parents} parent mechanisms, "
              f"{n_subs} sub-mechanisms")

    # Final verdict
    print(f"\n  Verdict:")
    print(f"    IC:  {len(by_parent.get('inward-collapse', set()))} sub-types "
          f"(direct/delayed/bounced)")
    print(f"    OP:  {len(by_parent.get('outward-PBC', set()))} sub-types "
          f"(single/glancing/winding)")
    print(f"    CR:  {len(by_parent.get('PBC-corner', set()))} sub-types "
          f"(single/multi/diagonal)")
    print(f"    OL:  {len(by_parent.get('other-late', set()))} sub-types "
          f"(stalled/wandering/periodic/multi)")
    print(f"    DE:  {len(by_parent.get('DECAY', set()))} sub-types "
          f"(tangent/frozen)")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 90)
    print("  MECHANISM ATLAS ANALYSIS")
    print("=" * 90)
    print()

    records = load_all_atlas()
    print(f"  Loaded {len(records)} atlas records\n")

    if not records:
        print("  ERROR: No atlas files found.")
        return

    print_taxonomy(records)
    print_sub_profiles(records)
    print_sub_by_n(records)
    print_opbc_structure(records)
    print_corner_structure(records)
    print_ol_structure(records)
    print_n8_fingerprint(records)
    print_large_n_reduction(records)
    print_cross_n_invariants(records)
    print_drift_effect(records)
    print_summary(records)

    print("\n" + "=" * 90)


if __name__ == "__main__":
    main()
