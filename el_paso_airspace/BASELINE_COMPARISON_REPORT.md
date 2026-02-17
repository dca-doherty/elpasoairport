# KEPZ Baseline Comparison Report

## Feb 10, 2026 (Control Night) vs Feb 11, 2026 (Event Night)

### Purpose

This analysis answers the critical question: **Do the KEPZ radar returns on the event night look different from a normal night?** Without this comparison, we cannot distinguish engagement debris from persistent ground clutter.

### Methodology

- Downloaded 29 KEPZ Level 2 scans from Feb 10, 2026 (03:00–05:04 UTC)
- Applied identical sector parameters: az 60–110°, range 15–55 km
- Used identical thresholds, dual-pol classification, and analysis pipeline
- Split into North cluster (az 60–90°, Franklin Mountains) and South cluster (az 90–110°)
- Ran Welch's t-test, Kolmogorov-Smirnov, and Mann-Whitney U tests

---

## Key Findings

### 1. The Baseline ALSO Has Massive Returns

| Metric | Baseline (Feb 10) | Event (Feb 11) | Delta |
|--------|-------------------|----------------|-------|
| Avg pixels >10 dBZ/scan | 764 | 869 | +105 (+14%) |
| Avg pixels >20 dBZ/scan | 225 | 334 | +110 (+49%) |
| Max reflectivity | 64.5 dBZ | 65.0 dBZ | +0.5 |
| Mean CC | 0.4033 | 0.3396 | -0.064 |
| Mean ZDR | -1.22 dB | -0.46 dB | +0.77 |

**The baseline night shows ~764 pixels above 10 dBZ per scan and peaks up to 64.5 dBZ.** This means the ~870 pixels and 65 dBZ on event night are NOT a dramatic departure — they represent a ~14% increase above an already-strong background.

### 2. The CC = 0.34 Finding is Overstated

The baseline CC is **0.4033**. The event CC is **0.3396**. Both are in the "chaff/debris" range using standard classification thresholds.

This means the low CC is a **characteristic of this sector's ground clutter environment**, not a unique signature of engagement debris. The difference (0.064) is statistically significant due to the large sample sizes (26,583 vs 33,207 returns), but the practical significance is modest — both values are in the same regime.

### 3. The "Debris Classification" Was Capturing Clutter

| Classification | Baseline | Event |
|---------------|----------|-------|
| Debris/chaff | **54.4%** | 47.9% |
| Military chaff | 5.1% | 5.6% |
| Metallic fragments | 2.1% | 2.1% |
| Ground clutter | 0.7% | 0.5% |

**The baseline has MORE returns classified as "debris" (54.4%) than the event night (47.9%).** The dual-pol classification thresholds (CC < 0.80, ZDR > 4 or ZDR < -2) are not discriminating in this terrain environment — they misclassify ground clutter as debris because the complex terrain produces anomalous dual-pol signatures.

### 4. Max Reflectivity is NOT Anomalous

The baseline peak of 64.5 dBZ (from a scan at 03:09 UTC) is statistically indistinguishable from the event peak of 65.0 dBZ (p = 0.11). The mean max-per-scan is 45.6 ± 9.4 dBZ (baseline) vs 49.2 ± 7.4 dBZ (event).

### 5. Returns Are Concentrated at 15–20 km (Terrain)

| Range Bin | Baseline (rets/scan) | Event (rets/scan) |
|-----------|---------------------|-------------------|
| 15–20 km | 831 | 927 |
| 20–25 km | 61 | 64 |
| 25–30 km | 7 | 8 |
| 30–55 km | <11 combined | <4 combined |

The overwhelming majority of returns are at 15–20 km — precisely the range of the Franklin Mountains. Very few returns exist at longer ranges on either night.

---

## What IS Different on Event Night

Despite the clutter dominance, there are real differences:

### A. South Cluster Shows Genuine Anomalies

| South Cluster Metric | Baseline | Event | Change |
|---------------------|----------|-------|--------|
| Avg returns/scan | 371 | 404 | +9% |
| Avg >20 dBZ/scan | 128 | 215 | **+67%** |
| Mean CC | 0.4030 | 0.3483 | -0.055 |
| Mean ZDR | -1.57 dB | +0.01 dB | **+1.58 dB** |

The south cluster (az 90–110°, away from Franklin Mountains) shows:
- **67% more returns above 20 dBZ** on event night
- A statistically significant CC drop (KS p = 2×10⁻⁹⁷)
- A ZDR shift of +1.58 dB (KS p = 5×10⁻¹⁵²)

The ZDR shift from -1.57 to +0.01 is the most notable change — it indicates a shift toward more symmetric or prolate scatterers, consistent with debris fragments rather than the oblate clutter signature.

### B. Total >20 dBZ Returns Increased 49%

While the >10 dBZ count increased only 14%, the >20 dBZ count increased 49%. This suggests the event night added stronger scatterers on top of the clutter baseline, not just more of the same.

### C. Statistical Tests Confirm Real Differences

| Test | Statistic | p-value | Significant? |
|------|-----------|---------|-------------|
| >10 dBZ count (t-test) | -12.40 | <10⁻⁶ | Yes |
| >20 dBZ count (t-test) | -16.66 | <10⁻⁶ | Yes |
| Overall CC (KS) | 0.157 | 7×10⁻³¹⁷ | Yes |
| South CC (KS) | 0.137 | 2×10⁻⁹⁷ | Yes |
| South ZDR (KS) | 0.171 | 5×10⁻¹⁵² | Yes |
| Max reflectivity (t-test) | -1.62 | 0.112 | **No** |

---

## Revised Assessment

### What We Got Wrong

1. **The CC = 0.34 is NOT a smoking gun.** The baseline CC is 0.40 — the sector always has anomalously low CC due to terrain effects.
2. **The "47.9% debris classification" is meaningless.** The baseline classifies 54.4% as debris using the same thresholds. The classification doesn't work in this terrain.
3. **The 65 dBZ peak is NOT anomalous.** The baseline hit 64.5 dBZ.
4. **The ZDR range of -13 to +20 dB is NOT unique.** The baseline has the identical range.

### What Holds Up

1. **There IS an incremental signal above baseline.** ~49% more >20 dBZ returns, concentrated in the south cluster.
2. **The south cluster ZDR shifted by +1.58 dB.** This is a real change in scatterer characteristics.
3. **The south cluster CC dropped by 0.055.** While both values are low, the event night is measurably lower.
4. **The KHDX cross-confirmation remains valid.** KHDX at 135 km doesn't have the KEPZ terrain problem and independently confirmed activity.

### Revised Confidence Levels

| Claim | Previous | Revised |
|-------|----------|---------|
| "CC proves metallic debris" | ~55% | **25–30%** — baseline shows this sector always has low CC |
| "Dual-pol classification proves debris" | ~60% | **15–20%** — classification doesn't discriminate in this terrain |
| "Something anomalous happened" | ~85% | **75–80%** — incremental signal exists but is modest |
| "KEPZ alone proves engagements" | ~70% | **35–40%** — most returns are ground clutter |
| "KHDX confirms activity" | ~70% | **70%** — unaffected by this analysis |
| "Govt not fully transparent" | ~85% | **80–85%** — news reports still contradict party balloon narrative |

### Bottom Line

The KEPZ Level 2 data shows a **modest incremental signal above persistent ground clutter**, not the dramatic metallic debris cloud previously described. The strongest evidence for the event comes from:

1. The KHDX independent cross-confirmation (unaffected by this finding)
2. The south cluster ZDR shift (+1.58 dB)
3. The 49% increase in >20 dBZ returns
4. The news reporting itself (ABC, CBS, Military Times)

The CC = 0.34 headline number should be retired as evidence. The correct framing is: "KEPZ shows a statistically significant but modest increase in returns above its persistent clutter baseline, with the strongest anomaly in the south cluster ZDR."
