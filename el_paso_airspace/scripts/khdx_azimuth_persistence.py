#!/usr/bin/env python3
"""
KHDX Azimuth Persistence Check
Checks if returns at az 193-195° appear on multiple dates.
If they do → ground clutter. If only on Feb 11 → transient event.

Analyzes 4 dates:
  Feb 8  (3 days before)
  Feb 10 (1 day before / baseline)
  Feb 11 (event night)
  Feb 12 (1 day after)
"""

import numpy as np
import pyart
import os
import glob
import json
import warnings
from collections import defaultdict
warnings.filterwarnings('ignore')

KHDX_LAT, KHDX_LON = 33.0803, -106.1225

# Fort Bliss sector
TARGET_AZ_MIN, TARGET_AZ_MAX = 185.0, 199.0
TARGET_RNG_MIN, TARGET_RNG_MAX = 125000, 148000

# Narrow azimuth bins for persistence check
AZ_BINS = [(185, 187), (187, 189), (189, 191), (191, 193), (193, 195), (195, 197), (197, 199)]

datasets = {
    'Feb 8': ('/home/user/uap-transient-research/el_paso_airspace/khdx_extra_20260208/', 'Level2_KHDX_20260208_*.ar2v'),
    'Feb 10': ('/home/user/uap-transient-research/el_paso_airspace/khdx_baseline_data/', 'Level2_KHDX_20260210_*.ar2v'),
    'Feb 11 (EVENT)': ('/home/user/uap-transient-research/el_paso_airspace/khdx_data/', 'Level2_KHDX_20260211_*.ar2v'),
    'Feb 12': ('/home/user/uap-transient-research/el_paso_airspace/khdx_extra_20260212/', 'Level2_KHDX_20260212_*.ar2v'),
}

def radar_to_latlon(rlat, rlon, range_m, az_deg):
    rkm = range_m / 1000.0
    lat = rlat + (rkm * np.cos(np.radians(az_deg))) / 111.0
    lon = rlon + (rkm * np.sin(np.radians(az_deg))) / (111.0 * np.cos(np.radians(rlat)))
    return lat, lon

def safe_val(val):
    if hasattr(val, 'item'):
        if np.ma.is_masked(val):
            return None
        return float(val.item())
    if val is None or np.isnan(val):
        return None
    return float(val)

print("=" * 80)
print("KHDX AZIMUTH PERSISTENCE CHECK")
print("Do returns at az 193-195° appear on multiple dates?")
print("=" * 80)

all_results = {}

for label, (data_dir, pattern) in datasets.items():
    files = sorted(glob.glob(os.path.join(data_dir, pattern)))
    print(f"\n  Processing {label}: {len(files)} files...")

    date_returns = []
    date_scans = []
    az_bin_counts = defaultdict(list)  # per scan, count in each az bin

    for filepath in files:
        fname = os.path.basename(filepath)
        time_str = fname.split('_')[-1].replace('.ar2v', '')

        try:
            radar = pyart.io.read_nexrad_archive(filepath)
        except Exception:
            continue

        ranges = radar.range['data']
        ref_field = None
        for c in ['reflectivity', 'REF', 'ref']:
            if c in radar.fields:
                ref_field = c; break
        if not ref_field:
            del radar; continue

        sweep_idx = 0
        start = radar.sweep_start_ray_index['data'][sweep_idx]
        end = radar.sweep_end_ray_index['data'][sweep_idx]
        azimuths = radar.azimuth['data'][start:end+1]
        ref_data = radar.fields[ref_field]['data'][start:end+1]

        az_mask = np.array([(TARGET_AZ_MIN <= az <= TARGET_AZ_MAX) for az in azimuths])
        rng_mask = (ranges >= TARGET_RNG_MIN) & (ranges <= TARGET_RNG_MAX)

        sector_az = azimuths[az_mask]
        sector_rng = ranges[rng_mask]
        sector_ref = ref_data[az_mask][:, rng_mask]

        scan_returns = []
        scan_az_bin = defaultdict(int)

        for i in range(sector_ref.shape[0]):
            for j in range(sector_ref.shape[1]):
                rv = safe_val(sector_ref[i, j])
                if rv is None or rv < 5.0:
                    continue
                az = float(sector_az[i])
                rng = float(sector_rng[j])
                lat, lon = radar_to_latlon(KHDX_LAT, KHDX_LON, rng, az)

                ret = {'time': time_str, 'az': az, 'range_km': rng/1000,
                       'lat': lat, 'lon': lon, 'ref_dBZ': rv}
                scan_returns.append(ret)
                date_returns.append(ret)

                # Bin by azimuth
                for az_lo, az_hi in AZ_BINS:
                    if az_lo <= az < az_hi:
                        scan_az_bin[(az_lo, az_hi)] += 1
                        break

        for ab in AZ_BINS:
            az_bin_counts[ab].append(scan_az_bin.get(ab, 0))

        n10 = sum(1 for r in scan_returns if r['ref_dBZ'] > 10)
        max_ref = max((r['ref_dBZ'] for r in scan_returns), default=-999)

        date_scans.append({
            'time': time_str, 'n_returns': len(scan_returns),
            'n_above_10': n10, 'max_ref': max_ref
        })
        del radar

    # Compute stats
    all_rets = len(date_returns)
    n10_total = sum(1 for r in date_returns if r['ref_dBZ'] > 10)
    n20_total = sum(1 for r in date_returns if r['ref_dBZ'] > 20)
    max_ref = max((r['ref_dBZ'] for r in date_returns), default=-999)
    avg_per_scan = all_rets / max(len(date_scans), 1)
    avg_n10 = n10_total / max(len(date_scans), 1)

    all_results[label] = {
        'n_scans': len(date_scans),
        'total_returns': all_rets,
        'avg_per_scan': avg_per_scan,
        'avg_n10': avg_n10,
        'max_ref': max_ref,
        'az_bin_avgs': {f"{lo}-{hi}": float(np.mean(counts)) if counts else 0
                        for (lo, hi), counts in az_bin_counts.items()},
        'az_bin_maxs': {f"{lo}-{hi}": int(np.max(counts)) if counts else 0
                        for (lo, hi), counts in az_bin_counts.items()},
    }

# ============================================================================
# COMPARISON TABLE
# ============================================================================
print("\n\n" + "=" * 80)
print("MULTI-DATE COMPARISON — KHDX Fort Bliss Sector")
print("=" * 80)

print(f"\n{'Date':<20} {'Scans':>6} {'Total':>8} {'Avg/scan':>10} {'Avg>10':>8} {'Max dBZ':>8}")
print(f"{'─'*20} {'─'*6} {'─'*8} {'─'*10} {'─'*8} {'─'*8}")
for label, stats in all_results.items():
    marker = " ***" if 'EVENT' in label else ""
    print(f"{label:<20} {stats['n_scans']:>6} {stats['total_returns']:>8} "
          f"{stats['avg_per_scan']:>10.1f} {stats['avg_n10']:>8.1f} {stats['max_ref']:>8.1f}{marker}")

# ============================================================================
# AZIMUTH BIN PERSISTENCE
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"AZIMUTH BIN ANALYSIS — Average returns per scan by 2° bin")
print(f"{'─'*80}")

print(f"\n{'Az Bin':>12}", end='')
for label in all_results:
    short = label.replace(' (EVENT)', '*')
    print(f" {short:>12}", end='')
print(f" {'Persistent?':>12}")

print(f"{'─'*12}", end='')
for _ in all_results:
    print(f" {'─'*12}", end='')
print(f" {'─'*12}")

persistence_results = {}
for az_lo, az_hi in AZ_BINS:
    key = f"{az_lo}-{az_hi}"
    print(f" {key:>11}°", end='')

    values = []
    for label, stats in all_results.items():
        v = stats['az_bin_avgs'].get(key, 0)
        values.append(v)
        marker = ""
        if 'EVENT' in label and v > 0:
            # Check if this is higher than other dates
            others = [stats2['az_bin_avgs'].get(key, 0) for l2, stats2 in all_results.items()
                      if 'EVENT' not in l2]
            if others and v > max(others) * 1.5:
                marker = " **"
        print(f" {v:>11.1f}{marker}", end='')

    # Determine persistence
    non_event = [v for v, l in zip(values, all_results.keys()) if 'EVENT' not in l]
    event_val = [v for v, l in zip(values, all_results.keys()) if 'EVENT' in l][0]

    if all(v > 0.5 for v in non_event):
        status = "PERSISTENT"
    elif event_val > max(non_event) * 2:
        status = "TRANSIENT"
    elif event_val > max(non_event) * 1.3:
        status = "ELEVATED"
    else:
        status = "normal"

    persistence_results[key] = status
    print(f" {status:>12}")

# ============================================================================
# KEY AZIMUTH 193-195° DEEP DIVE
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"FOCUSED CHECK: az 193-195° (primary KHDX detection azimuth)")
print(f"{'─'*80}")

key = "193-195"
for label, stats in all_results.items():
    avg = stats['az_bin_avgs'].get(key, 0)
    mx = stats['az_bin_maxs'].get(key, 0)
    print(f"  {label:<20} avg={avg:.1f}/scan, max={mx} in single scan")

non_event_193 = [stats['az_bin_avgs'].get(key, 0) for l, stats in all_results.items() if 'EVENT' not in l]
event_193 = [stats['az_bin_avgs'].get(key, 0) for l, stats in all_results.items() if 'EVENT' in l][0]

if all(v > 0.5 for v in non_event_193):
    print(f"\n  VERDICT: Returns at az 193-195° are PERSISTENT across all dates.")
    print(f"  This indicates GROUND CLUTTER, not transient event debris.")
    print(f"  Event night avg ({event_193:.1f}) vs other dates avg ({np.mean(non_event_193):.1f})")
    if event_193 > np.mean(non_event_193) * 1.3:
        print(f"  However, event night is {100*(event_193/np.mean(non_event_193)-1):.0f}% above the multi-date average.")
        print(f"  This suggests: CLUTTER BASELINE + TRANSIENT EXCESS")
elif event_193 > max(non_event_193) * 2:
    print(f"\n  VERDICT: Returns at az 193-195° are PRIMARILY TRANSIENT.")
    print(f"  Event night ({event_193:.1f}) is >2x any other date ({max(non_event_193):.1f}).")
    print(f"  This supports the engagement/debris hypothesis.")
else:
    print(f"\n  VERDICT: Returns at az 193-195° show MIXED evidence.")
    print(f"  Some baseline clutter exists, with possible event-night elevation.")

# ============================================================================
# OVERALL VERDICT
# ============================================================================
print(f"\n\n{'='*80}")
print(f"OVERALL KHDX PERSISTENCE VERDICT")
print(f"{'='*80}")

n_persistent = sum(1 for v in persistence_results.values() if v == 'PERSISTENT')
n_transient = sum(1 for v in persistence_results.values() if v == 'TRANSIENT')
n_elevated = sum(1 for v in persistence_results.values() if v == 'ELEVATED')

print(f"\n  Of {len(AZ_BINS)} azimuth bins:")
print(f"    PERSISTENT (ground clutter): {n_persistent}")
print(f"    TRANSIENT (event-only):      {n_transient}")
print(f"    ELEVATED (above normal):     {n_elevated}")
print(f"    Normal:                      {len(AZ_BINS) - n_persistent - n_transient - n_elevated}")

# Save results
output = {
    'dates_analyzed': list(all_results.keys()),
    'per_date_stats': all_results,
    'azimuth_persistence': persistence_results,
    'az_193_195_verdict': 'persistent_with_excess' if all(v > 0.5 for v in non_event_193) and event_193 > np.mean(non_event_193) * 1.3 else 'persistent' if all(v > 0.5 for v in non_event_193) else 'transient' if event_193 > max(non_event_193) * 2 else 'mixed',
}

output_path = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/khdx_azimuth_persistence.json'
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")
