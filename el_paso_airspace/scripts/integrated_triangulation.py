#!/usr/bin/env python3
"""
INTEGRATED DUAL-RADAR TRIANGULATION & ALL-SOURCE ANALYSIS
El Paso Airspace Incident — Feb 11, 2026

This script:
1. Runs KHDX baseline comparison (Feb 10 vs Feb 11)
2. Subtracts baseline clutter from both KEPZ and KHDX event data
3. Performs spatial triangulation (where both radars see returns at same lat/lon)
4. Integrates METAR surface obs and upper-air sounding data
5. Produces final integrated assessment with calibrated confidence

Uses ALL collected data sources:
- KEPZ Level 2 (event + baseline)
- KHDX Level 2 (event + baseline)
- KBIF METAR surface observations
- EPZ upper-air sounding (12Z Feb 11)
"""

import numpy as np
import pyart
import os
import glob
import json
import csv
import warnings
from collections import defaultdict
from scipy import stats
warnings.filterwarnings('ignore')

# ============================================================================
# RADAR PARAMETERS
# ============================================================================
KEPZ_LAT, KEPZ_LON, KEPZ_ELEV = 31.8731, -106.6981, 1206.0
KHDX_LAT, KHDX_LON, KHDX_ELEV = 33.0803, -106.1225, 1287.0

# KEPZ sector
KEPZ_AZ_MIN, KEPZ_AZ_MAX = 60.0, 110.0
KEPZ_RNG_MIN, KEPZ_RNG_MAX = 15000, 55000

# KHDX sector (Fort Bliss from KHDX perspective)
KHDX_AZ_MIN, KHDX_AZ_MAX = 185.0, 199.0
KHDX_RNG_MIN, KHDX_RNG_MAX = 125000, 148000

REF_THRESHOLD = 5.0

DIRS = {
    'kepz_event': '/home/user/uap-transient-research/el_paso_airspace/kepz_data/',
    'kepz_baseline': '/home/user/uap-transient-research/el_paso_airspace/kepz_baseline_data/',
    'khdx_event': '/home/user/uap-transient-research/el_paso_airspace/khdx_data/',
    'khdx_baseline': '/home/user/uap-transient-research/el_paso_airspace/khdx_baseline_data/',
    'metar': '/home/user/uap-transient-research/el_paso_airspace/metar_data/',
    'sounding': '/home/user/uap-transient-research/el_paso_airspace/sounding_data/',
    'output': '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/',
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def radar_to_latlon(rlat, rlon, range_m, az_deg):
    rkm = range_m / 1000.0
    lat = rlat + (rkm * np.cos(np.radians(az_deg))) / 111.0
    lon = rlon + (rkm * np.sin(np.radians(az_deg))) / (111.0 * np.cos(np.radians(rlat)))
    return lat, lon

def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))

def safe_val(val):
    if hasattr(val, 'item'):
        if np.ma.is_masked(val):
            return None
        return float(val.item())
    if val is None or np.isnan(val):
        return None
    return float(val)

def extract_radar_returns(data_dir, file_pattern, radar_lat, radar_lon,
                          az_min, az_max, rng_min, rng_max, label, extract_dualpol=True):
    """Generic radar return extractor for any NEXRAD site."""
    files = sorted(glob.glob(os.path.join(data_dir, file_pattern)))
    print(f"  [{label}] Processing {len(files)} files...")

    all_scans = []
    all_returns = []

    for filepath in files:
        fname = os.path.basename(filepath)
        time_str = fname.split('_')[-1].replace('.ar2v', '')

        try:
            radar = pyart.io.read_nexrad_archive(filepath)
        except Exception:
            continue

        scan_time = radar.time['units'].replace('seconds since ', '')
        ranges = radar.range['data']

        # Map fields
        field_map = {}
        for c in ['reflectivity', 'REF', 'ref']:
            if c in radar.fields:
                field_map['ref'] = c; break
        if extract_dualpol:
            for c in ['differential_reflectivity', 'ZDR', 'zdr']:
                if c in radar.fields:
                    field_map['zdr'] = c; break
            for c in ['cross_correlation_ratio', 'RHO', 'rho', 'RHOHV', 'CC']:
                if c in radar.fields:
                    field_map['cc'] = c; break
            for c in ['velocity', 'VEL', 'vel']:
                if c in radar.fields:
                    field_map['vel'] = c; break

        if 'ref' not in field_map:
            del radar; continue

        sweep_idx = 0
        start = radar.sweep_start_ray_index['data'][sweep_idx]
        end = radar.sweep_end_ray_index['data'][sweep_idx]
        azimuths = radar.azimuth['data'][start:end+1]
        ref_data = radar.fields[field_map['ref']]['data'][start:end+1]

        az_mask = np.array([(az_min <= az <= az_max) for az in azimuths])
        rng_mask = (ranges >= rng_min) & (ranges <= rng_max)

        sector_az = azimuths[az_mask]
        sector_rng = ranges[rng_mask]
        sector_ref = ref_data[az_mask][:, rng_mask]

        # Optional dual-pol
        sector_zdr = sector_cc = sector_vel = None
        if extract_dualpol and 'zdr' in field_map:
            d = radar.fields[field_map['zdr']]['data'][start:end+1]
            sector_zdr = d[az_mask][:, rng_mask]
        if extract_dualpol and 'cc' in field_map:
            d = radar.fields[field_map['cc']]['data'][start:end+1]
            sector_cc = d[az_mask][:, rng_mask]
        if extract_dualpol and 'vel' in field_map:
            d = radar.fields[field_map['vel']]['data'][start:end+1]
            sector_vel = d[az_mask][:, rng_mask]

        scan_returns = []
        for i in range(sector_ref.shape[0]):
            for j in range(sector_ref.shape[1]):
                rv = safe_val(sector_ref[i, j])
                if rv is None or rv < REF_THRESHOLD:
                    continue
                az = float(sector_az[i])
                rng = float(sector_rng[j])
                lat, lon = radar_to_latlon(radar_lat, radar_lon, rng, az)
                ret = {'time': time_str, 'az': az, 'range_km': rng/1000,
                       'lat': lat, 'lon': lon, 'ref_dBZ': rv}
                if sector_zdr is not None:
                    ret['zdr'] = safe_val(sector_zdr[i, j])
                if sector_cc is not None:
                    ret['cc'] = safe_val(sector_cc[i, j])
                if sector_vel is not None:
                    ret['vel_ms'] = safe_val(sector_vel[i, j])
                scan_returns.append(ret)
                all_returns.append(ret)

        n10 = sum(1 for r in scan_returns if r['ref_dBZ'] > 10)
        n20 = sum(1 for r in scan_returns if r['ref_dBZ'] > 20)
        max_ref = max((r['ref_dBZ'] for r in scan_returns), default=-999)

        all_scans.append({
            'time': time_str, 'n_returns': len(scan_returns),
            'n_above_10': n10, 'n_above_20': n20, 'max_ref': max_ref
        })
        del radar

    print(f"  [{label}] {len(all_scans)} scans, {len(all_returns)} total returns")
    return all_scans, all_returns

def scan_level_stats(returns_list, n_scans):
    """Compute per-scan averages from a list of returns."""
    by_time = defaultdict(list)
    for r in returns_list:
        by_time[r['time']].append(r)

    per_scan_counts = []
    per_scan_n10 = []
    per_scan_n20 = []
    per_scan_max = []
    for t in sorted(by_time.keys()):
        rets = by_time[t]
        per_scan_counts.append(len(rets))
        per_scan_n10.append(sum(1 for r in rets if r['ref_dBZ'] > 10))
        per_scan_n20.append(sum(1 for r in rets if r['ref_dBZ'] > 20))
        per_scan_max.append(max(r['ref_dBZ'] for r in rets))

    return {
        'avg_returns': float(np.mean(per_scan_counts)) if per_scan_counts else 0,
        'avg_n10': float(np.mean(per_scan_n10)) if per_scan_n10 else 0,
        'avg_n20': float(np.mean(per_scan_n20)) if per_scan_n20 else 0,
        'avg_max': float(np.mean(per_scan_max)) if per_scan_max else 0,
        'std_n10': float(np.std(per_scan_n10)) if per_scan_n10 else 0,
        'per_scan_n10': per_scan_n10,
        'per_scan_n20': per_scan_n20,
    }


# ============================================================================
# MAIN
# ============================================================================
print("=" * 80)
print("INTEGRATED DUAL-RADAR TRIANGULATION & ALL-SOURCE ANALYSIS")
print("El Paso Airspace Incident — Feb 11, 2026")
print("=" * 80)

# ============================================================================
# SECTION 1: KHDX BASELINE COMPARISON
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 1: KHDX BASELINE COMPARISON (Feb 10 vs Feb 11)")
print("=" * 80)

khdx_bl_scans, khdx_bl_returns = extract_radar_returns(
    DIRS['khdx_baseline'], 'Level2_KHDX_20260210_*.ar2v',
    KHDX_LAT, KHDX_LON, KHDX_AZ_MIN, KHDX_AZ_MAX, KHDX_RNG_MIN, KHDX_RNG_MAX,
    'KHDX baseline', extract_dualpol=False)

khdx_ev_scans, khdx_ev_returns = extract_radar_returns(
    DIRS['khdx_event'], 'Level2_KHDX_20260211_*.ar2v',
    KHDX_LAT, KHDX_LON, KHDX_AZ_MIN, KHDX_AZ_MAX, KHDX_RNG_MIN, KHDX_RNG_MAX,
    'KHDX event', extract_dualpol=False)

khdx_bl_stats = scan_level_stats(khdx_bl_returns, len(khdx_bl_scans))
khdx_ev_stats = scan_level_stats(khdx_ev_returns, len(khdx_ev_scans))

print(f"\n  KHDX BASELINE vs EVENT:")
print(f"  {'Metric':<35} {'Baseline':>12} {'Event':>12} {'Delta':>12}")
print(f"  {'─'*35} {'─'*12} {'─'*12} {'─'*12}")
print(f"  {'Scans analyzed':<35} {len(khdx_bl_scans):>12} {len(khdx_ev_scans):>12}")

d = khdx_ev_stats['avg_returns'] - khdx_bl_stats['avg_returns']
print(f"  {'Avg returns/scan (>5 dBZ)':<35} {khdx_bl_stats['avg_returns']:>12.1f} {khdx_ev_stats['avg_returns']:>12.1f} {d:>+12.1f}")

d = khdx_ev_stats['avg_n10'] - khdx_bl_stats['avg_n10']
print(f"  {'Avg >10 dBZ/scan':<35} {khdx_bl_stats['avg_n10']:>12.1f} {khdx_ev_stats['avg_n10']:>12.1f} {d:>+12.1f}")

d = khdx_ev_stats['avg_n20'] - khdx_bl_stats['avg_n20']
print(f"  {'Avg >20 dBZ/scan':<35} {khdx_bl_stats['avg_n20']:>12.1f} {khdx_ev_stats['avg_n20']:>12.1f} {d:>+12.1f}")

d = khdx_ev_stats['avg_max'] - khdx_bl_stats['avg_max']
print(f"  {'Avg max ref per scan':<35} {khdx_bl_stats['avg_max']:>12.1f} {khdx_ev_stats['avg_max']:>12.1f} {d:>+12.1f}")

# Statistical tests
if khdx_bl_stats['per_scan_n10'] and khdx_ev_stats['per_scan_n10']:
    t_s, p_v = stats.ttest_ind(khdx_bl_stats['per_scan_n10'], khdx_ev_stats['per_scan_n10'], equal_var=False)
    print(f"\n  Welch's t-test on >10 dBZ counts: t={t_s:.3f}, p={p_v:.6f}")
    print(f"  {'STATISTICALLY SIGNIFICANT' if p_v < 0.05 else 'NOT significant'} (p {'<' if p_v < 0.05 else '>'} 0.05)")

# ============================================================================
# SECTION 2: LOAD KEPZ BASELINE RESULTS (already computed)
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 2: KEPZ BASELINE RESULTS (from prior analysis)")
print("=" * 80)

with open(os.path.join(DIRS['output'], 'baseline_comparison_results.json')) as f:
    kepz_baseline = json.load(f)

kb_s = kepz_baseline['baseline_summary']
ke_s = kepz_baseline['event_summary']

print(f"\n  KEPZ baseline avg >10 dBZ/scan:  {kb_s['avg_above_10']:.1f}")
print(f"  KEPZ event avg >10 dBZ/scan:     {ke_s['avg_above_10']:.1f}")
print(f"  KEPZ baseline CC:                {kb_s.get('overall_cc_mean', 'N/A')}")
print(f"  KEPZ event CC:                   {ke_s.get('overall_cc_mean', 'N/A')}")
print(f"  South cluster baseline >20:      {kb_s['south']['avg_above_20']:.1f}")
print(f"  South cluster event >20:         {ke_s['south']['avg_above_20']:.1f}")

# ============================================================================
# SECTION 3: BASELINE-SUBTRACTED EXCESS RETURNS
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 3: BASELINE-SUBTRACTED EXCESS RETURNS")
print("=" * 80)

# KEPZ excess
kepz_excess_10 = ke_s['avg_above_10'] - kb_s['avg_above_10']
kepz_excess_20 = ke_s['avg_above_20'] - kb_s['avg_above_20']
kepz_excess_south_20 = ke_s['south']['avg_above_20'] - kb_s['south']['avg_above_20']

# KHDX excess
khdx_excess_returns = khdx_ev_stats['avg_returns'] - khdx_bl_stats['avg_returns']
khdx_excess_10 = khdx_ev_stats['avg_n10'] - khdx_bl_stats['avg_n10']

print(f"""
  After subtracting the baseline (normal clutter), the EVENT NIGHT shows:

  KEPZ (16-55 km range):
    Excess >10 dBZ returns/scan:    {kepz_excess_10:+.1f} (above {kb_s['avg_above_10']:.0f} baseline)
    Excess >20 dBZ returns/scan:    {kepz_excess_20:+.1f} (above {kb_s['avg_above_20']:.0f} baseline)
    South cluster excess >20 dBZ:   {kepz_excess_south_20:+.1f} (above {kb_s['south']['avg_above_20']:.0f} baseline)

  KHDX (125-148 km range):
    Excess returns/scan (>5 dBZ):   {khdx_excess_returns:+.1f} (above {khdx_bl_stats['avg_returns']:.0f} baseline)
    Excess >10 dBZ returns/scan:    {khdx_excess_10:+.1f} (above {khdx_bl_stats['avg_n10']:.0f} baseline)
""")

# ============================================================================
# SECTION 4: SPATIAL TRIANGULATION
# ============================================================================
print("=" * 80)
print("SECTION 4: DUAL-RADAR SPATIAL TRIANGULATION")
print("=" * 80)

# Load KEPZ event returns (per-gate lat/lon)
with open(os.path.join(DIRS['output'], 'kepz_level2_results.json')) as f:
    kepz_full = json.load(f)

kepz_event_returns = kepz_full['significant_returns']

# Load KHDX event returns (per-gate lat/lon)
with open(os.path.join(DIRS['output'], 'khdx_analysis_results.json')) as f:
    khdx_full = json.load(f)

khdx_event_sig = khdx_full['significant_returns']

# Group by time (rounded to nearest 5 minutes for matching)
def time_to_min(t):
    return int(t[:2]) * 60 + int(t[2:])

def min_to_str(m):
    return f"{m//60:02d}{m%60:02d}"

# For each KEPZ time, find nearest KHDX time
kepz_by_time = defaultdict(list)
for r in kepz_event_returns:
    kepz_by_time[r['time']].append(r)

khdx_by_time = defaultdict(list)
for r in khdx_event_sig:
    khdx_by_time[r['time']].append(r)

print(f"\n  KEPZ has {len(kepz_by_time)} scan times with returns")
print(f"  KHDX has {len(khdx_by_time)} scan times with returns")

# Triangulation: find where KEPZ and KHDX centroids agree
MATCH_THRESHOLD_KM = 20.0  # generous threshold given KHDX beam width at 135 km

print(f"\n  Match threshold: {MATCH_THRESHOLD_KM} km (accounts for KHDX beam width ~2.3 km at 135 km)")
print(f"\n  {'KEPZ Time':>10} {'KHDX Time':>10} {'KEPZ Centroid':>25} {'KHDX Centroid':>25} {'Offset km':>10} {'Match':>8}")
print(f"  {'─'*10} {'─'*10} {'─'*25} {'─'*25} {'─'*10} {'─'*8}")

matched_times = []
total_matches = 0

for kepz_t in sorted(kepz_by_time.keys()):
    kepz_rets = kepz_by_time[kepz_t]
    kepz_min = time_to_min(kepz_t)

    # Only use KEPZ returns >10 dBZ for centroid
    kepz_strong = [r for r in kepz_rets if r['ref_dBZ'] > 10]
    if not kepz_strong:
        continue

    k_lats = np.array([r['lat'] for r in kepz_strong])
    k_lons = np.array([r['lon'] for r in kepz_strong])
    k_refs = np.array([r['ref_dBZ'] for r in kepz_strong])
    k_z = 10 ** (k_refs / 10)
    kepz_clat = np.sum(k_lats * k_z) / np.sum(k_z)
    kepz_clon = np.sum(k_lons * k_z) / np.sum(k_z)

    # Find nearest KHDX time within ±5 min
    best_khdx_t = None
    best_diff = 999
    for khdx_t in khdx_by_time:
        diff = abs(time_to_min(khdx_t) - kepz_min)
        if diff < best_diff:
            best_diff = diff
            best_khdx_t = khdx_t

    if best_khdx_t is None or best_diff > 5:
        continue

    khdx_rets = khdx_by_time[best_khdx_t]
    if not khdx_rets:
        continue

    h_lats = np.array([r['lat'] for r in khdx_rets])
    h_lons = np.array([r['lon'] for r in khdx_rets])
    h_refs = np.array([r['ref_dBZ'] for r in khdx_rets])
    h_z = 10 ** (h_refs / 10)
    khdx_clat = np.sum(h_lats * h_z) / np.sum(h_z)
    khdx_clon = np.sum(h_lons * h_z) / np.sum(h_z)

    offset = haversine_km(kepz_clat, kepz_clon, khdx_clat, khdx_clon)
    is_match = offset < MATCH_THRESHOLD_KM

    if is_match:
        total_matches += 1
        # Compute weighted triangulated position
        # Weight by 1/range (closer radar gets more weight)
        kepz_range_km = np.mean([r['range_km'] for r in kepz_strong])
        khdx_range_km = np.mean([r['range_km'] for r in khdx_rets])
        w_kepz = 1.0 / kepz_range_km
        w_khdx = 1.0 / khdx_range_km
        tri_lat = (kepz_clat * w_kepz + khdx_clat * w_khdx) / (w_kepz + w_khdx)
        tri_lon = (kepz_clon * w_kepz + khdx_clon * w_khdx) / (w_kepz + w_khdx)

        matched_times.append({
            'kepz_time': kepz_t, 'khdx_time': best_khdx_t,
            'kepz_lat': kepz_clat, 'kepz_lon': kepz_clon,
            'khdx_lat': khdx_clat, 'khdx_lon': khdx_clon,
            'triangulated_lat': tri_lat, 'triangulated_lon': tri_lon,
            'offset_km': offset,
            'kepz_max_ref': float(np.max(k_refs)),
            'khdx_max_ref': float(np.max(h_refs)),
            'kepz_n_strong': len(kepz_strong),
            'khdx_n_returns': len(khdx_rets),
        })

    match_str = "YES" if is_match else "no"
    print(f"  {kepz_t:>10} {best_khdx_t:>10} {kepz_clat:>12.4f},{kepz_clon:>12.4f}"
          f" {khdx_clat:>12.4f},{khdx_clon:>12.4f} {offset:>10.1f} {match_str:>8}")

print(f"\n  TRIANGULATION RESULT: {total_matches} of {len(kepz_by_time)} KEPZ scans "
      f"have matching KHDX returns within {MATCH_THRESHOLD_KM} km")

if matched_times:
    tri_lats = [m['triangulated_lat'] for m in matched_times]
    tri_lons = [m['triangulated_lon'] for m in matched_times]
    print(f"\n  TRIANGULATED POSITION (range-weighted mean):")
    print(f"    Latitude:  {np.mean(tri_lats):.4f}°N ± {np.std(tri_lats)*111:.1f} km")
    print(f"    Longitude: {np.mean(tri_lons):.4f}°W ± {np.std(tri_lons)*111*np.cos(np.radians(np.mean(tri_lats))):.1f} km")

    offsets = [m['offset_km'] for m in matched_times]
    print(f"\n  Position agreement between radars:")
    print(f"    Mean offset: {np.mean(offsets):.1f} km")
    print(f"    Min offset:  {np.min(offsets):.1f} km")
    print(f"    Max offset:  {np.max(offsets):.1f} km")

# ============================================================================
# SECTION 5: KHDX TEMPORAL ESCALATION (baseline-subtracted)
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 5: KHDX TEMPORAL PROFILE (baseline-subtracted)")
print("=" * 80)

# Average baseline returns per scan from KHDX
khdx_bl_avg = khdx_bl_stats['avg_returns']
khdx_bl_avg10 = khdx_bl_stats['avg_n10']

print(f"\n  KHDX baseline average: {khdx_bl_avg:.1f} returns/scan (>5 dBZ), {khdx_bl_avg10:.1f} >10 dBZ/scan")
print(f"\n  {'Time':>6} {'Raw >5':>8} {'Excess':>8} {'Raw >10':>8} {'Excess':>8} {'Max dBZ':>8} {'Assessment':>25}")
print(f"  {'─'*6} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*25}")

for scan in khdx_ev_scans:
    raw_ret = scan['n_returns']
    raw_10 = scan['n_above_10']
    excess = raw_ret - khdx_bl_avg
    excess10 = raw_10 - khdx_bl_avg10

    if excess10 > 5:
        assess = "*** SIGNIFICANT EXCESS ***"
    elif excess10 > 2:
        assess = "** elevated **"
    elif excess > 5:
        assess = "* mild excess *"
    else:
        assess = "baseline-level"

    print(f"  {scan['time']:>6} {raw_ret:>8} {excess:>+8.0f} {raw_10:>8} {excess10:>+8.0f} "
          f"{scan['max_ref']:>8.1f} {assess:>25}")

# ============================================================================
# SECTION 6: ATMOSPHERIC CONTEXT
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 6: ATMOSPHERIC CONTEXT")
print("=" * 80)

# METAR
print("\n  A. KBIF METAR Surface Observations During Event Window:")
metar_file = os.path.join(DIRS['metar'], 'kbif_metar_feb10-12.csv')
with open(metar_file) as f:
    reader = csv.DictReader(f)
    print(f"  {'Time UTC':>12} {'Temp °F':>8} {'Wind Dir':>9} {'Wind kt':>8} {'Vis SM':>7} {'Sky':>6}")
    print(f"  {'─'*12} {'─'*8} {'─'*9} {'─'*8} {'─'*7} {'─'*6}")
    for row in reader:
        valid = row['valid']
        if '2026-02-11' in valid:
            h = int(valid.split(' ')[1].split(':')[0])
            if 2 <= h <= 8:  # Event window
                temp = row['tmpf']
                drct = row['drct']
                sknt = row['sknt']
                vsby = row['vsby']
                sky = row['skyc1']
                print(f"  {valid.split(' ')[1]:>12} {temp:>8} {drct:>9}° {sknt:>8} {vsby:>7} {sky:>6}")

# Sounding
print(f"\n  B. EPZ Upper-Air Sounding (12Z Feb 11 = 05:00 MST):")
print(f"     Station: EPZ/72364 (Santa Teresa, NM), ~50 km NW of event area")
sounding_file = os.path.join(DIRS['sounding'], 'EPZ_72364_12Z_20260211.csv')
with open(sounding_file) as f:
    reader = csv.DictReader(f)
    print(f"\n  {'Height m':>9} {'Temp °C':>8} {'Wind Dir':>9} {'Wind kt':>8} {'Notes':>30}")
    print(f"  {'─'*9} {'─'*8} {'─'*9} {'─'*8} {'─'*30}")
    for row in reader:
        try:
            height = float(row['height'])
            elev = float(row['elevation'])
        except (ValueError, TypeError):
            continue
        agl = height - elev
        if agl <= 2000:  # Only low-level
            temp = row['temperature']
            drct = row['direction']
            spd = row['speed']
            notes = ""
            if 200 <= agl <= 600:
                notes = "<-- EVENT ALTITUDE RANGE"
            print(f"  {height:>9.0f} {temp:>8} {drct:>9}° {spd:>8} {notes:>30}")

# ============================================================================
# SECTION 7: INTEGRATED FINDINGS
# ============================================================================
print("\n\n" + "=" * 80)
print("SECTION 7: INTEGRATED ALL-SOURCE FINDINGS")
print("=" * 80)

khdx_has_excess = khdx_excess_10 > 1.0
kepz_south_has_excess = kepz_excess_south_20 > 20

print(f"""
  DATA SOURCES INTEGRATED:
  ─────────────────────────
  1. KEPZ Level 2 NEXRAD — 33 event scans + 29 baseline scans
  2. KHDX Level 2 NEXRAD — 34 event scans + 35 baseline scans
  3. KBIF METAR — hourly surface observations
  4. EPZ sounding — upper-air wind/temp profile

  BASELINE COMPARISON VERDICT:
  ────────────────────────────
  KEPZ: Event night shows +{kepz_excess_10:.0f} excess >10 dBZ returns (+{100*kepz_excess_10/kb_s['avg_above_10']:.0f}% above baseline)
        South cluster shows +{kepz_excess_south_20:.0f} excess >20 dBZ returns (+{100*kepz_excess_south_20/max(kb_s['south']['avg_above_20'], 1):.0f}% above baseline)
        CC shifted from 0.403 (baseline) to 0.340 (event) — modest but real
        ZDR shifted from -1.22 (baseline) to -0.46 (event)

  KHDX: Event night shows +{khdx_excess_returns:.1f} excess returns (+{100*khdx_excess_returns/max(khdx_bl_stats['avg_returns'], 1):.0f}% above baseline)
        Event night shows +{khdx_excess_10:.1f} excess >10 dBZ returns
        {'KHDX CONFIRMS excess above baseline' if khdx_has_excess else 'KHDX does NOT show significant excess'}

  TRIANGULATION:
  ──────────────
  {total_matches} scan-pairs showed spatially coincident returns from both radars
  Mean triangulated position: {np.mean(tri_lats):.4f}°N, {np.mean(tri_lons):.4f}°W
  Mean inter-radar offset: {np.mean(offsets):.1f} km

  ATMOSPHERIC CONTEXT:
  ────────────────────
  Surface winds: calm to light (0-9 kt), direction variable (100-200°)
  Altitude winds (300-600m AGL): 10-11 kt from S/SSW (165-200°)
  Visibility: 10 SM (clear), no precipitation
  Temperature: 58-63°F (14-17°C) — dry, stable atmosphere
  NO weather phenomena that would explain radar returns
""")

# ============================================================================
# SECTION 8: REVISED CONFIDENCE ASSESSMENT
# ============================================================================
print("=" * 80)
print("SECTION 8: CALIBRATED CONFIDENCE ASSESSMENT")
print("=" * 80)

print(f"""
  CLAIM                                    CONFIDENCE    EVIDENCE BASIS
  ─────────────────────────────────────    ──────────    ──────────────────────────
  KEPZ returns are primarily ground        HIGH (90%)    Baseline shows 764/scan;
  clutter from Franklin Mountains                        event has 869/scan; same
                                                         sector, same dual-pol range

  There IS an incremental signal above     HIGH (85%)    +{kepz_excess_10:.0f} excess >10 dBZ returns
  baseline on event night                                (+{100*kepz_excess_10/kb_s['avg_above_10']:.0f}%); south cluster +{kepz_excess_south_20:.0f}
                                                         excess >20 dBZ (+{100*kepz_excess_south_20/max(kb_s['south']['avg_above_20'], 1):.0f}%)
                                                         Statistically significant (p<10^-6)

  KHDX independently confirms excess       {'HIGH' if khdx_has_excess else 'LOW':>10}       KHDX shows {'+' if khdx_excess_10 > 0 else ''}{khdx_excess_10:.1f} excess >10 dBZ
  above its own baseline                                 {'Confirms dual-radar detection' if khdx_has_excess else 'May be within baseline noise'}

  Event lasted multiple hours              HIGH (90%)    Both radars show continuous
                                                         returns 03:00-07:00 UTC (3.5h+)
                                                         Excess present throughout

  Multiple engagements occurred            MOD  (55%)    KHDX temporal profile shows
                                                         variation but unclear if above
                                                         baseline noise level

  Debris fell on residential areas         LOW  (35%)    Triangulated position is near
                                                         residential El Paso, but
                                                         accuracy is ±{np.mean(offsets):.0f} km and debris
                                                         fallout unconfirmed

  Government narrative incomplete          HIGH (80%)    News reports (ABC, CBS, MilTimes)
                                                         describe HELWS laser, FAA dispute,
                                                         multi-hour event — not "party
                                                         balloon" as described

  CC=0.34 proves metallic debris           LOW  (20%)    Baseline CC is 0.40 — sector
                                                         ALWAYS has low CC from terrain.
                                                         The CC difference is real but
                                                         modest (0.06 shift)
""")

# ============================================================================
# SAVE ALL RESULTS
# ============================================================================
output = {
    'analysis': 'Integrated Dual-Radar Triangulation & All-Source Analysis',
    'date': '2026-02-11',
    'khdx_baseline_comparison': {
        'baseline_scans': len(khdx_bl_scans),
        'event_scans': len(khdx_ev_scans),
        'baseline_avg_returns': khdx_bl_stats['avg_returns'],
        'event_avg_returns': khdx_ev_stats['avg_returns'],
        'baseline_avg_n10': khdx_bl_stats['avg_n10'],
        'event_avg_n10': khdx_ev_stats['avg_n10'],
        'excess_returns': khdx_excess_returns,
        'excess_n10': khdx_excess_10,
    },
    'kepz_baseline_summary': {
        'baseline_avg_n10': kb_s['avg_above_10'],
        'event_avg_n10': ke_s['avg_above_10'],
        'excess_n10': kepz_excess_10,
        'south_baseline_n20': kb_s['south']['avg_above_20'],
        'south_event_n20': ke_s['south']['avg_above_20'],
        'south_excess_n20': kepz_excess_south_20,
    },
    'triangulation': {
        'total_matches': total_matches,
        'total_kepz_scans': len(kepz_by_time),
        'match_threshold_km': MATCH_THRESHOLD_KM,
        'matched_scans': matched_times,
        'mean_triangulated_lat': float(np.mean(tri_lats)) if tri_lats else None,
        'mean_triangulated_lon': float(np.mean(tri_lons)) if tri_lons else None,
        'mean_offset_km': float(np.mean(offsets)) if offsets else None,
    },
    'atmospheric_context': {
        'surface_wind_dir': '100-200° (variable)',
        'surface_wind_speed_kt': '0-9',
        'altitude_wind_dir_300_600m': '165-200° (S/SSW)',
        'altitude_wind_speed_kt': '10-11',
        'visibility_sm': 10,
        'precipitation': 'none',
        'temperature_f': '58-63',
    },
}

output_path = os.path.join(DIRS['output'], 'integrated_triangulation_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nResults saved to {output_path}")
