#!/usr/bin/env python3
"""
KEPZ BASELINE COMPARISON
Compares Feb 10 (baseline/control night) vs Feb 11 (event night)
using IDENTICAL sector parameters and analysis methodology.

This is the critical test: if the baseline looks the same as the event night,
then the KEPZ returns are ground clutter, not engagement debris.
"""

import numpy as np
import pyart
import os
import glob
import json
import warnings
from collections import defaultdict
warnings.filterwarnings('ignore')

# ============================================================================
# PARAMETERS — IDENTICAL TO EVENT ANALYSIS
# ============================================================================
KEPZ_LAT = 31.8731
KEPZ_LON = -106.6981
KEPZ_ELEV = 1206.0

TARGET_AZ_MIN = 60.0
TARGET_AZ_MAX = 110.0
TARGET_RANGE_MIN = 15000   # 15 km
TARGET_RANGE_MAX = 55000   # 55 km
REF_THRESHOLD = 5.0        # dBZ

# Sub-sectors for cluster analysis
NORTH_AZ_MIN = 60.0
NORTH_AZ_MAX = 90.0
SOUTH_AZ_MIN = 90.0
SOUTH_AZ_MAX = 110.0

BASELINE_DIR = '/home/user/uap-transient-research/el_paso_airspace/kepz_baseline_data/'
EVENT_DIR = '/home/user/uap-transient-research/el_paso_airspace/kepz_data/'
OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def radar_to_latlon(radar_lat, radar_lon, range_m, az_deg):
    range_km = range_m / 1000.0
    lat = radar_lat + (range_km * np.cos(np.radians(az_deg))) / 111.0
    lon = radar_lon + (range_km * np.sin(np.radians(az_deg))) / (111.0 * np.cos(np.radians(radar_lat)))
    return lat, lon

def get_sweep_data(radar, sweep_idx, field_name):
    if field_name not in radar.fields:
        return None
    start = radar.sweep_start_ray_index['data'][sweep_idx]
    end = radar.sweep_end_ray_index['data'][sweep_idx]
    return radar.fields[field_name]['data'][start:end+1]

def get_sweep_azimuths(radar, sweep_idx):
    start = radar.sweep_start_ray_index['data'][sweep_idx]
    end = radar.sweep_end_ray_index['data'][sweep_idx]
    return radar.azimuth['data'][start:end+1]

def safe_val(val):
    if hasattr(val, 'item'):
        if np.ma.is_masked(val):
            return None
        return float(val.item())
    if np.isnan(val):
        return None
    return float(val)

def analyze_dataset(data_dir, file_pattern, label):
    """Run the full KEPZ sector analysis on a set of Level 2 files."""
    files = sorted(glob.glob(os.path.join(data_dir, file_pattern)))
    print(f"\n{'='*80}")
    print(f"  ANALYZING: {label}")
    print(f"  Directory: {data_dir}")
    print(f"  Files: {len(files)}")
    print(f"  Sector: az {TARGET_AZ_MIN}-{TARGET_AZ_MAX}°, range {TARGET_RANGE_MIN/1000:.0f}-{TARGET_RANGE_MAX/1000:.0f} km")
    print(f"{'='*80}")

    all_scan_results = []
    all_returns = []

    for filepath in files:
        fname = os.path.basename(filepath)
        time_str = fname.split('_')[-1].replace('.ar2v', '')

        try:
            radar = pyart.io.read_nexrad_archive(filepath)
        except Exception as e:
            print(f"  {fname}: ERROR - {e}")
            continue

        scan_time = radar.time['units'].replace('seconds since ', '')
        ranges = radar.range['data']

        # Map field names
        field_map = {}
        for candidate in ['reflectivity', 'REF', 'ref']:
            if candidate in radar.fields:
                field_map['ref'] = candidate
                break
        for candidate in ['differential_reflectivity', 'ZDR', 'zdr']:
            if candidate in radar.fields:
                field_map['zdr'] = candidate
                break
        for candidate in ['cross_correlation_ratio', 'RHO', 'rho', 'RHOHV', 'CC']:
            if candidate in radar.fields:
                field_map['cc'] = candidate
                break
        for candidate in ['velocity', 'VEL', 'vel']:
            if candidate in radar.fields:
                field_map['vel'] = candidate
                break
        for candidate in ['spectrum_width', 'SW', 'sw', 'WIDTH']:
            if candidate in radar.fields:
                field_map['sw'] = candidate
                break

        if 'ref' not in field_map:
            del radar
            continue

        # Process lowest elevation sweep
        sweep_idx = 0
        elev_angle = radar.fixed_angle['data'][sweep_idx]
        azimuths = get_sweep_azimuths(radar, sweep_idx)

        az_mask = np.array([(TARGET_AZ_MIN <= az <= TARGET_AZ_MAX) for az in azimuths])
        rng_mask = (ranges >= TARGET_RANGE_MIN) & (ranges <= TARGET_RANGE_MAX)

        sector_az = azimuths[az_mask]
        sector_rng = ranges[rng_mask]

        ref_data = get_sweep_data(radar, sweep_idx, field_map['ref'])
        sector_ref = ref_data[az_mask][:, rng_mask]

        # Get dual-pol fields
        sector_zdr = None
        sector_cc = None
        sector_vel = None
        sector_sw = None

        if 'zdr' in field_map:
            zdr_data = get_sweep_data(radar, sweep_idx, field_map['zdr'])
            if zdr_data is not None:
                sector_zdr = zdr_data[az_mask][:, rng_mask]
        if 'cc' in field_map:
            cc_data = get_sweep_data(radar, sweep_idx, field_map['cc'])
            if cc_data is not None:
                sector_cc = cc_data[az_mask][:, rng_mask]
        if 'vel' in field_map:
            vel_data = get_sweep_data(radar, sweep_idx, field_map['vel'])
            if vel_data is not None:
                sector_vel = vel_data[az_mask][:, rng_mask]
        if 'sw' in field_map:
            sw_data = get_sweep_data(radar, sweep_idx, field_map['sw'])
            if sw_data is not None:
                sector_sw = sw_data[az_mask][:, rng_mask]

        # Extract returns
        scan_returns = []
        for i in range(sector_ref.shape[0]):
            for j in range(sector_ref.shape[1]):
                ref_val = safe_val(sector_ref[i, j])
                if ref_val is None or ref_val < REF_THRESHOLD:
                    continue

                az = float(sector_az[i])
                rng = float(sector_rng[j])
                lat, lon = radar_to_latlon(KEPZ_LAT, KEPZ_LON, rng, az)

                ret = {
                    'time': time_str,
                    'az': az,
                    'range_m': rng,
                    'range_km': rng / 1000.0,
                    'lat': lat,
                    'lon': lon,
                    'ref_dBZ': ref_val,
                }

                if sector_zdr is not None:
                    ret['zdr'] = safe_val(sector_zdr[i, j])
                if sector_cc is not None:
                    ret['cc'] = safe_val(sector_cc[i, j])
                if sector_vel is not None:
                    ret['vel_ms'] = safe_val(sector_vel[i, j])
                if sector_sw is not None:
                    ret['spectrum_width'] = safe_val(sector_sw[i, j])

                # Classify into north/south cluster
                if az < SOUTH_AZ_MIN:
                    ret['cluster'] = 'north'
                else:
                    ret['cluster'] = 'south'

                scan_returns.append(ret)
                all_returns.append(ret)

        # Scan-level stats
        valid_ref = []
        for i in range(sector_ref.shape[0]):
            for j in range(sector_ref.shape[1]):
                v = safe_val(sector_ref[i, j])
                if v is not None:
                    valid_ref.append(v)

        above_5 = sum(1 for v in valid_ref if v > 5)
        above_10 = sum(1 for v in valid_ref if v > 10)
        above_15 = sum(1 for v in valid_ref if v > 15)
        above_20 = sum(1 for v in valid_ref if v > 20)
        max_ref = max(valid_ref) if valid_ref else -999

        # North/south cluster breakdown
        north_returns = [r for r in scan_returns if r['cluster'] == 'north']
        south_returns = [r for r in scan_returns if r['cluster'] == 'south']

        # Dual-pol for each cluster
        def cluster_dp_stats(rets):
            stats = {}
            zdr_vals = [r['zdr'] for r in rets if r.get('zdr') is not None]
            cc_vals = [r['cc'] for r in rets if r.get('cc') is not None]
            vel_vals = [r['vel_ms'] for r in rets if r.get('vel_ms') is not None]
            if zdr_vals:
                stats['zdr_mean'] = float(np.mean(zdr_vals))
                stats['zdr_min'] = float(np.min(zdr_vals))
                stats['zdr_max'] = float(np.max(zdr_vals))
                stats['zdr_std'] = float(np.std(zdr_vals))
            if cc_vals:
                stats['cc_mean'] = float(np.mean(cc_vals))
                stats['cc_min'] = float(np.min(cc_vals))
                stats['cc_max'] = float(np.max(cc_vals))
                stats['cc_std'] = float(np.std(cc_vals))
            if vel_vals:
                stats['vel_mean'] = float(np.mean(vel_vals))
                stats['vel_min'] = float(np.min(vel_vals))
                stats['vel_max'] = float(np.max(vel_vals))
            stats['n_returns'] = len(rets)
            stats['max_ref'] = max(r['ref_dBZ'] for r in rets) if rets else -999
            stats['n_above_10'] = sum(1 for r in rets if r['ref_dBZ'] > 10)
            stats['n_above_20'] = sum(1 for r in rets if r['ref_dBZ'] > 20)
            return stats

        north_dp = cluster_dp_stats(north_returns)
        south_dp = cluster_dp_stats(south_returns)

        # Overall sector dual-pol
        sector_dp = cluster_dp_stats(scan_returns)

        scan_result = {
            'file': fname,
            'time': time_str,
            'n_valid_pixels': len(valid_ref),
            'n_above_5dBZ': above_5,
            'n_above_10dBZ': above_10,
            'n_above_15dBZ': above_15,
            'n_above_20dBZ': above_20,
            'max_ref_dBZ': max_ref,
            'n_significant_returns': len(scan_returns),
            'north_cluster': north_dp,
            'south_cluster': south_dp,
            'sector_dual_pol': sector_dp,
        }
        all_scan_results.append(scan_result)

        # Print per-scan summary
        n_str = f"N:{north_dp['n_returns']:>3}"
        s_str = f"S:{south_dp['n_returns']:>3}"
        cc_n = f"{north_dp['cc_mean']:.3f}" if 'cc_mean' in north_dp else "---"
        cc_s = f"{south_dp['cc_mean']:.3f}" if 'cc_mean' in south_dp else "---"
        print(f"  {time_str}: max={max_ref:>5.1f}dBZ, >10={above_10:>4}, >20={above_20:>4} | "
              f"{n_str} CC={cc_n} | {s_str} CC={cc_s}")

        del radar

    return all_scan_results, all_returns


def compute_summary(scan_results, all_returns, label):
    """Compute aggregate statistics for a dataset."""
    summary = {'label': label, 'n_scans': len(scan_results)}

    if not scan_results:
        return summary

    # Aggregate scan-level stats
    summary['avg_returns_per_scan'] = float(np.mean([s['n_significant_returns'] for s in scan_results]))
    summary['avg_above_10'] = float(np.mean([s['n_above_10dBZ'] for s in scan_results]))
    summary['avg_above_20'] = float(np.mean([s['n_above_20dBZ'] for s in scan_results]))
    summary['max_ref_overall'] = float(max(s['max_ref_dBZ'] for s in scan_results))
    summary['min_max_ref'] = float(min(s['max_ref_dBZ'] for s in scan_results))
    summary['std_above_10'] = float(np.std([s['n_above_10dBZ'] for s in scan_results]))

    # North cluster aggregates
    north_cc_means = [s['north_cluster']['cc_mean'] for s in scan_results
                      if 'cc_mean' in s['north_cluster']]
    north_zdr_means = [s['north_cluster']['zdr_mean'] for s in scan_results
                       if 'zdr_mean' in s['north_cluster']]
    north_returns = [s['north_cluster']['n_returns'] for s in scan_results]
    north_above10 = [s['north_cluster']['n_above_10'] for s in scan_results]
    north_above20 = [s['north_cluster']['n_above_20'] for s in scan_results]
    north_max_ref = [s['north_cluster']['max_ref'] for s in scan_results
                     if s['north_cluster']['max_ref'] > -999]

    summary['north'] = {
        'avg_returns': float(np.mean(north_returns)),
        'std_returns': float(np.std(north_returns)),
        'avg_above_10': float(np.mean(north_above10)),
        'avg_above_20': float(np.mean(north_above20)),
        'avg_cc': float(np.mean(north_cc_means)) if north_cc_means else None,
        'std_cc': float(np.std(north_cc_means)) if north_cc_means else None,
        'avg_zdr': float(np.mean(north_zdr_means)) if north_zdr_means else None,
        'max_ref': float(max(north_max_ref)) if north_max_ref else None,
    }

    # South cluster aggregates
    south_cc_means = [s['south_cluster']['cc_mean'] for s in scan_results
                      if 'cc_mean' in s['south_cluster']]
    south_zdr_means = [s['south_cluster']['zdr_mean'] for s in scan_results
                       if 'zdr_mean' in s['south_cluster']]
    south_returns = [s['south_cluster']['n_returns'] for s in scan_results]
    south_above10 = [s['south_cluster']['n_above_10'] for s in scan_results]
    south_above20 = [s['south_cluster']['n_above_20'] for s in scan_results]
    south_max_ref = [s['south_cluster']['max_ref'] for s in scan_results
                     if s['south_cluster']['max_ref'] > -999]

    summary['south'] = {
        'avg_returns': float(np.mean(south_returns)),
        'std_returns': float(np.std(south_returns)),
        'avg_above_10': float(np.mean(south_above10)),
        'avg_above_20': float(np.mean(south_above20)),
        'avg_cc': float(np.mean(south_cc_means)) if south_cc_means else None,
        'std_cc': float(np.std(south_cc_means)) if south_cc_means else None,
        'avg_zdr': float(np.mean(south_zdr_means)) if south_zdr_means else None,
        'max_ref': float(max(south_max_ref)) if south_max_ref else None,
    }

    # Overall dual-pol from all returns
    dp_returns = [r for r in all_returns if r.get('cc') is not None and r.get('zdr') is not None]
    if dp_returns:
        all_cc = [r['cc'] for r in dp_returns]
        all_zdr = [r['zdr'] for r in dp_returns]
        summary['overall_cc_mean'] = float(np.mean(all_cc))
        summary['overall_cc_median'] = float(np.median(all_cc))
        summary['overall_cc_std'] = float(np.std(all_cc))
        summary['overall_zdr_mean'] = float(np.mean(all_zdr))
        summary['overall_zdr_std'] = float(np.std(all_zdr))
        summary['overall_zdr_min'] = float(np.min(all_zdr))
        summary['overall_zdr_max'] = float(np.max(all_zdr))

    # Dual-pol classification
    if dp_returns:
        n_rain = sum(1 for r in dp_returns if r['cc'] > 0.95 and -1 < r['zdr'] < 4)
        n_ground = sum(1 for r in dp_returns if r['cc'] > 0.85 and abs(r['zdr']) < 2)
        n_debris = sum(1 for r in dp_returns if r['cc'] < 0.80 and (r['zdr'] > 4 or r['zdr'] < -2))
        n_chaff = sum(1 for r in dp_returns if r['cc'] < 0.70 and r['zdr'] > 6)
        n_metallic = sum(1 for r in dp_returns if r['cc'] < 0.60 and r['zdr'] > 8)
        n_bio = sum(1 for r in dp_returns if 0.3 < r['cc'] < 0.7 and r['zdr'] > 3)

        total = len(dp_returns)
        summary['classification'] = {
            'rain_pct': 100 * n_rain / total,
            'ground_clutter_pct': 100 * n_ground / total,
            'debris_pct': 100 * n_debris / total,
            'chaff_pct': 100 * n_chaff / total,
            'metallic_pct': 100 * n_metallic / total,
            'biological_pct': 100 * n_bio / total,
            'total_dp_returns': total,
        }

    return summary


# ============================================================================
# MAIN ANALYSIS
# ============================================================================
print("=" * 80)
print("KEPZ BASELINE COMPARISON ANALYSIS")
print("=" * 80)
print("PURPOSE: Compare Feb 10 (control night) vs Feb 11 (event night)")
print("         to determine if KEPZ returns are ground clutter or engagement debris")
print()
print("METHODOLOGY: Identical sector parameters, thresholds, and analysis pipeline")
print(f"SECTOR: az {TARGET_AZ_MIN}-{TARGET_AZ_MAX}°, range {TARGET_RANGE_MIN/1000:.0f}-{TARGET_RANGE_MAX/1000:.0f} km")
print(f"CLUSTERS: North (az {NORTH_AZ_MIN}-{NORTH_AZ_MAX}°) vs South (az {SOUTH_AZ_MIN}-{SOUTH_AZ_MAX}°)")
print()

# Analyze baseline
baseline_scans, baseline_returns = analyze_dataset(
    BASELINE_DIR, 'Level2_KEPZ_20260210_*.ar2v', 'BASELINE — Feb 10, 2026 (night before)')

# Analyze event
event_scans, event_returns = analyze_dataset(
    EVENT_DIR, 'Level2_KEPZ_20260211_*.ar2v', 'EVENT — Feb 11, 2026 (incident night)')

# Compute summaries
baseline_summary = compute_summary(baseline_scans, baseline_returns, 'Baseline Feb 10')
event_summary = compute_summary(event_scans, event_returns, 'Event Feb 11')

# ============================================================================
# SIDE-BY-SIDE COMPARISON
# ============================================================================
print("\n\n" + "=" * 80)
print("SIDE-BY-SIDE COMPARISON: BASELINE vs EVENT")
print("=" * 80)

def fmt(val, decimals=1):
    if val is None:
        return "N/A"
    return f"{val:.{decimals}f}"

print(f"\n{'Metric':<40} {'Baseline (Feb 10)':>20} {'Event (Feb 11)':>20} {'Delta':>12}")
print(f"{'─'*40} {'─'*20} {'─'*20} {'─'*12}")

# Total scans
print(f"{'Scans analyzed':<40} {baseline_summary['n_scans']:>20} {event_summary['n_scans']:>20} {'':>12}")

# Returns per scan
b_rps = baseline_summary.get('avg_returns_per_scan', 0)
e_rps = event_summary.get('avg_returns_per_scan', 0)
delta = e_rps - b_rps
print(f"{'Avg returns/scan (>5 dBZ)':<40} {fmt(b_rps):>20} {fmt(e_rps):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

# Above 10 dBZ
b_10 = baseline_summary.get('avg_above_10', 0)
e_10 = event_summary.get('avg_above_10', 0)
delta = e_10 - b_10
print(f"{'Avg pixels > 10 dBZ/scan':<40} {fmt(b_10):>20} {fmt(e_10):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

# Above 20 dBZ
b_20 = baseline_summary.get('avg_above_20', 0)
e_20 = event_summary.get('avg_above_20', 0)
delta = e_20 - b_20
print(f"{'Avg pixels > 20 dBZ/scan':<40} {fmt(b_20):>20} {fmt(e_20):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

# Variability
b_std = baseline_summary.get('std_above_10', 0)
e_std = event_summary.get('std_above_10', 0)
print(f"{'Std dev of >10dBZ count (variability)':<40} {fmt(b_std):>20} {fmt(e_std):>20} {'':>12}")

# Max reflectivity
b_max = baseline_summary.get('max_ref_overall', 0)
e_max = event_summary.get('max_ref_overall', 0)
delta = e_max - b_max
print(f"{'Max reflectivity (dBZ)':<40} {fmt(b_max):>20} {fmt(e_max):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

# Overall CC
b_cc = baseline_summary.get('overall_cc_mean')
e_cc = event_summary.get('overall_cc_mean')
delta_cc = (e_cc - b_cc) if (b_cc and e_cc) else None
print(f"{'Overall mean CC':<40} {fmt(b_cc, 4):>20} {fmt(e_cc, 4):>20} {'+' if delta_cc and delta_cc > 0 else ''}{fmt(delta_cc, 4) if delta_cc else 'N/A':>11}")

# Overall ZDR
b_zdr = baseline_summary.get('overall_zdr_mean')
e_zdr = event_summary.get('overall_zdr_mean')
delta_zdr = (e_zdr - b_zdr) if (b_zdr and e_zdr) else None
print(f"{'Overall mean ZDR (dB)':<40} {fmt(b_zdr, 2):>20} {fmt(e_zdr, 2):>20} {'+' if delta_zdr and delta_zdr > 0 else ''}{fmt(delta_zdr, 2) if delta_zdr else 'N/A':>11}")

# ZDR range
b_zdr_min = baseline_summary.get('overall_zdr_min')
b_zdr_max = baseline_summary.get('overall_zdr_max')
e_zdr_min = event_summary.get('overall_zdr_min')
e_zdr_max = event_summary.get('overall_zdr_max')
b_range = f"{fmt(b_zdr_min)} to {fmt(b_zdr_max)}" if b_zdr_min is not None else "N/A"
e_range = f"{fmt(e_zdr_min)} to {fmt(e_zdr_max)}" if e_zdr_min is not None else "N/A"
print(f"{'ZDR range':<40} {b_range:>20} {e_range:>20}")

# ============================================================================
# CLUSTER COMPARISON
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"NORTH CLUSTER (az {NORTH_AZ_MIN}-{NORTH_AZ_MAX}° — Franklin Mountains sector)")
print(f"{'─'*80}")

bn = baseline_summary.get('north', {})
en = event_summary.get('north', {})

print(f"{'Metric':<40} {'Baseline':>20} {'Event':>20} {'Delta':>12}")
print(f"{'─'*40} {'─'*20} {'─'*20} {'─'*12}")

delta = (en.get('avg_returns', 0) - bn.get('avg_returns', 0))
print(f"{'Avg returns/scan':<40} {fmt(bn.get('avg_returns')):>20} {fmt(en.get('avg_returns')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

delta = (en.get('avg_above_10', 0) - bn.get('avg_above_10', 0))
print(f"{'Avg >10 dBZ/scan':<40} {fmt(bn.get('avg_above_10')):>20} {fmt(en.get('avg_above_10')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

delta = (en.get('avg_above_20', 0) - bn.get('avg_above_20', 0))
print(f"{'Avg >20 dBZ/scan':<40} {fmt(bn.get('avg_above_20')):>20} {fmt(en.get('avg_above_20')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

print(f"{'Max reflectivity':<40} {fmt(bn.get('max_ref')):>20} {fmt(en.get('max_ref')):>20}")
print(f"{'Mean CC':<40} {fmt(bn.get('avg_cc'), 4):>20} {fmt(en.get('avg_cc'), 4):>20}")
print(f"{'CC std dev':<40} {fmt(bn.get('std_cc'), 4):>20} {fmt(en.get('std_cc'), 4):>20}")
print(f"{'Mean ZDR':<40} {fmt(bn.get('avg_zdr'), 2):>20} {fmt(en.get('avg_zdr'), 2):>20}")

print(f"\n\n{'─'*80}")
print(f"SOUTH CLUSTER (az {SOUTH_AZ_MIN}-{SOUTH_AZ_MAX}° — away from Franklin Mountains)")
print(f"{'─'*80}")

bs = baseline_summary.get('south', {})
es = event_summary.get('south', {})

print(f"{'Metric':<40} {'Baseline':>20} {'Event':>20} {'Delta':>12}")
print(f"{'─'*40} {'─'*20} {'─'*20} {'─'*12}")

delta = (es.get('avg_returns', 0) - bs.get('avg_returns', 0))
print(f"{'Avg returns/scan':<40} {fmt(bs.get('avg_returns')):>20} {fmt(es.get('avg_returns')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

delta = (es.get('avg_above_10', 0) - bs.get('avg_above_10', 0))
print(f"{'Avg >10 dBZ/scan':<40} {fmt(bs.get('avg_above_10')):>20} {fmt(es.get('avg_above_10')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

delta = (es.get('avg_above_20', 0) - bs.get('avg_above_20', 0))
print(f"{'Avg >20 dBZ/scan':<40} {fmt(bs.get('avg_above_20')):>20} {fmt(es.get('avg_above_20')):>20} {'+' if delta > 0 else ''}{fmt(delta):>11}")

print(f"{'Max reflectivity':<40} {fmt(bs.get('max_ref')):>20} {fmt(es.get('max_ref')):>20}")
print(f"{'Mean CC':<40} {fmt(bs.get('avg_cc'), 4):>20} {fmt(es.get('avg_cc'), 4):>20}")
print(f"{'CC std dev':<40} {fmt(bs.get('std_cc'), 4):>20} {fmt(es.get('std_cc'), 4):>20}")
print(f"{'Mean ZDR':<40} {fmt(bs.get('avg_zdr'), 2):>20} {fmt(es.get('avg_zdr'), 2):>20}")


# ============================================================================
# DUAL-POL CLASSIFICATION COMPARISON
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"DUAL-POL TARGET CLASSIFICATION COMPARISON")
print(f"{'─'*80}")

bc = baseline_summary.get('classification', {})
ec = event_summary.get('classification', {})

print(f"\n{'Category':<40} {'Baseline %':>15} {'Event %':>15} {'Delta':>12}")
print(f"{'─'*40} {'─'*15} {'─'*15} {'─'*12}")

for key, label in [
    ('rain_pct', 'Rain-like'),
    ('ground_clutter_pct', 'Ground clutter'),
    ('biological_pct', 'Biological'),
    ('debris_pct', 'Debris/chaff'),
    ('chaff_pct', 'Military chaff'),
    ('metallic_pct', 'Metallic fragments'),
]:
    bv = bc.get(key, 0)
    ev = ec.get(key, 0)
    delta = ev - bv
    print(f"{'  ' + label:<40} {fmt(bv):>14}% {fmt(ev):>14}% {'+' if delta > 0 else ''}{fmt(delta):>10}%")


# ============================================================================
# AZIMUTH-BY-AZIMUTH BREAKDOWN
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"AZIMUTH BIN COMPARISON (5° bins)")
print(f"{'─'*80}")

az_bins = list(range(60, 115, 5))

print(f"\n{'Az Bin':<12} {'--- Baseline ---':>30} {'--- Event ---':>30}")
print(f"{'':>12} {'Rets':>8} {'MaxdBZ':>8} {'CC':>8} {'ZDR':>8} {'Rets':>8} {'MaxdBZ':>8} {'CC':>8} {'ZDR':>8}")
print(f"{'─'*12} {'─'*32} {'─'*32}")

for i in range(len(az_bins) - 1):
    az_lo = az_bins[i]
    az_hi = az_bins[i + 1]

    b_rets = [r for r in baseline_returns if az_lo <= r['az'] < az_hi]
    e_rets = [r for r in event_returns if az_lo <= r['az'] < az_hi]

    b_n = len(b_rets) / max(baseline_summary['n_scans'], 1)
    e_n = len(e_rets) / max(event_summary['n_scans'], 1)

    b_max = max((r['ref_dBZ'] for r in b_rets), default=-999)
    e_max = max((r['ref_dBZ'] for r in e_rets), default=-999)

    b_cc_vals = [r['cc'] for r in b_rets if r.get('cc') is not None]
    e_cc_vals = [r['cc'] for r in e_rets if r.get('cc') is not None]

    b_zdr_vals = [r['zdr'] for r in b_rets if r.get('zdr') is not None]
    e_zdr_vals = [r['zdr'] for r in e_rets if r.get('zdr') is not None]

    b_cc = f"{np.mean(b_cc_vals):.3f}" if b_cc_vals else "---"
    e_cc = f"{np.mean(e_cc_vals):.3f}" if e_cc_vals else "---"

    b_zdr = f"{np.mean(b_zdr_vals):.1f}" if b_zdr_vals else "---"
    e_zdr = f"{np.mean(e_zdr_vals):.1f}" if e_zdr_vals else "---"

    b_max_s = f"{b_max:.0f}" if b_max > -999 else "---"
    e_max_s = f"{e_max:.0f}" if e_max > -999 else "---"

    print(f"  {az_lo:>3}-{az_hi:<3}°   {b_n:>8.0f} {b_max_s:>8} {b_cc:>8} {b_zdr:>8}"
          f" {e_n:>8.0f} {e_max_s:>8} {e_cc:>8} {e_zdr:>8}")


# ============================================================================
# RANGE BIN COMPARISON
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"RANGE BIN COMPARISON (5 km bins)")
print(f"{'─'*80}")

range_bins = list(range(15, 60, 5))

print(f"\n{'Rng Bin':<12} {'--- Baseline ---':>30} {'--- Event ---':>30}")
print(f"{'':>12} {'Rets':>8} {'MaxdBZ':>8} {'CC':>8} {'ZDR':>8} {'Rets':>8} {'MaxdBZ':>8} {'CC':>8} {'ZDR':>8}")
print(f"{'─'*12} {'─'*32} {'─'*32}")

for i in range(len(range_bins) - 1):
    r_lo = range_bins[i]
    r_hi = range_bins[i + 1]

    b_rets = [r for r in baseline_returns if r_lo <= r['range_km'] < r_hi]
    e_rets = [r for r in event_returns if r_lo <= r['range_km'] < r_hi]

    b_n = len(b_rets) / max(baseline_summary['n_scans'], 1)
    e_n = len(e_rets) / max(event_summary['n_scans'], 1)

    b_max = max((r['ref_dBZ'] for r in b_rets), default=-999)
    e_max = max((r['ref_dBZ'] for r in e_rets), default=-999)

    b_cc_vals = [r['cc'] for r in b_rets if r.get('cc') is not None]
    e_cc_vals = [r['cc'] for r in e_rets if r.get('cc') is not None]

    b_zdr_vals = [r['zdr'] for r in b_rets if r.get('zdr') is not None]
    e_zdr_vals = [r['zdr'] for r in e_rets if r.get('zdr') is not None]

    b_cc = f"{np.mean(b_cc_vals):.3f}" if b_cc_vals else "---"
    e_cc = f"{np.mean(e_cc_vals):.3f}" if e_cc_vals else "---"

    b_zdr = f"{np.mean(b_zdr_vals):.1f}" if b_zdr_vals else "---"
    e_zdr = f"{np.mean(e_zdr_vals):.1f}" if e_zdr_vals else "---"

    b_max_s = f"{b_max:.0f}" if b_max > -999 else "---"
    e_max_s = f"{e_max:.0f}" if e_max > -999 else "---"

    print(f"  {r_lo:>3}-{r_hi:<3}km  {b_n:>8.0f} {b_max_s:>8} {b_cc:>8} {b_zdr:>8}"
          f" {e_n:>8.0f} {e_max_s:>8} {e_cc:>8} {e_zdr:>8}")


# ============================================================================
# STATISTICAL SIGNIFICANCE TEST
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"STATISTICAL SIGNIFICANCE TESTS")
print(f"{'─'*80}")

from scipy import stats

# Test: Are the per-scan return counts different?
b_counts = [s['n_above_10dBZ'] for s in baseline_scans]
e_counts = [s['n_above_10dBZ'] for s in event_scans]
t_stat, p_val = stats.ttest_ind(b_counts, e_counts, equal_var=False)
print(f"\n  Welch's t-test on per-scan >10 dBZ counts:")
print(f"    Baseline mean: {np.mean(b_counts):.1f} ± {np.std(b_counts):.1f}")
print(f"    Event mean:    {np.mean(e_counts):.1f} ± {np.std(e_counts):.1f}")
print(f"    t-statistic:   {t_stat:.3f}")
print(f"    p-value:       {p_val:.6f}")
print(f"    Significant?   {'YES (p < 0.05)' if p_val < 0.05 else 'NO (p >= 0.05)'}")

# Test: per-scan >20 dBZ counts
b_20counts = [s['n_above_20dBZ'] for s in baseline_scans]
e_20counts = [s['n_above_20dBZ'] for s in event_scans]
t_stat20, p_val20 = stats.ttest_ind(b_20counts, e_20counts, equal_var=False)
print(f"\n  Welch's t-test on per-scan >20 dBZ counts:")
print(f"    Baseline mean: {np.mean(b_20counts):.1f} ± {np.std(b_20counts):.1f}")
print(f"    Event mean:    {np.mean(e_20counts):.1f} ± {np.std(e_20counts):.1f}")
print(f"    t-statistic:   {t_stat20:.3f}")
print(f"    p-value:       {p_val20:.6f}")
print(f"    Significant?   {'YES (p < 0.05)' if p_val20 < 0.05 else 'NO (p >= 0.05)'}")

# Test: CC distribution comparison
b_cc_all = [r['cc'] for r in baseline_returns if r.get('cc') is not None]
e_cc_all = [r['cc'] for r in event_returns if r.get('cc') is not None]
if b_cc_all and e_cc_all:
    ks_stat, ks_p = stats.ks_2samp(b_cc_all, e_cc_all)
    print(f"\n  KS test on CC distributions:")
    print(f"    Baseline: n={len(b_cc_all)}, mean={np.mean(b_cc_all):.4f}, median={np.median(b_cc_all):.4f}")
    print(f"    Event:    n={len(e_cc_all)}, mean={np.mean(e_cc_all):.4f}, median={np.median(e_cc_all):.4f}")
    print(f"    KS stat:  {ks_stat:.4f}")
    print(f"    p-value:  {ks_p:.2e}")
    print(f"    Significant? {'YES' if ks_p < 0.05 else 'NO'}")

# Test: South cluster specifically
b_south_cc = [r['cc'] for r in baseline_returns if r.get('cc') is not None and r['cluster'] == 'south']
e_south_cc = [r['cc'] for r in event_returns if r.get('cc') is not None and r['cluster'] == 'south']
if b_south_cc and e_south_cc:
    ks_s, ks_p_s = stats.ks_2samp(b_south_cc, e_south_cc)
    mw_s, mw_p = stats.mannwhitneyu(b_south_cc, e_south_cc, alternative='two-sided')
    print(f"\n  SOUTH CLUSTER CC comparison:")
    print(f"    Baseline south: n={len(b_south_cc)}, mean={np.mean(b_south_cc):.4f}")
    print(f"    Event south:    n={len(e_south_cc)}, mean={np.mean(e_south_cc):.4f}")
    print(f"    KS stat:        {ks_s:.4f}, p={ks_p_s:.2e}")
    print(f"    Mann-Whitney U: stat={mw_s:.0f}, p={mw_p:.2e}")
    print(f"    Significant?    {'YES' if ks_p_s < 0.05 else 'NO'}")

# Test: South cluster ZDR
b_south_zdr = [r['zdr'] for r in baseline_returns if r.get('zdr') is not None and r['cluster'] == 'south']
e_south_zdr = [r['zdr'] for r in event_returns if r.get('zdr') is not None and r['cluster'] == 'south']
if b_south_zdr and e_south_zdr:
    ks_z, ks_p_z = stats.ks_2samp(b_south_zdr, e_south_zdr)
    print(f"\n  SOUTH CLUSTER ZDR comparison:")
    print(f"    Baseline south: n={len(b_south_zdr)}, mean={np.mean(b_south_zdr):.2f}, range=[{np.min(b_south_zdr):.1f}, {np.max(b_south_zdr):.1f}]")
    print(f"    Event south:    n={len(e_south_zdr)}, mean={np.mean(e_south_zdr):.2f}, range=[{np.min(e_south_zdr):.1f}, {np.max(e_south_zdr):.1f}]")
    print(f"    KS stat:        {ks_z:.4f}, p={ks_p_z:.2e}")
    print(f"    Significant?    {'YES' if ks_p_z < 0.05 else 'NO'}")

# Test: Max reflectivity per scan
b_max_refs = [s['max_ref_dBZ'] for s in baseline_scans]
e_max_refs = [s['max_ref_dBZ'] for s in event_scans]
t_max, p_max = stats.ttest_ind(b_max_refs, e_max_refs, equal_var=False)
print(f"\n  Max reflectivity per scan:")
print(f"    Baseline: mean={np.mean(b_max_refs):.1f} ± {np.std(b_max_refs):.1f}")
print(f"    Event:    mean={np.mean(e_max_refs):.1f} ± {np.std(e_max_refs):.1f}")
print(f"    t-stat:   {t_max:.3f}, p={p_max:.6f}")
print(f"    Significant? {'YES' if p_max < 0.05 else 'NO'}")


# ============================================================================
# VERDICT
# ============================================================================
print("\n\n" + "=" * 80)
print("VERDICT: BASELINE vs EVENT COMPARISON")
print("=" * 80)

# Determine verdict based on data
north_similar = True
south_different = False

if bn.get('avg_returns') and en.get('avg_returns'):
    n_ratio = en['avg_returns'] / max(bn['avg_returns'], 1)
    if abs(n_ratio - 1.0) < 0.3:
        north_similar = True
    else:
        north_similar = False

if bs.get('avg_returns') and es.get('avg_returns'):
    s_ratio = es['avg_returns'] / max(bs['avg_returns'], 1)
    if s_ratio > 1.5 or s_ratio < 0.5:
        south_different = True

# CC comparison
cc_different = False
if b_cc_all and e_cc_all:
    if abs(np.mean(b_cc_all) - np.mean(e_cc_all)) > 0.05:
        cc_different = True

print(f"""
FINDINGS:

1. NORTH CLUSTER (Franklin Mountains sector, az 60-90°):
   Baseline avg returns/scan:  {fmt(bn.get('avg_returns'))}
   Event avg returns/scan:     {fmt(en.get('avg_returns'))}
   Baseline CC:                {fmt(bn.get('avg_cc'), 4)}
   Event CC:                   {fmt(en.get('avg_cc'), 4)}
   INTERPRETATION: {'SIMILAR — confirms this is PERSISTENT GROUND CLUTTER' if north_similar else 'DIFFERENT — something changed in this sector'}

2. SOUTH CLUSTER (away from mountains, az 90-110°):
   Baseline avg returns/scan:  {fmt(bs.get('avg_returns'))}
   Event avg returns/scan:     {fmt(es.get('avg_returns'))}
   Baseline CC:                {fmt(bs.get('avg_cc'), 4)}
   Event CC:                   {fmt(es.get('avg_cc'), 4)}
   INTERPRETATION: {'SIGNIFICANTLY DIFFERENT — event night shows anomalous activity' if south_different else 'SIMILAR — suggests persistent returns in this sector too'}

3. OVERALL CC DISTRIBUTION:
   Baseline mean CC: {fmt(baseline_summary.get('overall_cc_mean'), 4)}
   Event mean CC:    {fmt(event_summary.get('overall_cc_mean'), 4)}
   INTERPRETATION: {'CC distributions are STATISTICALLY DIFFERENT' if cc_different else 'CC distributions are SIMILAR'}
""")

# Final assessment
if north_similar and south_different:
    print("""CONCLUSION: MIXED — Ground clutter confirmed in north cluster, but south cluster
shows significant differences between baseline and event night. The south cluster
anomaly supports the engagement debris hypothesis for az 90-110°.""")
elif north_similar and not south_different:
    print("""CONCLUSION: GROUND CLUTTER DOMINANT — Both clusters look similar on baseline
and event nights. The KEPZ dual-pol signatures are primarily terrain/ground clutter
artifacts, NOT engagement debris. The CC=0.34 finding from the event analysis is
likely a normal feature of this sector's clutter environment.""")
elif not north_similar and south_different:
    print("""CONCLUSION: SIGNIFICANT EVENT DETECTED — Both clusters show differences
between baseline and event night. This supports the engagement/debris hypothesis
for the entire sector.""")
else:
    print("""CONCLUSION: ANOMALOUS — North cluster changed but south didn't.
This is unexpected and requires further investigation.""")


# ============================================================================
# SAVE RESULTS
# ============================================================================
output = {
    'analysis': 'KEPZ Baseline Comparison',
    'baseline_date': '2026-02-10',
    'event_date': '2026-02-11',
    'parameters': {
        'sector_az': [TARGET_AZ_MIN, TARGET_AZ_MAX],
        'sector_range_km': [TARGET_RANGE_MIN/1000, TARGET_RANGE_MAX/1000],
        'ref_threshold': REF_THRESHOLD,
        'north_cluster_az': [NORTH_AZ_MIN, NORTH_AZ_MAX],
        'south_cluster_az': [SOUTH_AZ_MIN, SOUTH_AZ_MAX],
    },
    'baseline_summary': baseline_summary,
    'event_summary': event_summary,
    'baseline_scans': baseline_scans,
    'event_scans': event_scans,
    'statistical_tests': {
        'above_10_ttest': {'t': float(t_stat), 'p': float(p_val)},
        'above_20_ttest': {'t': float(t_stat20), 'p': float(p_val20)},
        'max_ref_ttest': {'t': float(t_max), 'p': float(p_max)},
    },
}

if b_cc_all and e_cc_all:
    output['statistical_tests']['cc_ks'] = {'stat': float(ks_stat), 'p': float(ks_p)}
if b_south_cc and e_south_cc:
    output['statistical_tests']['south_cc_ks'] = {'stat': float(ks_s), 'p': float(ks_p_s)}
if b_south_zdr and e_south_zdr:
    output['statistical_tests']['south_zdr_ks'] = {'stat': float(ks_z), 'p': float(ks_p_z)}

output_path = os.path.join(OUTPUT_DIR, 'baseline_comparison_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nResults saved to {output_path}")
