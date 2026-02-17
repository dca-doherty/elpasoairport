#!/usr/bin/env python3
"""
KEPZ Level 2 NEXRAD Full Analysis
El Paso Airspace Incident — Feb 11, 2026

Now with ACTUAL Level 2 data (not imagery analysis).
Extracts for each scan:
  - Reflectivity (dBZ) in target sector
  - Differential reflectivity (ZDR) — shape/orientation
  - Correlation coefficient (CC/RhoHV) — target homogeneity
  - Doppler velocity — radial motion
  - Spectrum width — turbulence/multi-scatterer
  - Precise azimuth, range, lat/lon of peak returns
"""

import numpy as np
import pyart
import os
import glob
import json
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PARAMETERS
# ============================================================================
KEPZ_LAT = 31.8731
KEPZ_LON = -106.6981
KEPZ_ELEV = 1206.0  # meters MSL

# Target sector from KEPZ (Fort Bliss / engagement area)
# From prior imagery analysis: spikes at az ~82-88°, range ~30 km
# Expand search window to catch everything
TARGET_AZ_MIN = 60.0
TARGET_AZ_MAX = 110.0
TARGET_RANGE_MIN = 15000   # 15 km (meters)
TARGET_RANGE_MAX = 55000   # 55 km

# Reflectivity threshold for "significant" return
REF_THRESHOLD = 5.0  # dBZ

DATA_DIR = '/home/user/uap-transient-research/el_paso_airspace/kepz_data/'
OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def radar_to_latlon(radar_lat, radar_lon, range_m, az_deg):
    """Convert radar polar coords to lat/lon (flat earth approx, fine at 30km)."""
    range_km = range_m / 1000.0
    lat = radar_lat + (range_km * np.cos(np.radians(az_deg))) / 111.0
    lon = radar_lon + (range_km * np.sin(np.radians(az_deg))) / (111.0 * np.cos(np.radians(radar_lat)))
    return lat, lon

def get_sweep_data(radar, sweep_idx, field_name):
    """Extract data for a given sweep."""
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
    """Extract scalar from masked array."""
    if hasattr(val, 'item'):
        if np.ma.is_masked(val):
            return None
        return float(val.item())
    if np.isnan(val):
        return None
    return float(val)

# ============================================================================
# SCAN ALL FILES
# ============================================================================
files = sorted(glob.glob(os.path.join(DATA_DIR, 'Level2_KEPZ_*.ar2v')))
print("=" * 80)
print("KEPZ LEVEL 2 NEXRAD FULL ANALYSIS")
print("El Paso Airspace Incident — Feb 11, 2026")
print(f"Target sector: az {TARGET_AZ_MIN}-{TARGET_AZ_MAX}°, range {TARGET_RANGE_MIN/1000:.0f}-{TARGET_RANGE_MAX/1000:.0f} km")
print(f"Files to analyze: {len(files)}")
print("=" * 80)

all_scan_results = []
all_significant_returns = []
all_dual_pol = []

for filepath in files:
    fname = os.path.basename(filepath)
    time_str = fname.split('_')[-1].replace('.ar2v', '')

    try:
        radar = pyart.io.read_nexrad_archive(filepath)
    except Exception as e:
        print(f"\n{fname}: ERROR - {e}")
        continue

    scan_time = radar.time['units'].replace('seconds since ', '')
    ranges = radar.range['data']  # meters

    # Identify available fields
    available_fields = list(radar.fields.keys())

    # Map field names (NEXRAD naming varies)
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
    for candidate in ['differential_phase', 'PHI', 'phi', 'PHIDP', 'KDP']:
        if candidate in radar.fields:
            field_map['phidp'] = candidate
            break

    # Process lowest elevation sweep (0.5°)
    sweep_idx = 0
    elev_angle = radar.fixed_angle['data'][sweep_idx]
    azimuths = get_sweep_azimuths(radar, sweep_idx)

    # Build sector masks
    az_mask = np.zeros(len(azimuths), dtype=bool)
    for i, az in enumerate(azimuths):
        if TARGET_AZ_MIN <= az <= TARGET_AZ_MAX:
            az_mask[i] = True

    rng_mask = (ranges >= TARGET_RANGE_MIN) & (ranges <= TARGET_RANGE_MAX)

    sector_az = azimuths[az_mask]
    sector_rng = ranges[rng_mask]

    # Get reflectivity in sector
    if 'ref' not in field_map:
        print(f"\n{time_str}: No reflectivity field found. Fields: {available_fields}")
        del radar
        continue

    ref_data = get_sweep_data(radar, sweep_idx, field_map['ref'])
    sector_ref = ref_data[az_mask][:, rng_mask]

    # Get dual-pol fields in sector
    sector_zdr = None
    sector_cc = None
    sector_vel = None
    sector_sw = None
    sector_phidp = None

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

    if 'phidp' in field_map:
        phidp_data = get_sweep_data(radar, sweep_idx, field_map['phidp'])
        if phidp_data is not None:
            sector_phidp = phidp_data[az_mask][:, rng_mask]

    # Find significant returns
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
                'scan_time': scan_time,
                'az': az,
                'range_m': rng,
                'range_km': rng / 1000.0,
                'lat': lat,
                'lon': lon,
                'ref_dBZ': ref_val,
                'elev_angle': float(elev_angle),
            }

            # Add dual-pol values at this gate
            if sector_zdr is not None:
                ret['zdr'] = safe_val(sector_zdr[i, j])
            if sector_cc is not None:
                ret['cc'] = safe_val(sector_cc[i, j])
            if sector_vel is not None:
                ret['vel_ms'] = safe_val(sector_vel[i, j])
            if sector_sw is not None:
                ret['spectrum_width'] = safe_val(sector_sw[i, j])
            if sector_phidp is not None:
                ret['phidp'] = safe_val(sector_phidp[i, j])

            scan_returns.append(ret)
            all_significant_returns.append(ret)

    # Compute scan-level statistics
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

    # Peak return for this scan
    peak_ret = max(scan_returns, key=lambda r: r['ref_dBZ']) if scan_returns else None

    # Dual-pol summary for returns > 5 dBZ
    dp_summary = {}
    if scan_returns:
        zdr_vals = [r['zdr'] for r in scan_returns if r.get('zdr') is not None]
        cc_vals = [r['cc'] for r in scan_returns if r.get('cc') is not None]
        vel_vals = [r['vel_ms'] for r in scan_returns if r.get('vel_ms') is not None]
        sw_vals = [r['spectrum_width'] for r in scan_returns if r.get('spectrum_width') is not None]

        if zdr_vals:
            dp_summary['zdr_mean'] = float(np.mean(zdr_vals))
            dp_summary['zdr_min'] = float(np.min(zdr_vals))
            dp_summary['zdr_max'] = float(np.max(zdr_vals))
            dp_summary['zdr_std'] = float(np.std(zdr_vals))
        if cc_vals:
            dp_summary['cc_mean'] = float(np.mean(cc_vals))
            dp_summary['cc_min'] = float(np.min(cc_vals))
            dp_summary['cc_max'] = float(np.max(cc_vals))
        if vel_vals:
            dp_summary['vel_mean'] = float(np.mean(vel_vals))
            dp_summary['vel_min'] = float(np.min(vel_vals))
            dp_summary['vel_max'] = float(np.max(vel_vals))
        if sw_vals:
            dp_summary['sw_mean'] = float(np.mean(sw_vals))
            dp_summary['sw_max'] = float(np.max(sw_vals))

    scan_result = {
        'file': fname,
        'time': time_str,
        'scan_time': scan_time,
        'elev_angle': float(elev_angle),
        'available_fields': available_fields,
        'n_valid_pixels': len(valid_ref),
        'n_above_5dBZ': above_5,
        'n_above_10dBZ': above_10,
        'n_above_15dBZ': above_15,
        'n_above_20dBZ': above_20,
        'max_ref_dBZ': max_ref,
        'n_significant_returns': len(scan_returns),
        'dual_pol': dp_summary,
    }

    if peak_ret:
        scan_result['peak'] = {
            'az': peak_ret['az'],
            'range_km': peak_ret['range_km'],
            'lat': peak_ret['lat'],
            'lon': peak_ret['lon'],
            'ref_dBZ': peak_ret['ref_dBZ'],
            'zdr': peak_ret.get('zdr'),
            'cc': peak_ret.get('cc'),
            'vel_ms': peak_ret.get('vel_ms'),
        }

    all_scan_results.append(scan_result)

    # Print summary
    status = ""
    if above_10 > 0:
        status = f" ** {above_10} pixels > 10 dBZ"
    if above_20 > 0:
        status = f" *** STRONG: {above_20} pixels > 20 dBZ, max {max_ref:.1f} dBZ ***"

    print(f"\n{time_str} UTC | elev {elev_angle:.1f}° | {len(valid_ref)} valid, "
          f"{above_5} >5dBZ, {above_10} >10dBZ | max {max_ref:.1f} dBZ{status}")

    if peak_ret:
        print(f"  Peak: {peak_ret['ref_dBZ']:.1f} dBZ at az={peak_ret['az']:.1f}°, "
              f"range={peak_ret['range_km']:.1f}km -> {peak_ret['lat']:.4f}°N, {abs(peak_ret['lon']):.4f}°W")
        if peak_ret.get('zdr') is not None:
            print(f"  Peak dual-pol: ZDR={peak_ret.get('zdr', 'N/A'):.1f}, "
                  f"CC={peak_ret.get('cc', 'N/A')}, "
                  f"Vel={peak_ret.get('vel_ms', 'N/A')} m/s")

    if dp_summary:
        parts = []
        if 'zdr_mean' in dp_summary:
            parts.append(f"ZDR={dp_summary['zdr_mean']:.1f}±{dp_summary.get('zdr_std', 0):.1f}")
        if 'cc_mean' in dp_summary:
            parts.append(f"CC={dp_summary['cc_mean']:.3f}")
        if 'vel_mean' in dp_summary:
            parts.append(f"Vel={dp_summary['vel_mean']:.1f}m/s")
        if 'sw_mean' in dp_summary:
            parts.append(f"SW={dp_summary['sw_mean']:.1f}")
        if parts:
            print(f"  Sector dual-pol mean: {', '.join(parts)}")

    del radar

# ============================================================================
# TIMELINE ANALYSIS
# ============================================================================
print("\n\n" + "=" * 80)
print("KEPZ LEVEL 2 — TIMELINE SUMMARY")
print("=" * 80)

print(f"\n{'Time':>6} {'Max dBZ':>8} {'>5':>5} {'>10':>5} {'>15':>5} {'>20':>5} "
      f"{'ZDR':>8} {'CC':>8} {'Vel m/s':>8} {'Peak Az':>8} {'Peak Rng':>9} {'Assessment':>20}")
print(f"{'─'*6} {'─'*8} {'─'*5} {'─'*5} {'─'*5} {'─'*5} "
      f"{'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*9} {'─'*20}")

for sr in all_scan_results:
    t = sr['time']
    dp = sr['dual_pol']

    zdr_str = f"{dp['zdr_mean']:.1f}" if 'zdr_mean' in dp else "---"
    cc_str = f"{dp['cc_mean']:.3f}" if 'cc_mean' in dp else "---"
    vel_str = f"{dp['vel_mean']:.1f}" if 'vel_mean' in dp else "---"

    peak_az = f"{sr['peak']['az']:.0f}°" if sr.get('peak') else "---"
    peak_rng = f"{sr['peak']['range_km']:.1f}km" if sr.get('peak') else "---"

    # Classify the scan
    if sr['n_above_20dBZ'] > 0:
        assessment = "*** ENGAGEMENT ***"
    elif sr['n_above_15dBZ'] > 0:
        assessment = "** ACTIVE **"
    elif sr['n_above_10dBZ'] > 0:
        assessment = "* elevated *"
    elif sr['n_above_5dBZ'] > 0:
        assessment = "returns present"
    else:
        assessment = "quiet"

    print(f"{t:>6} {sr['max_ref_dBZ']:>8.1f} {sr['n_above_5dBZ']:>5} "
          f"{sr['n_above_10dBZ']:>5} {sr['n_above_15dBZ']:>5} {sr['n_above_20dBZ']:>5} "
          f"{zdr_str:>8} {cc_str:>8} {vel_str:>8} {peak_az:>8} {peak_rng:>9} {assessment:>20}")

# ============================================================================
# DUAL-POL DEEP DIVE
# ============================================================================
print("\n\n" + "=" * 80)
print("DUAL-POL SIGNATURE ANALYSIS")
print("=" * 80)

# Collect all returns with dual-pol data
returns_with_dp = [r for r in all_significant_returns
                   if r.get('zdr') is not None and r.get('cc') is not None]

if returns_with_dp:
    all_zdr = [r['zdr'] for r in returns_with_dp]
    all_cc = [r['cc'] for r in returns_with_dp]
    all_ref = [r['ref_dBZ'] for r in returns_with_dp]

    print(f"\n  Total returns with dual-pol data: {len(returns_with_dp)}")
    print(f"\n  ZDR (Differential Reflectivity):")
    print(f"    Mean:   {np.mean(all_zdr):.2f} dB")
    print(f"    Median: {np.median(all_zdr):.2f} dB")
    print(f"    Min:    {np.min(all_zdr):.2f} dB")
    print(f"    Max:    {np.max(all_zdr):.2f} dB")
    print(f"    Std:    {np.std(all_zdr):.2f} dB")

    print(f"\n  CC (Correlation Coefficient / RhoHV):")
    print(f"    Mean:   {np.mean(all_cc):.4f}")
    print(f"    Median: {np.median(all_cc):.4f}")
    print(f"    Min:    {np.min(all_cc):.4f}")
    print(f"    Max:    {np.max(all_cc):.4f}")

    # Velocity analysis
    vel_returns = [r for r in returns_with_dp if r.get('vel_ms') is not None]
    if vel_returns:
        all_vel = [r['vel_ms'] for r in vel_returns]
        print(f"\n  Doppler Velocity (radial component):")
        print(f"    Mean:   {np.mean(all_vel):.2f} m/s")
        print(f"    Min:    {np.min(all_vel):.2f} m/s (toward radar)")
        print(f"    Max:    {np.max(all_vel):.2f} m/s (away from radar)")
        print(f"    Std:    {np.std(all_vel):.2f} m/s")

    # Spectrum width
    sw_returns = [r for r in returns_with_dp if r.get('spectrum_width') is not None]
    if sw_returns:
        all_sw = [r['spectrum_width'] for r in sw_returns]
        print(f"\n  Spectrum Width:")
        print(f"    Mean:   {np.mean(all_sw):.2f} m/s")
        print(f"    Max:    {np.max(all_sw):.2f} m/s")

    # Classify returns by dual-pol signature
    print("\n\n  DUAL-POL SIGNATURE CLASSIFICATION:")
    print("  " + "─" * 60)

    # Standard thresholds
    n_rain = sum(1 for r in returns_with_dp if r['cc'] > 0.95 and -1 < r['zdr'] < 4)
    n_debris = sum(1 for r in returns_with_dp if r['cc'] < 0.80 and (r['zdr'] > 4 or r['zdr'] < -2))
    n_chaff = sum(1 for r in returns_with_dp if r['cc'] < 0.70 and r['zdr'] > 6)
    n_bio = sum(1 for r in returns_with_dp if 0.3 < r['cc'] < 0.7 and r['zdr'] > 3)
    n_ground = sum(1 for r in returns_with_dp if r['cc'] > 0.85 and abs(r['zdr']) < 2)
    n_metallic = sum(1 for r in returns_with_dp if r['cc'] < 0.60 and r['zdr'] > 8)

    print(f"    Rain-like (CC>0.95, ZDR -1 to 4):          {n_rain:>5} ({100*n_rain/len(returns_with_dp):.1f}%)")
    print(f"    Ground clutter-like (CC>0.85, |ZDR|<2):    {n_ground:>5} ({100*n_ground/len(returns_with_dp):.1f}%)")
    print(f"    Biological (CC 0.3-0.7, ZDR>3):            {n_bio:>5} ({100*n_bio/len(returns_with_dp):.1f}%)")
    print(f"    Debris/chaff (CC<0.80, ZDR>4 or ZDR<-2):   {n_debris:>5} ({100*n_debris/len(returns_with_dp):.1f}%)")
    print(f"    Military chaff (CC<0.70, ZDR>6):            {n_chaff:>5} ({100*n_chaff/len(returns_with_dp):.1f}%)")
    print(f"    Metallic fragments (CC<0.60, ZDR>8):        {n_metallic:>5} ({100*n_metallic/len(returns_with_dp):.1f}%)")

    # Reference values
    print("""
  REFERENCE — Expected dual-pol values for different targets:
  ┌────────────────────┬───────────┬───────────┬──────────────────────────┐
  │ Target Type        │ ZDR (dB)  │ CC/RhoHV  │ Notes                    │
  ├────────────────────┼───────────┼───────────┼──────────────────────────┤
  │ Rain               │ 0 to 4   │ >0.97     │ Spherical/oblate drops   │
  │ Snow               │ 0 to 1.5 │ >0.95     │ Near-spherical crystals  │
  │ Hail               │ -1 to 2  │ 0.85-0.95 │ Tumbling spheroids       │
  │ Ground clutter     │ ~0       │ >0.90     │ Stationary hard targets  │
  │ Birds/insects      │ 2 to 8   │ 0.3-0.7   │ Elongated, variable      │
  │ Military chaff     │ 4 to 12+ │ 0.2-0.5   │ Thin metallic dipoles    │
  │ Metallic debris    │ 5 to 20+ │ 0.1-0.6   │ Irregular metal fragments│
  │ Mylar balloon      │ 3 to 8   │ 0.4-0.7   │ Thin metallic film       │
  │ Drone/UAS          │ -2 to 2  │ 0.7-0.95  │ Rigid structure          │
  │ Destroyed drone    │ 5 to 15+ │ 0.2-0.6   │ Debris cloud of parts    │
  └────────────────────┴───────────┴───────────┴──────────────────────────┘
""")

    # Time evolution of dual-pol
    print("  DUAL-POL TIME EVOLUTION (for returns > 10 dBZ):")
    print("  " + "─" * 60)

    # Group by scan time
    from collections import defaultdict
    by_time = defaultdict(list)
    for r in returns_with_dp:
        if r['ref_dBZ'] > 10:
            by_time[r['time']].append(r)

    for t in sorted(by_time.keys()):
        rets = by_time[t]
        if not rets:
            continue
        zdr_vals = [r['zdr'] for r in rets if r['zdr'] is not None]
        cc_vals = [r['cc'] for r in rets if r['cc'] is not None]
        vel_vals = [r['vel_ms'] for r in rets if r.get('vel_ms') is not None]
        ref_vals = [r['ref_dBZ'] for r in rets]

        zdr_str = f"{np.mean(zdr_vals):.1f}±{np.std(zdr_vals):.1f}" if zdr_vals else "---"
        cc_str = f"{np.mean(cc_vals):.3f}" if cc_vals else "---"
        vel_str = f"{np.mean(vel_vals):.1f}" if vel_vals else "---"

        print(f"    {t} UTC: n={len(rets):>3}, ref={np.max(ref_vals):.0f}dBZ, "
              f"ZDR={zdr_str}, CC={cc_str}, Vel={vel_str}m/s")

# ============================================================================
# VELOCITY ANALYSIS — TARGET MOTION
# ============================================================================
print("\n\n" + "=" * 80)
print("DOPPLER VELOCITY ANALYSIS — RADIAL MOTION OF TARGET/DEBRIS")
print("=" * 80)

vel_returns = [r for r in all_significant_returns
               if r.get('vel_ms') is not None and r['ref_dBZ'] > 10]

if vel_returns:
    # Group by time
    from collections import defaultdict
    vel_by_time = defaultdict(list)
    for r in vel_returns:
        vel_by_time[r['time']].append(r)

    print(f"\n  Radial velocity is the component of motion toward/away from KEPZ.")
    print(f"  KEPZ is to the WEST of the target area.")
    print(f"  Positive velocity = moving AWAY from radar (eastward component)")
    print(f"  Negative velocity = moving TOWARD radar (westward component)")
    print()

    for t in sorted(vel_by_time.keys()):
        rets = vel_by_time[t]
        vels = [r['vel_ms'] for r in rets]
        refs = [r['ref_dBZ'] for r in rets]

        # Reflectivity-weighted mean velocity
        z_linear = np.array([10**(r/10) for r in refs])
        weighted_vel = np.sum(np.array(vels) * z_linear) / np.sum(z_linear)

        print(f"  {t} UTC: n={len(rets):>3}, mean_vel={np.mean(vels):>+6.2f} m/s, "
              f"weighted_vel={weighted_vel:>+6.2f} m/s, "
              f"spread={np.max(vels)-np.min(vels):.1f} m/s")

    all_vels = [r['vel_ms'] for r in vel_returns]
    print(f"\n  Overall radial velocity:")
    print(f"    Mean:     {np.mean(all_vels):+.2f} m/s")
    print(f"    Median:   {np.median(all_vels):+.2f} m/s")
    print(f"    Range:    {np.min(all_vels):.2f} to {np.max(all_vels):.2f} m/s")

    if np.mean(all_vels) > 1:
        print(f"\n  >>> NET EASTWARD MOTION (away from KEPZ)")
        print(f"  >>> Consistent with wind-driven debris drift (winds from SSE)")
    elif np.mean(all_vels) < -1:
        print(f"\n  >>> NET WESTWARD MOTION (toward KEPZ)")
        print(f"  >>> ANOMALOUS — moving AGAINST prevailing winds")
        print(f"  >>> Consistent with powered flight or active propulsion")
    else:
        print(f"\n  >>> NEAR-ZERO net radial motion")
        print(f"  >>> Target approximately stationary relative to KEPZ line of sight")
        print(f"  >>> Could be hovering, orbiting, or moving tangentially")

# ============================================================================
# SPATIAL ANALYSIS — HIGH-PRECISION POSITIONS
# ============================================================================
print("\n\n" + "=" * 80)
print("HIGH-PRECISION POSITION EXTRACTION (KEPZ @ 30 km)")
print("=" * 80)

print(f"\n  KEPZ gate spacing: 250m range × ~0.5° azimuth")
print(f"  At 30 km range: cross-range resolution = ~260m")
print(f"  This is ~10× better than KHDX at 135 km\n")

# Group returns by time and find centroids
from collections import defaultdict
pos_by_time = defaultdict(list)
for r in all_significant_returns:
    if r['ref_dBZ'] > 10:
        pos_by_time[r['time']].append(r)

print(f"  {'Time':>6} {'N':>4} {'Peak dBZ':>9} {'Centroid Lat':>13} {'Centroid Lon':>13} "
      f"{'Peak Lat':>10} {'Peak Lon':>11} {'Spread km':>10}")
print(f"  {'─'*6} {'─'*4} {'─'*9} {'─'*13} {'─'*13} {'─'*10} {'─'*11} {'─'*10}")

kepz_positions = []

for t in sorted(pos_by_time.keys()):
    rets = pos_by_time[t]
    if not rets:
        continue

    lats = np.array([r['lat'] for r in rets])
    lons = np.array([r['lon'] for r in rets])
    refs = np.array([r['ref_dBZ'] for r in rets])

    # Z-weighted centroid
    z_lin = 10 ** (refs / 10.0)
    c_lat = np.sum(lats * z_lin) / np.sum(z_lin)
    c_lon = np.sum(lons * z_lin) / np.sum(z_lin)

    peak = max(rets, key=lambda r: r['ref_dBZ'])

    # Spread
    spread = np.sqrt(
        ((np.max(lats) - np.min(lats)) * 111.0)**2 +
        ((np.max(lons) - np.min(lons)) * 111.0 * np.cos(np.radians(np.mean(lats))))**2
    )

    print(f"  {t:>6} {len(rets):>4} {peak['ref_dBZ']:>9.1f} {c_lat:>13.4f} {c_lon:>13.4f} "
          f"{peak['lat']:>10.4f} {peak['lon']:>11.4f} {spread:>10.1f}")

    kepz_positions.append({
        'time': t,
        'n_returns': len(rets),
        'centroid_lat': float(c_lat),
        'centroid_lon': float(c_lon),
        'peak_lat': peak['lat'],
        'peak_lon': peak['lon'],
        'peak_ref': peak['ref_dBZ'],
        'peak_az': peak['az'],
        'peak_range_km': peak['range_km'],
        'spread_km': float(spread),
    })

# ============================================================================
# SAVE RESULTS
# ============================================================================
output = {
    'parameters': {
        'radar': 'KEPZ',
        'lat': KEPZ_LAT,
        'lon': KEPZ_LON,
        'target_sector': {
            'az_min': TARGET_AZ_MIN,
            'az_max': TARGET_AZ_MAX,
            'range_min_km': TARGET_RANGE_MIN / 1000,
            'range_max_km': TARGET_RANGE_MAX / 1000,
        },
    },
    'scan_results': all_scan_results,
    'significant_returns': all_significant_returns,
    'kepz_positions': kepz_positions,
    'summary': {
        'total_scans': len(all_scan_results),
        'total_significant_returns': len(all_significant_returns),
        'total_dual_pol_returns': len(returns_with_dp) if returns_with_dp else 0,
        'max_ref_overall': max(sr['max_ref_dBZ'] for sr in all_scan_results),
        'scans_with_returns_gt10': sum(1 for sr in all_scan_results if sr['n_above_10dBZ'] > 0),
    },
}

output_path = os.path.join(OUTPUT_DIR, 'kepz_level2_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\n\nResults saved to {output_path}")
print(f"Total significant returns: {len(all_significant_returns)}")
print(f"Total dual-pol returns: {len(returns_with_dp) if returns_with_dp else 0}")
