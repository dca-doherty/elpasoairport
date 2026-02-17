#!/usr/bin/env python3
"""
KHDX (Holloman AFB) NEXRAD Analysis for El Paso Airspace Incident.
Cross-references KHDX radar data with KEPZ-detected spikes.

KHDX: 33.0803°N, 106.1225°W, elevation 1287m MSL
Target area (from KHDX): az ~189-192°, range ~135-137 km
"""

import numpy as np
import pyart
import os
import glob
import json
import warnings
warnings.filterwarnings('ignore')

# Parameters
KHDX_LAT = 33.0803
KHDX_LON = -106.1225
KHDX_ELEV = 1287.0

# Spike locations from KEPZ analysis
SPIKE1_LAT = 31.8771
SPIKE1_LON = -106.3423
SPIKE2_LAT = 31.8783
SPIKE2_LON = -106.4164

# Fort Bliss general area bounding box
FB_LAT_MIN = 31.82
FB_LAT_MAX = 31.93
FB_LON_MIN = -106.50
FB_LON_MAX = -106.30

# Spike times from KEPZ
SPIKE1_TIME_UTC = "04:24"  # UTC
SPIKE2_TIME_UTC = "05:13"  # UTC

def haversine(lat1, lon1, lat2, lon2):
    """Distance in meters."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return 6371000 * 2 * np.arcsin(np.sqrt(a))

def azimuth(lat1, lon1, lat2, lon2):
    """Azimuth in degrees."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(dlon)
    return np.degrees(np.arctan2(x, y)) % 360

# Compute target sector from KHDX perspective
for name, lat, lon in [("Spike 1", SPIKE1_LAT, SPIKE1_LON),
                        ("Spike 2", SPIKE2_LAT, SPIKE2_LON)]:
    dist = haversine(KHDX_LAT, KHDX_LON, lat, lon)
    az = azimuth(KHDX_LAT, KHDX_LON, lat, lon)
    print(f"{name}: az={az:.1f}°, range={dist/1000:.1f} km from KHDX")

# Compute Fort Bliss area sector
corners = [(FB_LAT_MIN, FB_LON_MIN), (FB_LAT_MIN, FB_LON_MAX),
           (FB_LAT_MAX, FB_LON_MIN), (FB_LAT_MAX, FB_LON_MAX)]
az_list = [azimuth(KHDX_LAT, KHDX_LON, lat, lon) for lat, lon in corners]
rng_list = [haversine(KHDX_LAT, KHDX_LON, lat, lon)/1000 for lat, lon in corners]
print(f"\nFort Bliss sector from KHDX:")
print(f"  Azimuth range: {min(az_list):.1f}° - {max(az_list):.1f}°")
print(f"  Range: {min(rng_list):.1f} - {max(rng_list):.1f} km")

TARGET_AZ_MIN = min(az_list) - 2  # Add margin
TARGET_AZ_MAX = max(az_list) + 2
TARGET_RANGE_MIN = (min(rng_list) - 5) * 1000  # meters
TARGET_RANGE_MAX = (max(rng_list) + 5) * 1000

print(f"  Search sector: az {TARGET_AZ_MIN:.1f}-{TARGET_AZ_MAX:.1f}°, "
      f"range {TARGET_RANGE_MIN/1000:.0f}-{TARGET_RANGE_MAX/1000:.0f} km")

# Scan KHDX files
data_dir = '/home/user/uap-transient-research/el_paso_airspace/khdx_data/'
files = sorted(glob.glob(os.path.join(data_dir, 'Level2_KHDX_*.ar2v')))
print(f"\nTotal KHDX files to analyze: {len(files)}")

results = []
significant_returns = []

print("\n" + "=" * 70)
print("KHDX NEXRAD SCAN-BY-SCAN ANALYSIS")
print("Fort Bliss sector (az {:.0f}-{:.0f}°, range {:.0f}-{:.0f} km)".format(
    TARGET_AZ_MIN, TARGET_AZ_MAX, TARGET_RANGE_MIN/1000, TARGET_RANGE_MAX/1000))
print("=" * 70)

for filepath in files:
    fname = os.path.basename(filepath)
    # Extract time
    time_str = fname.split('_')[-1].replace('.ar2v', '')

    try:
        radar = pyart.io.read_nexrad_archive(filepath)
    except Exception as e:
        print(f"\n{fname}: ERROR reading file - {e}")
        continue

    scan_time = radar.time['units'].replace('seconds since ', '')
    n_sweeps = radar.nsweeps

    # Get the lowest elevation sweep (0.5°)
    sweep_idx = 0
    sweep_start = radar.sweep_start_ray_index['data'][sweep_idx]
    sweep_end = radar.sweep_end_ray_index['data'][sweep_idx]

    # Get azimuth and range arrays
    azimuths = radar.azimuth['data'][sweep_start:sweep_end+1]
    ranges = radar.range['data']  # meters

    # Get reflectivity for the sweep
    ref_field = None
    for field_name in ['reflectivity', 'REF', 'ref']:
        if field_name in radar.fields:
            ref_field = field_name
            break

    if ref_field is None:
        print(f"\n{fname} ({time_str}): No reflectivity field found")
        continue

    ref_data = radar.fields[ref_field]['data'][sweep_start:sweep_end+1]

    # Filter to target sector
    az_mask = np.zeros(len(azimuths), dtype=bool)
    for i, az in enumerate(azimuths):
        if TARGET_AZ_MIN <= az <= TARGET_AZ_MAX:
            az_mask[i] = True

    rng_mask = (ranges >= TARGET_RANGE_MIN) & (ranges <= TARGET_RANGE_MAX)

    # Extract sector
    sector_ref = ref_data[az_mask][:, rng_mask]
    sector_az = azimuths[az_mask]
    sector_rng = ranges[rng_mask]

    if sector_ref.size == 0:
        print(f"\n{fname} ({time_str}): No data in target sector")
        continue

    # Find non-masked (valid) returns
    if hasattr(sector_ref, 'compressed'):
        valid_data = sector_ref.compressed()
    else:
        valid_data = sector_ref[~np.isnan(sector_ref)]

    # Count returns above threshold
    above_5 = np.sum(valid_data > 5.0) if len(valid_data) > 0 else 0
    above_10 = np.sum(valid_data > 10.0) if len(valid_data) > 0 else 0
    above_15 = np.sum(valid_data > 15.0) if len(valid_data) > 0 else 0
    max_ref = np.max(valid_data) if len(valid_data) > 0 else -999

    scan_result = {
        'file': fname,
        'time': time_str,
        'scan_time': scan_time,
        'n_valid_pixels': len(valid_data),
        'n_above_5dBZ': int(above_5),
        'n_above_10dBZ': int(above_10),
        'n_above_15dBZ': int(above_15),
        'max_ref_dBZ': float(max_ref),
    }

    results.append(scan_result)

    # Report
    status = ""
    if above_5 > 0:
        status = f" *** {above_5} pixels > 5 dBZ, max {max_ref:.1f} dBZ ***"
    if above_15 > 0:
        status = f" *** SIGNIFICANT: {above_15} pixels > 15 dBZ, max {max_ref:.1f} dBZ ***"

    print(f"\n{fname} ({time_str} UTC):")
    print(f"  Sweep 0 (0.5° elev), {len(valid_data)} valid pixels in sector")
    print(f"  >5 dBZ: {above_5}, >10 dBZ: {above_10}, >15 dBZ: {above_15}, Max: {max_ref:.1f} dBZ{status}")

    # If significant returns found, get positions
    if above_5 > 0:
        # Find positions of returns > 5 dBZ
        for i in range(sector_ref.shape[0]):
            for j in range(sector_ref.shape[1]):
                val = sector_ref[i, j]
                if hasattr(val, 'item'):
                    val = val.item() if not np.ma.is_masked(val) else -999
                if val > 5.0:
                    r = sector_rng[j]
                    a = sector_az[i]
                    # Convert to lat/lon
                    lat = KHDX_LAT + (r/1000) * np.cos(np.radians(a)) / 111.0
                    lon = KHDX_LON + (r/1000) * np.sin(np.radians(a)) / (111.0 * np.cos(np.radians(KHDX_LAT)))
                    significant_returns.append({
                        'time': time_str,
                        'az': float(a),
                        'range_km': float(r/1000),
                        'ref_dBZ': float(val),
                        'lat': float(lat),
                        'lon': float(lon),
                    })
                    if val > 10:
                        print(f"    Return at az={a:.1f}°, range={r/1000:.1f}km, "
                              f"ref={val:.1f}dBZ -> ~{lat:.4f}°N, {lon:.4f}°W")

    # Also check other dual-pol fields if available
    for dpfield in ['differential_reflectivity', 'cross_correlation_ratio',
                     'ZDR', 'zdr', 'CC', 'RHO', 'rho']:
        if dpfield in radar.fields and above_5 > 0:
            dp_data = radar.fields[dpfield]['data'][sweep_start:sweep_end+1]
            dp_sector = dp_data[az_mask][:, rng_mask]
            if hasattr(dp_sector, 'compressed'):
                dp_valid = dp_sector.compressed()
            else:
                dp_valid = dp_sector[~np.isnan(dp_sector)]
            if len(dp_valid) > 0:
                print(f"    {dpfield}: mean={np.mean(dp_valid):.2f}, "
                      f"min={np.min(dp_valid):.2f}, max={np.max(dp_valid):.2f}")

    del radar  # Free memory

# ============= SUMMARY =============
print("\n" + "=" * 70)
print("KHDX ANALYSIS SUMMARY")
print("=" * 70)

# Check around spike times
print("\nScans around Spike 1 time (04:24 UTC):")
for r in results:
    t = r['time']
    if t.startswith('04'):
        h, m = int(t[:2]), int(t[2:])
        if 15 <= m <= 40:
            print(f"  {t} UTC: max={r['max_ref_dBZ']:.1f} dBZ, "
                  f"pixels>5={r['n_above_5dBZ']}, pixels>10={r['n_above_10dBZ']}")

print("\nScans around Spike 2 time (05:13 UTC):")
for r in results:
    t = r['time']
    if t.startswith('05'):
        h, m = int(t[:2]), int(t[2:])
        if 5 <= m <= 25:
            print(f"  {t} UTC: max={r['max_ref_dBZ']:.1f} dBZ, "
                  f"pixels>5={r['n_above_5dBZ']}, pixels>10={r['n_above_10dBZ']}")

# Overall statistics
total_sig = sum(1 for r in results if r['n_above_5dBZ'] > 0)
print(f"\nTotal scans analyzed: {len(results)}")
print(f"Scans with returns > 5 dBZ in target sector: {total_sig}")
print(f"Total significant returns (> 5 dBZ): {len(significant_returns)}")

if len(significant_returns) > 0:
    print("\nAll significant returns:")
    for sr in significant_returns:
        print(f"  {sr['time']} UTC: {sr['ref_dBZ']:.1f} dBZ at "
              f"az={sr['az']:.1f}°, range={sr['range_km']:.1f}km "
              f"-> ({sr['lat']:.4f}°N, {sr['lon']:.4f}°W)")
else:
    print("\n*** NO significant returns detected from KHDX in the Fort Bliss sector ***")
    print("This means either:")
    print("  1. The target was below KHDX 0.5° beam height (~1200-1300m AGL at 135km)")
    print("  2. The target had insufficient RCS to be detected at 135km range")
    print("  3. There was no airborne target (ground clutter, sidelobe, etc.)")
    print("  4. The target was present but too small/low for KHDX to detect")
    print()
    print("Given KHDX beam height of ~2500m MSL (~1300m AGL) at this range,")
    print("a low-altitude target (< 500m AGL) would be invisible to KHDX.")
    print("This is CONSISTENT with either a low-flying drone or a balloon below 500m AGL.")

# Save results
output = {
    'scan_results': results,
    'significant_returns': significant_returns,
    'target_sector': {
        'az_min': TARGET_AZ_MIN,
        'az_max': TARGET_AZ_MAX,
        'range_min_km': TARGET_RANGE_MIN/1000,
        'range_max_km': TARGET_RANGE_MAX/1000,
    }
}
with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/khdx_analysis_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to khdx_analysis_results.json")
