#!/usr/bin/env python3
"""
GLM Feb 10 Detection Drill-Down
================================
Focused re-analysis of the 103 GLM detections found on Feb 10, 2026
in the 02-08 UTC extended window near El Paso / Fort Bliss.

Questions to answer:
1. Where do the detections geolocate? Clustered at Fort Bliss or scattered?
2. What time window do they fall in?
3. Do any coincide with the 03:43 UTC SkySentinel observation?
4. What are their energies and characteristics?
"""

import s3fs
import xarray as xr
import numpy as np
import json
import os
from io import BytesIO
from datetime import datetime, timedelta
from collections import defaultdict

# Fort Bliss target area
TARGET_LAT, TARGET_LON = 31.87, -106.52
SEARCH_RADIUS_DEG = 1.0  # Same as original analysis

# GOES-19 S3
BUCKET = 'noaa-goes19'
PRODUCT = 'GLM-L2-LCFA'

# SkySentinel last observation
SKYSENTINEL_TIME = datetime(2026, 2, 10, 3, 43)

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'

def list_glm_files(fs, start_dt, end_dt):
    files = []
    current = start_dt
    while current <= end_dt:
        day_of_year = current.timetuple().tm_yday
        prefix = f"{BUCKET}/{PRODUCT}/{current.year}/{day_of_year:03d}/{current.hour:02d}/"
        try:
            hour_files = fs.ls(prefix)
            hour_files = [f for f in hour_files if f.endswith('.nc')]
            files.extend(hour_files)
        except FileNotFoundError:
            pass
        current += timedelta(hours=1)
    return sorted(files)

def extract_timestamp(filepath):
    basename = os.path.basename(filepath)
    try:
        parts = basename.split('_')
        for p in parts:
            if p.startswith('s'):
                year = int(p[1:5])
                doy = int(p[5:8])
                hour = int(p[8:10])
                minute = int(p[10:12])
                second = int(p[12:14])
                dt = datetime(year, 1, 1) + timedelta(days=doy-1, hours=hour, minutes=minute, seconds=second)
                return dt
    except (ValueError, IndexError):
        pass
    return None

def haversine_km(lat1, lon1, lat2, lon2):
    """Distance in km between two lat/lon points."""
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat/2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(a))

print("=" * 80)
print("GLM FEB 10 DETECTION DRILL-DOWN")
print("Analyzing all detections within 1° of Fort Bliss on Feb 10 02-08 UTC")
print("=" * 80)

fs = s3fs.S3FileSystem(anon=True)

# Feb 10, 02-08 UTC (the window where 103 detections were found)
start = datetime(2026, 2, 10, 2, 0)
end = datetime(2026, 2, 10, 8, 0)

files = list_glm_files(fs, start, end)
print(f"\nFound {len(files)} GLM files for Feb 10 02-08 UTC")

all_detections = []
errors = 0
files_with_detections = 0

for i, fpath in enumerate(files):
    if i % 100 == 0:
        print(f"  Processing file {i+1}/{len(files)}...")

    ts = extract_timestamp(fpath)

    try:
        with fs.open(fpath, 'rb') as f:
            raw = f.read()
        ds = xr.open_dataset(BytesIO(raw), engine='h5netcdf')
    except Exception as e:
        errors += 1
        continue

    file_has_detection = False

    # Check all three detection levels
    for det_type, lat_var, lon_var, extra_vars in [
        ('flash', 'flash_lat', 'flash_lon', ['flash_energy', 'flash_area', 'flash_time_offset_of_first_event', 'flash_time_offset_of_last_event']),
        ('group', 'group_lat', 'group_lon', ['group_energy', 'group_area', 'group_time_offset']),
        ('event', 'event_lat', 'event_lon', ['event_energy', 'event_time_offset']),
    ]:
        if lat_var in ds and len(ds[lat_var]) > 0:
            lats = ds[lat_var].values
            lons = ds[lon_var].values

            mask = ((np.abs(lats - TARGET_LAT) < SEARCH_RADIUS_DEG) &
                    (np.abs(lons - TARGET_LON) < SEARCH_RADIUS_DEG))

            if np.any(mask):
                file_has_detection = True
                indices = np.where(mask)[0]

                for idx in indices:
                    det = {
                        'type': det_type,
                        'lat': float(lats[idx]),
                        'lon': float(lons[idx]),
                        'file_time': ts.strftime('%Y-%m-%d %H:%M:%S') if ts else 'unknown',
                        'file': os.path.basename(fpath),
                        'distance_km': float(haversine_km(TARGET_LAT, TARGET_LON, float(lats[idx]), float(lons[idx]))),
                    }

                    # Extract extra variables
                    for var in extra_vars:
                        if var in ds:
                            try:
                                val = ds[var].values[idx]
                                det[var] = float(val) if not np.isnan(val) else None
                            except (IndexError, TypeError):
                                pass

                    # Calculate precise time from file timestamp + offset
                    offset_var = None
                    if det_type == 'event' and 'event_time_offset' in ds:
                        offset_var = 'event_time_offset'
                    elif det_type == 'group' and 'group_time_offset' in ds:
                        offset_var = 'group_time_offset'
                    elif det_type == 'flash' and 'flash_time_offset_of_first_event' in ds:
                        offset_var = 'flash_time_offset_of_first_event'

                    if offset_var and offset_var in ds:
                        try:
                            # Time offset is in seconds from product_time
                            if 'product_time' in ds:
                                prod_time = ds['product_time'].values
                                det['product_time'] = str(prod_time)
                            offset_val = ds[offset_var].values[idx]
                            det['time_offset_s'] = float(offset_val) if not np.isnan(offset_val) else None
                        except (IndexError, TypeError):
                            pass

                    all_detections.append(det)

    if file_has_detection:
        files_with_detections += 1

    ds.close()

print(f"\n{'='*80}")
print(f"RESULTS: {len(all_detections)} total detections in {files_with_detections} files ({errors} errors)")
print(f"{'='*80}")

# Break down by type
by_type = defaultdict(list)
for d in all_detections:
    by_type[d['type']].append(d)

print(f"\n--- Detection Type Breakdown ---")
for dtype in ['flash', 'group', 'event']:
    print(f"  {dtype}s: {len(by_type[dtype])}")

# Break down by hour
print(f"\n--- Temporal Distribution ---")
by_hour = defaultdict(int)
for d in all_detections:
    hour = d['file_time'][:13] if d['file_time'] != 'unknown' else 'unknown'
    by_hour[hour] += 1

for hour in sorted(by_hour.keys()):
    marker = " <-- SkySentinel 03:43 UTC" if "03" in hour else ""
    print(f"  {hour} UTC: {by_hour[hour]} detections{marker}")

# Distance from Fort Bliss
print(f"\n--- Distance from Fort Bliss (31.87°N, 106.52°W) ---")
distances = [d['distance_km'] for d in all_detections]
if distances:
    print(f"  Min: {min(distances):.1f} km")
    print(f"  Max: {max(distances):.1f} km")
    print(f"  Mean: {np.mean(distances):.1f} km")
    print(f"  Median: {np.median(distances):.1f} km")

    within_20km = sum(1 for d in distances if d < 20)
    within_50km = sum(1 for d in distances if d < 50)
    beyond_50km = sum(1 for d in distances if d >= 50)
    print(f"  Within 20km of Fort Bliss: {within_20km}")
    print(f"  Within 50km of Fort Bliss: {within_50km}")
    print(f"  Beyond 50km: {beyond_50km}")

# Geographic clustering
print(f"\n--- Geographic Distribution ---")
lats = [d['lat'] for d in all_detections]
lons = [d['lon'] for d in all_detections]
if lats:
    print(f"  Lat range: {min(lats):.3f} to {max(lats):.3f} (span: {max(lats)-min(lats):.3f}°)")
    print(f"  Lon range: {min(lons):.3f} to {max(lons):.3f} (span: {max(lons)-min(lons):.3f}°)")
    print(f"  Fort Bliss: 31.87°N, -106.52°W")

    # Check for clustering
    if max(lats) - min(lats) < 0.1 and max(lons) - min(lons) < 0.1:
        print(f"  >>> TIGHTLY CLUSTERED (within ~10km box)")
    elif max(lats) - min(lats) < 0.5 and max(lons) - min(lons) < 0.5:
        print(f"  >>> MODERATELY CLUSTERED (within ~50km box)")
    else:
        print(f"  >>> SCATTERED across search area")

# SkySentinel coincidence check
print(f"\n--- SkySentinel Coincidence Check ---")
print(f"  SkySentinel last observation: {SKYSENTINEL_TIME} UTC")
near_skysentinel = [d for d in all_detections
                    if d['file_time'] != 'unknown'
                    and abs((datetime.strptime(d['file_time'], '%Y-%m-%d %H:%M:%S') - SKYSENTINEL_TIME).total_seconds()) < 600]
print(f"  Detections within ±10 minutes of 03:43 UTC: {len(near_skysentinel)}")

near_skysentinel_30 = [d for d in all_detections
                       if d['file_time'] != 'unknown'
                       and abs((datetime.strptime(d['file_time'], '%Y-%m-%d %H:%M:%S') - SKYSENTINEL_TIME).total_seconds()) < 1800]
print(f"  Detections within ±30 minutes of 03:43 UTC: {len(near_skysentinel_30)}")

if near_skysentinel:
    print(f"\n  COINCIDENT DETECTIONS:")
    for d in near_skysentinel:
        print(f"    {d['type']} at {d['file_time']} UTC: ({d['lat']:.3f}°N, {d['lon']:.3f}°W) "
              f"dist={d['distance_km']:.1f}km")

# Energy analysis
print(f"\n--- Energy Analysis ---")
for dtype in ['flash', 'group', 'event']:
    dets = by_type[dtype]
    energies = [d.get(f'{dtype}_energy') for d in dets if d.get(f'{dtype}_energy') is not None]
    if energies:
        print(f"  {dtype} energies: {min(energies):.2e} to {max(energies):.2e} J "
              f"(mean: {np.mean(energies):.2e})")

# Print ALL detections for full detail
print(f"\n--- ALL DETECTIONS (sorted by time) ---")
for d in sorted(all_detections, key=lambda x: x.get('file_time', '')):
    energy_key = f"{d['type']}_energy"
    energy = d.get(energy_key, 'N/A')
    energy_str = f"{energy:.2e}" if isinstance(energy, float) else str(energy)
    print(f"  {d['file_time']} | {d['type']:>5} | ({d['lat']:8.3f}, {d['lon']:9.3f}) | "
          f"dist={d['distance_km']:6.1f}km | E={energy_str}")

# Save detailed results
output = {
    'analysis': 'GLM Feb 10 Drill-Down',
    'search': {
        'center': [TARGET_LAT, TARGET_LON],
        'radius_deg': SEARCH_RADIUS_DEG,
        'time_window': '2026-02-10 02:00-08:00 UTC',
    },
    'summary': {
        'total_detections': len(all_detections),
        'files_processed': len(files),
        'files_with_detections': files_with_detections,
        'by_type': {k: len(v) for k, v in by_type.items()},
    },
    'temporal': dict(by_hour),
    'spatial': {
        'lat_range': [min(lats), max(lats)] if lats else None,
        'lon_range': [min(lons), max(lons)] if lons else None,
        'distance_stats_km': {
            'min': float(min(distances)) if distances else None,
            'max': float(max(distances)) if distances else None,
            'mean': float(np.mean(distances)) if distances else None,
            'median': float(np.median(distances)) if distances else None,
        },
        'within_20km': sum(1 for d in distances if d < 20) if distances else 0,
        'within_50km': sum(1 for d in distances if d < 50) if distances else 0,
    },
    'skysentinel_coincidence': {
        'skysentinel_time': str(SKYSENTINEL_TIME),
        'within_10min': len(near_skysentinel),
        'within_30min': len(near_skysentinel_30),
        'coincident_detections': near_skysentinel,
    },
    'all_detections': all_detections,
}

output_path = os.path.join(OUTPUT_DIR, 'glm_feb10_drilldown.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nDetailed results saved to {output_path}")
