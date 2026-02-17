#!/usr/bin/env python3
"""
GOES-19 GLM (Geostationary Lightning Mapper) Analysis
=====================================================
Checks for optical flash detections near Fort Bliss / El Paso during the
Feb 11, 2026 incident window. GLM images at ~500fps in near-IR (777.4nm)
and has detected non-lightning events like bolides.

If a HELWS laser engagement produced any optical flash (plasma spark,
target combustion, scattered NIR), GLM could capture it.

Data source: NOAA GOES-19 on AWS S3 (noaa-goes19 bucket, GLM-L2-LCFA)
Note: GOES-16 was replaced by GOES-19 as GOES-East on April 7, 2025.
"""

import s3fs
import xarray as xr
import numpy as np
import json
import os
from io import BytesIO
from datetime import datetime, timedelta
from collections import defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================
# Fort Bliss engagement area
TARGET_LAT, TARGET_LON = 31.87, -106.52
SEARCH_RADIUS_DEG = 1.0  # ~111 km search box

# Time windows (UTC)
EVENT_START = datetime(2026, 2, 11, 3, 0)
EVENT_END = datetime(2026, 2, 11, 7, 0)
BASELINE_START = datetime(2026, 2, 10, 3, 0)
BASELINE_END = datetime(2026, 2, 10, 7, 0)

# Extended window for booming sound reports
EXTENDED_START = datetime(2026, 2, 8, 0, 0)
EXTENDED_END = datetime(2026, 2, 12, 12, 0)

# S3 bucket — GOES-19 replaced GOES-16 as GOES-East on April 7, 2025
BUCKET = 'noaa-goes19'
PRODUCT = 'GLM-L2-LCFA'

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/glm_data/'
ANALYSIS_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

def list_glm_files(fs, start_dt, end_dt):
    """List GLM LCFA files for a time range on S3."""
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

def analyze_glm_file(fs, filepath):
    """Extract flash/event/group detections near target area from a GLM file."""
    detections = {'flashes': [], 'events': [], 'groups': []}

    try:
        with fs.open(filepath, 'rb') as f:
            raw = f.read()
        ds = xr.open_dataset(BytesIO(raw), engine='h5netcdf')

        # Check for flashes (highest-level product)
        if 'flash_lat' in ds and len(ds.flash_lat) > 0:
            flash_lats = ds.flash_lat.values
            flash_lons = ds.flash_lon.values

            # Spatial filter
            mask = ((np.abs(flash_lats - TARGET_LAT) < SEARCH_RADIUS_DEG) &
                    (np.abs(flash_lons - TARGET_LON) < SEARCH_RADIUS_DEG))

            if np.any(mask):
                flash_energy = ds.flash_energy.values[mask] if 'flash_energy' in ds else [None] * mask.sum()
                flash_area = ds.flash_area.values[mask] if 'flash_area' in ds else [None] * mask.sum()

                for i in range(mask.sum()):
                    idx = np.where(mask)[0][i]
                    det = {
                        'lat': float(flash_lats[idx]),
                        'lon': float(flash_lons[idx]),
                        'energy_J': float(flash_energy[i]) if flash_energy[i] is not None else None,
                        'area_km2': float(flash_area[i]) if flash_area[i] is not None else None,
                    }
                    detections['flashes'].append(det)

        # Check for events (lowest-level, individual pixel detections)
        if 'event_lat' in ds and len(ds.event_lat) > 0:
            event_lats = ds.event_lat.values
            event_lons = ds.event_lon.values

            mask = ((np.abs(event_lats - TARGET_LAT) < SEARCH_RADIUS_DEG) &
                    (np.abs(event_lons - TARGET_LON) < SEARCH_RADIUS_DEG))

            if np.any(mask):
                event_energy = ds.event_energy.values[mask] if 'event_energy' in ds else [None] * mask.sum()
                for i in range(mask.sum()):
                    idx = np.where(mask)[0][i]
                    det = {
                        'lat': float(event_lats[idx]),
                        'lon': float(event_lons[idx]),
                        'energy_J': float(event_energy[i]) if event_energy[i] is not None else None,
                    }
                    detections['events'].append(det)

        # Check for groups (intermediate level)
        if 'group_lat' in ds and len(ds.group_lat) > 0:
            group_lats = ds.group_lat.values
            group_lons = ds.group_lon.values

            mask = ((np.abs(group_lats - TARGET_LAT) < SEARCH_RADIUS_DEG) &
                    (np.abs(group_lons - TARGET_LON) < SEARCH_RADIUS_DEG))

            if np.any(mask):
                group_energy = ds.group_energy.values[mask] if 'group_energy' in ds else [None] * mask.sum()
                group_area = ds.group_area.values[mask] if 'group_area' in ds else [None] * mask.sum()
                for i in range(mask.sum()):
                    idx = np.where(mask)[0][i]
                    det = {
                        'lat': float(group_lats[idx]),
                        'lon': float(group_lons[idx]),
                        'energy_J': float(group_energy[i]) if group_energy[i] is not None else None,
                        'area_km2': float(group_area[i]) if group_area[i] is not None else None,
                    }
                    detections['groups'].append(det)

        ds.close()

    except Exception as e:
        return None, str(e)

    return detections, None

def extract_timestamp_from_filename(filepath):
    """Extract timestamp from GLM filename like OR_GLM-L2-LCFA_G16_sYYYYDDDHHMMSSS..."""
    basename = os.path.basename(filepath)
    try:
        # Find the start time marker 's'
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

# ============================================================================
# MAIN ANALYSIS
# ============================================================================
print("=" * 80)
print("GOES-16 GLM OPTICAL FLASH ANALYSIS")
print("Searching for NIR transients near Fort Bliss / El Paso")
print("=" * 80)

# Connect to S3 (anonymous access)
fs = s3fs.S3FileSystem(anon=True)

# ============================================================================
# PHASE 1: Event night (Feb 11, 03-07 UTC)
# ============================================================================
print(f"\n--- PHASE 1: Event Night ---")
print(f"Window: {EVENT_START} to {EVENT_END} UTC")

event_files = list_glm_files(fs, EVENT_START, EVENT_END)
print(f"Found {len(event_files)} GLM files")

event_detections = {'flashes': [], 'events': [], 'groups': []}
event_file_count = 0
event_errors = 0
event_by_hour = defaultdict(lambda: {'flashes': 0, 'events': 0, 'groups': 0})

for i, fpath in enumerate(event_files):
    if i % 50 == 0:
        print(f"  Processing file {i+1}/{len(event_files)}...")

    dets, err = analyze_glm_file(fs, fpath)
    if err:
        event_errors += 1
        continue

    event_file_count += 1
    ts = extract_timestamp_from_filename(fpath)
    hour_key = ts.strftime('%H:%M') if ts else 'unknown'

    for det_type in ['flashes', 'events', 'groups']:
        if dets[det_type]:
            for d in dets[det_type]:
                d['time'] = hour_key
                d['file'] = os.path.basename(fpath)
            event_detections[det_type].extend(dets[det_type])
            event_by_hour[hour_key][det_type] += len(dets[det_type])

print(f"\n  Event Night Results:")
print(f"    Files processed: {event_file_count} ({event_errors} errors)")
print(f"    Flashes detected: {len(event_detections['flashes'])}")
print(f"    Groups detected:  {len(event_detections['groups'])}")
print(f"    Events detected:  {len(event_detections['events'])}")

if event_detections['flashes']:
    energies = [d['energy_J'] for d in event_detections['flashes'] if d.get('energy_J')]
    if energies:
        print(f"    Flash energy range: {min(energies):.2e} - {max(energies):.2e} J")
    for d in event_detections['flashes'][:10]:
        print(f"      Flash at {d.get('time','?')} UTC: ({d['lat']:.3f}, {d['lon']:.3f}) E={d.get('energy_J','?')} J")

# ============================================================================
# PHASE 2: Baseline night (Feb 10, 03-07 UTC)
# ============================================================================
print(f"\n--- PHASE 2: Baseline Night ---")
print(f"Window: {BASELINE_START} to {BASELINE_END} UTC")

baseline_files = list_glm_files(fs, BASELINE_START, BASELINE_END)
print(f"Found {len(baseline_files)} GLM files")

baseline_detections = {'flashes': [], 'events': [], 'groups': []}
baseline_file_count = 0
baseline_errors = 0

for i, fpath in enumerate(baseline_files):
    if i % 50 == 0:
        print(f"  Processing file {i+1}/{len(baseline_files)}...")

    dets, err = analyze_glm_file(fs, fpath)
    if err:
        baseline_errors += 1
        continue

    baseline_file_count += 1
    ts = extract_timestamp_from_filename(fpath)

    for det_type in ['flashes', 'events', 'groups']:
        if dets[det_type]:
            for d in dets[det_type]:
                d['time'] = ts.strftime('%H:%M') if ts else 'unknown'
            baseline_detections[det_type].extend(dets[det_type])

print(f"\n  Baseline Night Results:")
print(f"    Files processed: {baseline_file_count} ({baseline_errors} errors)")
print(f"    Flashes detected: {len(baseline_detections['flashes'])}")
print(f"    Groups detected:  {len(baseline_detections['groups'])}")
print(f"    Events detected:  {len(baseline_detections['events'])}")

# ============================================================================
# PHASE 3: Extended window for booming sounds (Feb 8-12)
# ============================================================================
print(f"\n--- PHASE 3: Extended Window (Feb 8-12, nighttime hours only) ---")

extended_detections = {'flashes': [], 'events': [], 'groups': []}
extended_by_date = defaultdict(lambda: {'flashes': 0, 'events': 0, 'groups': 0})

for day_offset in range(5):  # Feb 8, 9, 10, 11, 12
    date = datetime(2026, 2, 8) + timedelta(days=day_offset)
    night_start = date.replace(hour=2)  # 02 UTC = 8 PM MST
    night_end = date.replace(hour=8)    # 08 UTC = 2 AM MST

    date_label = date.strftime('%b %d')
    if day_offset <= 1:
        date_label += " (booming reports)"
    elif day_offset == 3:
        date_label += " (EVENT)"

    night_files = list_glm_files(fs, night_start, night_end)
    date_key = date.strftime('%Y-%m-%d')
    n_det = 0

    for fpath in night_files:
        dets, err = analyze_glm_file(fs, fpath)
        if err:
            continue
        for det_type in ['flashes', 'events', 'groups']:
            if dets[det_type]:
                n = len(dets[det_type])
                extended_by_date[date_key][det_type] += n
                n_det += n
                for d in dets[det_type]:
                    ts = extract_timestamp_from_filename(fpath)
                    d['time'] = ts.strftime('%Y-%m-%d %H:%M:%S') if ts else 'unknown'
                    d['date'] = date_key
                extended_detections[det_type].extend(dets[det_type])

    print(f"  {date_label}: {len(night_files)} files, {n_det} total detections")

# ============================================================================
# COMPARISON & VERDICT
# ============================================================================
print(f"\n\n{'='*80}")
print(f"GLM ANALYSIS RESULTS")
print(f"{'='*80}")

print(f"\n  EVENT NIGHT (Feb 11):")
print(f"    Flashes: {len(event_detections['flashes'])}")
print(f"    Events:  {len(event_detections['events'])}")
print(f"    Groups:  {len(event_detections['groups'])}")

print(f"\n  BASELINE NIGHT (Feb 10):")
print(f"    Flashes: {len(baseline_detections['flashes'])}")
print(f"    Events:  {len(baseline_detections['events'])}")
print(f"    Groups:  {len(baseline_detections['groups'])}")

# Ratio
event_total = sum(len(v) for v in event_detections.values())
baseline_total = sum(len(v) for v in baseline_detections.values())

if baseline_total > 0:
    ratio = event_total / baseline_total
    print(f"\n  Event/baseline ratio: {ratio:.1f}x")
    if ratio > 2.0:
        print(f"  *** SIGNIFICANT EXCESS on event night ***")
    elif ratio > 1.3:
        print(f"  ** Moderately elevated on event night **")
    else:
        print(f"  No significant difference")
elif event_total > 0:
    print(f"\n  *** Detections ONLY on event night, NONE on baseline ***")
else:
    print(f"\n  No detections on either night (GLM did not see optical flashes)")
    print(f"  This could mean:")
    print(f"    - HELWS engagement did not produce bright optical flashes in GLM band")
    print(f"    - Laser wavelength outside GLM passband (777.4 nm OI line)")
    print(f"    - Flash too dim/brief for GLM threshold")
    print(f"    - Cloud cover obscured the area (check METAR)")

print(f"\n  MULTI-DATE OVERVIEW:")
print(f"  {'Date':<16} {'Flashes':>8} {'Groups':>8} {'Events':>8} {'Total':>8}")
print(f"  {'─'*16} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")
for date_key in sorted(extended_by_date.keys()):
    d = extended_by_date[date_key]
    total = d['flashes'] + d['groups'] + d['events']
    marker = " ***" if '02-11' in date_key else ""
    print(f"  {date_key:<16} {d['flashes']:>8} {d['groups']:>8} {d['events']:>8} {total:>8}{marker}")

# If any flashes found, detail them
all_flashes = event_detections['flashes'] + baseline_detections['flashes']
if all_flashes:
    print(f"\n  FLASH DETAIL (all detections within {SEARCH_RADIUS_DEG}° of target):")
    for fl in sorted(all_flashes, key=lambda x: x.get('time', '')):
        print(f"    {fl.get('time','?')} UTC: ({fl['lat']:.3f}°N, {fl['lon']:.3f}°W) "
              f"E={fl.get('energy_J','?')} J, A={fl.get('area_km2','?')} km²")

# ============================================================================
# SAVE RESULTS
# ============================================================================
output = {
    'search_area': {
        'center_lat': TARGET_LAT, 'center_lon': TARGET_LON,
        'radius_deg': SEARCH_RADIUS_DEG
    },
    'event_window': {'start': str(EVENT_START), 'end': str(EVENT_END)},
    'baseline_window': {'start': str(BASELINE_START), 'end': str(BASELINE_END)},
    'event_night': {
        'files_processed': event_file_count,
        'flashes': len(event_detections['flashes']),
        'groups': len(event_detections['groups']),
        'events': len(event_detections['events']),
        'flash_details': event_detections['flashes'][:50],
        'hourly_breakdown': dict(event_by_hour),
    },
    'baseline_night': {
        'files_processed': baseline_file_count,
        'flashes': len(baseline_detections['flashes']),
        'groups': len(baseline_detections['groups']),
        'events': len(baseline_detections['events']),
    },
    'extended_by_date': dict(extended_by_date),
    'verdict': 'excess' if event_total > baseline_total * 2 else 'elevated' if event_total > baseline_total * 1.3 else 'event_only' if event_total > 0 and baseline_total == 0 else 'no_detection' if event_total == 0 else 'normal'
}

output_path = os.path.join(ANALYSIS_DIR, 'goes16_glm_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")
