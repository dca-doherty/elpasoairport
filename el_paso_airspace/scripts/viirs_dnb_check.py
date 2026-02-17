#!/usr/bin/env python3
"""
VIIRS Day/Night Band Check
============================
Checks Suomi NPP / NOAA-20 VIIRS Day/Night Band for anomalous bright
point sources at Fort Bliss during the incident period.

VIIRS DNB is extremely sensitive to nighttime lights and has detected
gas flares, fires, aurora, and other anomalous sources.

Key consideration: VIIRS is polar-orbiting with ~1:30 AM local overpass
(ascending node). For El Paso (MST = UTC-7), the ascending pass is at
~08:30 UTC — which is AFTER our 03-07 UTC event window. However, the
descending pass (~1:30 PM local / 20:30 UTC) won't show nighttime data.

So VIIRS DNB timing may NOT overlap the engagement window directly,
but could show residual thermal/luminous effects.

Data sources:
  - NASA LAADS DAAC (VNP46A1 daily nighttime product)
  - NASA Worldview (visual inspection)
  - NOAA CLASS
"""

import requests
import json
import os
from datetime import datetime

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/viirs_data/'
ANALYSIS_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Fort Bliss target
TARGET_LAT, TARGET_LON = 31.87, -106.52

print("=" * 80)
print("VIIRS DAY/NIGHT BAND CHECK")
print("=" * 80)

# ============================================================================
# TIMING ANALYSIS
# ============================================================================
print(f"\n--- VIIRS Overpass Timing ---")
print(f"  Suomi NPP ascending node: ~13:30 local solar time")
print(f"  El Paso (MST, UTC-7): ascending pass ~13:30 MST = 20:30 UTC")
print(f"  Nighttime descending: ~01:30 MST = 08:30 UTC")
print(f"")
print(f"  Event window: 03:00-07:00 UTC (8 PM - midnight MST)")
print(f"  Nearest VIIRS nighttime pass: ~08:30 UTC (1:30 AM MST)")
print(f"")
print(f"  TIMING GAP: VIIRS overpass is ~1.5-5.5 hours AFTER event window")
print(f"  The DNB would NOT directly image the engagement in real-time.")
print(f"  But could detect: residual fires, persistent luminous effects,")
print(f"  or anomalous ground-level light patterns.")

# ============================================================================
# Check NASA Worldview for visual reference
# ============================================================================
print(f"\n--- NASA Worldview URLs ---")
worldview_base = "https://worldview.earthdata.nasa.gov"
for date in ['2026-02-10', '2026-02-11', '2026-02-12']:
    url = (f"{worldview_base}/?v=-108,30.5,-105,33"
           f"&l=VIIRS_SNPP_DayNightBand_At_Sensor_Radiance"
           f"&t={date}")
    print(f"  {date}: {url}")

# ============================================================================
# Check LAADS DAAC for VNP46A1 (Black Marble nighttime)
# ============================================================================
print(f"\n--- LAADS DAAC VNP46A1 (Black Marble Daily) ---")

# VNP46A1 is tiled in sinusoidal grid. El Paso is approximately h08v05
# LAADS DAAC search API
laads_url = "https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/details"

for product, name in [
    ('VNP46A1', 'VIIRS/NPP Daily Gridded DNB'),
    ('VNP46A2', 'VIIRS/NPP Gap-Filled DNB'),
]:
    print(f"\n  Checking {product} ({name})...")

    # LAADS uses collection/product/year/day_of_year structure
    for date_str, doy in [('2026-02-10', 41), ('2026-02-11', 42), ('2026-02-12', 43)]:
        params = {
            'product': product,
            'collection': '5200',  # Collection 2
            'dateRanges': date_str,
            'areaOfInterest': f'x{TARGET_LON}y{TARGET_LAT}',
        }

        search_url = f"https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/details/allData/5200/{product}/2026/{doy:03d}"

        try:
            resp = requests.get(search_url, timeout=15,
                                headers={'User-Agent': 'Mozilla/5.0'})
            if resp.status_code == 200:
                try:
                    files = resp.json()
                    if isinstance(files, list):
                        # Filter for our tile (h08v05 for El Paso area)
                        our_files = [f for f in files
                                     if isinstance(f, dict) and 'h08v05' in f.get('name', '')]
                        if our_files:
                            print(f"    {date_str}: {len(our_files)} files for tile h08v05")
                            for finfo in our_files[:3]:
                                print(f"      {finfo.get('name', '?')}")
                        else:
                            all_files_with_names = [f for f in files if isinstance(f, dict) and 'name' in f]
                            print(f"    {date_str}: {len(all_files_with_names)} total files (h08v05 not found)")
                            # Show what tiles are available
                            tiles = set()
                            for f in all_files_with_names[:50]:
                                name = f.get('name', '')
                                # Extract tile from filename like VNP46A1.A2026042.h08v05.002.xxx.h5
                                parts = name.split('.')
                                for p in parts:
                                    if p.startswith('h') and 'v' in p and len(p) <= 7:
                                        tiles.add(p)
                            if tiles:
                                print(f"      Available tiles: {sorted(tiles)[:10]}")
                    else:
                        print(f"    {date_str}: unexpected response format")
                except json.JSONDecodeError:
                    print(f"    {date_str}: non-JSON response")
            else:
                print(f"    {date_str}: HTTP {resp.status_code}")
        except Exception as e:
            print(f"    {date_str}: Error - {e}")

# ============================================================================
# Check NOAA VIIRS Nighttime Lights (EOG)
# ============================================================================
print(f"\n--- EOG/Mines.edu VIIRS Nighttime Lights ---")
print(f"  The EOG group at Colorado School of Mines produces monthly/annual")
print(f"  composites of VIIRS nighttime lights, but daily data would be needed")
print(f"  for event detection.")
print(f"  URL: https://eogdata.mines.edu/products/vnl/")
print(f"  For daily anomaly detection, the VNP46A1 product is the right choice.")

# ============================================================================
# ASSESSMENT
# ============================================================================
print(f"\n\n{'='*80}")
print(f"VIIRS DNB ASSESSMENT")
print(f"{'='*80}")

print(f"""
  TIMING: The VIIRS descending (nighttime) pass at ~08:30 UTC arrives
  1.5-5.5 hours AFTER the 03:00-07:00 UTC event window. This means VIIRS
  would NOT have directly imaged the laser engagement in progress.

  WHAT IT COULD SHOW:
    - Residual fires or thermal sources left after engagement
    - Changes in nighttime light patterns at Fort Bliss (new lights,
      emergency vehicles, etc.)
    - Comparison of Feb 10 vs Feb 11 nighttime radiance at Fort Bliss

  WHAT IT CAN'T SHOW:
    - Real-time laser beam or engagement flashes
    - Events that occurred and ended before 08:30 UTC

  VERDICT: VIIRS DNB is of LIMITED use for this specific event due to
  timing mismatch. The ~08:30 UTC overpass misses the 03:00-07:00 UTC
  engagement window by 1.5+ hours.

  HOWEVER: If there were persistent fires or notable light changes,
  comparing the Feb 10 and Feb 11 nighttime passes would reveal them.
  The VNP46A1 product at LAADS DAAC is the right dataset to check.

  RECOMMENDED: Download h08v05 tiles for Feb 10 and Feb 11, extract
  radiance values at 31.87°N, 106.52°W, and compare.
""")

# ============================================================================
# Try to get actual data if VNP46A1 tiles are available
# ============================================================================
print(f"\n--- Attempting VNP46A1 Download ---")

# This requires an Earthdata token which we may not have
# Try direct LAADS download
for date_str, doy in [('2026-02-10', 41), ('2026-02-11', 42)]:
    search_url = f"https://ladsweb.modaps.eosdis.nasa.gov/api/v2/content/details/allData/5200/VNP46A1/2026/{doy:03d}"

    try:
        resp = requests.get(search_url, timeout=15, headers={'User-Agent': 'Mozilla/5.0'})
        if resp.status_code == 200:
            files = resp.json()
            if isinstance(files, list):
                our_files = [f for f in files if isinstance(f, dict) and 'h08v05' in f.get('name', '')]
                if our_files:
                    finfo = our_files[0]
                    download_url = f"https://ladsweb.modaps.eosdis.nasa.gov{finfo.get('downloadsLink', finfo.get('fileURL', ''))}"
                    print(f"  {date_str}: Found file, download URL requires Earthdata auth:")
                    print(f"    {download_url}")
                    print(f"    (Register at https://urs.earthdata.nasa.gov for free access)")
    except Exception as e:
        print(f"  {date_str}: {e}")

# Save results
output = {
    'timing_analysis': {
        'event_window': '03:00-07:00 UTC',
        'viirs_overpass': '~08:30 UTC (descending/nighttime)',
        'gap_hours': 1.5,
        'assessment': 'limited_utility_timing_mismatch',
    },
    'data_availability': {
        'VNP46A1': 'Available at LAADS DAAC, requires Earthdata token',
        'tile': 'h08v05',
    },
    'verdict': 'timing_mismatch_limited_utility',
}

output_path = os.path.join(ANALYSIS_DIR, 'viirs_dnb_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")
