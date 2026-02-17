#!/usr/bin/env python3
"""
NASA FIRMS Thermal Hotspot Check
=================================
Checks for active fire / thermal anomaly detections near Fort Bliss / El Paso
during the incident window. Long shot but trivially easy to check.

Uses FIRMS REST API (MODIS + VIIRS products).
"""

import requests
import json
import csv
import os
from io import StringIO

# ============================================================================
# CONFIGURATION
# ============================================================================
# FIRMS MAP_KEY — using the public demo key (limited queries)
# For production, register at https://firms.modaps.eosdis.nasa.gov/api/area/
MAP_KEY = 'DEMO_KEY'  # Will try with this, fall back to direct URL

# Fort Bliss area bounding box
# Roughly El Paso metro + Fort Bliss + White Sands edge
WEST, SOUTH, EAST, NORTH = -107.0, 31.5, -106.0, 32.5
AREA = f"{WEST},{SOUTH},{EAST},{NORTH}"

# Date range: Feb 8-13, 2026
DATE = '2026-02-10'
DAY_RANGE = 5  # days from date

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/firms_data/'
ANALYSIS_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Sources to check
SOURCES = {
    'VIIRS_SNPP': 'VIIRS Suomi NPP (375m)',
    'VIIRS_NOAA20': 'VIIRS NOAA-20 (375m)',
    'VIIRS_NOAA21': 'VIIRS NOAA-21 (375m)',
    'MODIS_NRT': 'MODIS Aqua/Terra (1km)',
}

print("=" * 80)
print("NASA FIRMS THERMAL HOTSPOT CHECK")
print(f"Area: {AREA}")
print(f"Date range: Feb 8-13, 2026")
print("=" * 80)

all_hotspots = []
results_by_source = {}

for source_key, source_name in SOURCES.items():
    print(f"\n  Checking {source_name}...")

    # Try FIRMS API
    url = f"https://firms.modaps.eosdis.nasa.gov/api/area/csv/{MAP_KEY}/{source_key}/{AREA}/{DAY_RANGE}/{DATE}"

    try:
        resp = requests.get(url, timeout=30, headers={'User-Agent': 'Mozilla/5.0'})

        if resp.status_code == 200 and 'latitude' in resp.text.lower():
            reader = csv.DictReader(StringIO(resp.text))
            hotspots = list(reader)
            results_by_source[source_key] = hotspots
            all_hotspots.extend(hotspots)
            print(f"    Found {len(hotspots)} hotspots")

            for hs in hotspots[:5]:
                lat = hs.get('latitude', '?')
                lon = hs.get('longitude', '?')
                date = hs.get('acq_date', '?')
                time = hs.get('acq_time', '?')
                conf = hs.get('confidence', '?')
                frp = hs.get('frp', '?')
                print(f"      {date} {time} UTC: ({lat}, {lon}) conf={conf} FRP={frp} MW")

            # Save raw CSV
            csv_path = os.path.join(OUTPUT_DIR, f'firms_{source_key}.csv')
            with open(csv_path, 'w') as f:
                f.write(resp.text)
            print(f"    Saved to {csv_path}")

        elif resp.status_code == 200:
            print(f"    No hotspots (empty response)")
            results_by_source[source_key] = []
        elif resp.status_code == 403:
            print(f"    API key rejected (403). Trying alternative URL...")
            # Try the NRT (near-real-time) direct download
            alt_url = f"https://firms.modaps.eosdis.nasa.gov/data/active_fire/{source_key}/csv/{source_key}_Global_24h.csv"
            print(f"    (Would need registered API key for historical queries)")
            results_by_source[source_key] = 'key_required'
        else:
            print(f"    HTTP {resp.status_code}: {resp.text[:200]}")
            results_by_source[source_key] = f'error_{resp.status_code}'

    except Exception as e:
        print(f"    Error: {e}")
        results_by_source[source_key] = f'error: {str(e)}'

# ============================================================================
# ALTERNATIVE: Try FIRMS MAP viewer URL for reference
# ============================================================================
print(f"\n\n{'─'*80}")
print(f"FIRMS MAP Viewer URLs (for manual verification):")
print(f"{'─'*80}")
map_url = (f"https://firms.modaps.eosdis.nasa.gov/map/#d:24hrs;@{(WEST+EAST)/2},{(SOUTH+NORTH)/2},10z")
print(f"  General: {map_url}")
print(f"  For Feb 11: Use date selector to navigate to 2026-02-11")

# ============================================================================
# ANALYSIS
# ============================================================================
print(f"\n\n{'='*80}")
print(f"FIRMS ANALYSIS RESULTS")
print(f"{'='*80}")

if all_hotspots:
    # Filter to Fort Bliss specific area (tighter box)
    bliss_hotspots = [h for h in all_hotspots
                      if float(h.get('latitude', 0)) > 31.7
                      and float(h.get('latitude', 0)) < 32.2
                      and float(h.get('longitude', 0)) > -106.8
                      and float(h.get('longitude', 0)) < -106.2]

    print(f"\n  Total hotspots in search area: {len(all_hotspots)}")
    print(f"  Hotspots near Fort Bliss (tight box): {len(bliss_hotspots)}")

    # Group by date
    by_date = {}
    for h in all_hotspots:
        d = h.get('acq_date', 'unknown')
        by_date.setdefault(d, []).append(h)

    print(f"\n  By date:")
    for d in sorted(by_date.keys()):
        marker = " ***" if '02-11' in d else ""
        print(f"    {d}: {len(by_date[d])} hotspots{marker}")

    if bliss_hotspots:
        print(f"\n  *** THERMAL ANOMALIES DETECTED NEAR FORT BLISS ***")
        for h in bliss_hotspots:
            print(f"    {h.get('acq_date','')} {h.get('acq_time','')} UTC: "
                  f"({h.get('latitude','')}, {h.get('longitude','')}) "
                  f"FRP={h.get('frp','')} MW conf={h.get('confidence','')}")
    else:
        print(f"\n  No thermal anomalies specifically at Fort Bliss")
else:
    print(f"\n  No FIRMS hotspots detected in the El Paso area for this period.")
    print(f"  This is not surprising — FIRMS is designed for fires/volcanoes,")
    print(f"  not laser engagements. The null result is expected but worth documenting.")

# ============================================================================
# SAVE RESULTS
# ============================================================================
output = {
    'search_area': {'west': WEST, 'south': SOUTH, 'east': EAST, 'north': NORTH},
    'date_range': f'{DATE} +/- {DAY_RANGE} days',
    'sources_checked': list(SOURCES.keys()),
    'total_hotspots': len(all_hotspots),
    'by_source': {k: len(v) if isinstance(v, list) else v for k, v in results_by_source.items()},
    'fort_bliss_hotspots': len([h for h in all_hotspots
                                if float(h.get('latitude', 0)) > 31.7
                                and float(h.get('latitude', 0)) < 32.2
                                and float(h.get('longitude', 0)) > -106.8
                                and float(h.get('longitude', 0)) < -106.2]),
    'verdict': 'hotspot_detected' if any(float(h.get('latitude', 0)) > 31.7
                                         and float(h.get('latitude', 0)) < 32.2
                                         and float(h.get('longitude', 0)) > -106.8
                                         and float(h.get('longitude', 0)) < -106.2
                                         for h in all_hotspots) else 'no_detection',
}

output_path = os.path.join(ANALYSIS_DIR, 'firms_hotspot_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")
