#!/usr/bin/env python3
"""
IRIS Seismic / Infrasound Station Check
========================================
Searches for seismic/infrasound signals near El Paso / Fort Bliss during:
  - Feb 8-9 (booming sound reports)
  - Feb 11 (LOCUST engagement window, 03-07 UTC)

If something exploded or produced acoustic energy, broadband seismic sensors
within 200 km would register it.

Data source: IRIS FDSN Web Services (https://service.iris.edu/fdsnws/)
"""

import requests
import json
import os
import numpy as np
from datetime import datetime

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/seismic_data/'
ANALYSIS_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# CONFIGURATION
# ============================================================================
TARGET_LAT, TARGET_LON = 31.87, -106.52
MAX_RADIUS_KM = 200

# IRIS FDSN web service endpoints
STATION_URL = "https://service.iris.edu/fdsnws/station/1/query"
DATASELECT_URL = "https://service.iris.edu/fdsnws/dataselect/1/query"
TIMESERIES_URL = "https://service.iris.edu/irisws/timeseries/1/query"

# USGS earthquake catalog
USGS_EARTHQUAKE_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"

print("=" * 80)
print("IRIS SEISMIC / INFRASOUND ANALYSIS")
print(f"Search: {MAX_RADIUS_KM}km radius from ({TARGET_LAT}°N, {TARGET_LON}°W)")
print("=" * 80)

# ============================================================================
# STEP 1: Find nearby seismic stations
# ============================================================================
print(f"\n--- STEP 1: Finding seismic stations within {MAX_RADIUS_KM}km ---")

station_params = {
    'latitude': TARGET_LAT,
    'longitude': TARGET_LON,
    'maxradius': MAX_RADIUS_KM / 111.0,  # Convert km to degrees
    'level': 'station',
    'format': 'text',
    'starttime': '2025-01-01',
    'endtime': '2026-12-31',
}

try:
    resp = requests.get(STATION_URL, params=station_params, timeout=30)
    if resp.status_code == 200:
        lines = resp.text.strip().split('\n')
        header = lines[0] if lines[0].startswith('#') else None
        stations = []
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('|')
            if len(parts) >= 6:
                stations.append({
                    'network': parts[0].strip(),
                    'station': parts[1].strip(),
                    'lat': float(parts[2]),
                    'lon': float(parts[3]),
                    'elevation_m': float(parts[4]) if parts[4].strip() else None,
                    'name': parts[5].strip() if len(parts) > 5 else '',
                })

        print(f"  Found {len(stations)} stations")
        for s in stations:
            dist_km = np.sqrt(((s['lat'] - TARGET_LAT) * 111)**2 +
                              ((s['lon'] - TARGET_LON) * 111 * np.cos(np.radians(TARGET_LAT)))**2)
            s['dist_km'] = dist_km
            print(f"    {s['network']}.{s['station']:<6} {s['name'][:40]:<40} "
                  f"({s['lat']:.3f}, {s['lon']:.3f}) {dist_km:.0f}km")

        # Save station list
        with open(os.path.join(OUTPUT_DIR, 'nearby_stations.json'), 'w') as f:
            json.dump(stations, f, indent=2)
    else:
        print(f"  Station query returned {resp.status_code}")
        print(f"  Response: {resp.text[:500]}")
        stations = []
except Exception as e:
    print(f"  Error querying stations: {e}")
    stations = []

# ============================================================================
# STEP 2: Check USGS earthquake catalog for reported events
# ============================================================================
print(f"\n--- STEP 2: USGS Earthquake Catalog ---")

for label, start, end in [
    ("Booming sounds (Feb 8-9)", "2026-02-08T00:00:00", "2026-02-10T00:00:00"),
    ("Event night (Feb 11)", "2026-02-11T00:00:00", "2026-02-12T00:00:00"),
    ("Extended (Feb 7-14)", "2026-02-07T00:00:00", "2026-02-14T00:00:00"),
]:
    print(f"\n  {label}:")
    eq_params = {
        'format': 'geojson',
        'starttime': start,
        'endtime': end,
        'latitude': TARGET_LAT,
        'longitude': TARGET_LON,
        'maxradiuskm': MAX_RADIUS_KM,
        'minmagnitude': 0.0,
    }

    try:
        resp = requests.get(USGS_EARTHQUAKE_URL, params=eq_params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            features = data.get('features', [])
            print(f"    Events found: {len(features)}")

            for feat in features:
                props = feat['properties']
                geom = feat['geometry']['coordinates']
                print(f"      M{props.get('mag', '?')} at {props.get('place', '?')}")
                print(f"        Time: {props.get('time', '?')} Type: {props.get('type', '?')}")
                print(f"        Location: ({geom[1]:.3f}, {geom[0]:.3f}), Depth: {geom[2]:.1f}km")

            if not features:
                print(f"    No cataloged seismic events")

            # Save
            with open(os.path.join(OUTPUT_DIR, f'usgs_events_{label[:10].replace(" ","_")}.json'), 'w') as f:
                json.dump(data, f, indent=2)

        elif resp.status_code == 204:
            print(f"    No events found (204)")
        else:
            print(f"    HTTP {resp.status_code}")
    except Exception as e:
        print(f"    Error: {e}")

# ============================================================================
# STEP 3: Download waveform data for closest stations during event window
# ============================================================================
print(f"\n--- STEP 3: Waveform Data from Closest Stations ---")

if stations:
    # Sort by distance, take closest 5
    stations_sorted = sorted(stations, key=lambda x: x['dist_km'])[:5]

    for s in stations_sorted:
        net = s['network']
        sta = s['station']
        print(f"\n  Station {net}.{sta} ({s['dist_km']:.0f}km):")

        # Try to get channel info first
        chan_params = {
            'network': net,
            'station': sta,
            'level': 'channel',
            'format': 'text',
            'starttime': '2026-02-01',
            'endtime': '2026-02-28',
        }

        channels = []
        try:
            resp = requests.get(STATION_URL, params=chan_params, timeout=15)
            if resp.status_code == 200:
                for line in resp.text.strip().split('\n'):
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split('|')
                    if len(parts) >= 4:
                        chan = parts[3].strip()
                        if chan not in channels:
                            channels.append(chan)
                print(f"    Available channels: {', '.join(channels[:10])}")
        except Exception:
            channels = ['BHZ', 'HHZ', 'BDF']  # Default broadband + infrasound

        # Try to download short waveform segments during key times
        for label, start, end in [
            ("Event peak (04:00-04:30 UTC)", "2026-02-11T04:00:00", "2026-02-11T04:30:00"),
            ("Booming (Feb 9 night)", "2026-02-09T04:00:00", "2026-02-09T04:30:00"),
        ]:
            # Pick best channel: prefer BHZ (broadband vertical) or HHZ (high-rate)
            target_chan = None
            for pref in ['BHZ', 'HHZ', 'BDF', 'LHZ', 'SHZ', 'EHZ']:
                if pref in channels:
                    target_chan = pref
                    break
            if not target_chan and channels:
                target_chan = channels[0]
            if not target_chan:
                target_chan = 'BHZ'

            print(f"    Trying {label} on {target_chan}...")

            # Use timeseries service for small data downloads
            ts_params = {
                'network': net,
                'station': sta,
                'channel': target_chan,
                'location': '--',
                'starttime': start,
                'endtime': end,
                'output': 'ascii',  # Human-readable
            }

            try:
                resp = requests.get(TIMESERIES_URL, params=ts_params, timeout=30)
                if resp.status_code == 200 and len(resp.text) > 100:
                    # Parse ASCII timeseries
                    lines = resp.text.strip().split('\n')
                    header_lines = [l for l in lines if l.startswith('TIMESERIES')]
                    data_lines = [l for l in lines if not l.startswith('TIMESERIES') and l.strip()]

                    if data_lines:
                        values = []
                        for dl in data_lines:
                            for v in dl.split():
                                try:
                                    values.append(float(v))
                                except ValueError:
                                    pass

                        if values:
                            vals = np.array(values)
                            print(f"      Got {len(vals)} samples")
                            print(f"      Mean: {np.mean(vals):.2f}, Std: {np.std(vals):.2f}")
                            print(f"      Peak: {np.max(np.abs(vals)):.2f}")
                            print(f"      RMS: {np.sqrt(np.mean(vals**2)):.2f}")

                            # Check for anomalous signals (>5 sigma spikes)
                            mean, std = np.mean(vals), np.std(vals)
                            if std > 0:
                                spikes = np.sum(np.abs(vals - mean) > 5 * std)
                                if spikes > 0:
                                    print(f"      *** {spikes} samples exceed 5-sigma ***")

                            # Save waveform
                            fname = f"{net}_{sta}_{target_chan}_{label[:10].replace(' ','_')}.txt"
                            with open(os.path.join(OUTPUT_DIR, fname), 'w') as f:
                                f.write(resp.text)
                            print(f"      Saved to {fname}")
                        else:
                            print(f"      No numeric data parsed")
                    else:
                        print(f"      Empty data response")
                elif resp.status_code == 204:
                    print(f"      No data available (204)")
                elif resp.status_code == 404:
                    print(f"      Data not found (404)")
                else:
                    print(f"      HTTP {resp.status_code}: {resp.text[:200]}")
            except Exception as e:
                print(f"      Error: {e}")

# ============================================================================
# STEP 4: Check for infrasound-specific stations
# ============================================================================
print(f"\n\n--- STEP 4: Infrasound Array Search ---")
print(f"  Searching for BDF/BDG channels (infrasound) near El Paso...")

infrasound_params = {
    'latitude': TARGET_LAT,
    'longitude': TARGET_LON,
    'maxradius': 5.0,  # ~555 km for infrasound (travels far)
    'channel': 'BDF,BDG,BDH',  # Infrasound channels
    'level': 'channel',
    'format': 'text',
    'starttime': '2026-01-01',
    'endtime': '2026-12-31',
}

try:
    resp = requests.get(STATION_URL, params=infrasound_params, timeout=30)
    if resp.status_code == 200 and len(resp.text) > 50:
        lines = resp.text.strip().split('\n')
        infra_stations = []
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('|')
            if len(parts) >= 4:
                infra_stations.append({
                    'network': parts[0].strip(),
                    'station': parts[1].strip(),
                    'channel': parts[3].strip() if len(parts) > 3 else '',
                })
        # Deduplicate
        seen = set()
        unique = []
        for s in infra_stations:
            key = f"{s['network']}.{s['station']}"
            if key not in seen:
                seen.add(key)
                unique.append(s)

        print(f"  Found {len(unique)} infrasound stations within 555km:")
        for s in unique:
            print(f"    {s['network']}.{s['station']} [{s['channel']}]")
    else:
        print(f"  No infrasound stations found (or query failed: {resp.status_code})")
except Exception as e:
    print(f"  Error: {e}")

# Also check for IMS (International Monitoring System) infrasound
print(f"\n  Note: The nearest IMS infrasound arrays are likely:")
print(f"    IS57 — Piñon Flat, CA (~800km, may be too far)")
print(f"    IS56 — Newport, WA (too far)")
print(f"    For closer coverage, check LANL or Sandia infrasound arrays")

# ============================================================================
# SUMMARY
# ============================================================================
print(f"\n\n{'='*80}")
print(f"SEISMIC/INFRASOUND SUMMARY")
print(f"{'='*80}")

print(f"\n  Stations within {MAX_RADIUS_KM}km: {len(stations)}")
print(f"  Key findings:")
print(f"    - USGS earthquake catalog checked for Feb 7-14, 2026")
print(f"    - Waveform data attempted for closest stations")
print(f"    - Infrasound array search completed")
print(f"\n  Note: If broadband seismic stations detected anomalous signals,")
print(f"  this would appear as elevated RMS or >5-sigma spikes in the")
print(f"  waveform data during the event window vs. baseline.")

# Save results
output = {
    'search_params': {
        'center_lat': TARGET_LAT, 'center_lon': TARGET_LON,
        'radius_km': MAX_RADIUS_KM,
    },
    'stations_found': len(stations),
    'station_list': stations[:20] if stations else [],
    'analysis_complete': True,
}

output_path = os.path.join(ANALYSIS_DIR, 'iris_seismic_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")
