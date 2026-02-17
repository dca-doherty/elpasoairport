#!/usr/bin/env python3
"""
Engagement Event Geolocation
El Paso Airspace Incident — Feb 11, 2026

For each identified engagement event, extracts:
  - Peak reflectivity position (lat/lon of brightest return)
  - Reflectivity-weighted centroid of all returns >10 dBZ
  - Debris field extent and bounding box
  - Position shift between events (target motion)

Data sources:
  - KHDX Level 2 significant returns (azimuth, range, lat, lon, dBZ)
  - KEPZ spike parameters (range, azimuth → lat, lon)
"""

import json
import numpy as np

# ============================================================================
# LOAD DATA
# ============================================================================
with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/khdx_analysis_results.json') as f:
    khdx_data = json.load(f)

returns = khdx_data['significant_returns']

# ============================================================================
# RADAR SITE POSITIONS
# ============================================================================
KEPZ_LAT = 31.8731
KEPZ_LON = -106.6981
KHDX_LAT = 33.0803
KHDX_LON = -106.1225

# KEPZ spike parameters (from imagery analysis)
# Spike 1: ~30 km range, ~88° azimuth from KEPZ
# Spike 2: ~30 km range, ~82° azimuth from KEPZ
def radar_to_latlon(radar_lat, radar_lon, range_km, azimuth_deg):
    """Convert radar polar coords (range, azimuth) to lat/lon."""
    lat = radar_lat + (range_km * np.cos(np.radians(azimuth_deg))) / 111.0
    lon = radar_lon + (range_km * np.sin(np.radians(azimuth_deg))) / (111.0 * np.cos(np.radians(radar_lat)))
    return lat, lon

# KEPZ-derived positions
kepz_s1_lat, kepz_s1_lon = radar_to_latlon(KEPZ_LAT, KEPZ_LON, 30.0, 88.0)
kepz_s2_lat, kepz_s2_lon = radar_to_latlon(KEPZ_LAT, KEPZ_LON, 30.0, 82.0)

print("=" * 78)
print("ENGAGEMENT EVENT GEOLOCATION")
print("El Paso Airspace Incident — Feb 11, 2026")
print("=" * 78)

print(f"\n  KEPZ-derived spike positions:")
print(f"    Spike 1 (04:24): {kepz_s1_lat:.4f}°N, {abs(kepz_s1_lon):.4f}°W")
print(f"    Spike 2 (05:13): {kepz_s2_lat:.4f}°N, {abs(kepz_s2_lon):.4f}°W")

# ============================================================================
# ENGAGEMENT EVENT TIME WINDOWS
# ============================================================================
# Each event defined by the scan time(s) at peak activity
engagement_events = [
    {
        'id': 'E1',
        'name': 'Initial Detection',
        'confidence': 'Initial',
        'scan_times': ['0317', '0325', '0332'],  # Pre-NOTAM ramp-up
        'kepz_visible': False,
    },
    {
        'id': 'E2',
        'name': 'Primary Engagement (KEPZ S1)',
        'confidence': 'Confirmed',
        'scan_times': ['0421', '0429', '0436'],  # Around KEPZ Spike 1
        'kepz_visible': True,
        'kepz_lat': kepz_s1_lat,
        'kepz_lon': kepz_s1_lon,
    },
    {
        'id': 'E3',
        'name': 'Hidden Engagement',
        'confidence': 'Likely',
        'scan_times': ['0457', '0504'],  # KHDX-only spike at 04:57
        'kepz_visible': False,
    },
    {
        'id': 'E4',
        'name': 'Second Engagement (KEPZ S2)',
        'confidence': 'Confirmed',
        'scan_times': ['0511', '0518', '0526'],  # Around KEPZ Spike 2
        'kepz_visible': True,
        'kepz_lat': kepz_s2_lat,
        'kepz_lon': kepz_s2_lon,
    },
    {
        'id': 'E5',
        'name': 'Escalation',
        'confidence': 'Likely',
        'scan_times': ['0533', '0540', '0547', '0554'],  # Post-S2 escalation
        'kepz_visible': False,
    },
    {
        'id': 'E6',
        'name': 'Maximum Return',
        'confidence': 'Probable',
        'scan_times': ['0601', '0608', '0615', '0622'],  # Strongest returns
        'kepz_visible': False,
    },
]

# ============================================================================
# EXTRACT POSITIONS FOR EACH EVENT
# ============================================================================
print("\n\n" + "=" * 78)
print("KHDX-DERIVED POSITIONS FOR EACH ENGAGEMENT EVENT")
print("=" * 78)

event_positions = []

for event in engagement_events:
    # Filter returns for this event's time window
    event_returns = [r for r in returns if r['time'] in event['scan_times']]

    if not event_returns:
        print(f"\n  {event['id']} ({event['name']}): No KHDX returns in time window")
        event_positions.append({
            'event_id': event['id'],
            'name': event['name'],
            'confidence': event['confidence'],
            'n_returns': 0,
        })
        continue

    # Separate by reflectivity threshold
    strong_returns = [r for r in event_returns if r['ref_dBZ'] > 10]
    all_returns = event_returns

    # Peak return (highest reflectivity)
    peak = max(event_returns, key=lambda r: r['ref_dBZ'])

    # Reflectivity-weighted centroid (for returns > 5 dBZ)
    lats = np.array([r['lat'] for r in event_returns])
    lons = np.array([r['lon'] for r in event_returns])
    refs = np.array([r['ref_dBZ'] for r in event_returns])

    # Convert reflectivity to linear scale for weighting
    z_linear = 10 ** (refs / 10.0)
    total_weight = np.sum(z_linear)

    centroid_lat = np.sum(lats * z_linear) / total_weight
    centroid_lon = np.sum(lons * z_linear) / total_weight

    # Bounding box
    lat_min, lat_max = np.min(lats), np.max(lats)
    lon_min, lon_max = np.min(lons), np.max(lons)
    extent_ns_km = (lat_max - lat_min) * 111.0
    extent_ew_km = (lon_max - lon_min) * 111.0 * np.cos(np.radians(np.mean(lats)))

    # Strong returns centroid (>10 dBZ only)
    strong_centroid_lat = None
    strong_centroid_lon = None
    if strong_returns:
        s_lats = np.array([r['lat'] for r in strong_returns])
        s_lons = np.array([r['lon'] for r in strong_returns])
        s_refs = np.array([r['ref_dBZ'] for r in strong_returns])
        s_z = 10 ** (s_refs / 10.0)
        strong_centroid_lat = np.sum(s_lats * s_z) / np.sum(s_z)
        strong_centroid_lon = np.sum(s_lons * s_z) / np.sum(s_z)

    pos = {
        'event_id': event['id'],
        'name': event['name'],
        'confidence': event['confidence'],
        'n_returns': len(event_returns),
        'n_strong_returns': len(strong_returns),
        'peak_lat': peak['lat'],
        'peak_lon': peak['lon'],
        'peak_ref_dBZ': peak['ref_dBZ'],
        'peak_az_from_khdx': peak['az'],
        'peak_range_from_khdx_km': peak['range_km'],
        'centroid_lat': float(centroid_lat),
        'centroid_lon': float(centroid_lon),
        'strong_centroid_lat': float(strong_centroid_lat) if strong_centroid_lat else None,
        'strong_centroid_lon': float(strong_centroid_lon) if strong_centroid_lon else None,
        'extent_ns_km': float(extent_ns_km),
        'extent_ew_km': float(extent_ew_km),
        'bbox': {
            'lat_min': float(lat_min),
            'lat_max': float(lat_max),
            'lon_min': float(lon_min),
            'lon_max': float(lon_max),
        },
        'kepz_visible': event['kepz_visible'],
    }

    if event.get('kepz_lat'):
        pos['kepz_lat'] = event['kepz_lat']
        pos['kepz_lon'] = event['kepz_lon']

    event_positions.append(pos)

    print(f"\n  {event['id']}: {event['name']} [{event['confidence']}]")
    print(f"  ─────────────────────────────────────────────────")
    print(f"    KHDX returns: {len(event_returns)} total, {len(strong_returns)} above 10 dBZ")
    print(f"    Peak return:  {peak['ref_dBZ']:.1f} dBZ at {peak['lat']:.4f}°N, {abs(peak['lon']):.4f}°W")
    print(f"                  (az={peak['az']:.1f}° range={peak['range_km']:.1f}km from KHDX)")
    print(f"    Z-weighted centroid: {centroid_lat:.4f}°N, {abs(centroid_lon):.4f}°W")
    if strong_centroid_lat:
        print(f"    Strong (>10dBZ) centroid: {strong_centroid_lat:.4f}°N, {abs(strong_centroid_lon):.4f}°W")
    print(f"    Debris field extent: {extent_ns_km:.1f} km N-S × {extent_ew_km:.1f} km E-W")
    print(f"    Bounding box: {lat_min:.4f}-{lat_max:.4f}°N, {abs(lon_max):.4f}-{abs(lon_min):.4f}°W")

    if event.get('kepz_lat'):
        dist_to_kepz = np.sqrt(
            ((pos['centroid_lat'] - event['kepz_lat']) * 111.0)**2 +
            ((pos['centroid_lon'] - event['kepz_lon']) * 111.0 * np.cos(np.radians(pos['centroid_lat'])))**2
        )
        print(f"    KEPZ position: {event['kepz_lat']:.4f}°N, {abs(event['kepz_lon']):.4f}°W")
        print(f"    KHDX-KEPZ offset: {dist_to_kepz:.1f} km")

# ============================================================================
# POSITION COMPARISON ACROSS EVENTS
# ============================================================================
print("\n\n" + "=" * 78)
print("POSITION COMPARISON ACROSS EVENTS — DID THE TARGET MOVE?")
print("=" * 78)

print(f"\n  {'Event':<6} {'Conf':>10} {'Centroid Lat':>13} {'Centroid Lon':>13} "
      f"{'Peak dBZ':>9} {'N Returns':>10}")
print(f"  {'─'*6} {'─'*10} {'─'*13} {'─'*13} {'─'*9} {'─'*10}")

prev_lat = None
prev_lon = None

for pos in event_positions:
    if pos['n_returns'] == 0:
        print(f"  {pos['event_id']:<6} {pos['confidence']:>10}   {'(no data)':>13}")
        continue

    lat = pos['centroid_lat']
    lon = pos['centroid_lon']

    motion_str = ""
    if prev_lat is not None:
        dist_km = np.sqrt(
            ((lat - prev_lat) * 111.0)**2 +
            ((lon - prev_lon) * 111.0 * np.cos(np.radians(lat)))**2
        )
        # Direction
        dlat = lat - prev_lat
        dlon = lon - prev_lon
        bearing = np.degrees(np.arctan2(dlon * np.cos(np.radians(lat)), dlat)) % 360
        dirs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
        dir_idx = int((bearing + 22.5) / 45) % 8
        motion_str = f"  → {dist_km:.1f} km {dirs[dir_idx]} from previous"

    print(f"  {pos['event_id']:<6} {pos['confidence']:>10} {lat:>13.4f} {lon:>13.4f} "
          f"{pos['peak_ref_dBZ']:>9.1f} {pos['n_returns']:>10}{motion_str}")

    prev_lat = lat
    prev_lon = lon

# ============================================================================
# BEST POSITION ESTIMATES
# ============================================================================
print("\n\n" + "=" * 78)
print("BEST POSITION ESTIMATES FOR EACH ENGAGEMENT")
print("=" * 78)

print("""
  METHOD: For KEPZ-visible events (E2, E4), we use the KEPZ-derived position
  as primary (better spatial resolution at 30 km) with KHDX as confirmation.
  For KHDX-only events, we use the reflectivity-weighted centroid of
  returns >10 dBZ (suppresses noise from peripheral weak returns).
""")

print(f"  {'Event':<6} {'Name':<30} {'Best Lat':>10} {'Best Lon':>12} {'Source':>12} {'Accuracy':>10}")
print(f"  {'─'*6} {'─'*30} {'─'*10} {'─'*12} {'─'*12} {'─'*10}")

best_positions = []

for pos in event_positions:
    if pos['n_returns'] == 0:
        continue

    if pos['kepz_visible'] and pos.get('kepz_lat'):
        # Use KEPZ position (30 km range → ~250m gate resolution)
        best_lat = pos['kepz_lat']
        best_lon = pos['kepz_lon']
        source = "KEPZ"
        accuracy = "~0.5 km"
    elif pos.get('strong_centroid_lat'):
        # Use strong-return centroid from KHDX
        best_lat = pos['strong_centroid_lat']
        best_lon = pos['strong_centroid_lon']
        source = "KHDX >10dBZ"
        accuracy = "~2-5 km"
    else:
        # Use all-return centroid
        best_lat = pos['centroid_lat']
        best_lon = pos['centroid_lon']
        source = "KHDX all"
        accuracy = "~5-10 km"

    best_positions.append({
        'event_id': pos['event_id'],
        'name': pos['name'],
        'lat': best_lat,
        'lon': best_lon,
        'source': source,
        'accuracy': accuracy,
    })

    print(f"  {pos['event_id']:<6} {pos['name']:<30} {best_lat:>10.4f} {best_lon:>12.4f} {source:>12} {accuracy:>10}")

# ============================================================================
# GEOGRAPHIC CONTEXT — What's at each location?
# ============================================================================
print("\n\n" + "=" * 78)
print("GEOGRAPHIC CONTEXT")
print("=" * 78)

# Key reference points
landmarks = [
    ("Biggs Army Airfield (KBIF)", 31.8495, -106.3800),
    ("Fort Bliss main post", 31.8115, -106.4222),
    ("HELWS test range (McGregor)", 32.05, -106.18),
    ("Doña Ana Range Complex", 32.38, -106.55),
    ("El Paso city center", 31.7619, -106.4850),
    ("KEPZ radar", KEPZ_LAT, KEPZ_LON),
    ("KHDX radar", KHDX_LAT, KHDX_LON),
    ("Franklin Mountains summit", 31.9113, -106.4954),
    ("US-Mexico border (nearest)", 31.755, -106.45),
]

print()
for bp in best_positions:
    print(f"\n  {bp['event_id']}: {bp['name']}")
    print(f"  Position: {bp['lat']:.4f}°N, {abs(bp['lon']):.4f}°W")
    print(f"  Distances to landmarks:")

    for lm_name, lm_lat, lm_lon in landmarks:
        dist = np.sqrt(
            ((bp['lat'] - lm_lat) * 111.0)**2 +
            ((bp['lon'] - lm_lon) * 111.0 * np.cos(np.radians(bp['lat'])))**2
        )
        bearing = np.degrees(np.arctan2(
            (lm_lon - bp['lon']) * np.cos(np.radians(bp['lat'])),
            lm_lat - bp['lat']
        )) % 360
        dirs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
        dir_idx = int((bearing + 22.5) / 45) % 8

        if dist < 50:
            marker = " ◄" if dist < 5 else ""
            print(f"    {dist:>6.1f} km {dirs[dir_idx]:>2} to {lm_name}{marker}")

# ============================================================================
# CROSS-RADAR TRIANGULATION
# ============================================================================
print("\n\n" + "=" * 78)
print("CROSS-RADAR TRIANGULATION (E2 and E4)")
print("=" * 78)

print("""
  For events E2 and E4 (visible to both KEPZ and KHDX), we can
  compare the independently-derived positions:
""")

for pos in event_positions:
    if pos['kepz_visible'] and pos.get('kepz_lat') and pos.get('strong_centroid_lat'):
        kepz_lat = pos['kepz_lat']
        kepz_lon = pos['kepz_lon']
        khdx_lat = pos['strong_centroid_lat']
        khdx_lon = pos['strong_centroid_lon']

        offset_km = np.sqrt(
            ((kepz_lat - khdx_lat) * 111.0)**2 +
            ((kepz_lon - khdx_lon) * 111.0 * np.cos(np.radians(kepz_lat)))**2
        )

        # Average position (weighted by range — closer radar gets more weight)
        # KEPZ at 30 km, KHDX at 135 km → KEPZ gets ~4.5× more weight
        w_kepz = 135.0 / (30.0 + 135.0)
        w_khdx = 30.0 / (30.0 + 135.0)
        tri_lat = w_kepz * kepz_lat + w_khdx * khdx_lat
        tri_lon = w_kepz * kepz_lon + w_khdx * khdx_lon

        print(f"  {pos['event_id']}: {pos['name']}")
        print(f"    KEPZ position:  {kepz_lat:.4f}°N, {abs(kepz_lon):.4f}°W  (30 km range)")
        print(f"    KHDX centroid:  {khdx_lat:.4f}°N, {abs(khdx_lon):.4f}°W  (135 km range)")
        print(f"    Offset: {offset_km:.1f} km")
        print(f"    Triangulated:   {tri_lat:.4f}°N, {abs(tri_lon):.4f}°W")
        print(f"    (KEPZ weighted {w_kepz*100:.0f}%, KHDX weighted {w_khdx*100:.0f}%)")
        print()

# ============================================================================
# COORDINATE SUMMARY TABLE
# ============================================================================
print("\n" + "=" * 78)
print("COORDINATE SUMMARY — ALL ENGAGEMENT EVENTS")
print("=" * 78)
print()
print(f"  {'Event':<5} {'Time':>8} {'Latitude':>12} {'Longitude':>12} {'dBZ':>5} "
      f"{'Accuracy':>10} {'From Biggs':>11}")
print(f"  {'─'*5} {'─'*8} {'─'*12} {'─'*12} {'─'*5} {'─'*10} {'─'*11}")

event_times = {
    'E1': '03:25',
    'E2': '04:24',
    'E3': '04:57',
    'E4': '05:13',
    'E5': '05:40',
    'E6': '06:08',
}

for bp in best_positions:
    eid = bp['event_id']
    time = event_times.get(eid, '??:??')

    # Distance from Biggs AAF
    dist_biggs = np.sqrt(
        ((bp['lat'] - 31.8495) * 111.0)**2 +
        ((bp['lon'] - (-106.3800)) * 111.0 * np.cos(np.radians(bp['lat'])))**2
    )

    # Find peak dBZ for this event
    peak_dbz = 0
    for pos in event_positions:
        if pos['event_id'] == eid and pos['n_returns'] > 0:
            peak_dbz = pos['peak_ref_dBZ']

    print(f"  {eid:<5} {time:>8} {bp['lat']:>12.4f} {bp['lon']:>12.4f} {peak_dbz:>5.1f} "
          f"{bp['accuracy']:>10} {dist_biggs:>8.1f} km")

print()
print("  Note: KHDX gate spacing is 250m × ~1° at 135 km range = ~2.4 km cross-range")
print("  KEPZ gate spacing is 250m × ~1° at 30 km range = ~0.5 km cross-range")
print("  Positions from KHDX-only events have inherently lower spatial precision.")

# ============================================================================
# SAVE RESULTS
# ============================================================================
output = {
    'event_positions': event_positions,
    'best_positions': best_positions,
    'kepz_derived': {
        'spike1': {'lat': float(kepz_s1_lat), 'lon': float(kepz_s1_lon)},
        'spike2': {'lat': float(kepz_s2_lat), 'lon': float(kepz_s2_lon)},
    },
    'radar_sites': {
        'KEPZ': {'lat': KEPZ_LAT, 'lon': KEPZ_LON},
        'KHDX': {'lat': KHDX_LAT, 'lon': KHDX_LON},
    },
}

output_path = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/engagement_geolocation_results.json'
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)

print(f"\nResults saved to {output_path}")
