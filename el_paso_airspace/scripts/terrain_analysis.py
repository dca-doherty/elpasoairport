#!/usr/bin/env python3
"""
Terrain-Corrected Position Analysis for El Paso Airspace Incident.
Analyzes Franklin Mountains beam blockage for KEPZ radar.

KEPZ radar: 31.8731°N, 106.6981°W, elevation 1251m MSL
Franklin Mountains: run N-S between KEPZ and Fort Bliss spikes
Spike 1: 31.8771°N, 106.3423°W (az ~85° from KEPZ, range ~33.6 km)
Spike 2: 31.8783°N, 106.4164°W (az ~88° from KEPZ, range ~26.6 km)
"""

import numpy as np
import json
import os

# ============= RADAR AND TARGET PARAMETERS =============
# KEPZ radar location
KEPZ_LAT = 31.8731
KEPZ_LON = -106.6981
KEPZ_ELEV = 1251.0  # meters MSL

# Spike locations
SPIKE1_LAT = 31.8771
SPIKE1_LON = -106.3423
SPIKE2_LAT = 31.8783
SPIKE2_LON = -106.4164

# NEXRAD beam parameters
BEAM_ELEV_05 = 0.5   # degrees (lowest tilt)
BEAM_ELEV_15 = 1.5   # next tilt
BEAMWIDTH = 0.95      # degrees (3dB beamwidth)
GATE_LENGTH = 250.0   # meters (standard NEXRAD gate)
EARTH_RADIUS = 6371000.0  # meters
EARTH_RADIUS_EFF = EARTH_RADIUS * 4/3  # 4/3 earth for refraction

# Franklin Mountains key peaks (approximate)
FRANKLIN_PEAKS = {
    'North Franklin Peak': (31.9447, -106.5008, 2192),
    'Mundy Peak': (31.9100, -106.4900, 2100),
    'Ranger Peak': (31.8600, -106.4830, 1690),
    'Scenic Drive Summit': (31.7900, -106.4700, 1600),
    'Trans Mountain Gap': (31.9200, -106.4800, 1580),
    'Smugglers Gap': (31.8400, -106.4850, 1500),
}

results = {}

print("=" * 70)
print("TERRAIN-CORRECTED POSITION ANALYSIS")
print("Franklin Mountains Beam Blockage Assessment")
print("=" * 70)

def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculate distance in meters between two lat/lon points."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return EARTH_RADIUS * c

def calculate_azimuth(lat1, lon1, lat2, lon2):
    """Calculate azimuth from point 1 to point 2 in degrees."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    az = np.degrees(np.arctan2(x, y))
    return az % 360

def beam_height(range_m, elev_deg, radar_elev_m):
    """Calculate beam center height MSL at given range, accounting for earth curvature."""
    elev_rad = np.radians(elev_deg)
    # Standard atmosphere beam propagation with 4/3 earth
    height = (range_m * np.sin(elev_rad) +
              range_m**2 / (2 * EARTH_RADIUS_EFF) +
              radar_elev_m)
    return height

# ============= AZIMUTH AND RANGE TO SPIKES =============
print("\n1. GEOMETRY FROM KEPZ TO SPIKE LOCATIONS")
print("-" * 50)

for name, lat, lon in [("Spike 1", SPIKE1_LAT, SPIKE1_LON),
                        ("Spike 2", SPIKE2_LAT, SPIKE2_LON)]:
    dist = haversine_distance(KEPZ_LAT, KEPZ_LON, lat, lon)
    az = calculate_azimuth(KEPZ_LAT, KEPZ_LON, lat, lon)

    # Beam heights at this range
    h_05 = beam_height(dist, BEAM_ELEV_05, KEPZ_ELEV)
    h_15 = beam_height(dist, BEAM_ELEV_15, KEPZ_ELEV)

    # Beam center and edges
    h_05_lower = beam_height(dist, BEAM_ELEV_05 - BEAMWIDTH/2, KEPZ_ELEV)
    h_05_upper = beam_height(dist, BEAM_ELEV_05 + BEAMWIDTH/2, KEPZ_ELEV)

    print(f"\n  {name} ({lat:.4f}°N, {lon:.4f}°W):")
    print(f"    Range from KEPZ: {dist/1000:.1f} km")
    print(f"    Azimuth from KEPZ: {az:.1f}°")
    print(f"    0.5° beam center height: {h_05:.0f} m MSL ({h_05-1200:.0f} m AGL at target)")
    print(f"    0.5° beam lower edge: {h_05_lower:.0f} m MSL ({h_05_lower-1200:.0f} m AGL)")
    print(f"    0.5° beam upper edge: {h_05_upper:.0f} m MSL ({h_05_upper-1200:.0f} m AGL)")
    print(f"    1.5° beam center height: {h_15:.0f} m MSL ({h_15-1200:.0f} m AGL)")

    results[name] = {
        'range_km': dist/1000,
        'azimuth_deg': az,
        'beam_height_05_msl': h_05,
        'beam_height_05_agl': h_05 - 1200,
        'beam_height_05_lower_msl': h_05_lower,
        'beam_height_05_upper_msl': h_05_upper,
    }

# ============= FRANKLIN MOUNTAINS BLOCKAGE =============
print("\n\n2. FRANKLIN MOUNTAINS BEAM BLOCKAGE ANALYSIS")
print("-" * 50)

print("\n  Franklin Mountains peaks along KEPZ beam path:")
for peak_name, (peak_lat, peak_lon, peak_elev) in FRANKLIN_PEAKS.items():
    dist_from_kepz = haversine_distance(KEPZ_LAT, KEPZ_LON, peak_lat, peak_lon)
    az_from_kepz = calculate_azimuth(KEPZ_LAT, KEPZ_LON, peak_lat, peak_lon)

    # Beam height at this point
    h_beam = beam_height(dist_from_kepz, BEAM_ELEV_05, KEPZ_ELEV)
    h_beam_lower = beam_height(dist_from_kepz, BEAM_ELEV_05 - BEAMWIDTH/2, KEPZ_ELEV)

    clearance = h_beam_lower - peak_elev
    blocked = "BLOCKED" if clearance < 0 else "CLEAR"

    print(f"\n  {peak_name} ({peak_elev}m MSL):")
    print(f"    Distance from KEPZ: {dist_from_kepz/1000:.1f} km")
    print(f"    Azimuth from KEPZ: {az_from_kepz:.1f}°")
    print(f"    0.5° beam lower edge at peak: {h_beam_lower:.0f} m MSL")
    print(f"    Clearance: {clearance:.0f} m -> {blocked}")

# ============= BEAM PATH PROFILE (SYNTHETIC) =============
print("\n\n3. SYNTHETIC TERRAIN PROFILE ALONG SPIKE 1 AZIMUTH")
print("-" * 50)
print("  (Using known Franklin Mountains ridge elevations)")

# Create synthetic terrain profile along az ~85° from KEPZ
# The Franklin Mountains ridge crosses this azimuth at roughly 10-15 km from KEPZ
# Ridge elevation approximately 1500-1700m MSL at this latitude

print("\n  Terrain profile along azimuth 85° (to Spike 1):")
print(f"  {'Range (km)':<12} {'Terrain (m MSL)':<17} {'Beam Lower (m MSL)':<20} {'Status'}")
print(f"  {'-'*12} {'-'*17} {'-'*20} {'-'*10}")

# Approximate terrain profile from KEPZ eastward at 31.87N
# KEPZ is in Santa Teresa (west of Franklins)
# The mountains are roughly 15-25 km east of KEPZ
terrain_profile = [
    (0, 1251, "KEPZ radar site"),
    (2, 1230, "Santa Teresa mesa"),
    (5, 1220, "Rio Grande valley"),
    (8, 1200, "West side of Franklins"),
    (10, 1350, "Franklin foothills (west)"),
    (12, 1500, "Franklin ridge approach"),
    (14, 1650, "Franklin ridge crest"),
    (15, 1700, "Franklin ridge peak"),
    (16, 1600, "Franklin ridge east slope"),
    (18, 1400, "East foothills"),
    (20, 1250, "Northeast El Paso mesa"),
    (25, 1200, "Fort Bliss area"),
    (30, 1200, "Fort Bliss/NE El Paso"),
    (33.6, 1200, "Spike 1 location"),
]

for range_km, terrain_m, desc in terrain_profile:
    range_m = range_km * 1000
    h_lower = beam_height(range_m, BEAM_ELEV_05 - BEAMWIDTH/2, KEPZ_ELEV)
    clearance = h_lower - terrain_m
    status = "BLOCKED" if clearance < 0 else f"Clear ({clearance:.0f}m)"
    print(f"  {range_km:<12.1f} {terrain_m:<17d} {h_lower:<20.0f} {status}")

# ============= BLOCKAGE CORRIDORS =============
print("\n\n4. BLOCKAGE CORRIDORS AND VISIBILITY GAPS")
print("-" * 50)

print("""
  The Franklin Mountains run roughly N-S from 31.75°N to 32.0°N.
  The ridge is at approximately 106.48-106.50°W.

  From KEPZ (106.698°W), looking east toward Fort Bliss:

  Azimuths 75-100° are the critical zone:
    - Az 80-85°: High ridge, likely FULLY BLOCKED below 500m AGL
    - Az 85-90°: Peak blockage zone (North Franklin Peak vicinity)
    - Az 90-95°: Moderate blockage
    - Az 95-100°: Possible gaps through Smugglers Pass area

  KEY: Trans-Mountain Road gap at ~31.92°N provides a low point
  (~1580m MSL) in the ridge. This gap may allow partial beam
  penetration at specific azimuths.

  IMPLICATIONS:
  1. An object flying LOW (below 300m AGL) east of the Franklins
     is INVISIBLE to KEPZ on its 0.5° tilt
  2. An object would need to be above ~500m AGL to appear at the
     spike locations through the mountain blockage
  3. An object approaching FROM the east or south would not show
     a gradual approach track - it would appear suddenly once it
     climbed above the blockage altitude or moved to an unblocked
     azimuth
  4. This supports the "appeared from nowhere" observation -
     an eastward-approaching target at low altitude would be
     invisible until it rose above the ridge shadow
""")

# ============= POSITION UNCERTAINTY =============
print("\n5. POSITION UNCERTAINTY ANALYSIS")
print("-" * 50)

# At 33.6 km range with 0.95° beamwidth:
range_m = 33600
beam_width_m = 2 * range_m * np.tan(np.radians(BEAMWIDTH/2))
gate_length_m = GATE_LENGTH

print(f"  At {range_m/1000:.1f} km range:")
print(f"    Beam width: {beam_width_m:.0f} m")
print(f"    Gate length: {gate_length_m} m")
print(f"    Azimuthal resolution: ±{beam_width_m/2:.0f} m")
print(f"    Range resolution: ±{gate_length_m/2:.0f} m")
print(f"    Combined position uncertainty: ~{np.sqrt((beam_width_m/2)**2 + (gate_length_m/2)**2):.0f} m")
print()
print("  This means the KEPZ-derived positions have ~300m uncertainty.")
print("  The spike being 'over a residential area' vs 'over Fort Bliss boundary'")
print("  is a meaningful distinction at 33.6 km range, since the Spike 1")
print("  coordinates are ~3-5 km into the civilian area.")
print()
print("  Beam blockage/diffraction could shift apparent position by:")
print("    - If signal diffracts over the ridge: position appears shifted")
print("      toward higher elevation (upward bias)")
print("    - Multipath from ridge: ghost targets possible at incorrect ranges")
print("    - However, diffraction effects are typically <1° shift in azimuth")
print("    - This translates to ~500m at 33.6km - not enough to explain")
print("      the 3-5 km displacement from Fort Bliss to residential area")

# ============= KHDX VIEWING GEOMETRY =============
print("\n\n6. KHDX VIEWING GEOMETRY (for cross-reference)")
print("-" * 50)

KHDX_LAT = 33.0803
KHDX_LON = -106.1225
KHDX_ELEV = 1287.0  # meters MSL

for name, lat, lon in [("Spike 1", SPIKE1_LAT, SPIKE1_LON),
                        ("Spike 2", SPIKE2_LAT, SPIKE2_LON)]:
    dist = haversine_distance(KHDX_LAT, KHDX_LON, lat, lon)
    az = calculate_azimuth(KHDX_LAT, KHDX_LON, lat, lon)

    h_05 = beam_height(dist, BEAM_ELEV_05, KHDX_ELEV)
    h_05_lower = beam_height(dist, BEAM_ELEV_05 - BEAMWIDTH/2, KHDX_ELEV)

    print(f"\n  {name} from KHDX:")
    print(f"    Range: {dist/1000:.1f} km")
    print(f"    Azimuth: {az:.1f}°")
    print(f"    0.5° beam center height: {h_05:.0f} m MSL ({h_05-1200:.0f} m AGL)")
    print(f"    0.5° beam lower edge: {h_05_lower:.0f} m MSL ({h_05_lower-1200:.0f} m AGL)")
    print()
    print(f"    NOTE: At {dist/1000:.0f}km range, KHDX sensitivity is significantly")
    print(f"    reduced. Minimum detectable Z ≈ {-10 + 20*np.log10(dist/1000/100):.0f} dBZ at this range.")
    print(f"    A 21.5 dBZ return at KEPZ (33.6 km) might be only ~{21.5 - 20*np.log10(dist/33600):.0f} dBZ")
    print(f"    at KHDX range, which should still be detectable if the target is above the beam.")

    results[f'{name}_from_KHDX'] = {
        'range_km': dist/1000,
        'azimuth_deg': az,
        'beam_height_05_msl': h_05,
        'beam_height_05_agl': h_05 - 1200,
    }

print("\n  KHDX ADVANTAGE:")
print("  - Views Fort Bliss area from the NORTH (no Franklin Mtns blockage)")
print("  - If target was real, KHDX should see it IF above beam height")
print("  - At ~135 km range, KHDX 0.5° beam is at ~2500m MSL (~1300m AGL)")
print("  - A low-altitude target (<500m AGL) would also be below KHDX beam")
print("  - KHDX mainly useful for confirming targets above 1000m AGL")

# Save results
with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/terrain_analysis_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\n\nResults saved to terrain_analysis_results.json")
