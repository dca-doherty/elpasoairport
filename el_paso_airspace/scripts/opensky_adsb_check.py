#!/usr/bin/env python3
"""
OpenSky Network ADS-B Aircraft Check
=====================================
Checks for military/government aircraft (USAF, CBP, DHS) in the El Paso
airspace during the Feb 11, 2026 incident window.

Notes on military ADS-B:
  - Many military aircraft operate with transponders OFF (Mode S only, no ADS-B Out)
  - However, support aircraft, tankers, and some spotters may broadcast ADS-B
  - CBP aircraft (Air and Marine Operations) often use ADS-B
  - Even partial data can show pattern of activity

Data source: OpenSky Network REST API and Trino historical database
"""

import requests
import json
import os
import time
from datetime import datetime, timezone

# ============================================================================
# CONFIGURATION
# ============================================================================
# El Paso / Fort Bliss bounding box
LAMIN, LAMAX = 31.0, 33.0
LOMIN, LOMAX = -107.5, -105.5

# Time windows (Unix timestamps)
EVENT_START = int(datetime(2026, 2, 11, 3, 0, tzinfo=timezone.utc).timestamp())
EVENT_END = int(datetime(2026, 2, 11, 7, 0, tzinfo=timezone.utc).timestamp())
BASELINE_START = int(datetime(2026, 2, 10, 3, 0, tzinfo=timezone.utc).timestamp())
BASELINE_END = int(datetime(2026, 2, 10, 7, 0, tzinfo=timezone.utc).timestamp())

# Known military/govt ICAO24 prefixes
# USAF: AE (common), other mil prefixes
# CBP: varies, often N-registered
MIL_CALLSIGN_PREFIXES = [
    'RCH',   # USAF airlift (Reach)
    'EVAC',  # USAF medevac
    'DUKE',  # USAF special ops
    'PHNX',  # CBP
    'CBP',   # CBP callsigns
    'OMAHA', # E-6B TACAMO
    'SNTRY', # AWACS
    'DRAGN', # USAF
    'VIPER', # USAF
    'COBRA', # Army
    'TOPCT', # Top Cat (various)
    'HAWK',  # Various military
    'IRON',  # Various military
    'TALON', # USAF
    'NIGHT', # Night operations
    'SKULL', # USAF
    'BATT',  # Battle management
    'DARKSTAR',
    'RAIDR', # B-1 / B-52
    'BONE',  # B-1
    'DEATH', # B-2
]

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/adsb_data/'
ANALYSIS_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 80)
print("OPENSKY ADS-B MILITARY AIRCRAFT CHECK")
print(f"Area: {LAMIN}-{LAMAX}°N, {LOMIN}-{LOMAX}°W")
print("=" * 80)

def query_opensky_states(begin, end, label):
    """Query OpenSky for aircraft states in our bounding box."""
    print(f"\n  Querying {label}...")
    all_aircraft = {}

    # OpenSky API limits: sample every 10 seconds for historical
    # We'll sample every 60 seconds to stay within rate limits
    current_time = begin
    step = 300  # 5-minute intervals
    queries_made = 0
    errors = 0

    while current_time <= end:
        url = "https://opensky-network.org/api/states/all"
        params = {
            'time': current_time,
            'lamin': LAMIN,
            'lamax': LAMAX,
            'lomin': LOMIN,
            'lomax': LOMAX,
        }

        try:
            resp = requests.get(url, params=params, timeout=15)
            queries_made += 1

            if resp.status_code == 200:
                data = resp.json()
                states = data.get('states', []) or []

                for state in states:
                    icao24 = state[0] if state[0] else ''
                    callsign = (state[1] or '').strip()
                    origin = state[2] if state[2] else ''
                    lat = state[6]
                    lon = state[5]
                    alt_m = state[7]  # barometric altitude
                    velocity = state[9]
                    heading = state[10]
                    on_ground = state[8]

                    aircraft_key = icao24 or callsign

                    if aircraft_key not in all_aircraft:
                        all_aircraft[aircraft_key] = {
                            'icao24': icao24,
                            'callsign': callsign,
                            'origin_country': origin,
                            'positions': [],
                            'is_military': False,
                        }

                    # Check if military
                    cs_upper = callsign.upper()
                    for prefix in MIL_CALLSIGN_PREFIXES:
                        if cs_upper.startswith(prefix):
                            all_aircraft[aircraft_key]['is_military'] = True
                            break

                    # US military origin countries
                    if origin in ['United States'] and (
                        icao24.startswith('ae') or icao24.startswith('af') or
                        icao24.startswith('a0') or icao24.startswith('~')
                    ):
                        all_aircraft[aircraft_key]['is_military'] = True

                    all_aircraft[aircraft_key]['positions'].append({
                        'time': current_time,
                        'time_utc': datetime.fromtimestamp(current_time, tz=timezone.utc).strftime('%H:%M:%S'),
                        'lat': lat,
                        'lon': lon,
                        'alt_m': alt_m,
                        'velocity_ms': velocity,
                        'heading': heading,
                        'on_ground': on_ground,
                    })

            elif resp.status_code == 429:
                print(f"    Rate limited at {queries_made} queries, waiting 10s...")
                time.sleep(10)
                continue
            elif resp.status_code == 404:
                # Historical data may not be available
                pass
            else:
                errors += 1

        except requests.exceptions.Timeout:
            errors += 1
        except Exception as e:
            errors += 1
            if queries_made < 3:
                print(f"    Error: {e}")

        current_time += step
        # Rate limiting: 1 query per second for anonymous users
        time.sleep(1.5)

    print(f"    Queries: {queries_made}, Errors: {errors}")
    print(f"    Unique aircraft seen: {len(all_aircraft)}")

    return all_aircraft

# ============================================================================
# QUERY EVENT AND BASELINE WINDOWS
# ============================================================================
event_aircraft = query_opensky_states(EVENT_START, EVENT_END, "Event Night (Feb 11, 03-07 UTC)")
baseline_aircraft = query_opensky_states(BASELINE_START, BASELINE_END, "Baseline Night (Feb 10, 03-07 UTC)")

# ============================================================================
# ANALYSIS
# ============================================================================
print(f"\n\n{'='*80}")
print(f"ADS-B ANALYSIS RESULTS")
print(f"{'='*80}")

# Event night
event_mil = {k: v for k, v in event_aircraft.items() if v['is_military']}
event_civ = {k: v for k, v in event_aircraft.items() if not v['is_military']}

print(f"\n  EVENT NIGHT (Feb 11, 03-07 UTC):")
print(f"    Total aircraft: {len(event_aircraft)}")
print(f"    Military/govt:  {len(event_mil)}")
print(f"    Civilian:       {len(event_civ)}")

if event_mil:
    print(f"\n    Military Aircraft Details:")
    for key, ac in event_mil.items():
        n_pos = len(ac['positions'])
        first_t = ac['positions'][0]['time_utc'] if ac['positions'] else '?'
        last_t = ac['positions'][-1]['time_utc'] if ac['positions'] else '?'
        alt_avg = sum(p['alt_m'] for p in ac['positions'] if p['alt_m']) / max(1, sum(1 for p in ac['positions'] if p['alt_m']))
        print(f"      {ac['callsign']:<12} ICAO:{ac['icao24']:<8} "
              f"Origin:{ac['origin_country']:<15} "
              f"Seen:{n_pos}x {first_t}-{last_t} Alt:{alt_avg:.0f}m")

# All aircraft detail
print(f"\n    All Aircraft on Event Night:")
for key, ac in sorted(event_aircraft.items(), key=lambda x: x[1].get('callsign', '')):
    n_pos = len(ac['positions'])
    mil = " [MIL]" if ac['is_military'] else ""
    print(f"      {ac['callsign']:<12} ICAO:{ac['icao24']:<8} "
          f"Origin:{ac['origin_country']:<15} Seen:{n_pos}x{mil}")

# Baseline comparison
baseline_mil = {k: v for k, v in baseline_aircraft.items() if v['is_military']}
baseline_civ = {k: v for k, v in baseline_aircraft.items() if not v['is_military']}

print(f"\n  BASELINE NIGHT (Feb 10, 03-07 UTC):")
print(f"    Total aircraft: {len(baseline_aircraft)}")
print(f"    Military/govt:  {len(baseline_mil)}")
print(f"    Civilian:       {len(baseline_civ)}")

# Compare
print(f"\n  COMPARISON:")
print(f"    Aircraft count: Event={len(event_aircraft)} vs Baseline={len(baseline_aircraft)}")
print(f"    Military count: Event={len(event_mil)} vs Baseline={len(baseline_mil)}")

if len(event_mil) > len(baseline_mil):
    print(f"    *** MORE military aircraft on event night ***")
elif len(event_mil) == len(baseline_mil) == 0:
    print(f"    No military ADS-B visible either night (expected — most mil ops run transponders off)")
else:
    print(f"    Normal or fewer military aircraft on event night")

# ============================================================================
# SAVE RESULTS
# ============================================================================
def serialize_aircraft(aircraft_dict):
    """Convert for JSON serialization."""
    out = {}
    for k, v in aircraft_dict.items():
        ac = dict(v)
        ac['n_positions'] = len(v['positions'])
        ac['positions'] = v['positions'][:10]  # Limit for JSON size
        out[k] = ac
    return out

output = {
    'search_area': {'lamin': LAMIN, 'lamax': LAMAX, 'lomin': LOMIN, 'lomax': LOMAX},
    'event_window': {
        'start_utc': datetime.fromtimestamp(EVENT_START, tz=timezone.utc).isoformat(),
        'end_utc': datetime.fromtimestamp(EVENT_END, tz=timezone.utc).isoformat(),
        'total_aircraft': len(event_aircraft),
        'military_aircraft': len(event_mil),
        'military_details': serialize_aircraft(event_mil),
        'all_callsigns': [v['callsign'] for v in event_aircraft.values() if v['callsign']],
    },
    'baseline_window': {
        'start_utc': datetime.fromtimestamp(BASELINE_START, tz=timezone.utc).isoformat(),
        'end_utc': datetime.fromtimestamp(BASELINE_END, tz=timezone.utc).isoformat(),
        'total_aircraft': len(baseline_aircraft),
        'military_aircraft': len(baseline_mil),
    },
    'verdict': 'military_excess' if len(event_mil) > len(baseline_mil) + 1 else 'no_military_adsb' if len(event_mil) == 0 else 'normal',
}

output_path = os.path.join(ANALYSIS_DIR, 'opensky_adsb_results.json')
with open(output_path, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"\nResults saved to {output_path}")

# Also save full aircraft list
full_path = os.path.join(OUTPUT_DIR, 'event_night_aircraft.json')
with open(full_path, 'w') as f:
    json.dump(serialize_aircraft(event_aircraft), f, indent=2, default=str)
print(f"Full aircraft data saved to {full_path}")
