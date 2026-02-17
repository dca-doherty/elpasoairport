#!/usr/bin/env python3
"""
IRIS Seismic/Infrasound Data Download via ObsPy
=================================================
Downloads waveform data from IRIS FDSN web services for seismic stations
near El Paso, TX / Fort Bliss (31.87N, 106.52W) during Feb 8-12, 2026.

Purpose: Detect potential explosions, acoustic energy events, or "booming sounds"
reported in the El Paso area on Feb 8-9, 2026.

Data source: IRIS/EarthScope FDSN Web Services
    - Station service: https://service.iris.edu/fdsnws/station/1/
    - Dataselect service: https://service.iris.edu/fdsnws/dataselect/1/
    - Timeseries service: https://service.iris.edu/irisws/timeseries/1/
    - MUSTANG metrics: https://service.iris.edu/mustang/

Requirements: pip install obspy numpy matplotlib
"""

import os
import sys
import json
import traceback
from datetime import datetime, timedelta

import numpy as np

try:
    from obspy import UTCDateTime, Stream, Inventory
    from obspy.clients.fdsn import Client
    from obspy.clients.fdsn.header import FDSNNoDataException, FDSNException
    from obspy.geodetics import gps2dist_azimuth
    HAS_OBSPY = True
except ImportError:
    HAS_OBSPY = False
    print("WARNING: ObsPy not installed. Install with: pip install obspy")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available; plots will be skipped.")

# ============================================================================
# CONFIGURATION
# ============================================================================
OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/seismic_data/'
PLOT_DIR = os.path.join(OUTPUT_DIR, 'plots')
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# Target location: Fort Bliss / El Paso, TX
TARGET_LAT = 31.87
TARGET_LON = -106.52
MAX_RADIUS_KM = 200

# Time windows of interest
EVENT_WINDOWS = {
    'booming_feb8_night': {
        'start': UTCDateTime("2026-02-08T00:00:00"),
        'end': UTCDateTime("2026-02-09T12:00:00"),
        'label': 'Booming sounds Feb 8-9',
    },
    'booming_peak_feb9': {
        'start': UTCDateTime("2026-02-09T02:00:00"),
        'end': UTCDateTime("2026-02-09T08:00:00"),
        'label': 'Booming peak Feb 9 night',
    },
    'event_window_feb11': {
        'start': UTCDateTime("2026-02-11T03:00:00"),
        'end': UTCDateTime("2026-02-11T07:00:00"),
        'label': 'Event window Feb 11 (03-07 UTC)',
    },
    'baseline_quiet': {
        'start': UTCDateTime("2026-02-10T12:00:00"),
        'end': UTCDateTime("2026-02-10T16:00:00"),
        'label': 'Baseline quiet period Feb 10',
    },
}

# Priority stations to try (real stations, not synthetics)
# Ordered by approximate distance from El Paso
PRIORITY_STATIONS = [
    # (network, station, description, approx_dist_km)
    ('EP', 'KIDD', 'UTEP Kidd Seismic Observatory, El Paso', 11),
    ('US', 'MNTX', 'Cornudas Mountains, TX (USGS)', 109),
    ('PN', 'PPDHS', 'Deming High School, NM', 125),
    ('SC', '121A', 'Cookes Peak, Deming, NM', 140),
    ('JW', 'BRR', 'BRR, NM', 140),
    ('TX', 'PB59', 'Pine Springs TxDot, TX (TexNet)', 161),
    ('TX', 'VHRN', 'Van Horn, TX (TexNet)', 188),
    ('TX', 'PB09', 'Culberson County, TX (TexNet)', 200),
    ('28', 'MVSS', 'Mountain Village, NM', 187),
    ('ZU', 'NMB23', 'Franklin Mountains, NM', 20),
    ('ZU', 'NMB22', 'West Potrillo Mountains, NM', 60),
    ('ZU', 'TXC24', 'San Felipe Arroyo, TX', 55),
]

# Preferred channel codes for different purposes
# BHZ = broadband seismometer vertical (20-40 sps)
# HHZ = high broadband vertical (80-100 sps)
# BDF = infrasound (microbarometer, 40 sps)
# BDO = infrasound (microbarometer, 40 sps, alternative)
# LHZ = long-period (1 sps)
PREFERRED_CHANNELS = ['BHZ', 'HHZ', 'BDF', 'BDO', 'EHZ', 'SHZ', 'LHZ']


def print_header(msg):
    print(f"\n{'='*80}")
    print(msg)
    print('='*80)


def compute_distance_km(lat1, lon1, lat2, lon2):
    """Compute distance in km between two lat/lon points."""
    if HAS_OBSPY:
        dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
        return dist_m / 1000.0
    else:
        # Approximate
        dlat = (lat2 - lat1) * 111.0
        dlon = (lon2 - lon1) * 111.0 * np.cos(np.radians(lat1))
        return np.sqrt(dlat**2 + dlon**2)


# ============================================================================
# STEP 1: DISCOVER STATIONS USING OBSPY FDSN CLIENT
# ============================================================================
def discover_stations(client, target_lat, target_lon, max_radius_km):
    """
    Query IRIS FDSN station service for all stations within radius.
    Returns an ObsPy Inventory object.
    """
    print_header("STEP 1: Station Discovery via IRIS FDSN")
    print(f"  Center: {target_lat}N, {target_lon}W")
    print(f"  Radius: {max_radius_km} km")

    max_radius_deg = max_radius_km / 111.0

    try:
        inventory = client.get_stations(
            latitude=target_lat,
            longitude=target_lon,
            maxradius=max_radius_deg,
            level='channel',
            starttime=UTCDateTime("2026-02-01"),
            endtime=UTCDateTime("2026-02-28"),
        )

        # Count and print
        total_nets = len(inventory)
        total_sta = sum(len(net) for net in inventory)
        total_cha = sum(
            len(sta.channels)
            for net in inventory
            for sta in net
        )

        print(f"\n  Results: {total_nets} networks, {total_sta} stations, {total_cha} channels")

        # List stations with distance
        station_list = []
        for net in inventory:
            for sta in net:
                dist = compute_distance_km(target_lat, target_lon,
                                           sta.latitude, sta.longitude)
                channels = [ch.code for ch in sta.channels]
                unique_channels = sorted(set(channels))

                station_list.append({
                    'network': net.code,
                    'station': sta.code,
                    'latitude': sta.latitude,
                    'longitude': sta.longitude,
                    'elevation_m': sta.elevation,
                    'name': sta.site.name if sta.site else '',
                    'channels': unique_channels,
                    'dist_km': dist,
                })

        # Sort by distance
        station_list.sort(key=lambda x: x['dist_km'])

        # Print top stations (skip SY=synthetic)
        print(f"\n  Closest real stations:")
        real_stations = [s for s in station_list if s['network'] != 'SY']
        for s in real_stations[:20]:
            chans = ','.join(s['channels'][:6])
            extra = f" +{len(s['channels'])-6}" if len(s['channels']) > 6 else ""
            print(f"    {s['network']:>3s}.{s['station']:<6s} {s['dist_km']:6.1f}km  "
                  f"{s['name'][:35]:<35s}  [{chans}{extra}]")

        # Save
        with open(os.path.join(OUTPUT_DIR, 'stations_obspy.json'), 'w') as f:
            json.dump(station_list, f, indent=2, default=str)
        print(f"\n  Saved station list to stations_obspy.json")

        return inventory, station_list

    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        return None, []


# ============================================================================
# STEP 2: DOWNLOAD WAVEFORM DATA
# ============================================================================
def download_waveforms(client, station_list, windows):
    """
    Download waveform data for priority stations during key time windows.
    """
    print_header("STEP 2: Waveform Data Download")

    results = {}

    # Build list of (net, sta) pairs to try, prioritizing our PRIORITY_STATIONS
    stations_to_try = []
    tried = set()

    # First: priority stations
    for net, sta, desc, dist in PRIORITY_STATIONS:
        key = f"{net}.{sta}"
        if key not in tried:
            tried.add(key)
            # Find matching entry in station_list for channels
            match = next((s for s in station_list
                          if s['network'] == net and s['station'] == sta), None)
            channels = match['channels'] if match else PREFERRED_CHANNELS
            stations_to_try.append((net, sta, desc, channels))

    # Second: closest non-synthetic stations from discovery
    for s in station_list:
        key = f"{s['network']}.{s['station']}"
        if key not in tried and s['network'] != 'SY' and s['dist_km'] <= MAX_RADIUS_KM:
            tried.add(key)
            stations_to_try.append((s['network'], s['station'], s['name'], s['channels']))
            if len(stations_to_try) >= 20:
                break

    print(f"  Will attempt {len(stations_to_try)} stations")

    for net, sta, desc, channels in stations_to_try:
        print(f"\n  --- {net}.{sta} ({desc}) ---")

        # Choose best channel
        target_channel = None
        for pref in PREFERRED_CHANNELS:
            if pref in channels:
                target_channel = pref
                break
        if not target_channel:
            # Try wildcard BH?, HH?
            bh_chans = [c for c in channels if c.startswith('BH')]
            hh_chans = [c for c in channels if c.startswith('HH')]
            if bh_chans:
                target_channel = bh_chans[0]
            elif hh_chans:
                target_channel = hh_chans[0]
            elif channels:
                target_channel = channels[0]
            else:
                target_channel = 'BHZ'

        print(f"    Available channels: {', '.join(channels[:10])}")
        print(f"    Selected: {target_channel}")

        for win_key, win_info in windows.items():
            start = win_info['start']
            end = win_info['end']
            label = win_info['label']

            # For long windows, just get a 30-minute sample
            duration = end - start
            if duration > 3600:
                # Get the midpoint 30 minutes
                midpoint = start + duration / 2
                dl_start = midpoint - 900  # 15 min before midpoint
                dl_end = midpoint + 900    # 15 min after midpoint
            else:
                dl_start = start
                dl_end = end

            print(f"    {label}: {dl_start} to {dl_end} ...")

            try:
                st = client.get_waveforms(
                    network=net,
                    station=sta,
                    location='*',
                    channel=target_channel,
                    starttime=dl_start,
                    endtime=dl_end,
                )

                if len(st) > 0:
                    tr = st[0]
                    npts = tr.stats.npts
                    sr = tr.stats.sampling_rate
                    data = tr.data

                    mean_val = np.mean(data)
                    std_val = np.std(data)
                    peak_val = np.max(np.abs(data))
                    rms_val = np.sqrt(np.mean(data.astype(float)**2))

                    # Check for anomalous signals
                    spikes_5sigma = np.sum(np.abs(data - mean_val) > 5 * std_val) if std_val > 0 else 0
                    spikes_10sigma = np.sum(np.abs(data - mean_val) > 10 * std_val) if std_val > 0 else 0

                    result_key = f"{net}.{sta}.{target_channel}.{win_key}"
                    results[result_key] = {
                        'network': net,
                        'station': sta,
                        'channel': target_channel,
                        'window': win_key,
                        'label': label,
                        'start': str(dl_start),
                        'end': str(dl_end),
                        'npts': npts,
                        'sampling_rate': sr,
                        'mean': float(mean_val),
                        'std': float(std_val),
                        'peak': float(peak_val),
                        'rms': float(rms_val),
                        'spikes_5sigma': int(spikes_5sigma),
                        'spikes_10sigma': int(spikes_10sigma),
                        'data_available': True,
                    }

                    print(f"      SUCCESS: {npts} samples at {sr} Hz ({npts/sr:.1f}s)")
                    print(f"      Mean={mean_val:.2f} Std={std_val:.2f} Peak={peak_val:.2f} RMS={rms_val:.2f}")
                    if spikes_5sigma > 0:
                        print(f"      *** {spikes_5sigma} samples > 5-sigma, {spikes_10sigma} > 10-sigma ***")

                    # Save miniSEED
                    fname = f"{net}_{sta}_{target_channel}_{win_key}.mseed"
                    fpath = os.path.join(OUTPUT_DIR, fname)
                    st.write(fpath, format='MSEED')
                    print(f"      Saved: {fname}")

                    # Plot if matplotlib available
                    if HAS_MATPLOTLIB:
                        try:
                            fig_path = os.path.join(PLOT_DIR, f"{net}_{sta}_{target_channel}_{win_key}.png")
                            st.plot(outfile=fig_path, size=(1200, 400))
                            print(f"      Plot: {os.path.basename(fig_path)}")
                        except Exception as pe:
                            print(f"      Plot error: {pe}")
                else:
                    print(f"      Empty stream returned")

            except FDSNNoDataException:
                print(f"      No data available")
                results[f"{net}.{sta}.{target_channel}.{win_key}"] = {
                    'network': net, 'station': sta, 'channel': target_channel,
                    'window': win_key, 'data_available': False,
                    'reason': 'no_data',
                }

            except FDSNException as e:
                print(f"      FDSN error: {e}")
                results[f"{net}.{sta}.{target_channel}.{win_key}"] = {
                    'network': net, 'station': sta, 'channel': target_channel,
                    'window': win_key, 'data_available': False,
                    'reason': str(e)[:200],
                }

            except Exception as e:
                print(f"      Error: {e}")

        # Also try infrasound channel if available
        infra_channels = [c for c in channels if c in ('BDF', 'BDO', 'BDG', 'HDF')]
        if infra_channels:
            infra_chan = infra_channels[0]
            print(f"    Also trying infrasound channel: {infra_chan}")

            for win_key, win_info in list(windows.items())[:2]:  # Just first two windows
                start = win_info['start']
                end = win_info['end']
                duration = end - start
                if duration > 3600:
                    midpoint = start + duration / 2
                    dl_start = midpoint - 900
                    dl_end = midpoint + 900
                else:
                    dl_start = start
                    dl_end = end

                try:
                    st = client.get_waveforms(
                        network=net, station=sta, location='*',
                        channel=infra_chan, starttime=dl_start, endtime=dl_end,
                    )
                    if len(st) > 0:
                        tr = st[0]
                        print(f"      {infra_chan}: {tr.stats.npts} samples, "
                              f"peak={np.max(np.abs(tr.data)):.2f}")

                        fname = f"{net}_{sta}_{infra_chan}_{win_key}.mseed"
                        st.write(os.path.join(OUTPUT_DIR, fname), format='MSEED')
                        print(f"      Saved: {fname}")
                except FDSNNoDataException:
                    print(f"      {infra_chan}: No data")
                except Exception as e:
                    print(f"      {infra_chan}: Error - {e}")

    return results


# ============================================================================
# STEP 3: CHECK USGS EARTHQUAKE CATALOG
# ============================================================================
def check_earthquake_catalog(client):
    """Query USGS earthquake catalog via FDSN event service."""
    print_header("STEP 3: USGS Earthquake Catalog (via FDSN)")

    try:
        usgs_client = Client("https://earthquake.usgs.gov")
    except Exception:
        try:
            usgs_client = Client("USGS")
        except Exception as e:
            print(f"  Cannot connect to USGS FDSN: {e}")
            # Fall back to requests-based query
            print("  Using requests fallback for USGS catalog...")
            import requests
            for label, start, end in [
                ("Feb 8-9 (booming sounds)", "2026-02-08", "2026-02-10"),
                ("Feb 11 (event night)", "2026-02-11", "2026-02-12"),
                ("Feb 7-14 (extended)", "2026-02-07", "2026-02-14"),
            ]:
                print(f"\n  {label}:")
                try:
                    resp = requests.get(
                        "https://earthquake.usgs.gov/fdsnws/event/1/query",
                        params={
                            'format': 'geojson', 'starttime': start, 'endtime': end,
                            'latitude': TARGET_LAT, 'longitude': TARGET_LON,
                            'maxradiuskm': MAX_RADIUS_KM, 'minmagnitude': 0.0,
                        }, timeout=30)
                    if resp.status_code == 200:
                        data = resp.json()
                        features = data.get('features', [])
                        print(f"    Found {len(features)} events (via HTTP)")
                        for feat in features:
                            p = feat['properties']
                            g = feat['geometry']['coordinates']
                            print(f"      M{p.get('mag','?')} | {p.get('place','?')}")
                except Exception as re:
                    print(f"    HTTP fallback error: {re}")
            return []

    events_found = []

    for label, start, end in [
        ("Feb 8-9 (booming sounds)", "2026-02-08", "2026-02-10"),
        ("Feb 11 (event night)", "2026-02-11", "2026-02-12"),
        ("Feb 7-14 (extended)", "2026-02-07", "2026-02-14"),
    ]:
        print(f"\n  {label}:")
        try:
            cat = usgs_client.get_events(
                starttime=UTCDateTime(start),
                endtime=UTCDateTime(end),
                latitude=TARGET_LAT,
                longitude=TARGET_LON,
                maxradius=MAX_RADIUS_KM / 111.0,
                minmagnitude=0.0,
            )
            print(f"    Found {len(cat)} events")
            for ev in cat:
                origin = ev.preferred_origin() or ev.origins[0]
                mag_obj = ev.preferred_magnitude() or (ev.magnitudes[0] if ev.magnitudes else None)
                mag = mag_obj.mag if mag_obj else '?'
                desc = ev.event_descriptions[0].text if ev.event_descriptions else '?'
                print(f"      M{mag} | {origin.time} | {desc}")
                print(f"        ({origin.latitude:.3f}, {origin.longitude:.3f}) depth={origin.depth/1000:.1f}km")

                events_found.append({
                    'time': str(origin.time),
                    'magnitude': float(mag) if isinstance(mag, (int, float)) else None,
                    'latitude': origin.latitude,
                    'longitude': origin.longitude,
                    'depth_km': origin.depth / 1000 if origin.depth else None,
                    'description': desc,
                    'period': label,
                })

        except FDSNNoDataException:
            print(f"    No events found")
        except Exception as e:
            print(f"    Error: {e}")

    # Save
    with open(os.path.join(OUTPUT_DIR, 'usgs_events_obspy.json'), 'w') as f:
        json.dump(events_found, f, indent=2, default=str)

    return events_found


# ============================================================================
# STEP 4: GENERATE COMPARISON ANALYSIS
# ============================================================================
def compare_windows(results):
    """Compare signal levels between event and baseline windows."""
    print_header("STEP 4: Event vs. Baseline Comparison")

    # Group by station
    by_station = {}
    for key, r in results.items():
        if not r.get('data_available', False):
            continue
        sta_key = f"{r['network']}.{r['station']}.{r['channel']}"
        if sta_key not in by_station:
            by_station[sta_key] = {}
        by_station[sta_key][r['window']] = r

    for sta_key, windows in by_station.items():
        print(f"\n  {sta_key}:")

        baseline = windows.get('baseline_quiet')
        if not baseline:
            print(f"    No baseline window data; cannot compare")
            for wk, wd in windows.items():
                print(f"    {wk}: RMS={wd.get('rms', '?'):.2f}, "
                      f"Peak={wd.get('peak', '?'):.2f}, "
                      f"5s-spikes={wd.get('spikes_5sigma', '?')}")
            continue

        base_rms = baseline['rms']
        base_std = baseline['std']

        for wk, wd in windows.items():
            if wk == 'baseline_quiet':
                continue
            rms_ratio = wd['rms'] / base_rms if base_rms > 0 else float('inf')
            std_ratio = wd['std'] / base_std if base_std > 0 else float('inf')
            print(f"    {wd['label']}:")
            print(f"      RMS ratio vs baseline: {rms_ratio:.2f}x")
            print(f"      Std ratio vs baseline: {std_ratio:.2f}x")
            print(f"      5-sigma spikes: {wd['spikes_5sigma']}")
            if rms_ratio > 2.0:
                print(f"      *** ELEVATED SIGNAL (>{rms_ratio:.1f}x baseline) ***")


# ============================================================================
# STEP 5: WEB SERVICE URL REFERENCE
# ============================================================================
def print_web_service_urls():
    """Print reference URLs for manual data access."""
    print_header("REFERENCE: IRIS/FDSN Web Service URLs")

    print("""
  FDSN Station Service (metadata):
    https://service.iris.edu/fdsnws/station/1/query
    https://service.earthscope.org/fdsnws/station/1/query

  FDSN Dataselect Service (miniSEED waveforms):
    https://service.iris.edu/fdsnws/dataselect/1/query
    https://service.earthscope.org/fdsnws/dataselect/1/query

  IRIS Timeseries Service (processed waveforms, plots):
    https://service.iris.edu/irisws/timeseries/1/query

  IRIS Timeseriesplot Service (image plots):
    https://service.iris.edu/irisws/timeseriesplot/1/query

  MUSTANG Quality Metrics:
    https://service.iris.edu/mustang/metrics/1/query
    https://service.iris.edu/mustang/measurements/1/query

  MUSTANG Noise Mode Timeseries:
    https://service.iris.edu/mustang/noise-mode-timeseries/1/query

  USGS Earthquake Catalog (FDSN Event):
    https://earthquake.usgs.gov/fdsnws/event/1/query

  SPUD Infrasound Events:
    https://ds.iris.edu/spud/infrasoundevent

  TexNet Earthquake Catalog:
    https://catalog.texnet.beg.utexas.edu/

  IRIS GMap (interactive station map):
    https://ds.iris.edu/gmap/

  Example URLs for El Paso area:

    Station query (200km radius):
      https://service.iris.edu/fdsnws/station/1/query?lat=31.87&lon=-106.52&maxradius=1.8&level=station&format=text

    KIDD waveform (miniSEED):
      https://service.iris.edu/fdsnws/dataselect/1/query?net=EP&sta=KIDD&loc=--&cha=BHZ&start=2026-02-09T04:00:00&end=2026-02-09T04:30:00

    KIDD waveform plot:
      https://service.iris.edu/irisws/timeseries/1/query?net=EP&sta=KIDD&loc=--&cha=BHZ&start=2026-02-09T04:00:00&end=2026-02-09T04:30:00&output=plot

    MNTX waveform:
      https://service.iris.edu/fdsnws/dataselect/1/query?net=US&sta=MNTX&loc=00&cha=BHZ&start=2026-02-09T04:00:00&end=2026-02-09T04:30:00

    USGS events near El Paso, Feb 8-12:
      https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&starttime=2026-02-08&endtime=2026-02-12&latitude=31.87&longitude=-106.52&maxradiuskm=200

    Infrasound station search (500km):
      https://service.iris.edu/fdsnws/station/1/query?lat=31.87&lon=-106.52&maxradius=5&channel=BDF,BDO&level=channel&format=text
""")


# ============================================================================
# MAIN
# ============================================================================
def main():
    if not HAS_OBSPY:
        print("ERROR: ObsPy is required. Install with: pip install obspy")
        print("Printing reference URLs only...\n")
        print_web_service_urls()
        return

    print_header("IRIS SEISMIC/INFRASOUND DATA DOWNLOAD")
    print(f"  Target: Fort Bliss / El Paso, TX ({TARGET_LAT}N, {abs(TARGET_LON)}W)")
    print(f"  Radius: {MAX_RADIUS_KM} km")
    print(f"  Time: Feb 8-12, 2026")
    print(f"  Output: {OUTPUT_DIR}")

    # Initialize FDSN clients
    # Note: IRIS station service works for metadata, but waveform downloads
    # require the EarthScope endpoint (service.earthscope.org) as IRIS
    # redirects dataselect there. We use separate clients for each.
    print("\n  Connecting to FDSN web services...")

    # Station metadata client (IRIS)
    try:
        station_client = Client("IRIS")
        print("  Station metadata: IRIS (service.iris.edu)")
    except Exception as e:
        print(f"  IRIS station client failed: {e}")
        station_client = None

    # Waveform data client (EarthScope)
    try:
        client = Client("https://service.earthscope.org")
        print("  Waveform data: EarthScope (service.earthscope.org)")
    except Exception as e:
        print(f"  EarthScope connection failed: {e}")
        print_web_service_urls()
        return

    # Use station_client for discovery if available, else fall back to EarthScope
    discovery_client = station_client if station_client else client

    # Step 1: Discover stations
    inventory, station_list = discover_stations(discovery_client, TARGET_LAT, TARGET_LON, MAX_RADIUS_KM)

    # Step 2: Download waveforms
    if station_list:
        results = download_waveforms(client, station_list, EVENT_WINDOWS)

        # Save results
        with open(os.path.join(OUTPUT_DIR, 'waveform_results.json'), 'w') as f:
            json.dump(results, f, indent=2, default=str)
    else:
        results = {}

    # Step 3: Check earthquake catalog
    events = check_earthquake_catalog(client)

    # Step 4: Compare windows
    if results:
        compare_windows(results)

    # Step 5: Print reference URLs
    print_web_service_urls()

    # Final summary
    print_header("SUMMARY")
    total_data = sum(1 for r in results.values() if r.get('data_available', False))
    total_nodata = sum(1 for r in results.values() if not r.get('data_available', False))
    print(f"  Stations discovered: {len(station_list)}")
    print(f"  Waveform downloads: {total_data} successful, {total_nodata} no-data")
    print(f"  USGS events found: {len(events)}")
    print(f"  Output directory: {OUTPUT_DIR}")

    # Save complete summary
    summary = {
        'search_params': {
            'center_lat': TARGET_LAT,
            'center_lon': TARGET_LON,
            'radius_km': MAX_RADIUS_KM,
        },
        'stations_found': len(station_list),
        'waveforms_downloaded': total_data,
        'waveforms_no_data': total_nodata,
        'usgs_events': len(events),
        'events': events,
        'waveform_results': results,
        'priority_stations': [
            {'network': n, 'station': s, 'description': d, 'dist_km': dist}
            for n, s, d, dist in PRIORITY_STATIONS
        ],
        'key_station_details': {
            'EP.KIDD': {
                'description': 'UTEP Kidd Seismic Observatory',
                'location': 'University of Texas at El Paso campus',
                'lat': 31.771774, 'lon': -106.506378,
                'dist_km': 11,
                'channels': ['BHZ', 'BHN', 'BHE', 'HHZ', 'HHN', 'HHE',
                              'HNZ', 'HNN', 'HNE', 'LHZ', 'LHN', 'LHE'],
                'sensor': 'CMG-3T 120s broadband + RT131A accelerometer',
                'sample_rates': {'BH?': 40, 'HH?': 100, 'HN?': 200, 'LH?': 1},
            },
            'US.MNTX': {
                'description': 'Cornudas Mountains, TX (USGS ANSS)',
                'lat': 31.6985, 'lon': -105.3821,
                'dist_km': 109,
                'channels': ['BHZ', 'BHN', 'BHE', 'HHZ', 'HHN', 'HHE',
                              'LHZ', 'LHN', 'LHE', 'HNZ'],
                'sensor': 'Broadband + accelerometer',
                'sample_rates': {'BH?': 40, 'HH?': 100, 'LH?': 1},
            },
        },
        'infrasound_notes': {
            'nearest_infrasound_on_iris': [
                'SC.121A BDF (Cookes Peak, Deming NM, ~140km) - Transportable Array legacy',
                'IU.ANMO BDF (Albuquerque NM, ~350km) - IRIS/USGS global station',
                'N4.MSTX BDF (Midland TX, ~500km) - Central/East TX, far from El Paso',
            ],
            'ims_arrays': [
                'IS57 Pinon Flat, CA (~800km) - too far for local event detection',
                'IS56 Newport, WA (~2000km) - too far',
            ],
            'lanl_arrays': (
                'LANL operates infrasound arrays in NM but data may not be on IRIS. '
                'Network code LN (Los Alamos Seismic Network). '
                'Contact LANL Seismoacoustics group for access.'
            ),
            'txar': (
                'TXAR/PS46 (Lajitas, TX, ~400km SE) is a primary seismic array operated '
                'by SMU. Check SMU Geophysics for infrasound data.'
            ),
        },
        'analysis_timestamp': str(datetime.utcnow()),
    }

    with open(os.path.join(OUTPUT_DIR, 'seismic_analysis_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n  Summary saved to seismic_analysis_summary.json")


if __name__ == '__main__':
    main()
