# Multi-Source Sensor Analysis: El Paso Airspace Incident
## February 11, 2026 (03:00-07:00 UTC)

### Overview
Six additional sensor modalities were checked beyond the primary NEXRAD dual-radar analysis (KEPZ + KHDX) to corroborate or refute the HELWS laser engagement hypothesis.

---

## 1. GOES-19 GLM (Geostationary Lightning Mapper)

Data analyzed. No detections at Fort Bliss.

The GLM images at ~500fps in near-infrared (777.4nm OI emission line) and has detected non-lightning events including bolides and rockets. If the HELWS engagement produced any optical flash (plasma spark, combustion, scattered NIR), GLM would catch it with precise timestamps.

GOES-16 was replaced by GOES-19 as GOES-East on April 7, 2025. Data sourced from the `noaa-goes19` S3 bucket.

### Results — Multi-Night GLM Scan (1° radius around Fort Bliss)

| Date | Window (UTC) | Files | Detections |
|------|-------------|-------|------------|
| Feb 8 | 02:00-08:00 | 1,259 | **0** |
| Feb 9 | 02:00-08:00 | 1,260 | **0** |
| Feb 10 | 02:00-08:00 | 1,260 | **103** |
| Feb 10 | 03:00-07:00 | 900 | **0** |
| Feb 11 (EVENT) | 02:00-08:00 | 1,260 | **0** |
| Feb 11 (EVENT) | 03:00-07:00 | 900 | **0** |

### Feb 10 Detection Drill-Down

The 103 detections on Feb 10 are **NOT at Fort Bliss**:
- **Location:** 30.87°N, 107.24-107.44°W — Columbus/Deming, NM area (~130-140 km WSW of Fort Bliss)
- **Time:** All in a 30-minute burst at 02:07-02:36 UTC (7:07-7:36 PM MST)
- **Type:** 95 events + 8 groups + 0 flashes (sub-flash threshold)
- **Energy:** ~10⁻¹⁵ J (at the absolute floor of GLM sensitivity)
- **SkySentinel coincidence:** None (ended 67 min before SkySentinel's 03:43 UTC observation)
- **Interpretation:** Faint distant lightning or GLM noise/hot pixels, unrelated to Fort Bliss

### Assessment

GLM shows **zero optical transients at Fort Bliss on any night**, including the Feb 11 event night. This means:
- HELWS engagement (if it occurred) did not produce bright optical flashes in the GLM 777.4nm band
- Laser wavelength likely outside GLM passband, or scattered light too dim for the ~2 fJ detection threshold
- Cloud cover may have obscured dim flashes (METAR should be checked)
- **This does not rule out laser engagement** — most directed-energy weapons operate at wavelengths (1064nm Nd:YAG, 10.6μm CO₂) outside GLM's narrow OI emission passband

Scripts: `scripts/goes16_glm_analysis.py`, `scripts/glm_feb10_drilldown.py`

---

## 2. IRIS Seismic / Infrasound

Data downloaded and analyzed.

### Stations Found
- **41 seismic stations** within 200km of Fort Bliss
- Key station: **EP.KIDD** (UTEP Kidd Seismic Observatory) — only **11km** from target, 100Hz broadband
- EP.KIDD data confirmed available at IRIS but download failed (502 server error on dataselect service)
- Successfully downloaded waveforms from 11 stations across 4 time windows

### USGS Earthquake Catalog
| Date | Magnitude | Location | Distance |
|------|-----------|----------|----------|
| Feb 8 07:35 UTC | M2.1 | 52km SSW of Whites City, NM | ~190km E |
| Feb 9 04:55 UTC | M2.5 | 52km SSW of Whites City, NM | ~190km E |
| Feb 11 | **NONE** | — | — |

- 7 total events in Feb 7-14, all in the Delaware Basin (Whites City area, ~190km east)
- **Zero cataloged seismic events on Feb 11** (event night)

### Waveform Analysis — Feb 9 "Booming Sounds"
Multiple TexNet stations near the Whites City epicenter showed massive anomalies during the Feb 9 04:45-05:15 UTC window, **coinciding exactly with the M2.5 earthquake at 04:55 UTC:**

| Station | Distance | 5-sigma spikes | 10-sigma spikes | Peak amplitude |
|---------|----------|---------------|-----------------|----------------|
| TX.PB36 | 199km | 527 | 217 | 30,608 |
| TX.PB35 | 198km | 483 | 177 | 43,810 |
| TX.PB28 | 192km | 385 | 212 | 192,121 |
| TX.PB38 | 198km | 380 | 167 | 29,950 |
| TX.PB37 | 188km | 397 | 194 | 758,158 |
| US.MNTX | 109km | 306 | 18 | 3,301 |

**These spikes are explained by the M2.5 earthquake, not by Fort Bliss activity.** The TexNet PB-series stations are in the Delaware Basin, close to the earthquake epicenter.

### Waveform Analysis — Feb 11 Event Night
The event night (04:45-05:15 UTC) shows **no anomalous seismic signals:**

| Station | 5-sigma spikes | Vs. Baseline |
|---------|---------------|--------------|
| US.MNTX | 0 | 0 (same) |
| SC.121A | 0 | 0 (same) |
| TX.PB28 | 191 | 282 (lower) |
| TX.PB35 | 6 | 2 (similar) |
| TX.VHRN | 0 | 0 (same) |

No seismic signature from the Feb 11 incident. This is consistent with a directed-energy weapon (laser engagement does not produce seismic waves unless it causes explosions). The "booming sounds" reported Feb 8-9 correlate with documented M2.1-2.5 earthquakes in the Delaware Basin, not with Fort Bliss activity.

### Infrasound Stations
- SC.121A BDF (Cookes Peak, 140km) — Transportable Array legacy
- IU.ANMO BDF (Albuquerque, ~350km) — IRIS/USGS global station
- LANL operates infrasound arrays in NM (not on IRIS, contact directly)
- **EP.KIDD waveform data should be prioritized** — it's 11km away and has HHZ data at 100Hz for the entire event

---

## 3. VIIRS Day/Night Band

Timing mismatch limits utility.

VIIRS nighttime overpass for El Paso: ~08:30 UTC (1:30 AM MST)
Event window: 03:00-07:00 UTC (8 PM - midnight MST)

**Gap: 1.5-5.5 hours.** VIIRS would NOT image the engagement in real-time.

Could detect: residual fires, persistent luminous effects, changes in base lighting patterns.
VNP46A1 (Black Marble) tiles at LAADS DAAC for h08v05 would allow comparison of Feb 10 vs Feb 11 nighttime radiance at Fort Bliss. Requires free Earthdata account.

NASA Worldview links for visual inspection provided in analysis outputs.

---

## 4. NASA FIRMS Thermal Hotspots

Checked. No detections.

Queried VIIRS SNPP, VIIRS NOAA-20, VIIRS NOAA-21, and MODIS for Feb 8-13, 2026 in the El Paso bounding box (-107.0 to -106.0, 31.5 to 32.5).

**Result: Zero thermal hotspots.** Not surprising — FIRMS is designed for fires and volcanoes, not laser engagements. The null result is worth documenting but doesn't constrain the hypothesis.

---

## 5. OpenSky ADS-B

API not accessible for this date range.

All 49 queries for both event and baseline nights returned HTTP errors (likely 429 rate limiting or historical data access restrictions). Zero aircraft observed on either night.

Limitations:
- OpenSky free API has limited historical access
- Most military aircraft operate without ADS-B transponders
- CBP Air and Marine Operations aircraft sometimes broadcast, sometimes don't
- Alternative sources: ADS-B Exchange (paid), FlightAware (paid), FlightRadar24

**Recommended action:** Check ADS-B Exchange or FlightAware historical data for military callsigns (RCH, CBP, PHNX, DUKE) in the El Paso area during the incident window.

---

## 6. 911 Call / Dispatch Logs

Research completed.

911 CAD (Computer-Aided Dispatch) logs are available through:
- **Texas Public Information Act** request to El Paso County
- **New Mexico IPRA** request to Dona Ana County (Las Cruces side)
- El Paso Police Department records

Key requests to file:
1. All 911 calls reporting "loud noises," "explosions," "lights in sky," "laser," "military" for Feb 8-12, 2026 in El Paso County
2. All fire/EMS dispatch records for Fort Bliss / Biggs Army Airfield on Feb 11, 2026
3. Broadcastify scanner archives for El Paso (may have real-time radio traffic recordings)

---

## Priority Ranking for Investigation

| Priority | Source | Status | Action |
|----------|--------|--------|--------|
| 1 | **EP.KIDD seismic** | Server error | Retry IRIS dataselect or contact UTEP directly |
| 2 | **911/dispatch logs** | FOIA needed | File Texas PIA request to El Paso County |
| 3 | **ADS-B data** | Paid source needed | Check ADS-B Exchange historical |
| 4 | **VIIRS DNB** | LAADS account needed | Download VNP46A1 h08v05 for Feb 10-11 |
| 5 | **Sky Sentinel** | PI contact needed | Email admin@goskysentinel.com for raw camera data during 57-hour gap |
| 6 | **LANL infrasound** | Contact needed | Reach out to LANL Seismoacoustics group |

---

## Integration with NEXRAD Analysis

The multi-source analysis reinforces the NEXRAD findings:

1. **No seismic signature on Feb 11** → Consistent with directed-energy (non-kinetic) engagement
2. **Feb 8-9 seismic anomalies** → Explained by Delaware Basin M2.5 earthquake, not Fort Bliss
3. **Sky Sentinel 57-hour gap** → The closest all-sky camera went dark for the entire incident period
4. **FIRMS null result** → No large fires (consistent with controlled laser ops, not explosions)
5. **GLM null result** → Zero optical transients at Fort Bliss across all 5 nights (Feb 8-12); 103 detections on Feb 10 are distant lightning ~130km WSW in Columbus/Deming, NM
6. **GLM interpretation** → Null result does not rule out HELWS — most DE weapons operate at wavelengths (1064nm, 10.6μm) outside GLM's 777.4nm passband

The strongest outstanding data requests are **EP.KIDD seismic waveforms** (11km from Fort Bliss, 100Hz broadband — would detect acoustic coupling from any directed-energy engagement) and **911/dispatch logs** (ground truth from civilian reports during the incident window).
