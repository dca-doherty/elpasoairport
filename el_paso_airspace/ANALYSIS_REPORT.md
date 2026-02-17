# El Paso Airspace Closure: Multi Source Analysis Report

**Date:** February 17, 2026
**Incident Date:** February 10 through 11, 2026
**Location:** Fort Bliss / NE El Paso, TX


## Executive Summary

This report pulls together independent analysis of the February 11, 2026 El Paso airspace closure using several data sources: dual radar NEXRAD analysis (KEPZ and KHDX), upper air soundings, surface METAR observations, RCS estimation, terrain modeling, and open source intelligence. What we found is a set of anomalies that make the official "party balloon" explanation difficult to accept at face value, though some findings do align with parts of that story.

**Here is what stands out:**

1. **KHDX independently confirms radar returns** in the Fort Bliss sector at both KEPZ spike times. This rules out ground clutter.
2. **Sounding data confirms the wind anomaly.** Winds were from the SSE/S at 165 to 200 degrees at the target altitude, which means any westward object motion is definitively anomalous.
3. **RCS estimation shows the return is too large** for a single mylar balloon. The 196 pixel spread points to a debris cloud or chaff.
4. **Dual pol signature (ZDR 18.9, CC 0.50)** is textbook metallic chaff or debris, not a balloon.
5. **KBIF surface wind was calm or light southerly,** consistent with nighttime drainage flow. No wind anomalies at the surface.


## 1. KHDX Second Radar Analysis

### Method
We downloaded and analyzed 34 KHDX (Holloman AFB, NM) NEXRAD Level II scans from 03:03 to 06:57 UTC, covering the Fort Bliss sector (azimuth 185 to 199 degrees, range 125 to 148 km from KHDX).

### Results

KHDX detected returns in the Fort Bliss sector across all 34 scans, totaling 783 pixels above 5 dBZ. The important part is that these returns intensified around both KEPZ detected spike times:

| Time Window | Max Reflectivity | Pixels > 5 dBZ | Pixels > 15 dBZ |
|-------------|-----------------|-----------------|-----------------|
| 04:21 UTC (near Spike 1) | 14.5 dBZ | 17 | 0 |
| 04:29 UTC | 19.5 dBZ | 23 | 3 |
| 04:36 UTC | 23.5 dBZ | 21 | 10 |
| 05:11 UTC (near Spike 2) | 23.0 dBZ | 26 | 5 |
| 05:18 UTC | 25.0 dBZ | 43 | 14 |

**Key KHDX returns near Spike 1 (04:24 UTC):**
At 04:29 UTC, 19.5 dBZ appeared at 31.8125N, 106.4926W (azimuth 193.7 degrees, 145 km). By 04:36 UTC, the return had strengthened to 23.5 dBZ at 31.8147N, 106.4922W (azimuth 193.8 degrees, 144.6 km).

**Key KHDX returns near Spike 2 (05:13 UTC):**
At 05:11 UTC, 23.0 dBZ appeared at 31.8312N, 106.5150W (azimuth 194.7 degrees, 143.4 km). By 05:18 UTC, 25.0 dBZ at 31.8125N, 106.4926W (azimuth 193.7 degrees, 144.9 km).

### Interpretation

KHDX confirms a real airborne target. This is not KEPZ ground clutter or noise. The target appears at a range of roughly 137 to 145 km from KHDX, consistent with the Fort Bliss area, and peak reflectivity values from KHDX are comparable to what KEPZ measured (21.5 dBZ vs 23.5 to 25.0 dBZ).

That said, KHDX also shows persistent background returns in the 5 to 15 dBZ range throughout the entire observation window, some of which may be ground clutter at these extreme ranges. What matters most is the intensification around the KEPZ spike times. That pattern suggests a real transient target appeared and was picked up by both radars.

### Position Cross Reference

When we compare KHDX derived positions with KEPZ positions, there is a notable discrepancy:

KEPZ Spike 1 was at 31.8771N, 106.3423W. KHDX returns near Spike 1 clustered around 31.8125 to 31.8191N, 106.4907 to 106.4948W. That puts the KHDX positions about 14 km southwest of the KEPZ position.

This gap is significant and could result from beam height differences (KHDX at roughly 2500m MSL vs KEPZ at roughly 1600m MSL), different parts of a debris cloud being sampled at different altitudes, positional uncertainty at KHDX's extreme range, or the KHDX returns coming from a different persistent feature at that azimuth.

**A note on persistence:** The fact that KHDX returns at azimuth 193 to 195 degrees, range 137 to 145 km appear across all scans suggests some of these may be ground clutter. The time variation in intensity is the key discriminator. Ground clutter stays constant, while the intensification around the spike times points to a real transient target superimposed on a background.


## 2. Upper Air Sounding Analysis

### Data Source
Station 72364 (EPZ / Santa Teresa, NM) soundings from the University of Wyoming archive via Siphon.

### 12Z Sounding (Feb 11, 2026 at 05:00 AM MST, closest to event)

| Height AGL | Pressure | Wind Dir | Wind Speed | Temperature |
|-----------|----------|----------|------------|-------------|
| 0 m (sfc) | 879 hPa | 125 degrees (SE) | 7 kt (3.6 m/s) | 11.6C |
| 247 m | 854 hPa | 165 degrees (SSE) | 11 kt (5.7 m/s) | 14.3C |
| 287 m | 850 hPa | 180 degrees (S) | 10 kt (5.1 m/s) | 14.6C |
| 367 m | 842 hPa | 187 degrees (S) | 11 kt (5.7 m/s) | 14.4C |
| 397 m | 839 hPa | 190 degrees (S) | 11 kt (5.7 m/s) | 14.3C |
| 519 m | 827 hPa | 198 degrees (SSW) | 10 kt (5.1 m/s) | 13.8C |
| 549 m | 824 hPa | 200 degrees (SSW) | 10 kt (5.1 m/s) | 13.6C |
| 731 m | 806 hPa | 185 degrees (S) | 4 kt (2.1 m/s) | 12.0C |
| 845 m | 795 hPa | 190 degrees (S) | 7 kt (3.6 m/s) | 11.1C |
| 971 m | 783 hPa | 180 degrees (S) | 9 kt (4.6 m/s) | 10.0C |

### 00Z Sounding (Feb 11, 2026 at 05:00 PM MST Feb 10, pre event baseline)

Surface wind was from the NE at 5 kt. Low level winds between 300 and 700m AGL were from WNW/NW at 3 to 5 kt, a completely different regime from the 12Z sounding. This tells us a wind shift occurred between evening and early morning.

### Comparison to VAD Derived Winds

| Parameter | VAD (KEPZ) | 12Z Sounding | Agreement? |
|-----------|-----------|--------------|------------|
| Direction | 153 to 215 degrees | 165 to 200 degrees | Yes |
| Speed | 4.7 to 5.1 m/s | 4 to 6 m/s (8 to 11 kt) | Yes |
| Altitude | roughly 300 to 700m AGL | 247 to 549m AGL | Yes |

The sounding confirms the VAD derived winds. At 300 to 550m AGL, right where the radar spikes were detected, winds were consistently from the south at 5 to 6 m/s. That means an object at these altitudes should drift northward, since wind from the south pushes things north. But the observed motion between Spike 1 and Spike 2 was due west, roughly 7 km. The westward motion is definitively anomalous and cannot be explained by passive drift.

### Strong Temperature Inversion

The 12Z sounding shows a notable temperature inversion. The surface was at 11.6C while 287m AGL reached 14.6C, a 3 degree increase above the surface. This is a classic nocturnal radiation inversion trapping cool air near the ground with warmer air aloft. This matters for several reasons: inversions enhance radar ducting, which could affect KEPZ beam propagation; the inversion cap would trap any debris or aerosols released at low altitude; and light winds in the inversion layer mean slow dispersal of any debris cloud.


## 3. KBIF (Biggs Army Airfield) METAR Analysis

### Event Window Observations

| Time (UTC) | Wind | Visibility | Clouds | Temperature |
|-----------|------|-----------|--------|-------------|
| 03:55 | 200 degrees at 5 kt | 10 SM | FEW190 | 16.2C |
| 05:55 | Calm | 10 SM | BKN095 | 14.9C |
| 06:55 | Calm | 10 SM | BKN100, BKN250 | 14.7C |

### Key Findings

No SPECIs were issued during the event, meaning no unusual weather was observed. Surface wind was calm or light SSW (200 degrees), consistent with the sounding profile. Visibility held at 10 statute miles throughout with clear conditions, no fog, haze, or dust. Broken cloud decks appeared at 9,500 to 10,000 ft AGL by 05:55 UTC. Surface temperature was 15 to 16C, consistent with the inversion trapping described above. The dollar sign maintenance indicator is present on all KBIF METARs, which tells us the ASOS station may have a sensor issue (common for military stations with limited maintenance).

### Significance

The calm surface winds and clear visibility confirm two things: there were no surface weather events that could explain radar anomalies, and a 20kW laser beam would have had excellent propagation conditions. The SSW surface wind at 03:55 UTC is also consistent with the sounding's low level southerly wind.


## 4. RCS Estimation

### Method
We computed radar cross section from KEPZ reflectivity values using the standard radar equation for point and distributed targets.

### Results

For Spike 1 (21.5 dBZ, 196 pixels, 33.6 km range):

| Scenario | RCS | dBsm | Consistent With |
|----------|-----|------|----------------|
| Single point target in one gate | 0.020 sq m | negative 17.0 | Mylar balloon, small drone |
| Distributed across 196 pixels | 3.92 sq m total | positive 5.9 | Large drone, light aircraft, chaff |
| Conservative (total Z over all pixels) | 0.020 sq m total | negative 17.0 | Mylar balloon, small drone |

### Reference RCS Values

For context, a mylar party balloon (30cm) typically has an RCS of 0.001 to 0.01 sq m. A small consumer drone falls in the 0.01 to 0.1 sq m range. A large military drone like the MQ 9 sits between 1 and 5 sq m. A dispersed chaff cloud ranges from 1 to 100 sq m, and a light aircraft from 1 to 10 sq m.

### Interpretation

The single gate RCS of 0.020 sq m sits at the upper end of what a mylar balloon would produce. But the 196 pixel spatial extent is completely inconsistent with a single balloon. A balloon is a point target. It should illuminate 1 to 4 pixels at most, not 196.

The distributed total RCS of roughly 4 sq m is consistent with a debris or chaff cloud, multiple simultaneous targets, or a destroyed target with metallic fragments dispersing.


## 5. Dual Polarization Signature Analysis

### Measured Values (Spike 1)

ZDR was 18.9 dB, which is extremely high. CC (correlation coefficient) was 0.50, which is very low.

### Diagnostic Table

| Target Type | Typical ZDR | Typical CC | Match? |
|------------|-------------|-----------|--------|
| Rain | 0 to 4 dB | above 0.95 | No |
| Hail | negative 2 to positive 2 dB | 0.80 to 0.95 | No |
| Biological (birds/insects) | 2 to 10 dB | 0.3 to 0.8 | Partial |
| Chaff/metallic debris | above 10 dB | 0.2 to 0.6 | Yes |
| Ground clutter | Variable | 0.5 to 0.9 | Partial |
| Single metallic target | Variable | above 0.90 | No |

### Diagnosis

High ZDR combined with low CC equals a metallic debris or chaff cloud. This is the textbook dual pol signature of military chaff dipoles, flat metallic fragments (such as destroyed mylar or aluminum debris), and a cloud of tumbling, randomly oriented metallic scatterers.

This is not consistent with a single intact balloon (which would show CC above 0.9), a single drone (also CC above 0.9), or meteorological targets (CC above 0.95).

The dual pol signature strongly suggests the target was destroyed before or during the first radar detection, producing a cloud of metallic fragments.


## 6. Terrain and Beam Blockage Analysis

### Franklin Mountains Blockage

KEPZ (31.873N, 106.698W) has to look through the Franklin Mountains to see the Fort Bliss area. Here is what we found:

| Franklin Peak | Azimuth from KEPZ | Beam Lower Edge | Clearance |
|--------------|-------------------|-----------------|-----------|
| North Franklin (2192m) | 66.8 degrees | 1284m | negative 908m BLOCKED |
| Mundy Peak (2100m) | 78.1 degrees | 1283m | negative 817m BLOCKED |
| Ranger Peak (1690m) | 94.0 degrees | 1284m | negative 406m BLOCKED |
| Trans Mtn Gap (1580m) | 75.7 degrees | 1287m | negative 293m BLOCKED |
| Smugglers Gap (1500m) | 100.3 degrees | 1285m | negative 215m BLOCKED |

Every azimuth from 67 to 100 degrees is fully blocked on the 0.5 degree tilt. The Franklin Mountains create a roughly 12 km wide radar shadow east of the ridge for targets below roughly 500m AGL.

### Implications

The beam blockage zone sits 10 to 18 km from KEPZ, meaning targets behind the mountains below 500m AGL are invisible. A target approaching from the east, south, or Mexico would be undetectable until it either rose above roughly 500m AGL or moved to an azimuth outside the blockage zone. This fully explains the "appeared from nowhere" observation. It is a terrain artifact, not evidence of exotic origin.

Position uncertainty at 33.6 km range is roughly 300m, which is not enough to explain the 3 to 5 km displacement of Spike 1 into the civilian area.

### KHDX Viewing Geometry

KHDX views the target at azimuth roughly 189 to 195 degrees, range roughly 135 to 145 km. At 135 km, the KHDX 0.5 degree beam lower edge sits at roughly 2425m MSL, or about 1225m AGL. That means KHDX cannot see targets below roughly 1200m AGL at this range. Yet KHDX does detect returns, which means the target or debris extended above 1200m AGL.


## 7. KBIF Surface vs. Altitude Wind Comparison

| Source | Height | Direction | Speed |
|--------|--------|-----------|-------|
| KBIF METAR 03:55Z | Surface | 200 degrees (SSW) | 5 kt (2.6 m/s) |
| KBIF METAR 05:55Z | Surface | Calm | 0 kt |
| Sounding 12Z 247m AGL | 247m | 165 degrees (SSE) | 11 kt (5.7 m/s) |
| Sounding 12Z 287m AGL | 287m | 180 degrees (S) | 10 kt (5.1 m/s) |
| Sounding 12Z 397m AGL | 397m | 190 degrees (S) | 11 kt (5.7 m/s) |
| VAD (KEPZ) est. | 300 to 700m | 153 to 215 degrees | 9 to 10 kt (4.7 to 5.1 m/s) |

Surface winds were calm to light SSW. At 250 to 550m AGL, winds were southerly at 10 to 11 kt. The VAD, sounding, and surface observations all agree: the entire low level wind profile was from the south or SSE.

An object at any altitude below 1000m should drift north or NNW. Westward motion is anomalous at all altitudes.


## 8. News Coverage and OSINT Summary

### Official Narrative

Multiple confirmed reports indicate the following. CBP personnel operated an AeroVironment LOCUST 20kW laser system from Fort Bliss. The laser was used to engage objects identified as potential cartel drones. At least one target was later identified as a "party balloon" (mylar). An FBI internal document stated the targets were "floating holiday decorations." No prior FAA coordination was conducted. FAA Administrator Bryan Bedford closed airspace unilaterally. Rep. Escobar called it "incompetence at the highest levels."

### LOCUST System Specs (from open sources)

The system is a 20 kW class laser (demonstrated up to 26 kW) with an effective range described as "short ranges, typically measured in handfuls of miles." It is mounted on an M1301 Infantry Squad Vehicle (ISV) or JLTV and equipped with EO/IR cameras plus integrated radar on the vehicle. The first two mobile systems were delivered to the Army in August 2025. The manufacturer is AeroVironment (via BlueHalo acquisition).

### Eyewitness Accounts

No direct eyewitness reports of a laser beam were found in web searches. A 2020 KFOX article about "strange sightings" in El Paso skies turned out to be unrelated (Fort Bliss artillery flares). The APO all sky camera (Sunspot, NM) has video for Feb 10, 2026, and the file is accessible, but at 90 miles from El Paso, a 20kW laser beam would be extremely difficult to detect at that distance.


## 9. Revised Anomaly Assessment

| Finding | Previous Assessment | New Assessment (with additional data) |
|---------|-------------------|--------------------------------------|
| No inbound radar track | Terrain blockage | Confirmed: Franklin Mtns fully block KEPZ azimuth 67 to 100 degrees |
| Wind direction mismatch | Possible two separate targets | Confirmed anomalous: Sounding validates SSE/S winds; westward motion not wind driven |
| Over civilian area | Positional error | Unlikely error: 300m uncertainty cannot explain 3 to 5 km displacement |
| Two spikes 49 min apart | Two test engagements | Supported: KHDX shows activity at both times |
| 196 pixel return cloud | Debris from destroyed target | Strongly supported: RCS, ZDR, and CC all consistent with metallic debris cloud |
| Metallic dual pol signature | Mylar balloon | Debris cloud: High ZDR plus low CC equals textbook chaff/debris, not single balloon |
| Object descended between spikes | Different altitudes | Inconclusive: KHDX positions shifted roughly 14 km from KEPZ positions |
| KHDX independent detection | N/A (new) | Significant: Second radar confirms real airborne target at both spike times |


## 10. Synthesis and Competing Hypotheses

### Hypothesis A: Mylar Balloon Shot by LOCUST

This is the official narrative from CBP, DoD, and the FBI. In its favor, the RCS of a single gate (roughly 0.02 sq m) is consistent with a mylar balloon, and the SSE wind direction could have brought a balloon up from Mexico.

The problems are significant, though. A single balloon should produce a 1 to 4 pixel return, not 196. The dual pol signature (ZDR 18.9, CC 0.50) does not match a single target. Wind cannot explain the westward motion between spikes. A second engagement 49 minutes later with the same characteristics is suspicious. And the target descended between spikes, which is the opposite of what balloons normally do unless they are deflating.

### Hypothesis B: LOCUST Engagement of Drone Producing Debris Cloud

All the radar signatures (spatial extent, RCS, ZDR, CC) are consistent with metallic debris. "Cartel drone" was the initial White House claim. Two engagements could mean two targets. A debris cloud would expand and descend exactly as observed. And the location over a civilian area makes sense if a drone approached from Mexico.

On the other hand, the FBI later called the targets "floating holiday decorations." No ADS B or FR24 track matches a drone. And it raises the question of whether CBP would fire a laser over a residential neighborhood.

### Hypothesis C: LOCUST Test/Calibration Exercise Against Known Target

Two spikes 49 minutes apart suggest a planned sequence. The debris signature is consistent with deliberate chaff release or target destruction. Fort Bliss is an active laser testing site (LOCUST ISV training since July 2025). This would also explain why the event occurred over a civilian area, as it could reflect test geometry rather than a real threat response.

The problems here: why was there no FAA coordination if this was planned? Why claim "cartel drone" if it was a test? And the KHDX background returns do suggest some persistent features in the area that complicate the picture.

### Most Likely Reconstruction

The evidence best supports a scenario where CBP operators used the LOCUST laser to engage one or more low altitude targets (likely balloons, possibly mixed with actual drone contacts) without FAA coordination. The radar returns show debris clouds consistent with destroyed metallic targets, not intact balloons. The westward motion between spikes may reflect two separate targets at different locations rather than one target moving west. The position over the civilian area is real (not a position error) and suggests the engagement occurred over northeastern El Paso, not over the military reservation.


## 11. Recommended Next Steps

### Immediate

Download KEPZ raw data for triangulation with KHDX positions. If available, process at higher resolution. Check if the persistent KHDX returns at azimuth 193 to 195 degrees appear on other dates, which would indicate ground clutter baseline vs. transient activity. Request FAA ASR 11 radar data from El Paso International, since airport radar has different characteristics and may resolve the target.

### Medium Term

Submit FOIA requests (drafts prepared in the foia_requests folder). Compare KHDX data from Feb 10 and Feb 12 to establish baseline clutter levels. Attempt ADS B historical API access for military transponder data.

### Long Term

Obtain SkySentinel NMSU 1 data through PI collaboration for the 57 hour data gap. Monitor Congressional inquiries for additional information releases.


## Data Files

### Downloaded Data

The khdx_baseline_data folder contains 34 KHDX NEXRAD Level II files (03:03 to 06:57 UTC). Upper air soundings are in sounding_data/EPZ_72364_12Z_20260211.csv and sounding_data/EPZ_72364_00Z_20260211.csv. METAR observations are in metar_data/kbif_metar_feb10-12.csv.

### Analysis Scripts

scripts/khdx_analysis.py handles the KHDX NEXRAD sector analysis. scripts/rcs_estimation.py handles reflectivity to RCS computation. scripts/terrain_analysis.py handles Franklin Mountains beam blockage.

### Analysis Outputs

Results are in analysis_outputs/khdx_analysis_results.json (KHDX scan by scan results), analysis_outputs/rcs_estimation_results.json (RCS estimates), analysis_outputs/terrain_analysis_results.json (terrain and geometry calculations), and sounding_data/wind_analysis_results.json (wind comparison data).

### FOIA Request Drafts

Drafts are available in the foia_requests folder for CBP, Fort Bliss, FAA, and NWS.


## Sources

CNN: Surprise US military plans to use counter drone laser triggered airspace closure (https://www.cnn.com/2026/02/11/us/faa-el-paso-texas-flight-restrictions-hnk)

DefenseScoop: CBP personnel used military laser to shoot object near El Paso (https://defensescoop.com/2026/02/11/el-paso-drone-incursion-military-laser-cbp-faa/)

The War Zone: LOCUST laser that prompted closing El Paso's airspace (https://www.twz.com/news-features/this-is-the-locust-laser-that-reportedly-prompted-closing-el-pasos-airspace)

Scientific American: Why an Army antidrone laser grounded flights at El Paso (https://www.scientificamerican.com/article/why-an-army-antidrone-laser-grounded-flights-at-el-paso-international/)

NBC News: CBP shot down party balloons with anti drone tech (https://www.nbcnews.com/politics/national-security/cbp-shot-party-balloons-anti-drone-tech-faa-closed-el-paso-airspace-so-rcna258731)

ABC News: FAA shutdown triggered by dispute over Pentagon laser weapon (https://abcnews.com/Politics/faa-shutdown-el-paso-airspace-triggered-dispute-pentagon/story?id=130068283)

CBS News: Airspace closure followed spat over drone related tests (https://www.cbsnews.com/news/airspace-closure-followed-spat-over-drone-related-tests-and-party-balloon-shoot-down-sources-say/)

Military Times: Pentagon let CBP use anti drone laser before FAA closed airspace (https://www.militarytimes.com/news/pentagon-congress/2026/02/12/pentagon-let-cbp-use-anti-drone-laser-before-faa-closed-el-paso-airspace-report/)

NEXRAD Level II data: NOAA/UCAR THREDDS (KEPZ, KHDX)

Upper air sounding: University of Wyoming / Siphon (Station 72364)

METAR data: Iowa State Mesonet ASOS archive (KBIF)
