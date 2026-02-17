# I analyzed 6,299 satellite files, dual-radar data, and seismic records from the El Paso airspace closure. Seven independent public data sources all point to the same conclusion.

On the night of Feb 10-11, 2026, CBP personnel at Fort Bliss reportedly engaged airborne objects using a directed energy system without coordinating with the FAA. The FAA closed El Paso's airspace. Over the following week, the official characterization of the targets shifted from "cartel drone incursion" to "party balloon" to "four Valentine's Day balloons" to "floating holiday decorations."

I downloaded every piece of publicly available sensor data I could find and spent the last week doing independent analysis. Everything here is reproducible from public sources, and I'll link the full methodology and code at the bottom.

**A second radar 135 km away shows nearly double its normal returns.**

Most analysis of this incident focuses on KEPZ, the primary NWS radar near El Paso. KEPZ has a problem: the Franklin Mountains sit between it and Fort Bliss, creating ground clutter that complicates interpretation. So I also pulled data from KHDX at Holloman AFB, 135 km away, which has a clean line of sight.

I downloaded KHDX scans from the event night (Feb 11) and four clear-sky baseline nights for comparison, filtering returns to the azimuth and range window toward Fort Bliss (not the entire scan volume). The difference:

KHDX baseline: 5.6 returns above 10 dBZ per scan on average. KHDX event night: 10.0 returns above 10 dBZ per scan. That is a 78% increase with a p-value of 0.003. Overall returns above 5 dBZ increased 45%. Surface observations confirmed clear skies and calm winds on all nights compared.

Twelve of 34 KHDX scans show significant excess above baseline, concentrated in the 04:30 to 06:30 UTC window, with peak returns reaching 37.0 dBZ. The temporal pattern escalates over roughly 90 minutes rather than fading, which is not consistent with a single target fragmenting and dispersing.

All 33 KEPZ scan-pairs have spatially coincident KHDX returns within 20 km. Mean offset between the two radars is 5.7 km, which is expected agreement at KHDX's extreme range. The triangulated position lands in the Fort Bliss area, east of the Franklin Mountains.

**Two engagements, not one, and the debris signature is consistent with metallic chaff.**

The KEPZ radar shows two clear reflectivity spikes 49 minutes apart (04:24 and 05:13 UTC, or 9:24 and 10:13 PM local). Both have nearly identical characteristics. The first spike produced 196 radar range-gate bins above threshold across a 7x13 km area in the full-volume scan. A single mylar balloon is a point target that should illuminate 1 to 4 bins at most.

The dual-polarization data is where it gets specific. At the peak spike: ZDR of 18.9 dB and a correlation coefficient (CC) of 0.50. Across the event bins, ZDR values ranged from approximately 8 to 19 dB with CC values between 0.45 and 0.65. For reference, a single intact metallic target produces CC above 0.90. Rain produces CC above 0.95. The combination of high ZDR with low CC is consistent with the known radar signature of metallic debris clouds or military chaff: a collection of tumbling, randomly oriented flat metallic scatterers.

**The target moved in a direction the wind cannot explain.**

Between the two spikes, the returns shifted about 7 km due west. The upper air sounding from Santa Teresa (Station 72364) shows winds at the target altitude (300 to 550m AGL) consistently from the south at 10 to 11 knots. Anything drifting passively should have moved north. The sounding, the radar-derived wind profile (VAD), and surface observations from Biggs Army Airfield all agree. Westward motion is anomalous at every altitude below 1000m. This likely means the two spikes represent two separate targets at different locations rather than one object drifting.

**The first spike was over a residential neighborhood.**

Spike 1 plots to coordinates over NE El Paso near Hondo Pass, not over the military reservation. Position uncertainty at 33.6 km range is about 300 meters, which is nowhere near enough to explain the 3 to 5 km displacement into a civilian area.

**Four additional sensor systems all tell a consistent story.**

**Seismic (IRIS/EarthScope station EP.KIDD, 11 km from Fort Bliss):** Zero anomalies on the event night. At 11 km, a broadband seismometer should detect any explosive or kinetic energy release above roughly a few kilograms of TNT equivalent. The absence of any seismic signal is consistent with a non-kinetic engagement: a continuous-wave laser that silently heats targets until structural failure rather than detonating them.

The seismic data also resolved an unrelated question. Reports of "booming sounds" on Feb 8-9 that some people connected to Fort Bliss activity correlate exactly with documented M2.5 and M2.1 earthquakes near Whites City, NM, roughly 190 km east. Delaware Basin induced seismicity, not weapons testing.

**GOES-19 Geostationary Lightning Mapper (6,299 files across 5 nights, Feb 8-12):** The GLM is a near-infrared optical transient detector operating at 777.4 nm with a detection threshold of approximately 1-2 femtojoules. I processed every file covering the Fort Bliss area across a five-day window. Result: zero detections at Fort Bliss on any night. The only detections in the entire dataset were 103 faint events near Columbus, NM (130 km away, consistent with distant lightning), completely unrelated to the investigation area.

The null result constrains the engagement. The GLM rules out any optical transient at Fort Bliss above its detection threshold across the entire five-day period. The LOCUST system likely operates at 1064 nm (typical for this class of high-energy laser), which falls outside GLM's 777.4 nm passband. A continuous-wave laser slowly heating a target would not produce the kind of sudden bright flash GLM detects. This is consistent with a CW laser engagement and inconsistent with explosive detonations, missile intercepts, or kinetic kill events, all of which would produce broadband optical emissions detectable by GLM.

**NASA FIRMS thermal hotspot detection:** Zero hotspots at Fort Bliss on the event night. No fires, no large-scale combustion.

**Surface weather (KBIF METAR):** Clear skies, 10 statute miles visibility, calm to light winds all night. No weather phenomena to explain any of the radar anomalies.

**What the combined sensor picture looks like.**

Seven independent data sources, operated by different agencies, archived on different public servers, all examined independently:

Radar (two stations) sees metallic debris from two engagements, with returns escalating over 90 minutes. Seismic (11 km away) confirms no explosions. Optical satellite (6,299 files) confirms no bright transients. Thermal satellite confirms no fires. Weather confirms clear skies. Upper air sounding confirms the wind anomaly. Two radars independently triangulate to Fort Bliss.

The combined picture is consistent with a quiet, non-explosive, non-optically-bright energy delivery that produced metallic debris detected by two independent radars. That description matches the expected signature of a continuous-wave laser engagement: silent heating until structural failure, target fragments, metallic debris disperses.

**Physical observables I cannot explain:**

Zero residents have reported finding metallic debris despite calm winds, low engagement altitude (roughly 250m), and a location over a residential neighborhood. If metallic targets were fragmented at that altitude in those conditions, debris should have reached the ground within a relatively tight footprint. Either it was recovered quickly, the fragments were too small to notice, or the debris field is not where the radar indicates.

The radar returns escalated for approximately 90 minutes after the initial spikes rather than fading. A single destroyed target should produce a debris cloud that disperses and weakens over time. Escalation suggests either additional engagements occurred (possibly hidden from KEPZ by the Franklin Mountains but visible to KHDX), or the activity was more sustained than the official account describes.

The westward displacement between the two KEPZ spikes cannot be explained by wind at any altitude in the sounding profile. Two separate targets at different locations is the most straightforward explanation, but it contradicts the single-engagement narrative.

**What still needs to be investigated:**

Police scanner audio from the event night exists on Broadcastify (El Paso feeds 42309, 29918, 28117, 32168) with 365-day retention. That audio could contain real-time dispatch traffic about civilian reports. Public records requests for 911 call logs can go through El Paso's GovQA portal for police, fire, and sheriff separately. The New Mexico side is handled by MVRDA under IPRA. A construction camera at El Paso International Airport may have captured footage (request pending). Historical ADS-B data could show military aircraft patterns during the engagement window.

**Everything is public and reproducible.** NEXRAD Level II files are on NOAA's AWS servers. Upper air soundings are on the University of Wyoming archive. METAR observations are on the Iowa State Mesonet. GOES GLM files are on NOAA's S3 bucket. Seismic waveforms are on the IRIS DMC. FIRMS data is on NASA Earthdata. Full technical report with methodology, data tables, analysis scripts, and source references is linked below.

[LINK TO REPORT/REPO]

**If you live in NE El Paso near Fort Bliss, particularly the Hondo Pass or Gateway area, and found shiny metallic fragments in your yard on the morning of Feb 11, or if you have security camera footage from roughly 9:00 to 10:30 PM on Feb 10, please reach out.**
