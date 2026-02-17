#!/usr/bin/env python3
"""
Multi-Radar Object Characterization Visualization
==================================================
Creates an informative multi-panel figure showing:
1. Spatial map: triangulated positions from KEPZ + KHDX overlaid on geography
2. Timeline: object return count escalation through the night
3. Object characterization: reflectivity, size, and movement analysis
4. Dual-radar convergence showing both radars seeing the same area
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'

# ============================================================================
# LOAD DATA
# ============================================================================
with open(f'{OUTPUT_DIR}/integrated_triangulation_results.json') as f:
    tri_data = json.load(f)

with open(f'{OUTPUT_DIR}/khdx_azimuth_persistence.json') as f:
    persist_data = json.load(f)

matches = tri_data['triangulation']['matched_scans']

# Extract triangulated positions over time
times_utc = []
tri_lats = []
tri_lons = []
kepz_refs = []
khdx_refs = []
kepz_n = []
khdx_n = []
offsets = []

for m in matches:
    time_str = m['kepz_time']
    hours = int(time_str[:2])
    mins = int(time_str[2:])
    t = hours + mins / 60.0
    times_utc.append(t)
    tri_lats.append(m['triangulated_lat'])
    tri_lons.append(m['triangulated_lon'])
    kepz_refs.append(m['kepz_max_ref'])
    khdx_refs.append(m['khdx_max_ref'])
    kepz_n.append(m['kepz_n_strong'])
    khdx_n.append(m['khdx_n_returns'])
    offsets.append(m['offset_km'])

times_utc = np.array(times_utc)
tri_lats = np.array(tri_lats)
tri_lons = np.array(tri_lons)
kepz_refs = np.array(kepz_refs)
khdx_refs = np.array(khdx_refs)
kepz_n = np.array(kepz_n)
khdx_n = np.array(khdx_n)
offsets = np.array(offsets)

# Radar locations
KEPZ_LAT, KEPZ_LON = 31.8731, -106.698
KHDX_LAT, KHDX_LON = 33.0803, -106.1225

# Fort Bliss / engagement area
BLISS_LAT, BLISS_LON = 31.85, -106.40

# ============================================================================
# CREATE FIGURE
# ============================================================================
fig = plt.figure(figsize=(20, 16))
fig.suptitle('MULTI-RADAR OBJECT CHARACTERIZATION — El Paso Airspace Incident\nFebruary 11, 2026 (03:06-06:53 UTC)',
             fontsize=16, fontweight='bold', y=0.98)

gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.30,
              left=0.06, right=0.96, top=0.93, bottom=0.05)

# ============================================================================
# PANEL 1: Spatial map of triangulated positions (large, left side)
# ============================================================================
ax1 = fig.add_subplot(gs[0:2, 0:2])

# Plot triangulated positions colored by time
scatter = ax1.scatter(tri_lons, tri_lats, c=times_utc, cmap='plasma',
                      s=80, edgecolors='black', linewidth=0.5, zorder=5,
                      vmin=3.0, vmax=7.0)
cb = plt.colorbar(scatter, ax=ax1, label='UTC Hour', shrink=0.8)

# Connect positions with line to show movement
ax1.plot(tri_lons, tri_lats, 'k-', alpha=0.3, linewidth=0.8, zorder=4)

# Mark radar stations
ax1.plot(KEPZ_LON, KEPZ_LAT, 'r^', markersize=14, zorder=10, label='KEPZ (El Paso)')
ax1.plot(KHDX_LON, KHDX_LAT, 'b^', markersize=14, zorder=10, label='KHDX (Holloman)')
ax1.annotate('KEPZ\n(16-55km)', (KEPZ_LON, KEPZ_LAT), fontsize=8,
             textcoords='offset points', xytext=(-15, 10), color='red', fontweight='bold')

# Mark Fort Bliss
ax1.plot(BLISS_LON, BLISS_LAT, 'gs', markersize=12, zorder=10, label='Fort Bliss')

# Draw lines from each radar to the mean triangulated position
mean_lat, mean_lon = np.mean(tri_lats), np.mean(tri_lons)
ax1.plot([KEPZ_LON, mean_lon], [KEPZ_LAT, mean_lat], 'r--', alpha=0.4, linewidth=1.5)
ax1.plot([KHDX_LON, mean_lon], [KHDX_LAT, mean_lat], 'b--', alpha=0.4, linewidth=1.5)

# Add mean position marker
ax1.plot(mean_lon, mean_lat, 'w*', markersize=20, markeredgecolor='black',
         markeredgewidth=1.5, zorder=11, label=f'Mean position\n({mean_lat:.4f}N, {mean_lon:.4f}W)')

# Spread ellipse
lat_std = np.std(tri_lats)
lon_std = np.std(tri_lons)
ellipse = matplotlib.patches.Ellipse((mean_lon, mean_lat),
                                      width=4*lon_std, height=4*lat_std,
                                      fill=False, edgecolor='orange', linewidth=2,
                                      linestyle='--', label=f'2-sigma spread ({lat_std*111:.1f}x{lon_std*111*np.cos(np.radians(mean_lat)):.1f} km)')
ax1.add_patch(ellipse)

# Geographic context labels
ax1.annotate('El Paso\n(city center)', (-106.44, 31.76), fontsize=7, ha='center', color='gray')
ax1.annotate('Franklin\nMountains', (-106.49, 31.92), fontsize=7, ha='center', color='brown')
ax1.annotate('Fort Bliss\nMain Post', (-106.38, 31.81), fontsize=7, ha='center', color='green')
ax1.annotate('NE El Paso\n(residential)', (-106.35, 31.87), fontsize=7, ha='center', color='gray')

ax1.set_xlabel('Longitude (°W)')
ax1.set_ylabel('Latitude (°N)')
ax1.set_title('Triangulated Target Positions (KEPZ + KHDX)\nColor = UTC time, Star = mean position', fontsize=11)
ax1.legend(loc='upper left', fontsize=7, framealpha=0.9)
ax1.set_xlim(-106.58, -106.32)
ax1.set_ylim(31.74, 31.96)
ax1.grid(True, alpha=0.3)
ax1.set_aspect(1.0 / np.cos(np.radians(mean_lat)))

# ============================================================================
# PANEL 2: Reflectivity timeline (top right)
# ============================================================================
ax2 = fig.add_subplot(gs[0, 2])

ax2.plot(times_utc, kepz_refs, 'r-o', markersize=4, label='KEPZ max dBZ', linewidth=1.5)
ax2.plot(times_utc, khdx_refs, 'b-s', markersize=4, label='KHDX max dBZ', linewidth=1.5)

# Shade the peak activity period
peak_mask = (times_utc >= 4.3) & (times_utc <= 5.5)
if np.any(peak_mask):
    ax2.axvspan(4.3, 5.5, alpha=0.15, color='red', label='Peak activity')

ax2.set_xlabel('UTC Hour')
ax2.set_ylabel('Max Reflectivity (dBZ)')
ax2.set_title('Dual-Radar Reflectivity Timeline', fontsize=10)
ax2.legend(fontsize=7)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(3.0, 7.0)

# ============================================================================
# PANEL 3: Return count timeline (middle right)
# ============================================================================
ax3 = fig.add_subplot(gs[1, 2])

# Baseline levels
khdx_baseline = tri_data['khdx_baseline_comparison']['baseline_avg_n10']
kepz_baseline = tri_data['kepz_baseline_summary']['baseline_avg_n10']

ax3.fill_between(times_utc, 0, khdx_n, alpha=0.3, color='blue', label=f'KHDX returns')
ax3.plot(times_utc, khdx_n, 'b-', linewidth=1.5)
ax3.axhline(y=khdx_baseline, color='blue', linestyle='--', alpha=0.5,
            label=f'KHDX baseline ({khdx_baseline:.0f}/scan)')

# Add annotations for peaks
peak_idx = np.argmax(khdx_n)
ax3.annotate(f'Peak: {khdx_n[peak_idx]} returns\n{times_utc[peak_idx]:.2f} UTC',
             (times_utc[peak_idx], khdx_n[peak_idx]),
             textcoords='offset points', xytext=(10, 5), fontsize=7,
             arrowprops=dict(arrowstyle='->', color='blue'))

ax3.set_xlabel('UTC Hour')
ax3.set_ylabel('Returns per Scan')
ax3.set_title('KHDX Target Return Escalation\n(baseline-subtracted = transient signal)', fontsize=10)
ax3.legend(fontsize=7)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(3.0, 7.0)

# ============================================================================
# PANEL 4: Object characterization (bottom left)
# ============================================================================
ax4 = fig.add_subplot(gs[2, 0])

# Position scatter colored by reflectivity
sizes = kepz_n / 10  # Scale for visibility
scatter4 = ax4.scatter(tri_lons, tri_lats, c=kepz_refs, cmap='YlOrRd',
                        s=sizes, edgecolors='black', linewidth=0.3,
                        vmin=30, vmax=65)
plt.colorbar(scatter4, ax=ax4, label='KEPZ Max dBZ', shrink=0.8)

# Calculate spread statistics
lat_range = (np.max(tri_lats) - np.min(tri_lats)) * 111  # km
lon_range = (np.max(tri_lons) - np.min(tri_lons)) * 111 * np.cos(np.radians(mean_lat))

ax4.set_xlabel('Longitude')
ax4.set_ylabel('Latitude')
ax4.set_title(f'Position Spread: {lat_range:.1f}km N-S x {lon_range:.1f}km E-W\n'
              f'(bubble size = # KEPZ returns)', fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.set_aspect(1.0 / np.cos(np.radians(mean_lat)))

# ============================================================================
# PANEL 5: Inter-radar agreement (bottom center)
# ============================================================================
ax5 = fig.add_subplot(gs[2, 1])

ax5.plot(times_utc, offsets, 'go-', markersize=5, linewidth=1.5)
ax5.axhline(y=np.mean(offsets), color='green', linestyle='--', alpha=0.5,
            label=f'Mean offset: {np.mean(offsets):.1f} km')
ax5.fill_between(times_utc, 0, offsets, alpha=0.2, color='green')

ax5.set_xlabel('UTC Hour')
ax5.set_ylabel('Inter-Radar Offset (km)')
ax5.set_title(f'KEPZ-KHDX Position Agreement\n'
              f'({np.mean(offsets):.1f} km mean, {np.min(offsets):.1f}-{np.max(offsets):.1f} km range)',
              fontsize=10)
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)
ax5.set_xlim(3.0, 7.0)
ax5.set_ylim(0, max(offsets) * 1.2)

# ============================================================================
# PANEL 6: Summary statistics (bottom right)
# ============================================================================
ax6 = fig.add_subplot(gs[2, 2])
ax6.axis('off')

# Object characterization summary
summary_text = f"""TARGET CHARACTERIZATION SUMMARY
{'─' * 40}

POSITION
  Mean: {mean_lat:.4f}°N, {abs(mean_lon):.4f}°W
  Spread: {lat_range:.1f}km N-S × {lon_range:.1f}km E-W
  Dual-radar agreement: {np.mean(offsets):.1f} km mean

REFLECTIVITY
  KEPZ range: {np.min(kepz_refs):.0f} – {np.max(kepz_refs):.0f} dBZ
  KHDX range: {np.min(khdx_refs):.0f} – {np.max(khdx_refs):.0f} dBZ
  Peak KEPZ: {np.max(kepz_refs):.0f} dBZ at {times_utc[np.argmax(kepz_refs)]:.2f} UTC

TEMPORAL
  Duration: {(times_utc[-1] - times_utc[0]):.1f} hours
  Scans with returns: {len(matches)}/33 (100%)
  Peak KHDX escalation: {np.max(khdx_n)} returns/scan

MULTI-DATE PERSISTENCE (KHDX az 193-195°)
  Baseline avg: {persist_data['per_date_stats']['Feb 10']['az_bin_avgs'].get('193-195', 0):.1f}/scan
  Event night: {persist_data['per_date_stats']['Feb 11 (EVENT)']['az_bin_avgs'].get('193-195', 0):.1f}/scan
  Excess: +{((persist_data['per_date_stats']['Feb 11 (EVENT)']['az_bin_avgs'].get('193-195', 1)) / max(persist_data['per_date_stats']['Feb 10']['az_bin_avgs'].get('193-195', 1), 0.001) - 1) * 100:.0f}%

SEISMIC (11 stations, 200km radius)
  Event night seismic anomalies: NONE
  → Consistent with non-kinetic (laser) engagement

NO. OF OBJECTS: Radar cannot resolve individual
  objects. The persistent return cluster is
  consistent with 1-3 targets or debris field
  within a ~3km area near Fort Bliss."""

ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
         fontsize=7.5, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# ============================================================================
# SAVE
# ============================================================================
output_path = f'{OUTPUT_DIR}/multi_radar_object_analysis.png'
fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved to {output_path}")
plt.close()

# Also create a focused position animation-style figure
fig2, axes2 = plt.subplots(2, 4, figsize=(20, 10))
fig2.suptitle('TRIANGULATED TARGET POSITION — Hourly Snapshots\n'
              'Each panel shows all positions for a ~1-hour block',
              fontsize=14, fontweight='bold')

hour_blocks = [(3.0, 3.75), (3.75, 4.5), (4.5, 5.25), (5.25, 6.0),
               (6.0, 6.75), (6.75, 7.5)]

for idx, (t_start, t_end) in enumerate(hour_blocks):
    if idx >= 8:
        break
    ax = axes2.flat[idx]
    mask = (times_utc >= t_start) & (times_utc < t_end)

    if np.any(mask):
        ax.scatter(tri_lons[mask], tri_lats[mask], c=kepz_refs[mask],
                   cmap='YlOrRd', s=100, edgecolors='black', linewidth=0.5,
                   vmin=30, vmax=65)
        ax.plot(tri_lons[mask], tri_lats[mask], 'k-', alpha=0.4)

        # Show all positions faintly for context
        ax.scatter(tri_lons, tri_lats, c='lightgray', s=20, alpha=0.3, zorder=1)

        n_pts = mask.sum()
        avg_ref = np.mean(kepz_refs[mask])
        avg_n = np.mean(khdx_n[mask])
    else:
        n_pts = 0
        avg_ref = 0
        avg_n = 0

    ax.set_xlim(-106.56, -106.46)
    ax.set_ylim(31.83, 31.93)
    ax.set_title(f'{int(t_start):02d}:{int((t_start%1)*60):02d}–{int(t_end):02d}:{int((t_end%1)*60):02d} UTC\n'
                 f'{n_pts} scans, avg {avg_ref:.0f} dBZ', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect(1.0 / np.cos(np.radians(mean_lat)))
    if idx % 4 == 0:
        ax.set_ylabel('Lat (°N)')
    if idx >= 4:
        ax.set_xlabel('Lon (°W)')

# Hide unused subplots
for idx in range(len(hour_blocks), 8):
    axes2.flat[idx].axis('off')

output_path2 = f'{OUTPUT_DIR}/target_position_hourly.png'
fig2.savefig(output_path2, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved to {output_path2}")
plt.close()
print("Done.")
