#!/usr/bin/env python3
"""
Generate comprehensive visualizations for the El Paso airspace incident.

Produces:
1. KEPZ PPI radar plot (event night, peak scan)
2. KHDX PPI radar plot (event night, peak scan)
3. Dual-radar spatial map with triangulated positions
4. Time-series comparison (baseline vs event for both radars)
5. KHDX temporal escalation profile
6. South cluster dual-pol scatter plot
"""

import numpy as np
import pyart
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import json
import os

OUTPUT_DIR = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/'

# Radar positions
KEPZ_LAT, KEPZ_LON = 31.8731, -106.6981
KHDX_LAT, KHDX_LON = 33.0803, -106.1225

# Load pre-computed results
with open(os.path.join(OUTPUT_DIR, 'baseline_comparison_results.json')) as f:
    baseline_data = json.load(f)

with open(os.path.join(OUTPUT_DIR, 'integrated_triangulation_results.json')) as f:
    tri_data = json.load(f)

with open(os.path.join(OUTPUT_DIR, 'kepz_level2_results.json')) as f:
    kepz_data = json.load(f)

with open(os.path.join(OUTPUT_DIR, 'khdx_analysis_results.json')) as f:
    khdx_data = json.load(f)


# ============================================================================
# FIGURE 1: KEPZ PPI Radar Plot (Peak Scan)
# ============================================================================
print("Generating Figure 1: KEPZ PPI Radar Plot...")

# Find peak scan file
peak_scan = max(kepz_data['scan_results'], key=lambda s: s['max_ref_dBZ'])
peak_file = os.path.join('/home/user/uap-transient-research/el_paso_airspace/kepz_data/', peak_scan['file'])

radar = pyart.io.read_nexrad_archive(peak_file)

fig, ax = plt.subplots(1, 1, figsize=(12, 10))
display = pyart.graph.RadarDisplay(radar)
display.plot_ppi('reflectivity', 0, ax=ax, vmin=-10, vmax=65,
                 title=f'KEPZ Reflectivity — {peak_scan["time"]} UTC Feb 11, 2026\n'
                       f'Peak: {peak_scan["max_ref_dBZ"]:.1f} dBZ | '
                       f'Event Night (Lowest Elevation Sweep)',
                 colorbar_label='Reflectivity (dBZ)',
                 cmap='NWSRef')
display.set_limits(xlim=(-60, 60), ylim=(-60, 60), ax=ax)

# Draw target sector
from matplotlib.patches import Wedge
wedge = Wedge((0, 0), 55, 60, 110, width=40, fill=False, edgecolor='red',
              linewidth=2, linestyle='--', label='Target sector (az 60-110°)')
ax.add_patch(wedge)
ax.legend(loc='upper right', fontsize=10)

fig.savefig(os.path.join(OUTPUT_DIR, 'fig1_kepz_ppi.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
del radar
print("  Saved fig1_kepz_ppi.png")


# ============================================================================
# FIGURE 2: KHDX PPI Radar Plot (Peak Scan)
# ============================================================================
print("Generating Figure 2: KHDX PPI Radar Plot...")

khdx_peak = max(khdx_data['scan_results'], key=lambda s: s['max_ref_dBZ'])
khdx_peak_file = os.path.join('/home/user/uap-transient-research/el_paso_airspace/khdx_data/', khdx_peak['file'])

radar = pyart.io.read_nexrad_archive(khdx_peak_file)

fig, ax = plt.subplots(1, 1, figsize=(12, 10))
display = pyart.graph.RadarDisplay(radar)
display.plot_ppi('reflectivity', 0, ax=ax, vmin=-10, vmax=40,
                 title=f'KHDX Reflectivity — {khdx_peak["time"]} UTC Feb 11, 2026\n'
                       f'Peak in Fort Bliss sector: {khdx_peak["max_ref_dBZ"]:.1f} dBZ | '
                       f'135 km range independent confirmation',
                 colorbar_label='Reflectivity (dBZ)',
                 cmap='NWSRef')
display.set_limits(xlim=(-160, 160), ylim=(-160, 160), ax=ax)

# Draw Fort Bliss sector
wedge = Wedge((0, 0), 148, 185, 199, width=23, fill=False, edgecolor='red',
              linewidth=2, linestyle='--', label='Fort Bliss sector (az 185-199°)')
ax.add_patch(wedge)
ax.legend(loc='upper right', fontsize=10)

fig.savefig(os.path.join(OUTPUT_DIR, 'fig2_khdx_ppi.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
del radar
print("  Saved fig2_khdx_ppi.png")


# ============================================================================
# FIGURE 3: Dual-Radar Spatial Map
# ============================================================================
print("Generating Figure 3: Dual-Radar Spatial Map...")

fig, ax = plt.subplots(1, 1, figsize=(14, 12))

# Plot KEPZ returns (south cluster only, event night)
kepz_south = [r for r in kepz_data['significant_returns']
              if r['az'] >= 90 and r['ref_dBZ'] > 15]
if kepz_south:
    ks_lats = [r['lat'] for r in kepz_south]
    ks_lons = [r['lon'] for r in kepz_south]
    ks_refs = [r['ref_dBZ'] for r in kepz_south]
    sc1 = ax.scatter(ks_lons, ks_lats, c=ks_refs, cmap='YlOrRd', s=3, alpha=0.3,
                     vmin=15, vmax=50, label=f'KEPZ south cluster (n={len(kepz_south)})')

# Plot KHDX returns
if khdx_data['significant_returns']:
    kh_lats = [r['lat'] for r in khdx_data['significant_returns']]
    kh_lons = [r['lon'] for r in khdx_data['significant_returns']]
    kh_refs = [r['ref_dBZ'] for r in khdx_data['significant_returns']]
    sc2 = ax.scatter(kh_lons, kh_lats, c=kh_refs, cmap='Blues', s=15, alpha=0.5,
                     vmin=5, vmax=30, marker='s', label=f'KHDX returns (n={len(khdx_data["significant_returns"])})')

# Plot triangulated positions
if tri_data['triangulation']['matched_scans']:
    tri_lats = [m['triangulated_lat'] for m in tri_data['triangulation']['matched_scans']]
    tri_lons = [m['triangulated_lon'] for m in tri_data['triangulation']['matched_scans']]
    ax.scatter(tri_lons, tri_lats, c='lime', s=80, marker='*', edgecolors='black',
               linewidths=0.5, zorder=10, label='Triangulated positions')

    # Mean position
    mean_lat = np.mean(tri_lats)
    mean_lon = np.mean(tri_lons)
    ax.scatter([mean_lon], [mean_lat], c='red', s=200, marker='*', edgecolors='black',
               linewidths=1.5, zorder=11, label=f'Mean position ({mean_lat:.4f}°N, {abs(mean_lon):.4f}°W)')

# Radar locations
ax.scatter([KEPZ_LON], [KEPZ_LAT], c='blue', s=150, marker='^', edgecolors='black',
           linewidths=1, zorder=12, label='KEPZ radar')
ax.scatter([KHDX_LON], [KHDX_LAT], c='green', s=150, marker='^', edgecolors='black',
           linewidths=1, zorder=12, label='KHDX radar (135 km N)')

# KBIF
ax.scatter([-106.38], [31.8495], c='orange', s=100, marker='D', edgecolors='black',
           linewidths=1, zorder=12, label='KBIF (Biggs AAF)')

# Annotations
ax.annotate('Franklin Mountains', xy=(-106.52, 31.91), fontsize=9,
            fontstyle='italic', color='brown', ha='center')
ax.annotate('Fort Bliss', xy=(-106.40, 31.84), fontsize=9,
            fontstyle='italic', color='darkred', ha='center')
ax.annotate('El Paso (residential)', xy=(-106.45, 31.78), fontsize=9,
            fontstyle='italic', color='purple', ha='center')

ax.set_xlabel('Longitude (°W)', fontsize=12)
ax.set_ylabel('Latitude (°N)', fontsize=12)
ax.set_title('Dual-Radar Spatial Map — KEPZ + KHDX Triangulation\n'
             'El Paso Airspace Incident, Feb 11, 2026, 03:00-07:00 UTC',
             fontsize=13, fontweight='bold')
ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
ax.set_xlim(-106.75, -106.20)
ax.set_ylim(31.70, 32.00)
ax.grid(True, alpha=0.3)
ax.set_aspect(1.0 / np.cos(np.radians(31.85)))

fig.savefig(os.path.join(OUTPUT_DIR, 'fig3_spatial_map.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved fig3_spatial_map.png")


# ============================================================================
# FIGURE 4: Baseline vs Event Time-Series (Both Radars)
# ============================================================================
print("Generating Figure 4: Baseline vs Event Comparison...")

fig, axes = plt.subplots(2, 2, figsize=(16, 10))
fig.suptitle('Baseline (Feb 10) vs Event (Feb 11) — Both Radars', fontsize=14, fontweight='bold')

# KEPZ >10 dBZ
ax = axes[0, 0]
bl_scans = baseline_data['baseline_scans']
ev_scans = baseline_data['event_scans']

bl_times = list(range(len(bl_scans)))
ev_times = list(range(len(ev_scans)))
bl_n10 = [s['n_above_10dBZ'] for s in bl_scans]
ev_n10 = [s['n_above_10dBZ'] for s in ev_scans]

ax.fill_between(bl_times, bl_n10, alpha=0.3, color='blue', label='Baseline')
ax.fill_between(ev_times, ev_n10, alpha=0.3, color='red', label='Event')
ax.plot(bl_times, bl_n10, 'b-', linewidth=1)
ax.plot(ev_times, ev_n10, 'r-', linewidth=1)
ax.axhline(np.mean(bl_n10), color='blue', linestyle='--', alpha=0.5)
ax.axhline(np.mean(ev_n10), color='red', linestyle='--', alpha=0.5)
ax.set_title('KEPZ: Pixels >10 dBZ per Scan', fontsize=11)
ax.set_ylabel('Count')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# KEPZ >20 dBZ
ax = axes[0, 1]
bl_n20 = [s['n_above_20dBZ'] for s in bl_scans]
ev_n20 = [s['n_above_20dBZ'] for s in ev_scans]

ax.fill_between(bl_times, bl_n20, alpha=0.3, color='blue', label='Baseline')
ax.fill_between(ev_times, ev_n20, alpha=0.3, color='red', label='Event')
ax.plot(bl_times, bl_n20, 'b-', linewidth=1)
ax.plot(ev_times, ev_n20, 'r-', linewidth=1)
ax.axhline(np.mean(bl_n20), color='blue', linestyle='--', alpha=0.5)
ax.axhline(np.mean(ev_n20), color='red', linestyle='--', alpha=0.5)
ax.set_title('KEPZ: Pixels >20 dBZ per Scan (+49% on event night)', fontsize=11)
ax.set_ylabel('Count')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# KHDX returns per scan
ax = axes[1, 0]
khdx_bl = tri_data['khdx_baseline_comparison']
# Reconstruct from stored data — use event scans from khdx_data
khdx_ev_scans = khdx_data['scan_results']
khdx_ev_n5 = [s['n_above_5dBZ'] for s in khdx_ev_scans]
khdx_ev_n10 = [s['n_above_10dBZ'] for s in khdx_ev_scans]

khdx_bl_avg = khdx_bl['baseline_avg_returns']
khdx_bl_avg10 = khdx_bl['baseline_avg_n10']

ax.bar(range(len(khdx_ev_n5)), khdx_ev_n5, color='red', alpha=0.4, label='Event >5 dBZ')
ax.bar(range(len(khdx_ev_n10)), khdx_ev_n10, color='darkred', alpha=0.7, label='Event >10 dBZ')
ax.axhline(khdx_bl_avg, color='blue', linestyle='--', linewidth=2,
           label=f'Baseline avg >5 dBZ ({khdx_bl_avg:.0f})')
ax.axhline(khdx_bl_avg10, color='navy', linestyle=':', linewidth=2,
           label=f'Baseline avg >10 dBZ ({khdx_bl_avg10:.0f})')
ax.set_title('KHDX: Returns per Scan (event vs baseline average)', fontsize=11)
ax.set_ylabel('Count')
ax.set_xlabel('Scan number')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# KHDX max reflectivity
ax = axes[1, 1]
khdx_ev_max = [s['max_ref_dBZ'] for s in khdx_ev_scans]
khdx_bl_avg_max = tri_data['khdx_baseline_comparison'].get('baseline_avg_returns', 17.6)

colors = ['red' if m > 20 else 'orange' if m > 15 else 'gray' for m in khdx_ev_max]
ax.bar(range(len(khdx_ev_max)), khdx_ev_max, color=colors, alpha=0.7)
ax.axhline(17.6, color='blue', linestyle='--', linewidth=2,
           label='Baseline avg max (17.6 dBZ)')
ax.set_title('KHDX: Peak Reflectivity per Scan', fontsize=11)
ax.set_ylabel('dBZ')
ax.set_xlabel('Scan number')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Annotate the escalation
for i, (m, t) in enumerate(zip(khdx_ev_max, [s['time'] for s in khdx_ev_scans])):
    if m > 25:
        ax.annotate(f'{t}\n{m:.0f}dBZ', xy=(i, m), fontsize=7,
                    ha='center', va='bottom', color='darkred')

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'fig4_baseline_comparison.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved fig4_baseline_comparison.png")


# ============================================================================
# FIGURE 5: KHDX Escalation Timeline (baseline-subtracted)
# ============================================================================
print("Generating Figure 5: KHDX Escalation Timeline...")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
fig.suptitle('KHDX Temporal Escalation Profile (baseline-subtracted)\n'
             'Fort Bliss sector, Feb 11, 2026', fontsize=13, fontweight='bold')

times = [s['time'] for s in khdx_ev_scans]
time_labels = [f"{t[:2]}:{t[2:]}" for t in times]

# Panel 1: Excess returns
excess_5 = [s['n_above_5dBZ'] - khdx_bl_avg for s in khdx_ev_scans]
excess_10 = [s['n_above_10dBZ'] - khdx_bl_avg10 for s in khdx_ev_scans]

colors_5 = ['red' if e > 10 else 'orange' if e > 5 else 'green' if e > 0 else 'gray' for e in excess_5]
ax1.bar(range(len(excess_5)), excess_5, color=colors_5, alpha=0.6, label='Excess >5 dBZ')
ax1.bar(range(len(excess_10)), excess_10, color='darkred', alpha=0.8, width=0.4,
        label='Excess >10 dBZ')
ax1.axhline(0, color='black', linewidth=1)
ax1.set_ylabel('Excess returns above baseline')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_title('Returns above baseline (positive = more than normal)')

# Panel 2: Max reflectivity
ax2.plot(range(len(khdx_ev_max)), khdx_ev_max, 'r-o', markersize=4, label='Event night max dBZ')
ax2.axhline(17.6, color='blue', linestyle='--', linewidth=2, label='Baseline avg max')
ax2.fill_between(range(len(khdx_ev_max)), 17.6, khdx_ev_max,
                 where=[m > 17.6 for m in khdx_ev_max],
                 alpha=0.3, color='red', label='Excess above baseline')
ax2.set_ylabel('Peak Reflectivity (dBZ)')
ax2.set_xlabel('UTC Time')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# Set x-tick labels
ax2.set_xticks(range(0, len(times), 2))
ax2.set_xticklabels([time_labels[i] for i in range(0, len(times), 2)], rotation=45, fontsize=8)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'fig5_khdx_escalation.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved fig5_khdx_escalation.png")


# ============================================================================
# FIGURE 6: South Cluster Dual-Pol Analysis
# ============================================================================
print("Generating Figure 6: South Cluster Dual-Pol...")

# Get south cluster returns with dual-pol
south_event = [r for r in kepz_data['significant_returns']
               if r['az'] >= 90 and r.get('zdr') is not None and r.get('cc') is not None
               and r['ref_dBZ'] > 10]

# Get equivalent from baseline
bl_south = [r for r in baseline_data.get('baseline_scans', [])
            if r.get('south_cluster', {}).get('cc_mean') is not None]

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('KEPZ South Cluster (az 90-110°) Dual-Pol Signatures\n'
             'Returns >10 dBZ only', fontsize=13, fontweight='bold')

if south_event:
    se_zdr = [r['zdr'] for r in south_event]
    se_cc = [r['cc'] for r in south_event]
    se_ref = [r['ref_dBZ'] for r in south_event]

    # ZDR vs CC scatter
    ax = axes[0]
    sc = ax.scatter(se_cc, se_zdr, c=se_ref, cmap='hot_r', s=3, alpha=0.3, vmin=10, vmax=50)
    plt.colorbar(sc, ax=ax, label='Reflectivity (dBZ)')
    ax.set_xlabel('CC (Correlation Coefficient)')
    ax.set_ylabel('ZDR (Differential Reflectivity, dB)')
    ax.set_title('Event Night: ZDR vs CC')
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0.85, color='green', linestyle='--', alpha=0.5, label='Ground clutter threshold')
    ax.axvline(0.50, color='red', linestyle='--', alpha=0.5, label='Chaff/debris threshold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(-15, 22)

    # CC histogram — event vs baseline
    ax = axes[1]

    # Get baseline south cluster returns
    bl_south_rets = []
    for scan in baseline_data.get('baseline_scans', []):
        sc_data = scan.get('south_cluster', {})
        if sc_data.get('cc_mean') is not None:
            bl_south_rets.append(sc_data['cc_mean'])

    ax.hist(se_cc, bins=50, alpha=0.5, color='red', density=True, label=f'Event (n={len(se_cc)})')
    ax.axvline(np.mean(se_cc), color='red', linestyle='--', linewidth=2)
    if bl_south_rets:
        ax.axvline(np.mean(bl_south_rets), color='blue', linestyle='--', linewidth=2,
                   label=f'Baseline mean CC ({np.mean(bl_south_rets):.3f})')
    ax.axvline(np.mean(se_cc), color='red', linestyle='--', linewidth=2,
               label=f'Event mean CC ({np.mean(se_cc):.3f})')
    ax.set_xlabel('CC')
    ax.set_ylabel('Density')
    ax.set_title('CC Distribution — South Cluster')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ZDR histogram
    ax = axes[2]
    ax.hist(se_zdr, bins=60, alpha=0.5, color='red', density=True, label=f'Event (n={len(se_zdr)})')
    ax.axvline(np.mean(se_zdr), color='red', linestyle='--', linewidth=2,
               label=f'Event mean ZDR ({np.mean(se_zdr):.2f} dB)')
    ax.axvline(0, color='gray', linewidth=0.5)
    ax.set_xlabel('ZDR (dB)')
    ax.set_ylabel('Density')
    ax.set_title('ZDR Distribution — South Cluster')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'fig6_south_cluster_dualpol.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved fig6_south_cluster_dualpol.png")


# ============================================================================
# FIGURE 7: Summary Dashboard
# ============================================================================
print("Generating Figure 7: Summary Dashboard...")

fig = plt.figure(figsize=(18, 14))
fig.suptitle('EL PASO AIRSPACE INCIDENT — INTEGRATED RADAR ANALYSIS\n'
             'February 11, 2026, 03:00-07:00 UTC',
             fontsize=16, fontweight='bold', y=0.98)

# Layout: 3 rows, uneven columns
gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)

# Panel A: KHDX escalation
ax = fig.add_subplot(gs[0, :2])
excess_10_vals = [s['n_above_10dBZ'] - khdx_bl_avg10 for s in khdx_ev_scans]
colors = ['#d32f2f' if e > 5 else '#ff9800' if e > 2 else '#4caf50' if e > 0 else '#9e9e9e'
          for e in excess_10_vals]
ax.bar(range(len(excess_10_vals)), excess_10_vals, color=colors, alpha=0.8)
ax.axhline(0, color='black', linewidth=1)
ax.set_title('A. KHDX Excess Returns Above Baseline (>10 dBZ)', fontsize=11, fontweight='bold')
ax.set_ylabel('Excess count')
ax.set_xticks(range(0, len(times), 3))
ax.set_xticklabels([time_labels[i] for i in range(0, len(times), 3)], fontsize=7, rotation=45)
ax.grid(True, alpha=0.3)

# Panel B: KEPZ baseline vs event
ax = fig.add_subplot(gs[0, 2:])
categories = ['KEPZ\n>10 dBZ', 'KEPZ\n>20 dBZ', 'KEPZ South\n>20 dBZ',
              'KHDX\n>5 dBZ', 'KHDX\n>10 dBZ']
bl_vals = [764, 225, 128, 17.2, 5.6]
ev_vals = [869, 334, 215, 24.9, 10.0]

x = np.arange(len(categories))
w = 0.35
bars1 = ax.bar(x - w/2, bl_vals, w, label='Baseline (Feb 10)', color='#42a5f5', alpha=0.8)
bars2 = ax.bar(x + w/2, ev_vals, w, label='Event (Feb 11)', color='#ef5350', alpha=0.8)
ax.set_title('B. Baseline vs Event — Returns per Scan', fontsize=11, fontweight='bold')
ax.set_ylabel('Avg returns/scan')
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=8)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

# Add percentage labels
for i, (b, e) in enumerate(zip(bl_vals, ev_vals)):
    pct = 100 * (e - b) / max(b, 0.1)
    ax.annotate(f'+{pct:.0f}%', xy=(x[i] + w/2, e), fontsize=8, ha='center', va='bottom',
                color='darkred', fontweight='bold')

# Panel C: Spatial map (simplified)
ax = fig.add_subplot(gs[1, :2])
if tri_data['triangulation']['matched_scans']:
    tri_lats = [m['triangulated_lat'] for m in tri_data['triangulation']['matched_scans']]
    tri_lons = [m['triangulated_lon'] for m in tri_data['triangulation']['matched_scans']]
    ax.scatter(tri_lons, tri_lats, c='red', s=40, alpha=0.5, label='Triangulated positions')
    mean_lat = np.mean(tri_lats)
    mean_lon = np.mean(tri_lons)
    ax.scatter([mean_lon], [mean_lat], c='yellow', s=200, marker='*', edgecolors='red',
               linewidths=2, zorder=10, label=f'Mean: {mean_lat:.4f}°N')

ax.scatter([KEPZ_LON], [KEPZ_LAT], c='blue', s=100, marker='^', edgecolors='black',
           zorder=10, label='KEPZ')
ax.scatter([-106.38], [31.8495], c='orange', s=80, marker='D', edgecolors='black',
           zorder=10, label='Biggs AAF')

ax.set_title('C. Triangulated Event Positions', fontsize=11, fontweight='bold')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.legend(fontsize=7, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlim(-106.65, -106.35)
ax.set_ylim(31.80, 31.95)
ax.set_aspect(1.0 / np.cos(np.radians(31.85)))

# Panel D: CC comparison
ax = fig.add_subplot(gs[1, 2:])
se_cc_all = [r['cc'] for r in kepz_data['significant_returns']
             if r.get('cc') is not None and r['ref_dBZ'] > 10]

ax.hist(se_cc_all, bins=50, alpha=0.6, color='red', density=True, label='Event night')
ax.axvline(0.34, color='red', linestyle='--', linewidth=2, label='Event mean (0.34)')
ax.axvline(0.40, color='blue', linestyle='--', linewidth=2, label='Baseline mean (0.40)')
ax.axvspan(0.85, 1.05, alpha=0.1, color='green', label='Ground clutter range')
ax.axvspan(0.20, 0.50, alpha=0.1, color='red', label='Chaff/debris range')
ax.set_title('D. CC Distribution — Sector Returns >10 dBZ', fontsize=11, fontweight='bold')
ax.set_xlabel('Correlation Coefficient')
ax.set_ylabel('Density')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1.05)

# Panel E: Text summary
ax = fig.add_subplot(gs[2, :])
ax.axis('off')
summary_text = """
KEY FINDINGS (ALL BASELINE-CORRECTED):

• KEPZ south cluster shows +67% excess strong returns (>20 dBZ) above baseline, with a +1.58 dB ZDR shift
• KHDX independently confirms +78% excess >10 dBZ returns above its own baseline (p = 0.003)
• All 33 KEPZ scan-pairs have spatially matching KHDX returns — mean offset 5.7 km
• KHDX shows clear temporal escalation: 12 of 34 scans have significant excess, peaking at 05:47 UTC (31 dBZ)
• Triangulated position: 31.87°N, 106.52°W (Fort Bliss area, near residential El Paso)
• No weather present: clear skies, calm winds (METAR), dry stable atmosphere (sounding)
• The CC = 0.34 finding is WEAKER than originally claimed (baseline CC = 0.40), but a real 0.06 shift exists

BOTTOM LINE: Two independent radars both show statistically significant excess returns above their baselines,
confirmed at the same location and time. The KHDX escalation pattern is the strongest single piece of evidence.
"""
ax.text(0.02, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

fig.savefig(os.path.join(OUTPUT_DIR, 'fig7_summary_dashboard.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved fig7_summary_dashboard.png")

print("\nAll visualizations complete.")
