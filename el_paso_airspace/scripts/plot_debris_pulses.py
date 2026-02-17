#!/usr/bin/env python3
"""
Debris Pulse Visualization — KHDX vs KEPZ Timeline
El Paso Airspace Incident — Feb 11, 2026

Generates annotated plots showing:
  1. KHDX max reflectivity + pixel counts with engagement markers
  2. KEPZ spike overlay for cross-radar comparison
  3. Estimated engagement events highlighted
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import json

# Load data
with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/khdx_analysis_results.json') as f:
    khdx_data = json.load(f)

with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/engagement_count_results.json') as f:
    engagement_data = json.load(f)

scans = khdx_data['scan_results']

# Extract time series
times_min = []
times_hr = []  # decimal hours for plotting
max_z = []
n15 = []
n10 = []
n5 = []

for scan in scans:
    t = scan['time']
    h, m = int(t[:2]), int(t[2:])
    t_minutes = h * 60 + m
    times_min.append(t_minutes)
    times_hr.append(h + m / 60)
    max_z.append(scan['max_ref_dBZ'])
    n15.append(scan['n_above_15dBZ'])
    n10.append(scan['n_above_10dBZ'])
    n5.append(scan['n_above_5dBZ'])

times_hr = np.array(times_hr)
max_z = np.array(max_z)
n15 = np.array(n15)
n10 = np.array(n10)
n5 = np.array(n5)

# Key event times (decimal hours)
NOTAM = 3 + 32/60
S1 = 4 + 24/60
S2 = 5 + 13/60
TFR = 6 + 30/60

# Engagement event times
engagement_times_hr = [3+25/60, 4+24/60, 4+57/60, 5+13/60, 5+40/60, 6+8/60]
engagement_labels = ['E1?', 'E2', 'E3', 'E4', 'E5', 'E6']
engagement_conf = ['Initial', 'Confirmed', 'Likely', 'Confirmed', 'Likely', 'Probable']

# ============================================================================
# FIGURE 1: Three-panel debris pulse analysis
# ============================================================================
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 14), sharex=True,
                                     gridspec_kw={'hspace': 0.08})

fig.suptitle('KHDX NEXRAD Debris Pulse Analysis — El Paso Airspace Incident\n'
             'Feb 11, 2026 | KHDX (Holloman AFB) — 135 km from Fort Bliss',
             fontsize=14, fontweight='bold', y=0.98)

# Colors
COLOR_Z = '#e63946'
COLOR_N15 = '#e76f51'
COLOR_N10 = '#f4a261'
COLOR_N5 = '#2a9d8f'
COLOR_NOTAM = '#457b9d'
COLOR_S1 = '#e63946'
COLOR_S2 = '#d62828'
COLOR_TFR = '#1d3557'
COLOR_ENGAGE = '#6a0dad'
COLOR_DEBRIS_RISE = '#ff6b35'

# Add vertical event lines to all panels
for ax in [ax1, ax2, ax3]:
    ax.axvline(NOTAM, color=COLOR_NOTAM, linestyle='--', alpha=0.7, linewidth=1.5)
    ax.axvline(S1, color=COLOR_S1, linestyle='-', alpha=0.8, linewidth=2)
    ax.axvline(S2, color=COLOR_S2, linestyle='-', alpha=0.8, linewidth=2)
    ax.axvline(TFR, color=COLOR_TFR, linestyle='--', alpha=0.7, linewidth=1.5)

    # Shade engagement events
    for et in engagement_times_hr:
        ax.axvline(et, color=COLOR_ENGAGE, linestyle=':', alpha=0.3, linewidth=1)

    ax.set_xlim(2.9, 7.1)
    ax.grid(True, alpha=0.2)
    ax.tick_params(labelsize=10)

# --- Panel 1: Max Reflectivity ---
ax1.fill_between(times_hr, max_z, alpha=0.15, color=COLOR_Z)
ax1.plot(times_hr, max_z, 'o-', color=COLOR_Z, markersize=6, linewidth=2, label='Max Z (dBZ)')

# Mark the sudden jumps
for jump in engagement_data.get('reflectivity_jumps', []):
    t_str = jump['time']
    h, m = int(t_str.split(':')[0]), int(t_str.split(':')[1])
    t_hr = h + m/60
    ax1.annotate(f'+{jump["jump_db"]:.0f} dB',
                xy=(t_hr, jump['new_z']), xytext=(t_hr + 0.1, jump['new_z'] + 3),
                fontsize=9, fontweight='bold', color=COLOR_DEBRIS_RISE,
                arrowprops=dict(arrowstyle='->', color=COLOR_DEBRIS_RISE, lw=1.5))

# Label the 37 dBZ peak
ax1.annotate('37.0 dBZ\nSTRONGEST', xy=(6+8/60, 37.0), xytext=(6.3, 38),
            fontsize=10, fontweight='bold', color='#d62828',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffe0e0', edgecolor='#d62828'),
            arrowprops=dict(arrowstyle='->', color='#d62828', lw=2))

ax1.set_ylabel('Max Reflectivity (dBZ)', fontsize=12, fontweight='bold')
ax1.set_ylim(5, 42)

# Add event labels at top
ax1.text(NOTAM, 41, 'NOTAM', fontsize=9, ha='center', color=COLOR_NOTAM, fontweight='bold')
ax1.text(S1, 41, 'KEPZ\nSpike 1', fontsize=9, ha='center', color=COLOR_S1, fontweight='bold')
ax1.text(S2, 41, 'KEPZ\nSpike 2', fontsize=9, ha='center', color=COLOR_S2, fontweight='bold')
ax1.text(TFR, 41, 'TFR', fontsize=9, ha='center', color=COLOR_TFR, fontweight='bold')

# Mark KEPZ spike reflectivities for comparison
ax1.axhspan(16, 20, alpha=0.08, color='blue')
ax1.text(3.0, 18, 'KEPZ spike\nrange\n(16-18 dBZ)', fontsize=8, va='center',
        color='blue', alpha=0.6, style='italic')

# --- Panel 2: Pixel counts (stacked) ---
ax2.fill_between(times_hr, n5, alpha=0.2, color=COLOR_N5, label='Pixels >5 dBZ')
ax2.fill_between(times_hr, n10, alpha=0.3, color=COLOR_N10, label='Pixels >10 dBZ')
ax2.fill_between(times_hr, n15, alpha=0.5, color=COLOR_N15, label='Pixels >15 dBZ')

ax2.plot(times_hr, n5, 'o-', color=COLOR_N5, markersize=4, linewidth=1.5)
ax2.plot(times_hr, n10, 's-', color=COLOR_N10, markersize=4, linewidth=1.5)
ax2.plot(times_hr, n15, '^-', color=COLOR_N15, markersize=5, linewidth=2)

ax2.set_ylabel('Pixel Count\n(Fort Bliss Sector)', fontsize=12, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9)
ax2.set_ylim(0, 50)

# Annotate the n15 peaks
peak_idx = np.argmax(n15)
ax2.annotate(f'{n15[peak_idx]} pixels >15 dBZ',
            xy=(times_hr[peak_idx], n15[peak_idx]),
            xytext=(times_hr[peak_idx] + 0.15, n15[peak_idx] + 5),
            fontsize=9, fontweight='bold', color=COLOR_N15,
            arrowprops=dict(arrowstyle='->', color=COLOR_N15, lw=1.5))

# --- Panel 3: Engagement events with debris pulse identification ---
# Plot n15 as the primary signal
ax3.bar(times_hr, n15, width=0.09, color=COLOR_N15, alpha=0.6, label='Debris pixels (>15 dBZ)')
ax3.plot(times_hr, n15, 'k-', linewidth=0.8, alpha=0.5)

# Highlight each engagement event
colors_engage = ['#ffd166', '#e63946', '#06d6a0', '#118ab2', '#ef476f', '#8338ec', '#ff6b35']
for i, (et, label, conf) in enumerate(zip(engagement_times_hr, engagement_labels, engagement_conf)):
    color = colors_engage[i % len(colors_engage)]
    # Shade a band around each event
    ax3.axvspan(et - 0.05, et + 0.2, alpha=0.15, color=color)

    # Arrow markers
    if conf == 'Confirmed':
        marker = '★'
        msize = 16
    elif conf == 'Likely':
        marker = '●'
        msize = 12
    else:
        marker = '○'
        msize = 12

    # Find nearest n15 value
    idx = np.argmin(np.abs(times_hr - et))
    y_val = n15[idx]

    ax3.annotate(f'{label}\n({conf})',
                xy=(et, max(y_val, 1)),
                xytext=(et, y_val + 8 + (i % 2) * 4),
                fontsize=9, fontweight='bold', ha='center',
                color=color if color != '#ffd166' else '#b8860b',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                         edgecolor=color, alpha=0.9),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

# Mark the quiet period at 04:50
quiet_hr = 4 + 50/60
ax3.annotate('QUIET\n(9 dBZ)', xy=(quiet_hr, 0), xytext=(quiet_hr - 0.2, 10),
            fontsize=9, color='green', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#e8f5e9', edgecolor='green'),
            arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

ax3.set_ylabel('Debris Pixels\n(>15 dBZ)', fontsize=12, fontweight='bold')
ax3.set_xlabel('Time (UTC) — Feb 11, 2026', fontsize=12, fontweight='bold')
ax3.set_ylim(0, 35)

# X-axis formatting
xticks = np.arange(3.0, 7.5, 0.5)
xtick_labels = [f'{int(h):02d}:{int((h % 1)*60):02d}' for h in xticks]
ax3.set_xticks(xticks)
ax3.set_xticklabels(xtick_labels, fontsize=10)

# Bottom annotation
ax3.text(0.5, -0.18, '5-7 distinct engagement events identified over 3.5 hours  |  '
         'KEPZ detected only 2 of these  |  KHDX reveals 3-5 additional KEPZ-invisible events\n'
         'Mean inter-engagement interval: ~33 minutes  |  '
         '37 dBZ peak at 06:08 UTC is 15+ dB above any KEPZ spike',
         transform=ax3.transAxes, fontsize=10, ha='center', va='top',
         bbox=dict(boxstyle='round', facecolor='#f0f0f0', edgecolor='gray', alpha=0.8))

plt.tight_layout(rect=[0, 0.05, 1, 0.95])

output_path = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/debris_pulse_timeline.png'
plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Figure 1 saved: {output_path}")

# ============================================================================
# FIGURE 2: KHDX vs KEPZ comparison — side by side
# ============================================================================
fig2, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 8))

fig2.suptitle('KHDX (Holloman, 135 km) vs KEPZ (El Paso, 30 km) — What Each Radar Sees',
              fontsize=14, fontweight='bold')

# Left panel: What KEPZ sees (2 spikes)
ax_left.set_title('KEPZ View (close radar, terrain-blocked)', fontsize=12, fontweight='bold',
                  color=COLOR_S1)

# Synthetic KEPZ profile (based on reported spike data)
kepz_times = np.linspace(3, 7, 100)
kepz_signal = np.zeros_like(kepz_times)

# Background noise ~3-5 dBZ
kepz_signal += 4 + np.random.normal(0, 0.5, len(kepz_times))

# Spike 1 at 04:24 — Gaussian pulse
s1_center = 4 + 24/60
kepz_signal += 14 * np.exp(-((kepz_times - s1_center) / 0.08)**2)

# Spike 2 at 05:13 — smaller Gaussian
s2_center = 5 + 13/60
kepz_signal += 7 * np.exp(-((kepz_times - s2_center) / 0.08)**2)

# Clutter suppression dip before Spike 1
dip_center = 4 + 10/60
kepz_signal -= 3 * np.exp(-((kepz_times - dip_center) / 0.12)**2)
kepz_signal = np.maximum(kepz_signal, 0)

ax_left.fill_between(kepz_times, kepz_signal, alpha=0.3, color='#2196F3')
ax_left.plot(kepz_times, kepz_signal, color='#1565C0', linewidth=2)

ax_left.axvline(S1, color=COLOR_S1, linewidth=2, linestyle='-')
ax_left.axvline(S2, color=COLOR_S2, linewidth=2, linestyle='-')
ax_left.axvline(NOTAM, color=COLOR_NOTAM, linewidth=1.5, linestyle='--')
ax_left.axvline(TFR, color=COLOR_TFR, linewidth=1.5, linestyle='--')

ax_left.annotate('Spike 1\n551 px, 18 dBZ\nZDR 8-12, CC 0.5',
                xy=(S1, 18), xytext=(S1 + 0.3, 20),
                fontsize=10, fontweight='bold', color=COLOR_S1,
                bbox=dict(boxstyle='round', facecolor='#ffe0e0'),
                arrowprops=dict(arrowstyle='->', color=COLOR_S1, lw=2))

ax_left.annotate('Spike 2\n456 px, 11 dBZ',
                xy=(S2, 11), xytext=(S2 + 0.3, 14),
                fontsize=10, fontweight='bold', color=COLOR_S2,
                bbox=dict(boxstyle='round', facecolor='#ffe0e0'),
                arrowprops=dict(arrowstyle='->', color=COLOR_S2, lw=2))

ax_left.text(4.0, 16, '← Clutter\n   suppression\n   dip', fontsize=9,
            color='#555', style='italic')

ax_left.text(5.8, 16, 'KEPZ sees\nNOTHING here\n(Franklin Mtns\nblocking)', fontsize=11,
            ha='center', color='#999', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#f5f5f5', edgecolor='#ccc'))

ax_left.set_ylabel('Reflectivity (dBZ)', fontsize=12)
ax_left.set_xlabel('Time (UTC)', fontsize=12)
ax_left.set_ylim(0, 25)
ax_left.set_xlim(3, 7)
xtick_labels_l = [f'{int(h):02d}:{int((h % 1)*60):02d}' for h in np.arange(3, 7.5, 0.5)]
ax_left.set_xticks(np.arange(3, 7.5, 0.5))
ax_left.set_xticklabels(xtick_labels_l, fontsize=9)
ax_left.grid(True, alpha=0.2)

# KEPZ count
ax_left.text(0.95, 0.95, 'KEPZ sees:\n2 events', transform=ax_left.transAxes,
            fontsize=14, fontweight='bold', ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#e3f2fd', edgecolor='#1565C0'),
            color='#1565C0')

# Right panel: What KHDX sees (5-7 events)
ax_right.set_title('KHDX View (far radar, unblocked, altitude >1000m AGL)', fontsize=12,
                   fontweight='bold', color='#8338ec')

# KHDX data
ax_right.fill_between(times_hr, max_z, alpha=0.15, color='#8338ec')
ax_right.plot(times_hr, max_z, 'o-', color='#8338ec', markersize=5, linewidth=2)

# Engagement event markers
for i, (et, label, conf) in enumerate(zip(engagement_times_hr, engagement_labels, engagement_conf)):
    color = colors_engage[i % len(colors_engage)]
    idx = np.argmin(np.abs(times_hr - et))
    y = max_z[idx]

    if conf == 'Confirmed':
        ax_right.plot(et, y, '*', markersize=18, color=color, markeredgecolor='black',
                     markeredgewidth=1, zorder=5)
    else:
        ax_right.plot(et, y, 'o', markersize=12, color=color, markeredgecolor='black',
                     markeredgewidth=1, zorder=5)

    ax_right.annotate(label, xy=(et, y), xytext=(et, y + 3 + (i % 2) * 2),
                     fontsize=10, fontweight='bold', ha='center',
                     color=color if color != '#ffd166' else '#b8860b')

ax_right.axvline(NOTAM, color=COLOR_NOTAM, linewidth=1.5, linestyle='--')
ax_right.axvline(S1, color=COLOR_S1, linewidth=2, linestyle='-', alpha=0.5)
ax_right.axvline(S2, color=COLOR_S2, linewidth=2, linestyle='-', alpha=0.5)
ax_right.axvline(TFR, color=COLOR_TFR, linewidth=1.5, linestyle='--')

ax_right.set_ylabel('Max Reflectivity (dBZ)', fontsize=12)
ax_right.set_xlabel('Time (UTC)', fontsize=12)
ax_right.set_ylim(5, 42)
ax_right.set_xlim(3, 7)
ax_right.set_xticks(np.arange(3, 7.5, 0.5))
ax_right.set_xticklabels(xtick_labels_l, fontsize=9)
ax_right.grid(True, alpha=0.2)

# KHDX count
ax_right.text(0.95, 0.95, 'KHDX sees:\n5-7 events', transform=ax_right.transAxes,
             fontsize=14, fontweight='bold', ha='right', va='top',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#f3e5f5', edgecolor='#8338ec'),
             color='#8338ec')

# Bottom legend
fig2.text(0.5, 0.02,
         '★ = Confirmed engagement (both radars)  |  ● = Likely/Probable engagement (KHDX only)\n'
         'The Franklin Mountains block KEPZ from seeing most of the activity.',
         fontsize=11, ha='center',
         bbox=dict(boxstyle='round', facecolor='#f8f8f8', edgecolor='gray'))

plt.tight_layout(rect=[0, 0.06, 1, 0.95])

output2 = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/kepz_vs_khdx_comparison.png'
plt.savefig(output2, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Figure 2 saved: {output2}")

# ============================================================================
# FIGURE 3: Engagement timeline bar chart
# ============================================================================
fig3, ax = plt.subplots(figsize=(16, 6))

ax.set_title('Estimated Weapon Engagement Timeline — 5-7 Directed-Energy Firings\n'
             'El Paso Airspace Incident, Feb 11, 2026',
             fontsize=14, fontweight='bold')

# Each engagement as a horizontal bar
eng_data = [
    ('E1? (Initial)', 3+17/60, 3+32/60, 'Initial', '#ffd166', 'KHDX only — pre-NOTAM detection\nMay be target RCS, not debris'),
    ('E2 (Primary)', 4+20/60, 4+40/60, 'Confirmed', '#e63946', 'KEPZ Spike 1 — 551 px, 18 dBZ\nZDR 8-12, CC 0.5-0.7'),
    ('E3 (Hidden)', 4+50/60, 5+5/60, 'Likely', '#06d6a0', 'KHDX-only — 26.5 dBZ after quiet\n+17.5 dB jump in one scan'),
    ('E4 (Second)', 5+10/60, 5+22/60, 'Confirmed', '#118ab2', 'KEPZ Spike 2 — 456 px, 11 dBZ\nWeaker: smaller/damaged target'),
    ('E5 (Escalation)', 5+30/60, 5+55/60, 'Likely', '#ef476f', 'KHDX 30-31 dBZ — STRONGER than S1/S2\nReturns increasing, not decaying'),
    ('E6 (Maximum)', 6+5/60, 6+25/60, 'Probable', '#8338ec', 'KHDX 37 dBZ — strongest return\n15+ dB above any KEPZ spike'),
]

y_positions = range(len(eng_data))

for i, (label, t_start, t_end, conf, color, notes) in enumerate(eng_data):
    # Main bar
    bar_height = 0.6
    ax.barh(i, t_end - t_start, left=t_start, height=bar_height,
           color=color, alpha=0.7, edgecolor='black', linewidth=1.5)

    # Label inside bar
    mid = (t_start + t_end) / 2
    ax.text(mid, i, label, ha='center', va='center', fontsize=10, fontweight='bold',
           color='black' if color in ['#ffd166', '#06d6a0'] else 'white')

    # Confidence badge
    badge_colors = {
        'Initial': '#f0f0f0',
        'Confirmed': '#c8e6c9',
        'Likely': '#fff9c4',
        'Probable': '#ffe0b2',
    }
    ax.text(t_end + 0.03, i + 0.15, conf, fontsize=8, fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.2', facecolor=badge_colors.get(conf, 'white'),
                    edgecolor='gray'), va='center')

    # Notes
    ax.text(t_end + 0.03, i - 0.15, notes, fontsize=7, va='center', color='#555',
           style='italic')

# Event lines
for t, label, color in [(NOTAM, 'NOTAM', COLOR_NOTAM), (S1, 'KEPZ S1', COLOR_S1),
                         (S2, 'KEPZ S2', COLOR_S2), (TFR, 'TFR', COLOR_TFR)]:
    ax.axvline(t, color=color, linestyle='--' if label in ['NOTAM', 'TFR'] else '-',
              linewidth=1.5, alpha=0.6, zorder=0)
    ax.text(t, len(eng_data) - 0.3, label, fontsize=9, ha='center', color=color,
           fontweight='bold', rotation=0)

ax.set_yticks([])
ax.set_xlabel('Time (UTC) — Feb 11, 2026', fontsize=12, fontweight='bold')
ax.set_xlim(3, 7)

xticks = np.arange(3.0, 7.5, 0.25)
xtick_labels = [f'{int(h):02d}:{int((h % 1)*60):02d}' for h in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(xtick_labels, fontsize=8, rotation=45)
ax.set_ylim(-0.5, len(eng_data) - 0.2)

ax.invert_yaxis()
ax.grid(True, axis='x', alpha=0.2)

# Summary box
summary = ('Mean interval between engagements: ~33 min\n'
           'KEPZ detected: 2 of 6 events (Franklin Mtns blockage)\n'
           'KHDX reveals: escalating response over 3.5 hours')
ax.text(0.01, 0.99, summary, transform=ax.transAxes, fontsize=9,
       va='top', ha='left',
       bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9))

plt.tight_layout()

output3 = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/engagement_timeline_bars.png'
plt.savefig(output3, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Figure 3 saved: {output3}")

print("\nAll 3 figures generated successfully.")
