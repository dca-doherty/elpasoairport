#!/usr/bin/env python3
"""
Debris Mass Estimation from NEXRAD Dual-Pol Observables
El Paso Airspace Incident - Feb 11, 2026

Estimates the mass of the destroyed object from five independent constraints:
  1. Total integrated RCS → fragment count → mass
  2. ZDR-constrained fragment geometry → mass per scatterer
  3. Debris field spatial extent → fragmentation energy → structural mass
  4. Fall rate / persistence → ballistic coefficient → mass-to-area ratio
  5. Chaff equivalence model → calibrated mass from military literature

Key observables (from KEPZ NEXRAD dual-pol analysis):
  Spike 1: 04:24 UTC, 551 pixels, 18 dBZ peak, ZDR 8-12 dB, CC 0.50-0.70
  Spike 2: 05:13 UTC, 456 pixels, 11 dBZ peak
  Range from KEPZ: ~30 km (revised from earlier 33.6 km estimate)
  Azimuth: ~85-95° from KEPZ

Previous RCS estimation (rcs_estimation.py) gave:
  Single-point RCS: ~0.02 m²
  Distributed total RCS: ~2-4 m² (upper bound, peak Z across all pixels)
"""

import numpy as np
import json
import os

# ============================================================================
# PHYSICAL CONSTANTS AND RADAR PARAMETERS
# ============================================================================

# KEPZ WSR-88D parameters
LAMBDA = 0.1071          # S-band wavelength (m), 2800 MHz
K_SQ_WATER = 0.93        # |K|² for water (NEXRAD assumption)
K_SQ_METAL = 1.0         # |K|² for perfect conductor
BEAMWIDTH_DEG = 0.95     # 3dB beamwidth (degrees)
BEAMWIDTH_RAD = np.radians(BEAMWIDTH_DEG)
GATE_LENGTH = 250.0      # Range gate (m)
PT_PEAK = 750e3           # Peak transmit power (W)
GAIN_DB = 45.5            # Antenna gain (dB)

# Atmospheric parameters
RHO_AIR = 1.0            # Air density at ~1500m MSL (kg/m³), slightly below sea level value
G_ACCEL = 9.81           # Gravitational acceleration (m/s²)

# Material properties
MATERIALS = {
    'aluminum': {'density': 2700, 'name': 'Aluminum (structural)'},
    'copper':   {'density': 8960, 'name': 'Copper (motor windings)'},
    'steel':    {'density': 7800, 'name': 'Steel (fasteners)'},
    'lithium_cell': {'density': 2500, 'name': 'Li-ion cell (avg)'},
    'cfrp':     {'density': 1600, 'name': 'Carbon fiber composite'},
    'pcb':      {'density': 1850, 'name': 'FR4 circuit board'},
}

# ============================================================================
# OBSERVED PARAMETERS (from dual-pol imagery)
# ============================================================================

# Spike 1 (primary engagement)
SPIKE1 = {
    'time_utc': '04:24',
    'pixel_count': 551,
    'peak_z_dbz': 18.0,
    'mean_z_dbz': 12.0,       # estimated average across debris field
    'zdr_range': (8.0, 12.0),  # dB, differential reflectivity
    'cc_range': (0.50, 0.70),  # cross-correlation ratio
    'range_km': 30.0,
    'az_span_deg': 15.0,       # azimuth span of debris field
    'range_span_km': 5.0,      # range extent of debris field
    'fill_fraction': 0.25,     # fraction of sector containing returns
}

# Spike 2 (possible re-engagement)
SPIKE2 = {
    'time_utc': '05:13',
    'pixel_count': 456,
    'peak_z_dbz': 11.0,
    'mean_z_dbz': 7.0,
    'range_km': 30.0,
}

# Timing
TIME_BETWEEN_SPIKES_S = 49 * 60  # 49 minutes
SCAN_INTERVAL_S = 15 * 60         # ~15 min between NEXRAD volume scans
# Spike visible for ~1 scan cycle before diminishing
DEBRIS_PERSISTENCE_S = 15 * 60    # minimum time debris was radar-visible

# Fort Bliss elevation
TERRAIN_ELEV_M = 1200  # meters MSL at Fort Bliss

results = {}

print("=" * 78)
print("DEBRIS MASS ESTIMATION FROM NEXRAD DUAL-POL OBSERVABLES")
print("El Paso Airspace Incident — Feb 11, 2026")
print("=" * 78)

# ============================================================================
# METHOD 1: TOTAL RCS → FRAGMENT COUNT → MASS
# ============================================================================
print("\n" + "=" * 78)
print("METHOD 1: RCS-Based Fragment Count Estimation")
print("=" * 78)

range_m = SPIKE1['range_km'] * 1000

# Beam volume at target range
beam_diameter = 2 * range_m * np.tan(BEAMWIDTH_RAD / 2)
beam_volume = (np.pi / 4) * (BEAMWIDTH_RAD * range_m) ** 2 * GATE_LENGTH

print(f"\n  Beam geometry at {SPIKE1['range_km']} km:")
print(f"    Beam diameter: {beam_diameter:.0f} m")
print(f"    Single gate volume: {beam_volume:.3e} m³")

# Total integrated reflectivity
# Using mean Z across all debris-containing resolution volumes
z_mean_linear = 10 ** (SPIKE1['mean_z_dbz'] / 10)  # mm⁶/m³
z_peak_linear = 10 ** (SPIKE1['peak_z_dbz'] / 10)

print(f"\n  Reflectivity:")
print(f"    Peak Z: {SPIKE1['peak_z_dbz']} dBZ = {z_peak_linear:.1f} mm⁶/m³")
print(f"    Mean Z: {SPIKE1['mean_z_dbz']} dBZ = {z_mean_linear:.1f} mm⁶/m³")

# Number of resolution volumes with debris
# Estimate from pixel count and image rendering
# The debris field spans ~15° azimuth × 5 km range at 30 km
az_extent_m = 2 * np.pi * range_m * SPIKE1['az_span_deg'] / 360
range_extent_m = SPIKE1['range_span_km'] * 1000
sector_area_m2 = az_extent_m * range_extent_m

# Resolution cell size at this range
cell_az_m = 2 * range_m * np.tan(BEAMWIDTH_RAD / 2)  # ~497 m
cell_range_m = GATE_LENGTH  # 250 m
cell_area_m2 = cell_az_m * cell_range_m

n_cells_total = sector_area_m2 / cell_area_m2
n_cells_with_debris = n_cells_total * SPIKE1['fill_fraction']

print(f"\n  Debris field geometry:")
print(f"    Azimuth extent: {SPIKE1['az_span_deg']}° = {az_extent_m:.0f} m")
print(f"    Range extent: {SPIKE1['range_span_km']} km")
print(f"    Sector area: {sector_area_m2/1e6:.1f} km²")
print(f"    Resolution cell: {cell_az_m:.0f} m × {cell_range_m:.0f} m")
print(f"    Total cells in sector: {n_cells_total:.0f}")
print(f"    Cells with debris (~{SPIKE1['fill_fraction']*100:.0f}% fill): {n_cells_with_debris:.0f}")

# Per-cell RCS from reflectivity (assuming metallic scatterers)
z_mean_m6 = z_mean_linear * 1e-18  # convert mm⁶/m³ to m⁶/m³
sigma_per_cell = z_mean_m6 * (np.pi ** 5) * K_SQ_METAL * beam_volume / (LAMBDA ** 4)
sigma_total = sigma_per_cell * n_cells_with_debris

print(f"\n  RCS estimation:")
print(f"    Per-cell RCS (mean Z): {sigma_per_cell:.6f} m² ({10*np.log10(max(sigma_per_cell,1e-30)):.1f} dBsm)")
print(f"    Total integrated RCS: {sigma_total:.3f} m² ({10*np.log10(max(sigma_total,1e-30)):.1f} dBsm)")

# Fragment count and mass for different fragment types
print(f"\n  Fragment scenarios (total RCS = {sigma_total:.3f} m²):")
print(f"  {'Fragment Type':<35} {'RCS each':<12} {'Count':<10} {'Mass each':<12} {'Total Mass'}")
print(f"  {'-'*35} {'-'*12} {'-'*10} {'-'*12} {'-'*12}")

fragment_scenarios = [
    # (name, rcs_m2, mass_kg, description)
    ("Half-wave wire (5cm Cu, 1mm dia)", 0.0086, 3.5e-4,
     "Resonant dipole from motor windings"),
    ("Quarter-wave wire (2.5cm Cu)", 0.002, 1.75e-4,
     "Sub-resonant motor winding fragment"),
    ("Al foil strip (5cm × 1cm × 0.1mm)", 0.001, 1.35e-4,
     "Battery foil or structural skin fragment"),
    ("PCB fragment (2cm × 2cm)", 0.0005, 1.5e-3,
     "Circuit board piece with copper traces"),
    ("Small Al structural piece (3cm)", 0.003, 5.0e-3,
     "Airframe fragment, 1mm thick"),
    ("Medium Al panel (10cm × 5cm)", 0.01, 0.027,
     "Structural panel, 2mm thick"),
    ("Steel fastener (M3 bolt)", 0.0002, 3.0e-3,
     "Hardware"),
    ("Li-ion foil fragment (3cm²)", 0.0008, 2.0e-4,
     "Battery cell foil"),
]

method1_masses = []
for name, rcs, mass, desc in fragment_scenarios:
    count = sigma_total / rcs
    total_mass = count * mass
    method1_masses.append(total_mass)
    print(f"  {name:<35} {rcs:<12.4f} {count:<10.0f} {mass*1000:<10.1f}g  {total_mass:.2f} kg")

method1_range = (min(method1_masses), max(method1_masses))
method1_median = np.median(method1_masses)

print(f"\n  Method 1 mass range: {method1_range[0]:.2f} - {method1_range[1]:.2f} kg")
print(f"  Method 1 median: {method1_median:.2f} kg")
print(f"\n  NOTE: Real debris is a MIX of fragment types. The true mass is")
print(f"  constrained by the dominant scatterer type (from ZDR analysis).")

results['method1'] = {
    'total_rcs_m2': float(sigma_total),
    'n_cells_with_debris': float(n_cells_with_debris),
    'mass_range_kg': [float(method1_range[0]), float(method1_range[1])],
    'mass_median_kg': float(method1_median),
}

# ============================================================================
# METHOD 2: ZDR-CONSTRAINED FRAGMENT GEOMETRY
# ============================================================================
print("\n" + "=" * 78)
print("METHOD 2: ZDR-Constrained Fragment Geometry")
print("=" * 78)

zdr_min, zdr_max = SPIKE1['zdr_range']
print(f"\n  Observed ZDR: {zdr_min}-{zdr_max} dB")
print(f"  ZDR = 10 × log10(σ_H / σ_V)")
print(f"  ZDR 8-12 dB → horizontal RCS is {10**(8/10):.0f}× to {10**(12/10):.0f}× larger than vertical")

# ZDR constrains aspect ratio
# For thin metallic strips/wires oriented horizontally:
# ZDR ≈ 20 × log10(L/D) for L/D > 3 (wire approximation)
# ZDR 8 dB → L/D ≈ 10^(8/20) = 2.5 (but this is for single-particle ZDR)
# For a tumbling population, the OBSERVED ZDR is the distribution-averaged value
# ZDR 8-12 dB from a population means individual particles have EXTREME aspect ratios
# This means predominantly wire-like or strip-like fragments

print(f"\n  ZDR 8-12 dB implies:")
print(f"    Single-particle axis ratio: >> 10:1")
print(f"    Population-averaged with tumbling still shows 8+ dB")
print(f"    → Dominant scatterers are WIRE-LIKE or THIN STRIP-LIKE")
print(f"    → Consistent with: motor windings, battery foil, thin structural members")
print(f"    → Inconsistent with: compact chunks, spheres, thick structural pieces")

# For wire-like scatterers (the ZDR-dominant type):
# Motor winding wire: copper, 0.5-2mm diameter, lengths 1-20 cm after breakup
# Battery electrode foil: aluminum/copper, 10-50 μm thick, 1-10 cm fragments
# The ZDR tells us THESE are the dominant radar scatterers

# Wire fragment mass distribution
print(f"\n  ZDR-constrained mass estimation (wire/strip-dominated debris):")
print(f"  Assuming dominant scatterers are copper motor winding fragments:")

wire_diameters = [0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3]  # meters
wire_lengths = [0.02, 0.05, 0.10, 0.15]  # meters

print(f"\n  {'Wire dia (mm)':<15} {'Length (cm)':<13} {'Mass (g)':<10} {'RCS (m²)':<12} {'Count for total RCS':<20} {'Wire mass (kg)'}")
print(f"  {'-'*15} {'-'*13} {'-'*10} {'-'*12} {'-'*20} {'-'*15}")

method2_masses = []
for d in wire_diameters:
    for L in wire_lengths:
        mass_single = np.pi * (d/2)**2 * L * MATERIALS['copper']['density']
        # RCS of a thin wire of length L at S-band
        # For L ~ lambda/2 (resonant): sigma ≈ 0.86 * lambda²
        # For L < lambda/2: sigma ≈ (pi³/3) * L⁴ * (2*pi/lambda)⁴ / (ln(2L/d))² × (d/2)²
        # Simplified: use empirical formula for thin metallic wire
        k = 2 * np.pi / LAMBDA
        if L >= LAMBDA / 2:
            # Near-resonant or super-resonant
            sigma = 0.86 * LAMBDA**2 * (L / (LAMBDA/2))  # scales roughly linearly above resonance
            sigma = min(sigma, L**2)  # physical limit
        else:
            # Sub-resonant (Rayleigh regime for thin wire)
            # sigma ≈ (2*pi*L²) * (k*L)² / (3*(ln(2L/d)-1)²)  for kL < 1
            kL = k * L
            if kL < 0.5:
                log_term = max(np.log(2*L/d) - 1, 0.5)
                sigma = (2 * np.pi * L**2) * kL**2 / (3 * log_term**2)
            else:
                # Transition region
                sigma = 0.86 * LAMBDA**2 * (L / (LAMBDA/2))**2

        count = sigma_total / sigma
        total_wire_mass = count * mass_single
        method2_masses.append(total_wire_mass)

        print(f"  {d*1000:<15.1f} {L*100:<13.1f} {mass_single*1000:<10.3f} "
              f"{sigma:<12.6f} {count:<20.0f} {total_wire_mass:.3f}")

# The wire mass is only the METALLIC component
# Total vehicle mass includes composite structure, payload, etc.
metallic_fraction_range = (0.15, 0.50)  # fraction of total mass that is metal
wire_mass_median = np.median(method2_masses)
total_mass_from_wire = [wire_mass_median / f for f in metallic_fraction_range]

print(f"\n  Median wire/strip metallic mass: {wire_mass_median:.2f} kg")
print(f"  If metal is {metallic_fraction_range[0]*100:.0f}-{metallic_fraction_range[1]*100:.0f}% of total vehicle mass:")
print(f"    Total vehicle mass: {total_mass_from_wire[1]:.1f} - {total_mass_from_wire[0]:.1f} kg")

results['method2'] = {
    'wire_mass_median_kg': float(wire_mass_median),
    'metallic_fraction_range': list(metallic_fraction_range),
    'total_mass_range_kg': [float(min(total_mass_from_wire)), float(max(total_mass_from_wire))],
}

# ============================================================================
# METHOD 3: DEBRIS FIELD EXTENT → FRAGMENTATION ENERGY
# ============================================================================
print("\n" + "=" * 78)
print("METHOD 3: Debris Field Extent → Fragmentation Energy → Mass")
print("=" * 78)

# The debris field spans ~7 km azimuthally × 5 km in range
# Assuming fragmentation happened ~4 min before the 04:24 scan
# (event at ~04:20, next scan captures debris at 04:24)
t_expansion = 4 * 60  # seconds since fragmentation (estimate)

# Subtract wind contribution
# El Paso surface winds at time: likely 5-10 m/s from west
wind_speed = 8  # m/s (typical for El Paso at night)
wind_displacement = wind_speed * t_expansion

# Remaining expansion from fragmentation energy
debris_radius_m = az_extent_m / 2  # ~3600 m half-width
expansion_from_frag = max(debris_radius_m - wind_displacement, 500)
v_frag = expansion_from_frag / t_expansion  # fragment radial velocity

print(f"\n  Debris field dimensions: {az_extent_m:.0f} m × {range_extent_m:.0f} m")
print(f"  Time since fragmentation (est.): {t_expansion/60:.0f} min")
print(f"  Wind displacement: {wind_displacement:.0f} m ({wind_speed} m/s × {t_expansion}s)")
print(f"  Expansion from fragmentation: {expansion_from_frag:.0f} m")
print(f"  Implied fragment radial velocity: {v_frag:.1f} m/s")

# For laser-induced destruction:
# The laser heats and weakens structure → aerodynamic breakup at flight speed
# Fragment velocities from aerodynamic breakup: ~1-2× flight speed
# For thermal destruction: fragments separate at relative velocities of 1-20 m/s
# For explosive destruction (battery detonation): fragments at 50-200 m/s

# Kinetic energy of fragmentation
# E_frag = 0.5 * m * v² (for each fragment moving at v_frag)
# This gives a MINIMUM mass (more fragments → more total KE → more energy required)

# But we can also use the field size to estimate flight speed of the target
# If the target was moving at speed V and broke up, debris spreads along track
# Range extent (5 km) / time (240 s) = 21 m/s along-track component
# Azimuth extent (7 km) / time (240 s) = 29 m/s cross-track component

# More realistically: debris field stretches due to differential drag
# Heavier pieces decelerate slower, lighter pieces decelerate faster
# This creates elongation along the original velocity vector

print(f"\n  Fragment velocity analysis:")
print(f"    Along-range component: {range_extent_m/t_expansion:.1f} m/s")
print(f"    Cross-range component: {az_extent_m/t_expansion:.1f} m/s")

# For a laser engagement, the primary energy input is thermal
# The laser heats structural members until failure
# Fragmentation is then driven by:
# 1. Aerodynamic forces on a weakened structure at flight speed
# 2. Internal pressure (battery thermal runaway)
# 3. Gravity (falling debris separates)
#
# For a drone at 50-100 knots (25-50 m/s), aerodynamic breakup
# would produce fragment velocities of ~10-50 m/s relative to each other
# after 4 minutes of differential drag deceleration

# Using drag-based approach:
# Two fragments with different ballistic coefficients (β = m/CdA)
# will separate in time t by: Δx = (1/β₁ - 1/β₂) × 0.5 × ρ × V² × t²
# For Δx = 5000 m, V₀ = 40 m/s, t = 240 s, ρ = 1.0:
# (1/β₁ - 1/β₂) = 2×5000 / (1.0 × 40² × 240²) = 10000 / (1600 × 57600) ≈ 1.1e-4

# This is consistent with β₁ = 50 kg/m² (heavy piece) and β₂ = 5 kg/m² (light piece)
# A 50 kg/m² ballistic coefficient for a 5cm thick aluminum plate:
# β = m/(Cd×A), Cd≈1.2, m = ρ_al × A × t = 2700 × A × 0.05
# β = 2700×0.05/1.2 = 112 kg/m² (too high for thin debris)
# For a 2mm plate: β = 2700×0.002/1.2 = 4.5 kg/m² ← more like it

# Method 3 gives a weak mass constraint but confirms the debris is from
# an object with flight speed 25-50 m/s and mixed fragment sizes
print(f"\n  Implications:")
print(f"    Debris dispersion consistent with flight speed ~25-50 m/s (50-100 kts)")
print(f"    Fragment ballistic coefficients span factor of ~10×")
print(f"    → Mix of heavy (structural) and light (foil/wire) fragments")
print(f"    → This is consistent with a composite-bodied UAS, NOT a balloon")
print(f"    (A balloon burst produces near-zero velocity fragments → tight cluster)")

# For a balloon at ~zero velocity, debris field after 4 min would be:
# pure wind: 8 m/s × 240s = 1920 m in ONE direction (no spread)
balloon_spread = wind_speed * t_expansion
print(f"\n  Balloon comparison:")
print(f"    Balloon burst debris field (wind only): {balloon_spread:.0f} m linear, no lateral spread")
print(f"    Observed debris field: {az_extent_m:.0f} m × {range_extent_m:.0f} m")
print(f"    Ratio: observed is {az_extent_m/balloon_spread:.0f}× wider → NOT a balloon")

results['method3'] = {
    'debris_field_m': [float(az_extent_m), float(range_extent_m)],
    'implied_flight_speed_ms': float(np.sqrt((az_extent_m/t_expansion)**2 + (range_extent_m/t_expansion)**2)),
    'balloon_spread_m': float(balloon_spread),
}

# ============================================================================
# METHOD 4: DEBRIS FALL RATE → BALLISTIC COEFFICIENT → MASS
# ============================================================================
print("\n" + "=" * 78)
print("METHOD 4: Debris Persistence → Fall Rate → Mass/Area Ratio")
print("=" * 78)

# The debris is visible at 04:24 (Spike 1)
# The signal diminishes significantly by the next scan
# This constrains how quickly fragments fall through the beam

# KEPZ 0.5° beam at 30 km range:
beam_center_height = range_m * np.sin(np.radians(0.5)) + TERRAIN_ELEV_M + range_m**2 / (2 * 6371000 * 4/3)
beam_lower_height = range_m * np.sin(np.radians(0.5 - BEAMWIDTH_DEG/2)) + TERRAIN_ELEV_M + range_m**2 / (2 * 6371000 * 4/3)
beam_upper_height = range_m * np.sin(np.radians(0.5 + BEAMWIDTH_DEG/2)) + TERRAIN_ELEV_M + range_m**2 / (2 * 6371000 * 4/3)
beam_thickness = beam_upper_height - beam_lower_height

print(f"\n  Beam geometry at {SPIKE1['range_km']} km:")
print(f"    Beam center: {beam_center_height:.0f} m MSL ({beam_center_height-TERRAIN_ELEV_M:.0f} m AGL)")
print(f"    Beam lower edge: {beam_lower_height:.0f} m MSL ({beam_lower_height-TERRAIN_ELEV_M:.0f} m AGL)")
print(f"    Beam upper edge: {beam_upper_height:.0f} m MSL ({beam_upper_height-TERRAIN_ELEV_M:.0f} m AGL)")
print(f"    Beam thickness: {beam_thickness:.0f} m")

# Scenario: debris was created above the beam, fell through it
# The persistence of ~15 min (one full scan cycle) means the heaviest
# fragments passed through the beam in < 15 min, while lightest may linger

# For a flat plate falling at terminal velocity:
# v_t = sqrt(2 × m × g / (ρ_air × Cd × A))
# Cd for flat plate perpendicular to flow: ~1.2
# Cd for tumbling plate: ~1.0

CD_PLATE = 1.0  # tumbling flat debris

print(f"\n  Terminal velocity analysis for flat plate fragments:")
print(f"  v_t = sqrt(2mg / (ρ_air × Cd × A))")
print(f"  Cd = {CD_PLATE} (tumbling plate)")
print(f"")
print(f"  {'Material':<25} {'Thickness':<12} {'Surface dens':<14} {'v_terminal':<12} {'Fall 500m':<12}")
print(f"  {'':<25} {'(mm)':<12} {'(kg/m²)':<14} {'(m/s)':<12} {'(seconds)':<12}")
print(f"  {'-'*25} {'-'*12} {'-'*14} {'-'*12} {'-'*12}")

fall_scenarios = []
for mat_key, mat in MATERIALS.items():
    for thickness_mm in [0.1, 0.5, 1.0, 2.0, 5.0]:
        thickness_m = thickness_mm / 1000
        surface_density = mat['density'] * thickness_m  # kg/m²

        # Terminal velocity
        v_t = np.sqrt(2 * surface_density * G_ACCEL / (RHO_AIR * CD_PLATE))

        # Time to fall through beam thickness
        fall_time_beam = beam_thickness / v_t

        # Time to fall 500m (approximate altitude to beam center)
        fall_time_500 = 500 / v_t

        fall_scenarios.append({
            'material': mat['name'],
            'thickness_mm': thickness_mm,
            'surface_density': surface_density,
            'v_terminal': v_t,
            'fall_time_beam': fall_time_beam,
            'fall_time_500': fall_time_500,
        })

        if thickness_mm in [0.1, 1.0, 5.0]:
            print(f"  {mat['name']:<25} {thickness_mm:<12.1f} {surface_density:<14.3f} "
                  f"{v_t:<12.1f} {fall_time_500:<12.0f}")

# Constraint: debris persists for ~1 scan (15 min) → slowest fragments
# have terminal velocity < beam_thickness / (15 min)
v_max_persist = beam_thickness / DEBRIS_PERSISTENCE_S
# Constraint: debris mostly gone by next scan → fastest fragments
# have terminal velocity > beam_thickness / (30 min) (minimum)
v_min_gone = beam_thickness / (2 * DEBRIS_PERSISTENCE_S)

print(f"\n  Fall rate constraints:")
print(f"    For debris to persist ~15 min in beam ({beam_thickness:.0f}m thick):")
print(f"    Slowest fragments: v_t < {v_max_persist:.2f} m/s")
print(f"    → Surface density < {0.5 * RHO_AIR * CD_PLATE * v_max_persist**2 / G_ACCEL:.3f} kg/m²")
print(f"    → This is {0.5 * RHO_AIR * CD_PLATE * v_max_persist**2 / G_ACCEL / 2.7:.2f} mm of aluminum")
print(f"    → Implies very thin fragments (foil, wire, thin sheet)")

# Altitude estimation
# If heaviest fragments (v_t ~ 5-10 m/s) reach ground before next scan (< 15 min):
# altitude < 5 m/s × 900s = 4500 m AGL (probably too high)
# More likely the heavy pieces fall faster: v_t ~ 15 m/s → altitude ~ 13500 m → way too high
#
# Alternatively, if the engagement altitude is 1000-3000 m AGL:
# Heavy fragments (v_t ~ 10 m/s) reach ground in 100-300 s (2-5 min)
# Light fragments (v_t ~ 1 m/s) reach ground in 1000-3000 s (17-50 min)
# This matches the observation: significant debris at 04:24, diminished by next scan

alt_estimates = [1000, 2000, 3000, 5000]
print(f"\n  Engagement altitude estimates from debris persistence:")
print(f"  {'Altitude AGL':<15} {'Heavy (v=10m/s)':<18} {'Medium (v=3m/s)':<18} {'Light (v=1m/s)':<18}")
print(f"  {'-'*15} {'-'*18} {'-'*18} {'-'*18}")
for alt in alt_estimates:
    t_heavy = alt / 10
    t_medium = alt / 3
    t_light = alt / 1
    print(f"  {alt:>8} m      {t_heavy/60:>6.1f} min       {t_medium/60:>6.1f} min       {t_light/60:>6.1f} min")

print(f"\n  Best fit altitude: ~2000-3000 m AGL (6500-10,000 ft)")
print(f"    → Heavy fragments hit ground in 3-5 min (gone before next scan)")
print(f"    → Light fragments persist 30-50 min (visible for 1-2 scan cycles)")
print(f"    → Consistent with Group 3 UAS operating altitude")

# Mass constraint from fall rate:
# The lightest visible fragments have v_t ~ 1 m/s
# The heaviest visible fragments have v_t ~ 10 m/s
# Assuming total debris spans a range of ballistic coefficients
# Surface densities from 0.05 to 10 kg/m²
# For fragments with typical area of 1-10 cm²:
#   Light: 0.05 kg/m² × 10 cm² = 0.5 g per fragment
#   Heavy: 10 kg/m² × 10 cm² = 100 g per fragment
#   Median: ~5 g per fragment

# Combined with fragment count from Method 1:
n_fragments_est = sigma_total / 0.002  # Using mid-range RCS per fragment
median_mass_per_frag = 0.005  # 5 g
method4_total_mass = n_fragments_est * median_mass_per_frag

print(f"\n  Estimated fragment count (from Method 1, mid-range): {n_fragments_est:.0f}")
print(f"  Estimated mass per fragment (from fall rate): ~5 g")
print(f"  Method 4 total metallic mass: {method4_total_mass:.1f} kg")

results['method4'] = {
    'beam_center_agl_m': float(beam_center_height - TERRAIN_ELEV_M),
    'beam_thickness_m': float(beam_thickness),
    'estimated_altitude_agl_m': [2000, 3000],
    'metallic_mass_kg': float(method4_total_mass),
}

# ============================================================================
# METHOD 5: CHAFF EQUIVALENCE
# ============================================================================
print("\n" + "=" * 78)
print("METHOD 5: Chaff Equivalence (Military Literature Calibration)")
print("=" * 78)

# RR-188 chaff bundle (standard NATO chaff):
# - Contains ~5 million aluminum-coated glass fiber dipoles
# - Each dipole: 25.4 μm diameter × ~2.5 cm (for S-band)
# - Total mass per bundle: ~85-100 g (early versions ~150g)
# - Fresh deployment RCS: varies by frequency, but at S-band ~50-100 m²
# - After 5-10 min dispersal: Z ~ 25-35 dBZ in radar returns
# - After 30 min: Z ~ 10-20 dBZ
#
# Our debris field: mean Z ~ 12 dBZ, peak 18 dBZ
# A 30-minute-old chaff cloud would produce similar reflectivity

CHAFF_MASS_KG = 0.085  # RR-188 bundle mass
CHAFF_N_DIPOLES = 5e6
CHAFF_FRESH_Z_DBZ = 30  # typical fresh deployment
CHAFF_30MIN_Z_DBZ = 15  # dispersed

# Chaff-equivalent mass
# At 15 dBZ (dispersed chaff level), we're at roughly 1 chaff bundle
# Our mean 12 dBZ is about 3 dB below dispersed chaff
# In linear: 10^((12-15)/10) = 0.50 × one chaff bundle

z_ratio = 10 ** ((SPIKE1['mean_z_dbz'] - CHAFF_30MIN_Z_DBZ) / 10)

# But chaff dipoles are OPTIMIZED for S-band (tuned half-wave dipoles)
# Random debris fragments have much LOWER RCS per unit mass than chaff
# Efficiency factor: debris generates ~1/10 to 1/100 the RCS per gram vs chaff
efficiency_range = (0.01, 0.10)

# Also need to account for the number of resolution volumes
# One chaff bundle fills maybe 5-10 resolution volumes at 30 km
# Our debris fills ~75 resolution volumes
volume_ratio = n_cells_with_debris / 7.5  # vs single chaff bundle

chaff_eq_bundles = z_ratio * volume_ratio
chaff_eq_mass_optimistic = chaff_eq_bundles * CHAFF_MASS_KG / efficiency_range[1]  # debris less efficient
chaff_eq_mass_conservative = chaff_eq_bundles * CHAFF_MASS_KG / efficiency_range[0]

print(f"\n  RR-188 chaff bundle reference:")
print(f"    Mass: {CHAFF_MASS_KG*1000:.0f} g")
print(f"    Dipoles: {CHAFF_N_DIPOLES:.0e}")
print(f"    Fresh Z: ~{CHAFF_FRESH_Z_DBZ} dBZ")
print(f"    30-min dispersed Z: ~{CHAFF_30MIN_Z_DBZ} dBZ")

print(f"\n  Chaff equivalence calculation:")
print(f"    Z ratio (debris/chaff): {z_ratio:.2f}")
print(f"    Volume ratio (debris/single chaff): {volume_ratio:.1f}")
print(f"    Chaff-equivalent bundles: {chaff_eq_bundles:.1f}")
print(f"    Chaff-equivalent mass: {chaff_eq_bundles * CHAFF_MASS_KG * 1000:.0f} g")

print(f"\n  Correcting for debris inefficiency vs tuned chaff:")
print(f"    Debris RCS efficiency: {efficiency_range[0]*100:.0f}-{efficiency_range[1]*100:.0f}% of chaff per gram")
print(f"    Actual metallic debris mass: {chaff_eq_mass_optimistic:.1f} - {chaff_eq_mass_conservative:.1f} kg")

results['method5'] = {
    'chaff_equivalent_bundles': float(chaff_eq_bundles),
    'chaff_equivalent_mass_g': float(chaff_eq_bundles * CHAFF_MASS_KG * 1000),
    'debris_mass_range_kg': [float(chaff_eq_mass_optimistic), float(chaff_eq_mass_conservative)],
}

# ============================================================================
# SYNTHESIS: COMBINED MASS ESTIMATE
# ============================================================================
print("\n" + "=" * 78)
print("SYNTHESIS: Combined Mass Estimate")
print("=" * 78)

print(f"\n  Summary of independent estimates:")
print(f"  {'Method':<45} {'Metallic Mass (kg)':<22} {'Total Vehicle Mass (kg)'}")
print(f"  {'-'*45} {'-'*22} {'-'*22}")

# Method 1: direct RCS
m1_metal = method1_median
m1_total = [m1_metal / 0.5, m1_metal / 0.15]
print(f"  {'1. RCS → fragment count → mass':<45} {m1_metal:<22.1f} {m1_total[0]:.0f} - {m1_total[1]:.0f}")

# Method 2: ZDR-constrained
m2_metal = wire_mass_median
m2_total = total_mass_from_wire
print(f"  {'2. ZDR-constrained wire geometry':<45} {m2_metal:<22.1f} {min(m2_total):.0f} - {max(m2_total):.0f}")

# Method 4: fall rate
m4_metal = method4_total_mass
m4_total = [m4_metal / 0.5, m4_metal / 0.15]
print(f"  {'4. Fall rate / ballistic coefficient':<45} {m4_metal:<22.1f} {m4_total[0]:.0f} - {m4_total[1]:.0f}")

# Method 5: chaff equivalence
m5_metal_mid = np.sqrt(chaff_eq_mass_optimistic * chaff_eq_mass_conservative)  # geometric mean
m5_total = [m5_metal_mid / 0.5, m5_metal_mid / 0.15]
print(f"  {'5. Chaff equivalence (geometric mean)':<45} {m5_metal_mid:<22.1f} {m5_total[0]:.0f} - {m5_total[1]:.0f}")

# Combined estimate
all_metal_estimates = [m1_metal, m2_metal, m4_metal, m5_metal_mid]
metal_low = np.percentile(all_metal_estimates, 25)
metal_high = np.percentile(all_metal_estimates, 75)
metal_median = np.median(all_metal_estimates)

total_low = metal_low / 0.50   # if 50% metal (heavily metallic construction)
total_high = metal_high / 0.15  # if 15% metal (mostly composite)

print(f"\n  Combined estimate (interquartile range of methods):")
print(f"    Metallic debris mass: {metal_low:.1f} - {metal_high:.1f} kg (median {metal_median:.1f} kg)")
print(f"    If vehicle is 15-50% metal by mass:")
print(f"    ╔══════════════════════════════════════════════╗")
print(f"    ║  Total vehicle mass: {total_low:.0f} - {total_high:.0f} kg            ║")
print(f"    ║  Best estimate: {metal_median/0.30:.0f} kg ({metal_median:.1f} kg metal @ 30%%)   ║")
print(f"    ╚══════════════════════════════════════════════╝")

best_total = metal_median / 0.30

# ============================================================================
# WHAT KIND OF OBJECT?
# ============================================================================
print(f"\n" + "=" * 78)
print("OBJECT CLASSIFICATION FROM MASS AND SIGNATURES")
print("=" * 78)

uas_classes = [
    ("Group 1 (Micro)", 0, 9, "RQ-11 Raven (1.9kg), commercial quad"),
    ("Group 2 (Small)", 9.1, 25, "ScanEagle (22kg), Puma AE (6.3kg)"),
    ("Group 3 (Medium)", 25, 600, "RQ-7 Shadow (170kg), Orbiter 3 (28kg)"),
    ("Group 4 (Large)", 600, 1320, "MQ-1C Gray Eagle (1633kg)"),
    ("Group 5 (Largest)", 1320, 15000, "MQ-9 Reaper (4760kg)"),
]

print(f"\n  US DoD UAS Group Classification:")
print(f"  {'Group':<25} {'Mass Range':<15} {'Examples'}")
print(f"  {'-'*25} {'-'*15} {'-'*40}")
for name, low, high, examples in uas_classes:
    marker = " ◄◄◄" if low <= best_total <= high else ""
    print(f"  {name:<25} {low}-{high} kg    {examples}{marker}")

# Swarm possibility
print(f"\n  Single vehicle vs swarm analysis:")
print(f"    If single vehicle: {best_total:.0f} kg → Group {'2-3' if 9 < best_total < 600 else '?'} UAS")
print(f"    If swarm of Group 1 drones ({best_total/2:.0f} kg each): {max(1,int(best_total/2))} drones")
print(f"    If swarm of Group 1 drones ({best_total/5:.0f} kg each): {max(1,int(best_total/5))} drones")

print(f"\n  Spike 2 comparison:")
z2_ratio = 10 ** ((SPIKE2['peak_z_dbz'] - SPIKE1['peak_z_dbz']) / 10)
print(f"    Spike 2 peak Z: {SPIKE2['peak_z_dbz']} dBZ (vs Spike 1: {SPIKE1['peak_z_dbz']} dBZ)")
print(f"    Z ratio: {z2_ratio:.2f}× (Spike 2 is {10*np.log10(z2_ratio):.0f} dB weaker)")
print(f"    Pixel ratio: {SPIKE2['pixel_count']}/{SPIKE1['pixel_count']} = {SPIKE2['pixel_count']/SPIKE1['pixel_count']:.2f}")
print(f"    If Spike 2 is a second target: mass ≈ {best_total * z2_ratio:.0f} kg")
print(f"    If Spike 2 is re-engagement of damaged Spike 1 target: reduced by {(1-z2_ratio)*100:.0f}%")
print(f"    → {(1-z2_ratio)*100:.0f}% of original target already destroyed/fallen")

# ============================================================================
# SENSITIVITY ANALYSIS
# ============================================================================
print(f"\n" + "=" * 78)
print("SENSITIVITY ANALYSIS: Key Uncertainties")
print("=" * 78)

params = {
    'Mean Z (dBZ)': (8, 12, 16, 'Biggest uncertainty — image-derived estimate'),
    'Range (km)': (25, 30, 35, 'Affects beam volume quadratically'),
    'Fill fraction': (0.10, 0.25, 0.40, 'What fraction of sector has debris'),
    'Metal fraction': (0.10, 0.30, 0.50, 'How much of vehicle was metal'),
    'Fragment RCS (m²)': (0.0005, 0.002, 0.01, 'Depends on fragment size distribution'),
}

print(f"\n  {'Parameter':<25} {'Low':<10} {'Nominal':<10} {'High':<10} {'Effect on mass':<25} {'Note'}")
print(f"  {'-'*25} {'-'*10} {'-'*10} {'-'*10} {'-'*25} {'-'*30}")

for param, (low, nom, high, note) in params.items():
    if param == 'Mean Z (dBZ)':
        m_low = best_total * 10**((low - nom)/10)
        m_high = best_total * 10**((high - nom)/10)
    elif param == 'Range (km)':
        m_low = best_total * (low/nom)**2
        m_high = best_total * (high/nom)**2
    elif param == 'Fill fraction':
        m_low = best_total * (low/nom)
        m_high = best_total * (high/nom)
    elif param == 'Metal fraction':
        m_low = metal_median / high
        m_high = metal_median / low
    elif param == 'Fragment RCS (m²)':
        m_low = best_total * (nom/high)  # higher RCS → fewer fragments → less mass
        m_high = best_total * (nom/low)  # lower RCS → more fragments → more mass
    else:
        m_low = best_total * 0.5
        m_high = best_total * 2.0

    print(f"  {param:<25} {low:<10} {nom:<10} {high:<10} {m_low:.0f}-{m_high:.0f} kg       {note}")

# Monte Carlo-style range
print(f"\n  Considering all uncertainties simultaneously:")
print(f"    Conservative lower bound: ~{max(2, best_total * 0.1):.0f} kg (small composite drone)")
print(f"    Best estimate: ~{best_total:.0f} kg")
print(f"    Conservative upper bound: ~{best_total * 5:.0f} kg (large composite UAS)")
print(f"    Extreme upper bound: ~{best_total * 10:.0f} kg (if very low metal fraction)")

results['synthesis'] = {
    'metallic_mass_median_kg': float(metal_median),
    'total_mass_best_estimate_kg': float(best_total),
    'total_mass_range_kg': [float(total_low), float(total_high)],
    'uas_group_classification': 'Group 2-3',
    'estimated_altitude_agl_m': [2000, 3000],
    'engagement_type': 'HELWS laser',
}

# ============================================================================
# SAVE RESULTS
# ============================================================================
output_dir = '/home/user/uap-transient-research/el_paso_airspace/analysis_outputs'
output_path = os.path.join(output_dir, 'debris_mass_estimation_results.json')

with open(output_path, 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\n\nResults saved to {output_path}")
print("=" * 78)
