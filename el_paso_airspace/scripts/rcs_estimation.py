#!/usr/bin/env python3
"""
Reflectivity-to-RCS Estimation for El Paso Airspace Incident
Estimates radar cross section from KEPZ NEXRAD reflectivity values.

Key spike parameters from previous analysis:
  Spike 1: 21.5 dBZ, 196 pixels, range ~33.6 km, ZDR 18.9 dB, CC 0.50
  Spike 2: 20.0 dBZ, 144 pixels, range ~33.6 km
"""

import numpy as np
import json

# KEPZ Radar Parameters
LAMBDA = 0.1071  # S-band wavelength in meters (2800 MHz)
K_SQ_WATER = 0.93  # |K|^2 for water
K_SQ_METAL = 1.0   # |K|^2 for perfect conductor (approximation)
BEAMWIDTH = 0.95   # degrees (NEXRAD standard)
GATE_LENGTH = 250.0  # meters (standard NEXRAD gate)

# Spike parameters
RANGE_KM = 33.6
RANGE_M = RANGE_KM * 1000.0
Z_DBZ_SPIKE1 = 21.5
Z_DBZ_SPIKE2 = 20.0
N_PIXELS_SPIKE1 = 196
N_PIXELS_SPIKE2 = 144

results = {}

print("=" * 70)
print("REFLECTIVITY-TO-RCS ESTIMATION")
print("El Paso Airspace Incident - Feb 11, 2026")
print("=" * 70)

# ============= BEAM VOLUME CALCULATION =============
print("\n1. BEAM VOLUME AT TARGET RANGE")
print("-" * 40)

beamwidth_rad = np.radians(BEAMWIDTH)
# Beam diameter at range
beam_diameter = 2 * RANGE_M * np.tan(beamwidth_rad / 2)
print(f"  Range: {RANGE_KM} km")
print(f"  Beam diameter at range: {beam_diameter:.1f} m")

# Cross-sectional area of beam
beam_area = np.pi * (beam_diameter / 2) ** 2
print(f"  Beam cross-section area: {beam_area:.0f} m²")

# Volume of single resolution volume (gate)
# V = pi/4 * (beamwidth_rad * range)^2 * gate_length
beam_volume = (np.pi / 4) * (beamwidth_rad * RANGE_M) ** 2 * GATE_LENGTH
print(f"  Single gate volume: {beam_volume:.3e} m³")
print(f"  Gate length: {GATE_LENGTH} m")

# ============= POINT TARGET RCS FROM Z =============
print("\n2. POINT TARGET RCS ESTIMATION")
print("-" * 40)

# For a point target in the beam center:
# Z = (sigma * lambda^4) / (pi^5 * |K|^2 * V_eff)
# Where V_eff is the effective resolution volume
# Rearranging: sigma = Z * pi^5 * |K|^2 * V_eff / lambda^4
#
# But Z is expressed in mm^6/m^3, so:
# Z_linear = 10^(Z_dBZ/10) in mm^6/m^3
# Convert to m^6/m^3: Z_m = Z_linear * 1e-18

# Single pixel case
for spike_name, z_dbz, n_pixels in [("Spike 1", Z_DBZ_SPIKE1, N_PIXELS_SPIKE1),
                                      ("Spike 2", Z_DBZ_SPIKE2, N_PIXELS_SPIKE2)]:
    print(f"\n  {spike_name}: {z_dbz} dBZ, {n_pixels} pixels")

    Z_linear = 10 ** (z_dbz / 10)  # mm^6/m^3
    Z_m6 = Z_linear * 1e-18  # m^6/m^3 = m^3

    print(f"  Z_linear = {Z_linear:.2f} mm⁶/m³")

    # Case A: Single point target in one gate
    # The full Z is concentrated in one scatterer
    # sigma = Z * V * pi^5 * |K|^2 / lambda^4
    # Actually for point target in NEXRAD:
    # Z = (sigma * lambda^4) / (pi^5 * |K|^2 * V)
    # sigma = Z * pi^5 * |K|^2 * V / lambda^4

    sigma_single_water = Z_m6 * (np.pi ** 5) * K_SQ_WATER * beam_volume / (LAMBDA ** 4)
    sigma_single_metal = Z_m6 * (np.pi ** 5) * K_SQ_METAL * beam_volume / (LAMBDA ** 4)

    print(f"\n  Case A: Single point target in one gate")
    print(f"    RCS (|K|²=water): {sigma_single_water:.4f} m² ({10*np.log10(max(sigma_single_water,1e-20)):.1f} dBsm)")
    print(f"    RCS (|K|²=metal): {sigma_single_metal:.4f} m² ({10*np.log10(max(sigma_single_metal,1e-20)):.1f} dBsm)")

    # Case B: Signal distributed across N pixels
    # If the 21.5 dBZ is the PEAK value and it extends across N pixels,
    # the total reflectivity is roughly N * Z_peak (assuming uniform)
    # Total RCS = N * sigma_single
    sigma_total_water = n_pixels * sigma_single_water
    sigma_total_metal = n_pixels * sigma_single_metal

    print(f"\n  Case B: Distributed across {n_pixels} pixels (total)")
    print(f"    Total RCS (water): {sigma_total_water:.2f} m² ({10*np.log10(max(sigma_total_water,1e-20)):.1f} dBsm)")
    print(f"    Total RCS (metal): {sigma_total_metal:.2f} m² ({10*np.log10(max(sigma_total_metal,1e-20)):.1f} dBsm)")

    # Case C: The 21.5 dBZ is summed across the cloud (conservative)
    # This means the per-pixel Z is much lower
    z_per_pixel_dbz = z_dbz - 10 * np.log10(n_pixels)
    z_per_pixel_linear = 10 ** (z_per_pixel_dbz / 10)
    z_per_pixel_m6 = z_per_pixel_linear * 1e-18
    sigma_perpixel = z_per_pixel_m6 * (np.pi ** 5) * K_SQ_METAL * beam_volume / (LAMBDA ** 4)

    print(f"\n  Case C: {z_dbz} dBZ total distributed over {n_pixels} gates")
    print(f"    Per-pixel: {z_per_pixel_dbz:.1f} dBZ")
    print(f"    Total RCS: {sigma_perpixel * n_pixels:.4f} m² = {sigma_perpixel:.6f} m² per pixel")

    results[spike_name] = {
        'z_dbz': z_dbz,
        'n_pixels': n_pixels,
        'rcs_single_point_m2': float(sigma_single_metal),
        'rcs_distributed_total_m2': float(sigma_total_metal),
        'rcs_conservative_total_m2': float(sigma_perpixel * n_pixels),
    }

# ============= REFERENCE RCS VALUES =============
print("\n" + "=" * 70)
print("3. REFERENCE RCS VALUES FOR COMPARISON")
print("-" * 40)

reference_objects = {
    "Mylar party balloon (30cm diameter)": (0.001, 0.01),
    "Mylar birthday balloon (45cm)": (0.005, 0.05),
    "Weather balloon (1-2m, latex)": (0.01, 0.1),
    "Bird (sparrow to goose)": (0.0001, 0.01),
    "Small consumer drone (DJI Phantom)": (0.01, 0.1),
    "Medium drone (DJI Matrice 300)": (0.05, 0.5),
    "Large military drone (MQ-1 Predator)": (0.5, 2.0),
    "Large military drone (MQ-9 Reaper)": (1.0, 5.0),
    "Chaff bundle (fresh)": (10.0, 1000.0),
    "Chaff cloud (dispersed)": (1.0, 100.0),
    "Light aircraft (Cessna 172)": (1.0, 10.0),
    "Military helicopter (UH-60)": (5.0, 20.0),
    "Fighter aircraft (F-16)": (1.0, 5.0),
    "Metal debris cloud (1m² equivalent)": (0.5, 5.0),
}

print(f"  {'Object':<45} {'RCS Range (m²)':<20} {'dBsm Range'}")
print(f"  {'-'*45} {'-'*20} {'-'*20}")
for obj, (rcs_min, rcs_max) in reference_objects.items():
    dbsm_min = 10 * np.log10(rcs_min)
    dbsm_max = 10 * np.log10(rcs_max)
    print(f"  {obj:<45} {rcs_min:.4f}-{rcs_max:.1f}    {dbsm_min:.1f} to {dbsm_max:.1f}")

# ============= INTERPRETATION =============
print("\n" + "=" * 70)
print("4. INTERPRETATION")
print("-" * 40)

spike1_rcs = results['Spike 1']['rcs_single_point_m2']
print(f"\n  Spike 1 estimated RCS (single point target): {spike1_rcs:.4f} m²")
print(f"  ({10*np.log10(spike1_rcs):.1f} dBsm)")
print()

if spike1_rcs < 0.01:
    print("  -> Consistent with: small bird, insect swarm")
elif spike1_rcs < 0.1:
    print("  -> Consistent with: mylar balloon, small drone, large bird")
elif spike1_rcs < 1.0:
    print("  -> Consistent with: medium drone, weather balloon, small debris")
elif spike1_rcs < 10.0:
    print("  -> Consistent with: large drone, light aircraft, chaff cloud")
elif spike1_rcs < 100.0:
    print("  -> Consistent with: aircraft, large chaff cloud")
else:
    print("  -> Consistent with: very large target or chaff cloud")

print(f"\n  Spike 1 distributed over {N_PIXELS_SPIKE1} pixels total RCS: "
      f"{results['Spike 1']['rcs_distributed_total_m2']:.2f} m²")
print(f"  ({10*np.log10(results['Spike 1']['rcs_distributed_total_m2']):.1f} dBsm)")
print()

distributed_rcs = results['Spike 1']['rcs_distributed_total_m2']
if distributed_rcs > 10:
    print("  -> Large distributed return strongly suggests chaff, debris cloud,")
    print("     or metallic fragmentation rather than a single compact target.")
    print("     A single mylar balloon would NOT produce 196 pixels of return.")
    print("     Even a large drone would not spread across this many gates.")

print("\n  Key finding:")
print("  The 196-pixel spread at 21.5 dBZ is inconsistent with a single")
print("  mylar balloon or even a single drone. The spatial extent suggests")
print("  either a debris/chaff cloud or a VERY strong point target with")
print("  significant sidelobe contamination.")

# ZDR analysis
print("\n" + "=" * 70)
print("5. ZDR SIGNATURE ANALYSIS")
print("-" * 40)
print(f"  Measured ZDR: 18.9 dB (Spike 1)")
print(f"  Correlation coefficient: 0.50")
print()
print("  ZDR = 18.9 dB is EXTREMELY high. Typical values:")
print("    Rain: 0-4 dB")
print("    Hail: -2 to +2 dB")
print("    Biological: 2-10 dB")
print("    Chaff: >10 dB (metallic dipoles, highly elongated)")
print("    Ground clutter: variable, can be >10 dB")
print("    Metallic debris: >10 dB (irregular shapes)")
print()
print("  ZDR 18.9 dB indicates:")
print("    - Targets are MUCH larger in horizontal than vertical dimension")
print("    - Consistent with flat metallic surfaces (mylar, foil, chaff)")
print("    - NOT consistent with spherical targets (raindrops, round balloons)")
print("    - Could be: flat debris fragments, chaff dipoles, or tumbling foil")
print()
print("  CC = 0.50 indicates:")
print("    - Poor correlation between H and V polarization returns")
print("    - Targets have diverse shapes/orientations within beam volume")
print("    - Consistent with debris cloud or chaff (many scatterers)")
print("    - NOT consistent with single compact target (would have CC > 0.9)")
print("    - NOT consistent with meteorological targets (CC > 0.95)")
print()
print("  COMBINED ZDR + CC DIAGNOSIS:")
print("    High ZDR + Low CC = METALLIC DEBRIS/CHAFF CLOUD")
print("    This is the textbook dual-pol signature of metallic chaff or")
print("    a cloud of flat metallic fragments.")

# Save results
with open('/home/user/uap-transient-research/el_paso_airspace/analysis_outputs/rcs_estimation_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n\nResults saved to rcs_estimation_results.json")
