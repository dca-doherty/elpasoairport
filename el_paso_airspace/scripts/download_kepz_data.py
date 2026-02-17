#!/usr/bin/env python3
"""
Download KEPZ Level 2 NEXRAD data from UCAR THREDDS for the incident window.
Feb 11, 2026, 03:00 - 07:00 UTC
"""

import urllib.request
import xml.etree.ElementTree as ET
import os
import time
import re

THREDDS_CATALOG = "https://thredds.ucar.edu/thredds/catalog/nexrad/level2/KEPZ/20260211/catalog.xml"
THREDDS_BASE = "https://thredds.ucar.edu/thredds/fileServer/nexrad/level2/KEPZ/20260211/"
OUTPUT_DIR = "/home/user/uap-transient-research/el_paso_airspace/kepz_data/"

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Fetching THREDDS catalog...")
try:
    req = urllib.request.Request(THREDDS_CATALOG, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        catalog_xml = resp.read().decode()
except Exception as e:
    print(f"Error fetching catalog: {e}")
    raise

# Parse XML to get filenames
# THREDDS catalog uses namespace
ns = {'thredds': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'}
root = ET.fromstring(catalog_xml)

# Find all dataset elements
all_files = []
for dataset in root.iter('{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset'):
    name = dataset.get('name', '')
    if name.startswith('Level2_KEPZ_20260211_') and name.endswith('.ar2v'):
        all_files.append(name)

all_files.sort()
print(f"Found {len(all_files)} total KEPZ files for Feb 11, 2026")

# Filter to incident window: 0300 - 0700 UTC
# Filename format: Level2_KEPZ_20260211_HHMM.ar2v
incident_files = []
for f in all_files:
    match = re.search(r'_(\d{4})\.ar2v', f)
    if match:
        hhmm = match.group(1)
        hour = int(hhmm[:2])
        if 3 <= hour <= 6:  # 03xx through 06xx
            incident_files.append(f)

print(f"Files in incident window (03:00-07:00 UTC): {len(incident_files)}")
for f in incident_files:
    print(f"  {f}")

# Download each file with retry
print(f"\nDownloading to {OUTPUT_DIR}...")
downloaded = 0
failed = 0

for fname in incident_files:
    outpath = os.path.join(OUTPUT_DIR, fname)
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100000:
        print(f"  {fname}: already exists, skipping")
        downloaded += 1
        continue

    url = THREDDS_BASE + fname

    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = resp.read()

            with open(outpath, 'wb') as f:
                f.write(data)

            size_mb = len(data) / 1024 / 1024
            print(f"  {fname}: {size_mb:.1f} MB OK")
            downloaded += 1
            break
        except Exception as e:
            wait = 2 ** (attempt + 1)
            print(f"  {fname}: attempt {attempt+1} failed ({e}), retrying in {wait}s...")
            time.sleep(wait)
    else:
        print(f"  {fname}: FAILED after 4 attempts")
        failed += 1

print(f"\nDone: {downloaded} downloaded, {failed} failed")
print(f"Files in {OUTPUT_DIR}:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    size = os.path.getsize(os.path.join(OUTPUT_DIR, f))
    print(f"  {f}: {size/1024/1024:.1f} MB")
