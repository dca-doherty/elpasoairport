#!/usr/bin/env python3
"""
Download KHDX Level 2 data from additional dates to check if
returns at az 193-195° are persistent ground clutter.

Downloads Feb 8 and Feb 12 (2 days before and after event).
"""

import urllib.request
import xml.etree.ElementTree as ET
import os
import time
import re

dates = ['20260208', '20260212']
labels = {'20260208': 'Feb 8 (3 days before)', '20260212': 'Feb 12 (1 day after)'}

for date in dates:
    catalog_url = f"https://thredds.ucar.edu/thredds/catalog/nexrad/level2/KHDX/{date}/catalog.xml"
    base_url = f"https://thredds.ucar.edu/thredds/fileServer/nexrad/level2/KHDX/{date}/"
    output_dir = f"/home/user/uap-transient-research/el_paso_airspace/khdx_extra_{date}/"

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"KHDX {labels[date]} — {date}")
    print(f"{'='*70}")

    print("Fetching catalog...")
    try:
        req = urllib.request.Request(catalog_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            catalog_xml = resp.read().decode()
    except Exception as e:
        print(f"Error: {e}")
        continue

    root = ET.fromstring(catalog_xml)
    all_files = []
    for dataset in root.iter('{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset'):
        name = dataset.get('name', '')
        if name.startswith(f'Level2_KHDX_{date}_') and name.endswith('.ar2v'):
            all_files.append(name)

    all_files.sort()

    # Filter to 03:00-07:00 UTC
    target_files = []
    for f in all_files:
        match = re.search(r'_(\d{4})\.ar2v', f)
        if match:
            hour = int(match.group(1)[:2])
            if 3 <= hour <= 6:
                target_files.append(f)

    print(f"Found {len(target_files)} files in 03:00-07:00 window")

    downloaded = 0
    for fname in target_files:
        outpath = os.path.join(output_dir, fname)
        if os.path.exists(outpath) and os.path.getsize(outpath) > 100000:
            downloaded += 1
            continue

        url = base_url + fname
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=60) as resp:
                    data = resp.read()
                with open(outpath, 'wb') as f:
                    f.write(data)
                size_mb = len(data) / 1024 / 1024
                print(f"  {fname}: {size_mb:.1f} MB")
                downloaded += 1
                break
            except Exception as e:
                wait = 2 ** (attempt + 1)
                print(f"  {fname}: attempt {attempt+1} failed ({e}), retry in {wait}s...")
                time.sleep(wait)
        else:
            print(f"  {fname}: FAILED")

    print(f"Downloaded: {downloaded}/{len(target_files)}")
