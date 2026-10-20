# Najran_0Carbon — data and code

Data, Earth Engine scripts and analysis code for:

> Elkhrachy, I. *Guiding low-carbon growth in an arid oasis city: satellite evidence from Najran, Saudi Arabia.* Manuscript submitted to *Sustainable Cities and Society*.

**Interactive atlas:** `docs/index.html` (published with GitHub Pages).

## Contents

| Folder | What it holds |
|---|---|
| `gee/` | Google Earth Engine scripts. 01: Landsat harmonisation, Random Forest land cover, NDVI and LST for 2000–2026. 02: rainfall (CHIRPS), GRACE and GLDAS. 03: building footprints and terrain. 04: MODIS day and night LST. 05: viewer for the blind validation of the 2000 map. |
| `python/` | Analysis code: stratified accuracy and area estimation (Olofsson et al. 2014), trend statistics (Mann–Kendall, Sen, Pettitt), growth scenarios (ML–CA), Monte Carlo carbon model, MODIS contrasts and figures. `shpread.py`, `shpwrite.py` and `tiflzw.py` are small GDAL-free readers and writers. |
| `data/` | Tables behind the results: zone indicators (GeoJSON and CSV), blind validation points with reference labels, epoch statistics, annual vegetation, rainfall, GRACE, GLDAS, MODIS contrasts, scenario metrics and carbon results. |
| `rasters/` | GeoTIFFs (EPSG:32638, 30 m): land cover 2000–2026, summer LST and relative LST, NDVI trend, growth epoch, fate of 2000 vegetation, 2040 scenarios, SRTM elevation and slope. `.clr` files hold the class colours. |
| `docs/` | The bilingual (Arabic/English) zone atlas. |

## Class codes
Land cover: 0 built-up, 1 vegetation, 2 sand / bare soil, 4 rock / gravel, 255 no data.
Validation reference classes: B built-up, V vacant subdivided land, G vegetation, S sand, R rock.

## Reproducing the results
1. Run `gee/Najran_01_...js` to `gee/Najran_04_...js` in the Earth Engine Code Editor; the exports go to a Google Drive folder named `Najran_LowCarbon`.
2. Put the exported files next to the Python scripts and edit the path constants at the top of each script.
3. Run `accuracy_full.py`, `series_stats.py`, `scenarios.py`, `carbon.py`, `night_stats.py` and `day_stats.py`. Python 3.10+ with numpy, pandas, scipy, scikit-learn, matplotlib and tifffile.

## Data not included
Municipal plan and parcel layers were provided by Najran Municipality and cannot be redistributed. Zone indicators derived from them are given only as zone totals.

## Licence
Code: MIT. Data and rasters: CC BY 4.0. Landsat, MODIS, CHIRPS, GRACE, GLDAS and SRTM are public data from USGS, NASA and UCSB-CHG; Microsoft building footprints are licensed under ODbL.

## Citation
Please cite the article and this archive (see `CITATION.cff`).
