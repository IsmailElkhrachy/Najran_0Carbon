/**** Najran low-carbon city – script 04: MODIS day and NIGHT land surface temperature
 *
 * Purpose: the paper's heat results use Landsat (overpass ~10:30, daytime only).
 * Arid-city heat islands are often negative by day and positive at night, so a
 * reviewer will ask about night-time LST. MODIS Terra (~22:30) and Aqua (~01:30)
 * provide night LST at 1 km since 2000 / 2002.
 *
 * Outputs (Google Drive folder Najran_LowCarbon):
 *   1. Najran_MODIS_LST_summer_night_2003_2025.tif     - night, Terra 22:30 + Aqua 01:30
 *   2. Najran_MODIS_LST_summer_day1330_2003_2025.tif   - Aqua 13:30, near the afternoon peak
 *   3. Najran_MODIS_LST_summer_day1030_2003_2025.tif   - Terra 10:30, same time as Landsat
 *   4. Najran_MODIS_LST_counts.csv                      - valid 8-day composites per year
 * Each .tif has 23 bands, one per year (deg C). Summer = June-August 8-day composites
 * (MOD11A2/MYD11A2 v061),
 * good-quality pixels only (QC bits 0-1 = 00 or 01).
 * After export: upload the three .tif and the .csv to the chat; the class contrasts
 * (built-up vs sand vs rock vs farmland) are computed with the 30 m land-cover maps.
 ****/

var AOI = ee.Geometry.Rectangle([43.95, 17.38, 44.60, 17.70], null, false);
var CRS = 'EPSG:32638';
var DRIVE_FOLDER = 'Najran_LowCarbon';
var Y0 = 2003, Y1 = 2025;          // first full Aqua summer to last full summer

function prep(band, qcBand) {
  return function (img) {
    var qc = img.select(qcBand);
    var good = qc.bitwiseAnd(3).lte(1);                  // mandatory QA: produced, good or other quality
    var lst = img.select(band).multiply(0.02).subtract(273.15);
    return lst.updateMask(good).rename('lst')
              .copyProperties(img, ['system:time_start']);
  };
}

function summer(colId, band, qcBand, y) {
  return ee.ImageCollection(colId)
    .filterDate(ee.Date.fromYMD(y, 6, 1), ee.Date.fromYMD(y, 9, 1))
    .filterBounds(AOI)
    .map(prep(band, qcBand));
}

// kind: 'night' (Terra + Aqua), 'day1330' (Aqua only), 'day1030' (Terra only)
function yearImage(kind, y) {
  var col;
  if (kind === 'night') {
    col = summer('MODIS/061/MOD11A2', 'LST_Night_1km', 'QC_Night', y)
      .merge(summer('MODIS/061/MYD11A2', 'LST_Night_1km', 'QC_Night', y));
  } else if (kind === 'day1330') {
    col = summer('MODIS/061/MYD11A2', 'LST_Day_1km', 'QC_Day', y);
  } else {
    col = summer('MODIS/061/MOD11A2', 'LST_Day_1km', 'QC_Day', y);
  }
  return col.mean().rename(kind + '_' + y).clip(AOI);
}

var years = ee.List.sequence(Y0, Y1).getInfo();

['night', 'day1330', 'day1030'].forEach(function (kind) {
  var stack = ee.Image.cat(years.map(function (y) { return yearImage(kind, y); })).toFloat();
  Export.image.toDrive({
    image: stack,
    description: 'Najran_MODIS_LST_summer_' + kind + '_' + Y0 + '_' + Y1,
    folder: DRIVE_FOLDER, region: AOI, crs: CRS, scale: 1000, maxPixels: 1e9
  });
  if (kind === 'night') {
    Map.addLayer(stack.select('night_' + Y1), {min: 25, max: 35,
      palette: ['313695', '74add1', 'ffffbf', 'f46d43', 'a50026']}, 'Night LST summer ' + Y1);
  }
});

// number of valid composites per year (data availability check)
var counts = ee.FeatureCollection(years.map(function (y) {
  var n = summer('MODIS/061/MOD11A2', 'LST_Night_1km', 'QC_Night', y).size()
    .add(summer('MODIS/061/MYD11A2', 'LST_Night_1km', 'QC_Night', y).size());
  return ee.Feature(null, {year: y, n_composites: n});
}));
Export.table.toDrive({collection: counts, description: 'Najran_MODIS_LST_counts',
                      folder: DRIVE_FOLDER, fileFormat: 'CSV'});

Map.centerObject(AOI, 10);
Map.addLayer(ee.Image().paint(AOI, 0, 2), {palette: ['ffff00']}, 'Study area');
print('Run the 4 tasks in the Tasks tab, then upload the files to the chat.');
