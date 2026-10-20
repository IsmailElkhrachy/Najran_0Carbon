/*****************************************************************
 * Najran Low-Carbon Smart City Study — Script 02
 * Drivers of vegetation loss in Wadi Najran, 2000–2026:
 *   (1) rainfall (CHIRPS), (2) terrestrial water storage (GRACE /
 *   GRACE-FO mascons), (3) modelled groundwater storage (GLDAS CLSM).
 *
 * Light-weight: only exports small CSV tables (a few hundred rows),
 * so it runs even with a restricted Earth Engine quota.
 *
 * Outputs (Google Drive folder DRIVE_FOLDER)
 *   Najran_rain_annual.csv       year, rain_mm (CHIRPS, AOI mean)
 *   Najran_rain_monthly.csv      date, rain_mm
 *   Najran_GRACE_monthly.csv     date, lwe_cm (equivalent water thickness anomaly)
 *   Najran_GLDAS_GWS_monthly.csv date, gws_mm (groundwater storage)
 *****************************************************************/

var AOI = ee.Geometry.Rectangle([43.95, 17.38, 44.60, 17.70], null, false); // city + wadi
var REGION = ee.Geometry.Rectangle([43.5, 17.0, 45.0, 18.2], null, false);   // coarse grids
var DRIVE_FOLDER = 'Najran_LowCarbon_GEE_v2';
var Y0 = 2000, Y1 = 2026;

function regionMean(img, geom, scale, band) {
  return ee.Number(img.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, scale: scale,
    bestEffort: true, maxPixels: 1e9
  }).get(band));
}

// ---------------------------------------------------------- 1. CHIRPS rainfall
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation');

var annual = ee.FeatureCollection(ee.List.sequence(Y0, Y1).map(function (y) {
  y = ee.Number(y);
  var col = chirps.filter(ee.Filter.calendarRange(y, y, 'year'));
  var tot = col.sum().rename('rain');
  return ee.Feature(null, {
    year: y,
    rain_mm: ee.Algorithms.If(col.size().gt(0), regionMean(tot, AOI, 5566, 'rain'), null),
    n_days: col.size()
  });
}));

var months = ee.List.sequence(0, (Y1 - Y0 + 1) * 12 - 1);
var monthly = ee.FeatureCollection(months.map(function (i) {
  var start = ee.Date.fromYMD(Y0, 1, 1).advance(i, 'month');
  var col = chirps.filterDate(start, start.advance(1, 'month'));
  return ee.Feature(null, {
    date: start.format('YYYY-MM'),
    rain_mm: ee.Algorithms.If(col.size().gt(0),
      regionMean(col.sum().rename('rain'), AOI, 5566, 'rain'), null)
  });
}));

// ---------------------------------------------------------- 2. GRACE mascons
// JPL mascon solution (RL06.x, CRI-filtered). The dataset ends where the
// catalogue ends; GRACE-FO months are included if present.
var grace = ee.ImageCollection('NASA/GRACE/MASS_GRIDS_V04/MASCON_CRI').select('lwe_thickness');
var graceFC = grace.filterBounds(REGION).map(function (img) {
  return ee.Feature(null, {
    date: img.date().format('YYYY-MM'),
    lwe_cm: regionMean(img, REGION, 55660, 'lwe_thickness')
  });
});

// ---------------------------------------------------------- 3. GLDAS groundwater
// GLDAS-2.2 CLSM with GRACE data assimilation, daily, 0.25 deg, from 2003.
var gldas = ee.ImageCollection('NASA/GLDAS/V022/CLSM/G025/DA1D').select('GWS_tavg');
var gldasMonths = ee.List.sequence(0, (Y1 - 2003 + 1) * 12 - 1);
var gldasFC = ee.FeatureCollection(gldasMonths.map(function (i) {
  var start = ee.Date.fromYMD(2003, 1, 1).advance(i, 'month');
  var col = gldas.filterDate(start, start.advance(1, 'month'));
  return ee.Feature(null, {
    date: start.format('YYYY-MM'),
    gws_mm: ee.Algorithms.If(col.size().gt(0),
      regionMean(col.mean().rename('gws'), REGION, 27830, 'gws'), null)
  });
}));

// ---------------------------------------------------------- charts and exports
print(ui.Chart.feature.byFeature(annual, 'year', ['rain_mm'])
  .setChartType('ColumnChart')
  .setOptions({title: 'Najran annual rainfall (CHIRPS)', hAxis: {format: '####'},
               vAxis: {title: 'mm'}, legend: {position: 'none'}}));
print('GRACE images in region:', grace.filterBounds(REGION).size());

Export.table.toDrive({collection: annual, description: 'Najran_rain_annual',
  folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: ['year', 'rain_mm', 'n_days']});
Export.table.toDrive({collection: monthly, description: 'Najran_rain_monthly',
  folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: ['date', 'rain_mm']});
Export.table.toDrive({collection: graceFC, description: 'Najran_GRACE_monthly',
  folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: ['date', 'lwe_cm']});
Export.table.toDrive({collection: gldasFC, description: 'Najran_GLDAS_GWS_monthly',
  folder: DRIVE_FOLDER, fileFormat: 'CSV', selectors: ['date', 'gws_mm']});
