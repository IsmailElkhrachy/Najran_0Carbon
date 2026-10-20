/*****************************************************************
 * Najran Low-Carbon Smart City Study — Script 01 (version 2)
 * Urban growth, vegetation decline and land surface temperature
 * 2000–2026 | Landsat Collection 2 Level-2 | Google Earth Engine
 *
 * Author : Ismail Elkhrachy, Department of Civil Engineering,
 *          Najran University
 *
 * Changes from version 1 (after checking the first results)
 *   A. New class "Rock/gravel": dark basement rock and gravel plains
 *      were being mapped as built-up (built-up 2000 > 2010).
 *   B. Temporal consistency: a pixel cannot be built-up in an epoch
 *      and non-built in every later epoch (cities do not shrink).
 *   C. Annual vegetation uses the maximum of MONTHLY medians, so the
 *      result no longer rises with the number of images (Landsat 9
 *      doubled the image count from 2022).
 *   D. Brightness band added to the classification features.
 *
 * What it does
 *   1. Builds harmonised Landsat 5/7/8/9 composites for six epochs
 *      (each epoch = 3-year window centred on the epoch year).
 *   2. Classifies 5 land-cover classes with Random Forest, using
 *      temporally-constrained training labels (ESA WorldCover 2021
 *      + GHSL built-up surface of the same epoch + spectral rules).
 *   3. Applies the "built-up cannot revert" rule across epochs.
 *   4. Derives summer (Jun–Sep) LST, the urban core of each epoch
 *      and the surface urban heat island (SUHI).
 *   5. Computes the annual vegetation series and NDVI trend 2000–2026.
 *   6. Exports maps, tables and BLIND validation points.
 *
 * Class codes: 0 Built-up | 1 Vegetation | 2 Sand/bare soil | 3 Water | 4 Rock/gravel
 *****************************************************************/

// =============================== 0. SETTINGS ===============================
var AOI = ee.Geometry.Rectangle([43.95, 17.38, 44.60, 17.70], null, false);

var CRS = 'EPSG:32638';              // WGS 84 / UTM zone 38N
var SCALE = 30;
var EPOCHS = [2000, 2005, 2010, 2015, 2020, 2026];
var TREND_FIRST = 2000, TREND_LAST = 2026;
var LST_MONTHS = [6, 9];             // summer: June–September
var SEED = 42;
var TRAIN_PER_CLASS = 600;
var VAL_PER_CLASS = 60;              // 5 classes x 60 = 300 blind points per epoch
var URBAN_RADIUS_M = 500;
var URBAN_DENSITY = 0.25;
var RURAL_RING_M = [2000, 6000];
var VEG_NDVI = 0.30;
var ROCK_SLOPE = 5;                  // degrees: bare pixels on slopes >= 5° -> rock
var ROCK_BRIGHT = 0.18;              // mean VIS-NIR reflectance below this -> dark rock/gravel
var SAND_BRIGHT = 0.22;              // flat bare pixels brighter than this -> sand/soil
var DRIVE_FOLDER = 'Najran_LowCarbon_GEE_v2';
var RUN_EXPORTS = true;
var PRINT_STATS = false;             // console printing is slow; tables are exported

var CLASS_KEYS = ['built', 'veg', 'sand', 'water', 'rock'];
var CLASS_NAMES = ['Built-up', 'Vegetation', 'Sand / bare soil', 'Water', 'Rock / gravel'];
var PALETTE = ['d7191c', '1a9641', 'f2e6c2', '2b83ba', '7a6a5d'];
var CLASS_IDS = [0, 1, 2, 3, 4];
var SR_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
var FEATURES = SR_NAMES.concat(['NDVI', 'NDBI', 'MNDWI', 'BSI', 'bright',
                                'NDVI_p10', 'NDVI_p90', 'elevation', 'slope']);
var PROJ = ee.Projection(CRS);

// ========================= 1. LANDSAT PRE-PROCESSING ========================
function maskScale(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(58).eq(0);           // bits 1,3,4,5
  var unsat = img.select('QA_RADSAT').eq(0);
  var sr = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  var st = img.select('ST_B.*').multiply(0.00341802).add(149.0).subtract(273.15); // °C
  return img.addBands(sr, null, true).addBands(st, null, true)
            .updateMask(clear).updateMask(unsat);
}

// TM/ETM+ -> OLI surface-reflectance harmonisation (Roy et al., 2016, OLS).
var TM_SLOPE = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]);
var TM_ITCP = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]);

function fromTM(img) {
  var sr = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'], SR_NAMES)
              .multiply(TM_SLOPE).add(TM_ITCP);
  return ee.Image(sr.addBands(img.select(['ST_B6'], ['lst']))
                    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']));
}

function fromOLI(img) {
  var sr = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'], SR_NAMES);
  return ee.Image(sr.addBands(img.select(['ST_B10'], ['lst']))
                    .copyProperties(img, ['system:time_start', 'SPACECRAFT_ID']));
}

// Landsat 7 only before 2013 (before Landsat 8 and before its orbit drift).
function landsat(start, end, useL7) {
  var f = function (id) {
    return ee.ImageCollection(id).filterBounds(AOI).filterDate(start, end).map(maskScale);
  };
  var tm = f('LANDSAT/LT05/C02/T1_L2');
  if (useL7) tm = tm.merge(f('LANDSAT/LE07/C02/T1_L2'));
  var oli = f('LANDSAT/LC08/C02/T1_L2').merge(f('LANDSAT/LC09/C02/T1_L2'));
  return tm.map(fromTM).merge(oli.map(fromOLI));
}

function addIndices(img) {
  var ndvi = img.normalizedDifference(['nir', 'red']).rename('NDVI');
  var ndbi = img.normalizedDifference(['swir1', 'nir']).rename('NDBI');
  var mndwi = img.normalizedDifference(['green', 'swir1']).rename('MNDWI');
  var bsi = img.expression('((S + R) - (N + B)) / ((S + R) + (N + B))', {
    S: img.select('swir1'), R: img.select('red'), N: img.select('nir'), B: img.select('blue')
  }).rename('BSI');
  var bright = img.select(['blue', 'green', 'red', 'nir']).reduce(ee.Reducer.mean()).rename('bright');
  return img.addBands([ndvi, ndbi, mndwi, bsi, bright]);
}

var dem = ee.Image('USGS/SRTMGL1_003').select('elevation');
var terrain = dem.addBands(ee.Terrain.slope(dem)).float();   // elevation, slope

function epochData(year) {
  var start = ee.Date.fromYMD(year - 1, 1, 1);
  var end = ee.Date.fromYMD(year + 2, 1, 1);
  var col = landsat(start, end, year < 2013).map(addIndices);
  var med = col.select(SR_NAMES.concat(['NDVI', 'NDBI', 'MNDWI', 'BSI', 'bright'])).median();
  var pct = col.select('NDVI').reduce(ee.Reducer.percentile([10, 90]));  // NDVI_p10, NDVI_p90
  var summer = col.filter(ee.Filter.calendarRange(LST_MONTHS[0], LST_MONTHS[1], 'month'));
  var lst = summer.select('lst').median().rename('LST');
  return med.addBands(pct).addBands(terrain).addBands(lst).clip(AOI)
            .set({year: year, n_images: col.size(), n_summer: summer.size()});
}

// ====================== 2. TEMPORALLY-CONSTRAINED LABELS ====================
// ESA WorldCover 2021 (10 m) -> pure 30 m pixels (>= 90 % one class).
var wcClass = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map')
  .remap([50, 10, 40, 60, 80], [0, 1, 1, 2, 3])   // built | tree, crop | bare | water
  .clip(AOI);
var fractions = ee.Image.cat([0, 1, 2, 3].map(function (k) {
  return wcClass.eq(k).unmask(0)
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 100})
    .reproject({crs: CRS, scale: SCALE});
}));
var wcPure = fractions.toArray().arrayArgmax().arrayGet([0]).toInt()
  .updateMask(fractions.reduce(ee.Reducer.max()).gte(0.9)).rename('wc');

function ghslBuilt(year) {
  var y = Math.min(Math.floor(year / 5) * 5, 2020);
  var img = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S')
    .filter(ee.Filter.calendarRange(y, y, 'year')).first();
  return ee.Image(img).select('built_surface').divide(10000).rename('ghsl'); // 100 m cell
}

// WorldCover "bare" is split into sand/soil and rock/gravel using slope and
// brightness; ambiguous pixels in between are not used for training.
function epochLabels(img, year) {
  var gh = ghslBuilt(year).unmask(0);
  var p90 = img.select('NDVI_p90');
  var bright = img.select('bright');
  var slope = img.select('slope');
  var wc = wcPure.unmask(-1);
  var bareBase = wc.eq(2).and(gh.lt(0.02)).and(p90.lt(0.15));
  var built = wc.eq(0).and(gh.gte(0.25));
  var veg = wc.eq(1).and(gh.lt(0.02)).and(p90.gte(VEG_NDVI));
  var sand = bareBase.and(slope.lt(ROCK_SLOPE)).and(bright.gte(SAND_BRIGHT));
  var rock = bareBase.and(slope.gte(ROCK_SLOPE).or(bright.lt(ROCK_BRIGHT)));
  var water = wc.eq(3).and(gh.lt(0.02)).and(img.select('MNDWI').gt(0));
  var lab = ee.Image(-1)
    .where(built.unmask(0), 0).where(veg.unmask(0), 1)
    .where(sand.unmask(0), 2).where(water.unmask(0), 3)
    .where(rock.unmask(0), 4)
    .rename('label').toInt().clip(AOI);
  return lab.updateMask(lab.gte(0));
}

// ========================== 3. CLASSIFICATION ===============================
function classifyEpoch(year) {
  var img = epochData(year);
  var lab = epochLabels(img, year);
  var samples = img.select(FEATURES).addBands(lab).stratifiedSample({
    numPoints: TRAIN_PER_CLASS, classBand: 'label', region: AOI,
    scale: SCALE, projection: PROJ, seed: SEED + year,
    dropNulls: true, tileScale: 8, geometries: true
  }).randomColumn('rand', SEED);
  var train = samples.filter(ee.Filter.lt('rand', 0.7));
  var test = samples.filter(ee.Filter.gte('rand', 0.7));
  var rf = ee.Classifier.smileRandomForest({numberOfTrees: 300, seed: SEED})
    .train({features: train, classProperty: 'label', inputProperties: FEATURES});
  var lulc = img.select(FEATURES).classify(rf).rename('lulc')
    .reproject({crs: CRS, scale: SCALE})
    .focalMode(1)                                   // 3x3 majority filter
    .reproject({crs: CRS, scale: SCALE})
    .toByte().clip(AOI);
  // Holdout accuracy uses automatic (pure) labels and is optimistic;
  // the paper reports area-adjusted accuracy from the blind points.
  var em = test.classify(rf).errorMatrix('label', 'classification', CLASS_IDS);
  return {year: year, img: img, raw: lulc, em: em};
}

// Built-up cannot revert: processed from the newest epoch backwards.
// A pixel stays built-up in epoch t only if it is also built-up in the
// (already corrected) epoch t+1; otherwise it takes its class at t+1.
// The newest epoch (2026) is kept as classified.
function applyTemporalRule(results) {
  var n = results.length;
  results[n - 1].lulc = results[n - 1].raw;
  for (var i = n - 2; i >= 0; i--) {
    var next = results[i + 1].lulc;
    var cur = results[i].raw;
    var falseBuilt = cur.eq(0).and(next.neq(0));
    results[i].lulc = cur.where(falseBuilt, next).toByte();
    results[i].removedBuilt = falseBuilt;
  }
  return results;
}

// ===================== 4. URBAN CORE, LST AND SUHI ==========================
function urbanCore(lulc) {
  return lulc.eq(0).reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.circle(URBAN_RADIUS_M, 'meters')
  }).reproject({crs: CRS, scale: SCALE}).gte(URBAN_DENSITY).rename('urban').toByte();
}

// Rural reference: natural desert surfaces (sand or rock) 2–6 km from the urban core.
function ruralRing(urban, lulc) {
  var distM = urban.fastDistanceTransform(256, 'pixels', 'squared_euclidean')
    .sqrt().multiply(SCALE).reproject({crs: CRS, scale: SCALE});
  return distM.gte(RURAL_RING_M[0]).and(distM.lte(RURAL_RING_M[1]))
    .and(lulc.eq(2).or(lulc.eq(4))).rename('rural');
}

function meanIn(img, mask, band) {
  return img.select(band).updateMask(mask).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: AOI, scale: SCALE, crs: CRS,
    maxPixels: 1e10, tileScale: 4
  }).get(band);
}

function areaIn(mask) {
  return ee.Image.pixelArea().divide(1e6).rename('area').updateMask(mask).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: AOI, scale: SCALE, crs: CRS,
    maxPixels: 1e10, tileScale: 4
  }).get('area');
}

function groupsToDict(groups, key, prefix) {
  var keys = ee.List(CLASS_KEYS);
  return ee.Dictionary(ee.List(groups).iterate(function (g, acc) {
    g = ee.Dictionary(g);
    var name = ee.String(prefix).cat(ee.String(keys.get(ee.Number(g.get('class')).int())));
    return ee.Dictionary(acc).set(name, g.get(key));
  }, ee.Dictionary({})));
}

function epochStats(r) {
  var lulc = r.lulc, lst = r.img.select('LST');
  var urban = urbanCore(lulc);
  var rural = ruralRing(urban, lulc);
  var area = ee.Image.pixelArea().divide(1e6);
  var areaGroups = area.addBands(lulc).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: AOI, scale: SCALE, crs: CRS, maxPixels: 1e10, tileScale: 4
  }).get('groups');
  var lstGroups = lst.addBands(lulc).reduceRegion({
    reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'class'}),
    geometry: AOI, scale: SCALE, crs: CRS, maxPixels: 1e10, tileScale: 4
  }).get('groups');
  var lstU = meanIn(r.img, urban, 'LST');
  var lstR = meanIn(r.img, rural, 'LST');
  var removed = r.removedBuilt ? areaIn(r.removedBuilt) : 0;
  var props = ee.Dictionary({
    year: r.year,
    n_images: r.img.get('n_images'),
    n_summer_images: r.img.get('n_summer'),
    urban_core_km2: areaIn(urban),
    built_raw_km2: areaIn(r.raw.eq(0)),
    built_removed_by_rule_km2: removed,
    lst_urban_C: lstU,
    lst_rural_desert_C: lstR,
    suhi_C: ee.Number(lstU).subtract(ee.Number(lstR)),
    oa_holdout: r.em.accuracy(),
    kappa_holdout: r.em.kappa()
  }).combine(groupsToDict(areaGroups, 'sum', 'area_km2_'))
    .combine(groupsToDict(lstGroups, 'mean', 'lst_C_'));
  return {feature: ee.Feature(null, props), urban: urban};
}

// ====================== 5. BLIND VALIDATION POINTS ==========================
function validationPoints(lulc, year) {
  return lulc.rename('map_class').stratifiedSample({
    numPoints: VAL_PER_CLASS, classBand: 'map_class', region: AOI,
    scale: SCALE, projection: PROJ, seed: SEED + 1000 + year,
    dropNulls: true, tileScale: 8, geometries: true
  }).randomColumn('shuffle', SEED + year).sort('shuffle')
    .map(function (f) {
      return f.set({pid: ee.String(String(year) + '_').cat(f.id()), epoch: year});
    });
}

// ======================== 6. RUN ALL EPOCHS =================================
Map.centerObject(AOI, 11);
Map.setOptions('HYBRID');

var results = applyTemporalRule(EPOCHS.map(classifyEpoch));

var statFeatures = [];
results.forEach(function (r) {
  var year = r.year;
  var s = epochStats(r);
  statFeatures.push(s.feature);

  if (year === EPOCHS[0] || year === EPOCHS[EPOCHS.length - 1]) {
    Map.addLayer(r.img, {bands: ['swir1', 'nir', 'red'], min: 0.05, max: 0.45},
                 'Landsat ' + year + ' (SWIR1-NIR-Red)', false);
    Map.addLayer(r.raw, {min: 0, max: 4, palette: PALETTE}, 'LULC ' + year + ' (raw)', false);
    Map.addLayer(r.lulc, {min: 0, max: 4, palette: PALETTE}, 'LULC ' + year, year !== EPOCHS[0]);
    Map.addLayer(r.img.select('LST'), {min: 30, max: 60,
                 palette: ['313695', '74add1', 'ffffbf', 'f46d43', 'a50026']},
                 'Summer LST ' + year + ' (°C)', false);
    Map.addLayer(s.urban.selfMask(), {palette: ['000000']}, 'Urban core ' + year, false);
  }

  if (PRINT_STATS) {
    print('Epoch ' + year + ' | holdout confusion matrix (auto labels):', r.em);
  }

  if (RUN_EXPORTS) {
    Export.image.toDrive({
      image: r.lulc.addBands(s.urban), description: 'Najran_LULC_' + year,
      folder: DRIVE_FOLDER, region: AOI, crs: CRS, scale: SCALE, maxPixels: 1e10
    });
    Export.image.toDrive({
      image: r.img.select('LST').toFloat(), description: 'Najran_LST_' + year,
      folder: DRIVE_FOLDER, region: AOI, crs: CRS, scale: SCALE, maxPixels: 1e10
    });
    var val = validationPoints(r.lulc, year);
    Export.table.toDrive({
      collection: val.select(['pid']), description: 'Najran_VAL_blind_' + year,
      folder: DRIVE_FOLDER, fileFormat: 'KML'
    });
    Export.table.toDrive({
      collection: val.select(['pid', 'epoch', 'map_class']),
      description: 'Najran_VAL_key_' + year, folder: DRIVE_FOLDER, fileFormat: 'CSV'
    });
  }
});

var epochFC = ee.FeatureCollection(statFeatures);
if (PRINT_STATS) print('Epoch statistics:', epochFC);

// ===================== 7. ANNUAL VEGETATION 2000–2026 =======================
// Annual greenness = maximum of the 12 monthly median NDVI images.
// Monthly medians make the value independent of how many scenes a year has.
var annualImages = [];
var vegFeatures = [];
for (var y = TREND_FIRST; y <= TREND_LAST; y++) {
  var col = landsat(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y + 1, 1, 1), y < 2013)
    .map(addIndices).select('NDVI');
  var empty = ee.Image.constant(0).float().rename('NDVI').updateMask(0);
  var monthly = ee.ImageCollection(ee.List.sequence(1, 12).map(function (m) {
    var mCol = col.filter(ee.Filter.calendarRange(m, m, 'month'));
    return ee.ImageCollection([empty]).merge(mCol).median().rename('NDVI');
  }));
  var peak = monthly.max().rename('NDVI');
  var im = ee.Image.constant(y).float().rename('year').addBands(peak).clip(AOI).set('year', y);
  annualImages.push(im);

  var nMonths = ee.List.sequence(1, 12).map(function (m) {
    return ee.Number(col.filter(ee.Filter.calendarRange(m, m, 'month')).size()).min(1);
  }).reduce(ee.Reducer.sum());
  vegFeatures.push(ee.Feature(null, {
    year: y,
    veg_km2: areaIn(peak.gte(VEG_NDVI)),
    valid_km2: areaIn(peak.mask()),
    n_images: col.size(),
    n_months_with_data: nMonths
  }));
}
var annual = ee.ImageCollection(annualImages);
var vegFC = ee.FeatureCollection(vegFeatures);

var sen = annual.select(['year', 'NDVI']).reduce(ee.Reducer.sensSlope());           // slope, offset
var mk = annual.select(['year', 'NDVI']).reduce(ee.Reducer.kendallsCorrelation(2)); // tau, p-value
var trend = sen.select('slope').multiply(10).rename('NDVI_slope_per_decade').addBands(mk);

Map.addLayer(trend.select('NDVI_slope_per_decade'),
  {min: -0.1, max: 0.1, palette: ['8c510a', 'd8b365', 'f5f5f5', '5ab4ac', '01665e']},
  'NDVI trend per decade 2000–2026', false);

print(ui.Chart.feature.byFeature(vegFC, 'year', ['veg_km2'])
  .setChartType('LineChart')
  .setOptions({title: 'Najran: vegetated area (max monthly NDVI ≥ ' + VEG_NDVI + ')',
               hAxis: {title: 'Year', format: '####'}, vAxis: {title: 'km²'},
               pointSize: 4, legend: {position: 'none'}}));

if (RUN_EXPORTS) {
  Export.image.toDrive({
    image: trend.toFloat(), description: 'Najran_NDVI_trend', folder: DRIVE_FOLDER,
    region: AOI, crs: CRS, scale: SCALE, maxPixels: 1e10
  });
  Export.table.toDrive({
    collection: epochFC, description: 'Najran_epoch_stats',
    folder: DRIVE_FOLDER, fileFormat: 'CSV'
  });
  Export.table.toDrive({
    collection: vegFC, description: 'Najran_NDVI_annual',
    folder: DRIVE_FOLDER, fileFormat: 'CSV'
  });
  Export.table.toDrive({
    collection: ee.FeatureCollection([ee.Feature(AOI, {name: 'Najran_AOI'})]),
    description: 'Najran_AOI', folder: DRIVE_FOLDER, fileFormat: 'KML'
  });
}

// ============================== 8. LEGEND ===================================
var legend = ui.Panel({style: {position: 'bottom-left', padding: '6px 10px'}});
legend.add(ui.Label('Land cover', {fontWeight: 'bold'}));
CLASS_NAMES.forEach(function (n, i) {
  legend.add(ui.Panel([
    ui.Label('', {backgroundColor: '#' + PALETTE[i], padding: '8px', margin: '2px 6px 2px 0'}),
    ui.Label(n, {margin: '2px 0'})
  ], ui.Panel.Layout.Flow('horizontal')));
});
Map.add(legend);
Map.addLayer(ee.Image().paint(AOI, 0, 2), {palette: ['ffff00']}, 'Study area (AOI)');
