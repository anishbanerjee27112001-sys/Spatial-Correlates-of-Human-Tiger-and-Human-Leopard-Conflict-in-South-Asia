// BiasCorrected background sampling for HLC and HTC.
var assetRoot = "projects/hwc-esa/assets/";
var exportFolder = "HCC_Background_Exports";
function runBackground(label) {
var assetName = assetRoot + label + "_Final_Conflict_Data";
// ======================================================
// ROBUST BIAS-CORRECTED BACKGROUND POINTS IN GEE
// ======================================================

// -----------------------------
// 0. SETTINGS
// -----------------------------

var nBackground = 40000;
var exclusionBufferM = 5000;
var settlementDecay = 5000;
var roadDecay = 5000;
var minBias = 0.05;        // important: prevents near-zero acceptance
var sampleScale = 250;     // coarser = lighter + faster
var randomSeed = 123;

// -----------------------------
// 1. LOAD PRESENCE POINTS
// -----------------------------
var raw = ee.FeatureCollection(assetName);

var pts = raw
  .filter(ee.Filter.notNull(['Long', 'Lat']))
  .map(function(f) {
    var lon = ee.Number.parse(ee.String(f.get('Long')));
    var lat = ee.Number.parse(ee.String(f.get('Lat')));
    return ee.Feature(ee.Geometry.Point([lon, lat]), f.toDictionary());
  });

var studyRegion = pts.geometry().buffer(50000);

print(label + ' count', pts.size());
Map.centerObject(studyRegion, 6);
Map.addLayer(pts, {color: 'red'}, label + ' points');

// -----------------------------
// 2. GHSL SETTLEMENT DATA
// -----------------------------
var ghsl = ee.ImageCollection("JRC/GHSL/P2023A/GHS_SMOD")
  .mosaic()
  .clip(studyRegion);

var smod = ghsl.select('smod_code');

// broad settlement definition
var settlementBinary = smod.gte(13).selfMask();

Map.addLayer(settlementBinary, {min: 0, max: 1}, 'Settlement Binary');

// -----------------------------
// 3. DISTANCE TO SETTLEMENT
// -----------------------------
var pixelSizeSettlement = smod.projection().nominalScale();

var distSettlement = settlementBinary
  .unmask(0)
  .fastDistanceTransform(256, 'pixels')
  .sqrt()
  .multiply(pixelSizeSettlement)
  .rename('Dist_to_Settlement_m')
  .clip(studyRegion);

Map.addLayer(distSettlement, {min: 0, max: 50000}, 'Distance to Settlement');

// -----------------------------
// 4. LOAD ROADS
// -----------------------------
var roads = ee.FeatureCollection("projects/hwc-esa/assets/GRIP4")
  .filterBounds(studyRegion)
  .map(function(f) {
    return f.set('road', 1);
  });

Map.addLayer(roads, {color: 'yellow'}, 'Roads');

// Rasterize roads
var roadRaster = roads.reduceToImage({
  properties: ['road'],
  reducer: ee.Reducer.first()
}).unmask(0).rename('roads');

// -----------------------------
// 5. DISTANCE TO ROAD
// -----------------------------
var distRoad = roadRaster
  .fastDistanceTransform(256, 'pixels')
  .sqrt()
  .multiply(30)
  .rename('Dist_to_Road_m')
  .clip(studyRegion);

Map.addLayer(distRoad, {min: 0, max: 50000}, 'Distance to Road');

// -----------------------------
// 6. BUILD BIAS SURFACE
// -----------------------------
var wSettlement = distSettlement.multiply(-1).divide(settlementDecay).exp();
var wRoad = distRoad.multiply(-1).divide(roadDecay).exp();

var biasSurfaceRaw = wSettlement.add(wRoad).rename('bias_raw');

var maxBias = ee.Number(
  biasSurfaceRaw.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: studyRegion,
    scale: sampleScale,
    maxPixels: 1e12,
    bestEffort: true
  }).get('bias_raw')
);

print('maxBias', maxBias);

var biasSurface = biasSurfaceRaw
  .divide(maxBias)
  .max(minBias)   // keeps some acceptance everywhere
  .rename('bias')
  .reproject({
    crs: distSettlement.projection(),
    scale: sampleScale
  })
  .clip(studyRegion);

Map.addLayer(biasSurface, {min: 0, max: 1}, 'Bias surface');

// -----------------------------
// 7. EXCLUSION AREA
// -----------------------------
var exclusionArea = pts.geometry().buffer(exclusionBufferM);

var allowedMask = ee.Image.constant(1)
  .clip(studyRegion)
  .paint(ee.FeatureCollection([ee.Feature(exclusionArea)]), 0)
  .eq(1)
  .rename('allowed');

Map.addLayer(allowedMask.selfMask(), {palette: ['green']}, 'Allowed area');

// bias only in allowed area
var biasMasked = biasSurface.updateMask(allowedMask);

Map.addLayer(biasMasked, {min: 0, max: 1}, 'Bias masked');

// -----------------------------
// 8. RASTER-BASED ACCEPTANCE
// -----------------------------
// Random raster
var randImg = ee.Image.random(randomSeed).rename('rand');

// Accept where random <= bias
var acceptedMask = randImg.lte(biasMasked).selfMask();

Map.addLayer(acceptedMask, {palette: ['blue']}, 'Accepted mask');

// -----------------------------
// 9. SAMPLE ACCEPTED PIXELS AS BACKGROUND POINTS
// -----------------------------
var samplingImage = ee.Image.pixelLonLat()
  .addBands(biasMasked)
  .addBands(distSettlement)
  .addBands(distRoad)
  .updateMask(acceptedMask);

var sampledBg = samplingImage.sample({
  region: studyRegion,
  scale: sampleScale,
  numPixels: nBackground,
  seed: randomSeed,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('Bias-corrected sampled points', sampledBg.size());

// -----------------------------
// 10. FALLBACK UNIFORM SAMPLE IF TOO FEW / ZERO
// -----------------------------
// This prevents empty exports.
var uniformFallbackImg = ee.Image.pixelLonLat()
  .addBands(ee.Image.constant(minBias).rename('bias'))
  .addBands(distSettlement)
  .addBands(distRoad)
  .updateMask(allowedMask);

var fallbackBg = uniformFallbackImg.sample({
  region: studyRegion,
  scale: sampleScale,
  numPixels: nBackground,
  seed: randomSeed + 999,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('Fallback sampled points', fallbackBg.size());

// If bias-corrected sample is zero, use fallback.
// If bias-corrected sample is non-zero but smaller than requested,
// merge with fallback and limit to nBackground.
var bgCombined = ee.FeatureCollection(
  ee.Algorithms.If(
    sampledBg.size().gt(0),
    sampledBg.merge(fallbackBg).limit(nBackground),
    fallbackBg.limit(nBackground)
  )
);

// -----------------------------
// 11. FINAL OUTPUT FEATURES
// -----------------------------
var backgroundPoints = bgCombined.map(function(f) {
  var lon = ee.Number(f.get('longitude'));
  var lat = ee.Number(f.get('latitude'));

  return ee.Feature(ee.Geometry.Point([lon, lat]), {
    Long: lon,
    Lat: lat,
    conflict: 0,
    conflict_type: label,
    bg_type: 'bias_corrected',
    bias: f.get('bias'),
    Dist_to_Settlement_m: f.get('Dist_to_Settlement_m'),
    Dist_to_Road_m: f.get('Dist_to_Road_m')
  });
});

print('Final background count', backgroundPoints.size());
Map.addLayer(backgroundPoints, {color: 'cyan'}, label + ' background');

// -----------------------------
// 12. EXPORT BACKGROUND POINTS
// -----------------------------
Export.table.toDrive({

    folder: exportFolder,
collection: backgroundPoints,
  description: label + "_BiasCorrected_Background_Raw",
    fileNamePrefix: label + "_BiasCorrected_Background_Raw",
    fileFormat: 'CSV'
});

// -----------------------------
// 13. OPTIONAL: EXPORT BIAS SURFACE
// -----------------------------
Export.image.toDrive({
  image: biasMasked,
  description: label + '_Bias_Surface_GEE_Robust',
  region: studyRegion,
  scale: sampleScale,
  maxPixels: 1e13
});

}
["HLC", "HTC"].forEach(runBackground);
