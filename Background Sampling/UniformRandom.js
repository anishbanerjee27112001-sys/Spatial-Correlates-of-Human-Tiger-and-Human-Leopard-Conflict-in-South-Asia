// UniformRandom background sampling for HLC and HTC.
var assetRoot = "projects/hwc-esa/assets/";
var exportFolder = "HCC_Background_Exports";
function runBackground(label) {
var assetName = assetRoot + label + "_Final_Conflict_Data";
// ======================================================
// ROBUST UNIFORM RANDOM BACKGROUND POINTS IN GEE
// ======================================================

// -----------------------------
// 0. SETTINGS
// -----------------------------

var nBackground = 40000;
var exclusionBufferM = 5000;
var sampleScale = 250;
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
// 6. EXCLUSION AREA
// -----------------------------
var exclusionArea = pts.geometry().buffer(exclusionBufferM);

var allowedMask = ee.Image.constant(1)
  .clip(studyRegion)
  .paint(ee.FeatureCollection([ee.Feature(exclusionArea)]), 0)
  .eq(1)
  .rename('allowed');

Map.addLayer(allowedMask.selfMask(), {palette: ['green']}, 'Allowed area');

// -----------------------------
// 7. SAMPLE UNIFORM BACKGROUND
// -----------------------------
var samplingImage = ee.Image.pixelLonLat()
  .addBands(ee.Image.constant(1).rename('bias'))
  .addBands(distSettlement)
  .addBands(distRoad)
  .updateMask(allowedMask);

var sampledBg = samplingImage.sample({
  region: studyRegion,
  scale: sampleScale,
  numPixels: nBackground,
  seed: randomSeed,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('Uniform sampled points', sampledBg.size());

// -----------------------------
// 8. FALLBACK SAMPLE IF TOO FEW / ZERO
// -----------------------------
var fallbackBg = samplingImage.sample({
  region: studyRegion,
  scale: sampleScale,
  numPixels: nBackground,
  seed: randomSeed + 999,
  geometries: true,
  dropNulls: true,
  tileScale: 4
});

print('Fallback sampled points', fallbackBg.size());

var bgCombined = ee.FeatureCollection(
  ee.Algorithms.If(
    sampledBg.size().gt(0),
    sampledBg.merge(fallbackBg).limit(nBackground),
    fallbackBg.limit(nBackground)
  )
);

// -----------------------------
// 9. FINAL OUTPUT FEATURES
// -----------------------------
var backgroundPoints = bgCombined.map(function(f) {
  var lon = ee.Number(f.get('longitude'));
  var lat = ee.Number(f.get('latitude'));

  return ee.Feature(ee.Geometry.Point([lon, lat]), {
    Long: lon,
    Lat: lat,
    conflict: 0,
    conflict_type: label,
    bg_type: 'uniform_random',
    bias: f.get('bias'),
    Dist_to_Settlement_m: f.get('Dist_to_Settlement_m'),
    Dist_to_Road_m: f.get('Dist_to_Road_m')
  });
});

print('Final background count', backgroundPoints.size());
Map.addLayer(backgroundPoints, {color: 'cyan'}, label + ' uniform background');

// -----------------------------
// 10. EXPORT BACKGROUND POINTS
// -----------------------------
Export.table.toDrive({

    folder: exportFolder,
collection: backgroundPoints,
  description: label + "_UniformRandom_Background_Raw",
    fileNamePrefix: label + "_UniformRandom_Background_Raw",
    fileFormat: 'CSV'
});

// -----------------------------
// 11. OPTIONAL: EXPORT ALLOWED MASK
// -----------------------------
Export.image.toDrive({
  image: allowedMask,
  description: label + '_UniformRandom_AllowedMask_GEE',
  region: studyRegion,
  scale: sampleScale,
  maxPixels: 1e13
});

}
["HLC", "HTC"].forEach(runBackground);
