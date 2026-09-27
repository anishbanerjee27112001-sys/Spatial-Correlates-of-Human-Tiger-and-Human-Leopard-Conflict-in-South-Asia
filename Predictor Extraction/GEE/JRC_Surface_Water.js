var assetRoot = "projects/hwc-esa/assets/";
var exportFolder = "HCC_Predictor_Exports";
var assetConfigs = [
  {assetName: assetRoot + "HLC_AccessibleArea_Merged_Final", label: "HLC_AccessibleArea_Merged_Final", yearField: "Year"},
  {assetName: assetRoot + "HLC_BiasCorrected_Merged_Final", label: "HLC_BiasCorrected_Merged_Final", yearField: "Year"},
  {assetName: assetRoot + "HLC_UniformRandom_Merged_Final", label: "HLC_UniformRandom_Merged_Final", yearField: "Year"},
  {assetName: assetRoot + "HTC_AccessibleArea_Merged_Final", label: "HTC_AccessibleArea_Merged_Final", yearField: "Year"},
  {assetName: assetRoot + "HTC_BiasCorrected_Merged_Final", label: "HTC_BiasCorrected_Merged_Final", yearField: "Year"},
  {assetName: assetRoot + "HTC_UniformRandom_Merged_Final", label: "HTC_UniformRandom_Merged_Final", yearField: "Year"}
];

function processDataset(assetName, label) {
// ==========================================================
// JRC WATER VARIABLES — POINT DISTANCE + BUFFER MEANS
// Mean occurrence + seasonality at 1, 2, 5 km buffers
// ==========================================================

var lonField = "Long";
var latField = "Lat";
var scaleUse = 30;

// 2048 pixels * 30 m ≈ 61 km
var maxDistPixels = 2048;

var region = ee.Geometry.Rectangle([60, 5, 100, 38], null, false);

// ==========================================================
// JRC GLOBAL SURFACE WATER
// occurrence: % frequency of water presence, 0–100
// seasonality: number of months water is present, 0–12
// ==========================================================

var jrc = ee.Image("JRC/GSW1_4/GlobalSurfaceWater");

var waterOccurrence = jrc
  .select("occurrence")
  .unmask(0)
  .rename("water_occurrence");

var waterSeasonality = jrc
  .select("seasonality")
  .unmask(0)
  .rename("water_seasonality");

// Any detected surface water for distance calculation
var waterBinary = jrc
  .select("occurrence")
  .gt(0)
  .rename("water");

// ==========================================================
// HELPER FUNCTIONS
// ==========================================================

function sanitizeDescription(txt) {
  return String(txt).replace(/[^A-Za-z0-9_]/g, "_");
}

function makeDistanceKm(binaryImg, outName) {
  var source = binaryImg
    .clip(region)
    .unmask(0)
    .selfMask();

  var distM = source
    .fastDistanceTransform(maxDistPixels, "pixels", "squared_euclidean")
    .sqrt()
    .multiply(scaleUse);

  return distM.divide(1000).rename(outName).clip(region);
}

function samplePoint(img, geom, bandName) {
  var v = img.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: geom,
    scale: scaleUse,
    maxPixels: 1e8,
    tileScale: 8
  }).get(bandName);

  return ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), null, v);
}

function meanInBuffer(img, geom, radiusM, bandName) {
  var v = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom.buffer(radiusM),
    scale: scaleUse,
    maxPixels: 1e9,
    tileScale: 8
  }).get(bandName);

  return ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), null, v);
}

// ==========================================================
// LOAD POINTS
// ==========================================================

var raw = ee.FeatureCollection(assetName);

var pts = raw
  .filter(ee.Filter.notNull([lonField, latField]))
  .filter(ee.Filter.neq(lonField, "NA"))
  .filter(ee.Filter.neq(latField, "NA"))
  .filter(ee.Filter.neq(lonField, ""))
  .filter(ee.Filter.neq(latField, ""))
  .map(function(f) {
    var lon = ee.Number.parse(ee.String(f.get(lonField)));
    var lat = ee.Number.parse(ee.String(f.get(latField)));

    return ee.Feature(ee.Geometry.Point([lon, lat]), f.toDictionary());
  });

// ==========================================================
// WATER DISTANCE IMAGE
// ==========================================================

var distWater = makeDistanceKm(waterBinary, "dist_water_km");

// ==========================================================
// EXTRACT WATER VARIABLES
// ==========================================================

var out = pts.map(function(f) {
  var geom = f.geometry();

  return f
    // Point distance to nearest mapped surface water
    .set("JRC_dist_to_water_km_Final",
      samplePoint(distWater, geom, "dist_water_km"))

    // Mean water occurrence, 0–100%, within buffers
    .set("mean_water_occurrence_1km",
      meanInBuffer(waterOccurrence, geom, 1000, "water_occurrence"))
    .set("mean_water_occurrence_2km",
      meanInBuffer(waterOccurrence, geom, 2000, "water_occurrence"))
    .set("mean_water_occurrence_5km",
      meanInBuffer(waterOccurrence, geom, 5000, "water_occurrence"))

    // Mean water seasonality, 0–12 months, within buffers
    .set("mean_water_seasonality_1km",
      meanInBuffer(waterSeasonality, geom, 1000, "water_seasonality"))
    .set("mean_water_seasonality_2km",
      meanInBuffer(waterSeasonality, geom, 2000, "water_seasonality"))
    .set("mean_water_seasonality_5km",
      meanInBuffer(waterSeasonality, geom, 5000, "water_seasonality"));
});

// ==========================================================
// EXPORT
// ==========================================================

Export.table.toDrive({

    folder: exportFolder,
collection: out,
  description: sanitizeDescription(label + "_JRC_Water_Distance_Occurrence_Seasonality"),
    fileNamePrefix: sanitizeDescription(label + "_JRC_Water_Distance_Occurrence_Seasonality"),
    fileFormat: "CSV"
});

print("Prepared JRC water variables:", label);
print("Input rows:", raw.size());
print("Output rows:", out.size());

}
assetConfigs.forEach(function(config) {
  processDataset(config.assetName, config.label);
});
