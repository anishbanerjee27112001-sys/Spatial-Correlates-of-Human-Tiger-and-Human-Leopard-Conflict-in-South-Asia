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
// GHSL SETTLEMENT DISTANCE — ALL POINTS, ONE OUTPUT FILE
// Native GHSL resolution: 100 m
// Temporal matching:
//   2000–2007  -> GHSL 2000
//   2008 onward -> GHSL 2015
// Pre-2000 rows retained with GHSL variables = null
// ==========================================================

// ==========================================================
// 0. SETTINGS
// ==========================================================

var lonField = "Long";
var latField = "Lat";
var yearField = "Year";

var scaleUse = 100;
var searchRadiusKm = 50;
var searchRadiusPixels = 500;

var crsUse = "EPSG:3857";

// ==========================================================
// 1. LOAD AND PREPARE GHSL SETTLEMENT LAYERS
// ==========================================================

var ghsl2000 = ee.Image("JRC/GHSL/P2023A/GHS_BUILT_S/2000")
  .select(0)
  .gt(0)
  .rename("settlement");

var ghsl2015 = ee.Image("JRC/GHSL/P2023A/GHS_BUILT_S/2015")
  .select(0)
  .gt(0)
  .rename("settlement");

// ==========================================================
// 2. HELPER FUNCTIONS
// ==========================================================

function sanitizeDescription(txt) {
  return String(txt).replace(/[^A-Za-z0-9_]/g, "_");
}

function getGhslYearUsed(yearNum) {
  yearNum = ee.Number(yearNum);

  return ee.Number(
    ee.Algorithms.If(
      yearNum.lte(2007),
      2000,
      2015
    )
  );
}

function getGhslByYear(yearNum) {
  var ghslYear = getGhslYearUsed(yearNum);

  return ee.Image(
    ee.Algorithms.If(
      ghslYear.eq(2000),
      ghsl2000,
      ghsl2015
    )
  );
}

function prepFeature(f) {
  var lon = ee.Number.parse(ee.String(f.get(lonField)));
  var lat = ee.Number.parse(ee.String(f.get(latField)));
  var year = ee.Number.parse(
    ee.String(f.get(yearField)).slice(0, 4)
  );

  return ee.Feature(
    ee.Geometry.Point([lon, lat]),
    f.toDictionary()
  ).set("parsed_year", year);
}

// ==========================================================
// 3. LOAD AND CLEAN INPUT POINTS
// ==========================================================

var raw = ee.FeatureCollection(assetName);

var pts = raw
  .filter(ee.Filter.notNull([lonField, latField, yearField]))
  .filter(ee.Filter.neq(lonField, "NA"))
  .filter(ee.Filter.neq(latField, "NA"))
  .filter(ee.Filter.neq(yearField, "NA"))
  .filter(ee.Filter.neq(lonField, ""))
  .filter(ee.Filter.neq(latField, ""))
  .filter(ee.Filter.neq(yearField, ""))
  .map(prepFeature);

print("Total valid points:", pts.size());

// ==========================================================
// 4. CALCULATE DISTANCE TO SETTLEMENT
// ==========================================================

var out = pts.map(function(f) {
  var geom = f.geometry();
  var yearNum = ee.Number(f.get("parsed_year"));
  var isPre2000 = yearNum.lt(2000);

  var ghslYear = getGhslYearUsed(yearNum);

  var settlement = getGhslByYear(yearNum)
    .unmask(0);

  var localRegion = geom.buffer(
    searchRadiusKm * 1000 + 10000
  );

  var settlementLocal = settlement
    .clip(localRegion)
    .reproject({
      crs: crsUse,
      scale: scaleUse
    });

  var distPixels = settlementLocal.distance(
    ee.Kernel.euclidean(
      searchRadiusPixels,
      "pixels"
    ),
    false
  );

  var distKm = distPixels
    .multiply(scaleUse)
    .divide(1000)
    .rename("dist_to_settlement_km");

  var pointSettlement = settlement
    .reproject({
      crs: crsUse,
      scale: scaleUse
    })
    .rename("point_settlement");

  var values = ee.Image.cat([
    distKm,
    pointSettlement
  ]).reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: geom,
    scale: scaleUse,
    crs: crsUse,
    maxPixels: 1e9
  });

  return f
    .set(
      "ghsl_year_used_FINAL",
      ee.Algorithms.If(
        isPre2000,
        null,
        ghslYear
      )
    )
    .set(
      "dist_to_settlement_km_FINAL",
      ee.Algorithms.If(
        isPre2000,
        null,
        values.get("dist_to_settlement_km")
      )
    )
    .set(
      "point_settlement_FINAL",
      ee.Algorithms.If(
        isPre2000,
        null,
        values.get("point_settlement")
      )
    )
    .set(
      "settlement_search_radius_km_FINAL",
      ee.Algorithms.If(
        isPre2000,
        null,
        searchRadiusKm
      )
    );
});

// ==========================================================
// 5. EXPORT
// ==========================================================

Export.table.toDrive({

    folder: exportFolder,
collection: out,
  description: sanitizeDescription(
    label + "_SETTLEMENT_DISTANCE_FINAL_ALL_POINTS"
  ),
    fileNamePrefix: sanitizeDescription(
    label + "_SETTLEMENT_DISTANCE_FINAL_ALL_POINTS"
  ),
    fileFormat: "CSV"
});

print("Prepared settlement-distance export:", label);

}
assetConfigs.forEach(function(config) {
  processDataset(config.assetName, config.label);
});
