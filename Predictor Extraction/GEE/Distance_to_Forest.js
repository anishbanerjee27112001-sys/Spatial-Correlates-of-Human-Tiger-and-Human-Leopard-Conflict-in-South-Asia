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
// GLAD POINT DISTANCE TO FOREST ONLY — MEMORY-SAFE VERSION
// ==========================================================

var lonField = "Long";
var latField = "Lat";
var yearField = "Year";
var scaleUse = 30;

// 2048 pixels * 30 m ≈ 61 km. Increase to 4096 only if needed.
var maxDistPixels = 2048;

var landmask = ee.Image("projects/glad/OceanMask").lte(1);

var glad2000 = ee.Image("projects/glad/GLCLU2020/v2/LCLUC_2000").updateMask(landmask).rename("GLAD_LC");
var glad2005 = ee.Image("projects/glad/GLCLU2020/v2/LCLUC_2005").updateMask(landmask).rename("GLAD_LC");
var glad2010 = ee.Image("projects/glad/GLCLU2020/v2/LCLUC_2010").updateMask(landmask).rename("GLAD_LC");
var glad2015 = ee.Image("projects/glad/GLCLU2020/v2/LCLUC_2015").updateMask(landmask).rename("GLAD_LC");
var glad2020 = ee.Image("projects/glad/GLCLU2020/v2/LCLUC_2020").updateMask(landmask).rename("GLAD_LC");

var forestCodes = [5,6,7,8,9,10,11,12,13];

// Simple South Asia box. Much cheaper than pts.geometry().buffer().
var region = ee.Geometry.Rectangle([60, 5, 100, 38], null, false);

function sanitizeDescription(txt) {
  return String(txt).replace(/[^A-Za-z0-9_]/g, "_");
}

function assignGladYear(y) {
  y = ee.Number(y);
  return ee.Number(
    ee.Algorithms.If(y.lt(2000), -999,
    ee.Algorithms.If(y.lte(2004), 2000,
    ee.Algorithms.If(y.lte(2009), 2005,
    ee.Algorithms.If(y.lte(2014), 2010,
    ee.Algorithms.If(y.lte(2019), 2015, 2020)))))
  );
}

function getGladImage(y) {
  y = ee.Number(y);
  return ee.Image(
    ee.Algorithms.If(y.eq(2000), glad2000,
    ee.Algorithms.If(y.eq(2005), glad2005,
    ee.Algorithms.If(y.eq(2010), glad2010,
    ee.Algorithms.If(y.eq(2015), glad2015, glad2020))))
  );
}

function maskFromCodes(img, codes, name) {
  return img.remap(codes, ee.List.repeat(1, codes.length), 0).rename(name);
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

var raw = ee.FeatureCollection(assetName);

var pts = raw
  .filter(ee.Filter.notNull([lonField, latField, yearField]))
  .filter(ee.Filter.neq(lonField, "NA"))
  .filter(ee.Filter.neq(latField, "NA"))
  .filter(ee.Filter.neq(yearField, "NA"))
  .filter(ee.Filter.neq(lonField, ""))
  .filter(ee.Filter.neq(latField, ""))
  .filter(ee.Filter.neq(yearField, ""))
  .map(function(f) {
    var lon = ee.Number.parse(ee.String(f.get(lonField)));
    var lat = ee.Number.parse(ee.String(f.get(latField)));
    var yr = ee.Number.parse(ee.String(f.get(yearField)).slice(0, 4));

    return ee.Feature(ee.Geometry.Point([lon, lat]), f.toDictionary())
      .set("parsed_year", yr)
      .set("GLAD_year_used_internal", assignGladYear(yr));
  });

var years = ee.List([2000, 2005, 2010, 2015, 2020]);

var outputs = years.map(function(y) {
  y = ee.Number(y);
  var ptsYear = pts.filter(ee.Filter.eq("GLAD_year_used_internal", y));

  return ee.FeatureCollection(
    ee.Algorithms.If(
      ptsYear.size().gt(0),
      (function() {
        var lc = getGladImage(y);
        var forest = maskFromCodes(lc, forestCodes, "forest");
        var distForest = makeDistanceKm(forest, "dist_forest_km");

        return ptsYear.map(function(f) {
          return f
            .set("GLAD_year_used_Final", y)
            .set(
              "GLAD_dist_to_forest_km_Final",
              samplePoint(distForest, f.geometry(), "dist_forest_km")
            );
        });
      })(),
      ee.FeatureCollection([])
    )
  );
});

var pre2000 = pts
  .filter(ee.Filter.eq("GLAD_year_used_internal", -999))
  .map(function(f) {
    return f
      .set("GLAD_year_used_Final", null)
      .set("GLAD_dist_to_forest_km_Final", null);
  });

var out = pre2000.merge(ee.FeatureCollection(outputs).flatten());

Export.table.toDrive({

    folder: exportFolder,
collection: out,
  description: sanitizeDescription(label + "_GLAD_Forest_PointDistance_SAFE"),
    fileNamePrefix: sanitizeDescription(label + "_GLAD_Forest_PointDistance_SAFE"),
    fileFormat: "CSV"
});

print("Prepared forest distance:", label);
print("Input rows:", raw.size());
print("Output rows:", out.size());

}
assetConfigs.forEach(function(config) {
  processDataset(config.assetName, config.label);
});
