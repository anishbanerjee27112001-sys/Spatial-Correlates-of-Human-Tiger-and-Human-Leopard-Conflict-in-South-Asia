// ==========================================================
// HANSEN / UMD GLOBAL FOREST CHANGE v1.12
// Uses only GFC-provided variables:
// treecover2000 and lossyear
//
// Outputs:
// 1. pct_hansen_baseline_treecover2000_*
//    = mean tree canopy cover (%) in 2000 within buffer
//
// 2. pct_hansen_treecover_loss_area_*
//    = % of total buffer area with lossyear > 0 and <= observation year
//
// 3. pct_hansen_baseline_treecover_lost_*
//    = cumulative lost tree cover relative to baseline tree cover 2000
// ==========================================================

var lonField = "Long";
var latField = "Lat";
var scaleUse = 30;

var buf1 = 1000;
var buf2 = 2000;
var buf5 = 5000;

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

var gfc = ee.Image("UMD/hansen/global_forest_change_2024_v1_12");

var treecover2000 = gfc.select("treecover2000").rename("treecover2000");
var lossyear = gfc.select("lossyear");
var pixelArea = ee.Image.pixelArea().rename("pixel_area");

// Tree-cover-equivalent area in 2000.
// Example: a 900 m2 pixel with 60% tree cover contributes 540 m2.
var baselineTreecoverArea2000 = treecover2000
  .divide(100)
  .multiply(pixelArea)
  .rename("baseline_treecover_area_2000");

function sanitizeDescription(txt) {
  return String(txt).replace(/[^A-Za-z0-9_]/g, "_");
}

function safeNum(val, defaultVal) {
  return ee.Number(
    ee.Algorithms.If(
      ee.Algorithms.IsEqual(val, null),
      defaultVal,
      val
    )
  );
}

function meanInBuffer(img, geom, radiusM) {
  return img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom.buffer(radiusM),
    scale: scaleUse,
    maxPixels: 1e9,
    tileScale: 4
  });
}

function sumInBuffer(img, geom, radiusM) {
  return img.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom.buffer(radiusM),
    scale: scaleUse,
    maxPixels: 1e9,
    tileScale: 4
  });
}

function setNullHansen(f) {
  return f.set({
    Hansen_Status: "pre_2001_no_loss_data",
    Hansen_End_Year_used: null,

    pct_hansen_baseline_treecover2000_1km: null,
    pct_hansen_baseline_treecover2000_2km: null,
    pct_hansen_baseline_treecover2000_5km: null,

    pct_hansen_treecover_loss_area_1km: null,
    pct_hansen_treecover_loss_area_2km: null,
    pct_hansen_treecover_loss_area_5km: null,

    pct_hansen_baseline_treecover_lost_1km: null,
    pct_hansen_baseline_treecover_lost_2km: null,
    pct_hansen_baseline_treecover_lost_5km: null
  });
}

function processHansenFeature(f) {
  var geom = f.geometry();
  var yr = ee.Number(f.get("Input_Year"));

  // Hansen v1.12 ends in 2024.
  var endYear = ee.Number(ee.Algorithms.If(yr.gt(2024), 2024, yr));
  var lossCode = endYear.subtract(2000); // 2001 = 1, ..., 2024 = 24

  // GFC cumulative loss up to observation year.
  var lossToYear = lossyear
    .gt(0)
    .and(lossyear.lte(lossCode))
    .rename("loss_to_year");

  // Area of pixels with GFC loss.
  var lossArea = lossToYear
    .multiply(pixelArea)
    .rename("loss_area");

  // Baseline tree-cover-equivalent area that experienced loss.
  var baselineTreecoverLostArea = baselineTreecoverArea2000
    .updateMask(lossToYear)
    .rename("baseline_treecover_lost_area");

  // Mean baseline tree cover 2000 (%)
  var tc1 = safeNum(meanInBuffer(treecover2000, geom, buf1).get("treecover2000"), null);
  var tc2 = safeNum(meanInBuffer(treecover2000, geom, buf2).get("treecover2000"), null);
  var tc5 = safeNum(meanInBuffer(treecover2000, geom, buf5).get("treecover2000"), null);

  // Total buffer area
  var area1 = safeNum(sumInBuffer(pixelArea, geom, buf1).get("pixel_area"), 0);
  var area2 = safeNum(sumInBuffer(pixelArea, geom, buf2).get("pixel_area"), 0);
  var area5 = safeNum(sumInBuffer(pixelArea, geom, buf5).get("pixel_area"), 0);

  // Lost area
  var lostArea1 = safeNum(sumInBuffer(lossArea, geom, buf1).get("loss_area"), 0);
  var lostArea2 = safeNum(sumInBuffer(lossArea, geom, buf2).get("loss_area"), 0);
  var lostArea5 = safeNum(sumInBuffer(lossArea, geom, buf5).get("loss_area"), 0);

  var pctLossArea1 = ee.Number(ee.Algorithms.If(area1.gt(0), lostArea1.divide(area1).multiply(100), 0));
  var pctLossArea2 = ee.Number(ee.Algorithms.If(area2.gt(0), lostArea2.divide(area2).multiply(100), 0));
  var pctLossArea5 = ee.Number(ee.Algorithms.If(area5.gt(0), lostArea5.divide(area5).multiply(100), 0));

  // Baseline tree-cover-equivalent denominator
  var baseTC1 = safeNum(sumInBuffer(baselineTreecoverArea2000, geom, buf1).get("baseline_treecover_area_2000"), 0);
  var baseTC2 = safeNum(sumInBuffer(baselineTreecoverArea2000, geom, buf2).get("baseline_treecover_area_2000"), 0);
  var baseTC5 = safeNum(sumInBuffer(baselineTreecoverArea2000, geom, buf5).get("baseline_treecover_area_2000"), 0);

  var baseLost1 = safeNum(sumInBuffer(baselineTreecoverLostArea, geom, buf1).get("baseline_treecover_lost_area"), 0);
  var baseLost2 = safeNum(sumInBuffer(baselineTreecoverLostArea, geom, buf2).get("baseline_treecover_lost_area"), 0);
  var baseLost5 = safeNum(sumInBuffer(baselineTreecoverLostArea, geom, buf5).get("baseline_treecover_lost_area"), 0);

  var pctBaseLost1 = ee.Number(ee.Algorithms.If(baseTC1.gt(0), baseLost1.divide(baseTC1).multiply(100), 0));
  var pctBaseLost2 = ee.Number(ee.Algorithms.If(baseTC2.gt(0), baseLost2.divide(baseTC2).multiply(100), 0));
  var pctBaseLost5 = ee.Number(ee.Algorithms.If(baseTC5.gt(0), baseLost5.divide(baseTC5).multiply(100), 0));

  return f.set({
    Hansen_Status: "extracted",
    Hansen_End_Year_used: endYear,

    pct_hansen_baseline_treecover2000_1km: tc1,
    pct_hansen_baseline_treecover2000_2km: tc2,
    pct_hansen_baseline_treecover2000_5km: tc5,

    pct_hansen_treecover_loss_area_1km: pctLossArea1,
    pct_hansen_treecover_loss_area_2km: pctLossArea2,
    pct_hansen_treecover_loss_area_5km: pctLossArea5,

    pct_hansen_baseline_treecover_lost_1km: pctBaseLost1,
    pct_hansen_baseline_treecover_lost_2km: pctBaseLost2,
    pct_hansen_baseline_treecover_lost_5km: pctBaseLost5
  });
}

function processAsset(assetName, label, yearField) {
  var raw = ee.FeatureCollection(assetName);

  var pts = raw
    .filter(ee.Filter.notNull([lonField, latField, yearField]))
    .map(function(f) {
      var lon = ee.Number.parse(ee.String(f.get(lonField)));
      var lat = ee.Number.parse(ee.String(f.get(latField)));
      var yr = ee.Number.parse(ee.String(f.get(yearField)).slice(0, 4));

      return ee.Feature(ee.Geometry.Point([lon, lat]), f.toDictionary())
        .set("Input_Year", yr);
    });

  var out = pts.map(function(f) {
    var yr = ee.Number(f.get("Input_Year"));
    return ee.Feature(
      ee.Algorithms.If(
        yr.lt(2001),
        setNullHansen(f),
        processHansenFeature(f)
      )
    );
  });

  Export.table.toDrive({

    folder: exportFolder,
collection: out,
    description: sanitizeDescription(label + "_HansenGFC_TreeCover_Loss"),
    fileNamePrefix: sanitizeDescription(label + "_HansenGFC_TreeCover_Loss"),
    fileFormat: "CSV"
  });

  print("Prepared Hansen export:", label, out.size());
}

for (var i = 0; i < assetConfigs.length; i++) {
  processAsset(
    assetConfigs[i].assetName,
    assetConfigs[i].label,
    assetConfigs[i].yearField
  );
}
