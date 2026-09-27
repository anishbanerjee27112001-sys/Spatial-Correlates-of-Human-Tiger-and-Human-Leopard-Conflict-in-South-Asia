// ==========================================================
// LANDSAT ANNUAL NDVI EXTRACTION
// NDVI derived from:
// Landsat 5 TM SR: 1984–2012
// Landsat 7 ETM+ SR: 1999–present
// Landsat 8 OLI SR: 2013–present
//
// Outputs:
// 1) NDVI
// 2) NDVI_Year
//
// Notes:
// - Uses the standardized Year field for temporal matching
// - Filters imagery Jan 1–Dec 31 of each point year
// - Years < 1985 return -999
// - Cloud/shadow masking using QA_PIXEL
// ==========================================================

// ===============================
// 0. SETTINGS
// ===============================

var lonField = "Long";
var latField = "Lat";
var scaleUse = 30;

// ===============================
// 1. INPUT FILES
// ===============================
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

// ===============================
// 2. HELPERS
// ===============================
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

function parseYear(f, yearField) {
  return ee.Number.parse(ee.String(f.get(yearField)).slice(0, 4));
}

// ===============================
// 3. LANDSAT CLOUD MASKS
// ===============================

// QA_PIXEL bits:
// bit 1 = dilated cloud
// bit 2 = cirrus
// bit 3 = cloud
// bit 4 = cloud shadow
// bit 5 = snow

function maskLandsatSR(img) {
  var qa = img.select("QA_PIXEL");

  var mask = qa.bitwiseAnd(1 << 1).eq(0)   // dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))      // cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))      // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));     // snow

  return img.updateMask(mask);
}

// ===============================
// 4. SENSOR-SPECIFIC NDVI FUNCTIONS
// ===============================

// Landsat Collection 2 Level 2 scale factor:
// Reflectance = DN * 0.0000275 - 0.2

function addNdviL57(img) {
  var sr = img.select(["SR_B3", "SR_B4"])
    .multiply(0.0000275)
    .add(-0.2);

  var red = sr.select("SR_B3");
  var nir = sr.select("SR_B4");

  var ndvi = nir.subtract(red)
    .divide(nir.add(red))
    .rename("NDVI");

  return img.addBands(ndvi);
}

function addNdviL8(img) {
  var sr = img.select(["SR_B4", "SR_B5"])
    .multiply(0.0000275)
    .add(-0.2);

  var red = sr.select("SR_B4");
  var nir = sr.select("SR_B5");

  var ndvi = nir.subtract(red)
    .divide(nir.add(red))
    .rename("NDVI");

  return img.addBands(ndvi);
}

// ===============================
// 5. BUILD ANNUAL NDVI IMAGE
// ===============================
function getAnnualNdvi(year, geom) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, 1, 1);
  var end   = ee.Date.fromYMD(year, 12, 31).advance(1, "day");

  var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
    .filterDate(start, end)
    .filterBounds(geom)
    .filter(ee.Filter.lt("CLOUD_COVER", 60))
    .map(maskLandsatSR)
    .map(addNdviL57)
    .select("NDVI");

  var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
    .filterDate(start, end)
    .filterBounds(geom)
    .filter(ee.Filter.lt("CLOUD_COVER", 60))
    .map(maskLandsatSR)
    .map(addNdviL57)
    .select("NDVI");

  var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    .filterDate(start, end)
    .filterBounds(geom)
    .filter(ee.Filter.lt("CLOUD_COVER", 60))
    .map(maskLandsatSR)
    .map(addNdviL8)
    .select("NDVI");

  // Use available sensors depending on year
  var collection = ee.ImageCollection(
    ee.Algorithms.If(
      year.lt(1999),
      l5,
      ee.Algorithms.If(
        year.lt(2013),
        l5.merge(l7),
        l7.merge(l8)
      )
    )
  );

  return collection.median().rename("NDVI");
}

// ===============================
// 6. CORE PROCESS FUNCTION
// ===============================
function processAsset(assetName, label, yearField) {

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
      var yr = parseYear(f, yearField);

      return ee.Feature(
        ee.Geometry.Point([lon, lat]),
        f.toDictionary()
      )
        .set("NDVI_Year", yr);
    });

  var out = pts.map(function(f) {
    var geom = f.geometry();
    var yr = ee.Number(f.get("NDVI_Year"));

    var ndviValue = ee.Algorithms.If(
      yr.lt(1985),
      -999,
      ee.Algorithms.If(
        yr.gt(2023),
        -999,
        (function() {
          var ndviImg = getAnnualNdvi(yr, geom);

          var stats = ndviImg.reduceRegion({
            reducer: ee.Reducer.first(),
            geometry: geom,
            scale: scaleUse,
            maxPixels: 1e8
          });

          return safeNum(stats.get("NDVI"), -999);
        })()
      )
    );

    return f
      .set("NDVI", ndviValue)
      .set("NDVI_Year", yr);
  });

  Export.table.toDrive({

    folder: exportFolder,
collection: out,
    description: sanitizeDescription(label + "_Landsat_Annual_NDVI"),
    fileNamePrefix: sanitizeDescription(label + "_Landsat_Annual_NDVI"),
    fileFormat: "CSV"
  });

  print("Prepared NDVI export for:", label, "Count:", out.size());
}

// ===============================
// 7. RUN ALL FILES
// ===============================
for (var i = 0; i < assetConfigs.length; i++) {
  processAsset(
    assetConfigs[i].assetName,
    assetConfigs[i].label,
    assetConfigs[i].yearField
  );
}
