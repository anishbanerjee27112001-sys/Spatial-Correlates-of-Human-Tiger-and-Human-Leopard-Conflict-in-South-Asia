// ======================================================
// ROAD DENSITY ONLY
// Road length + road density at 1 km, 2 km, and 5 km
// ======================================================

// =============================================
// 0. SETTINGS
// =============================================

var lonField = "Long";
var latField = "Lat";

var densityBuffer1M = 1000;
var densityBuffer2M = 2000;
var densityBuffer5M = 5000;

var roadAsset = "projects/hwc-esa/assets/GRIP4";

// =============================================
// 1. INPUT FILES
// =============================================
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

// =============================================
// 2. HELPERS
// =============================================
function sanitizeDescription(txt) {
  return String(txt).replace(/[^A-Za-z0-9_]/g, "_");
}

function prepFeature(f, yearField) {
  var lon = ee.Number.parse(ee.String(f.get(lonField)));
  var lat = ee.Number.parse(ee.String(f.get(latField)));
  var geom = ee.Geometry.Point([lon, lat]);

  var year = ee.Number.parse(
    ee.String(f.get(yearField)).slice(0, 4)
  );

  return f
    .setGeometry(geom)
    .set("Conflict_Year", year);
}

function getRoadLengthInBuffer(roads, buf) {
  var roadLenM = ee.Number(
    roads.filterBounds(buf).map(function(r) {
      var clipped = r.geometry().intersection(buf, ee.ErrorMargin(1));
      return ee.Feature(null, {
        len_m: clipped.length(ee.ErrorMargin(1))
      });
    }).aggregate_sum("len_m")
  );

  return ee.Number(
    ee.Algorithms.If(roadLenM, roadLenM, 0)
  );
}

function getRoadDensity(roadLenM, buf) {
  var bufAreaKm2 = buf.area(1).divide(1e6);

  return ee.Number(
    ee.Algorithms.If(
      bufAreaKm2.gt(0),
      roadLenM.divide(1000).divide(bufAreaKm2),
      0
    )
  );
}

// =============================================
// 3. LOAD ROADS
// =============================================
var roadsRaw = ee.FeatureCollection(roadAsset);

// =============================================
// 4. PROCESS EACH ASSET
// =============================================
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
      return prepFeature(f, yearField);
    });

  var studyRegion = pts.geometry().buffer(50000);
  var roads = roadsRaw.filterBounds(studyRegion);

  var out = pts.map(function(f) {
    var geom = f.geometry();

    var buf1 = geom.buffer(densityBuffer1M);
    var buf2 = geom.buffer(densityBuffer2M);
    var buf5 = geom.buffer(densityBuffer5M);

    var roadLen1M = getRoadLengthInBuffer(roads, buf1);
    var roadDensity1 = getRoadDensity(roadLen1M, buf1);

    var roadLen2M = getRoadLengthInBuffer(roads, buf2);
    var roadDensity2 = getRoadDensity(roadLen2M, buf2);

    var roadLen5M = getRoadLengthInBuffer(roads, buf5);
    var roadDensity5 = getRoadDensity(roadLen5M, buf5);

    return f
      .set("road_len_1000m_buf_m", roadLen1M)
      .set("road_density_1000m_buf_km_per_km2", roadDensity1)
      .set("road_len_2000m_buf_m", roadLen2M)
      .set("road_density_2000m_buf_km_per_km2", roadDensity2)
      .set("road_len_5000m_buf_m", roadLen5M)
      .set("road_density_5000m_buf_km_per_km2", roadDensity5);
  });

  Export.table.toDrive({

    folder: exportFolder,
collection: out,
    description: sanitizeDescription(
      label + "_Road_Density_GRIP4_1km_2km_5km"
    ),
    fileNamePrefix: sanitizeDescription(
      label + "_Road_Density_GRIP4_1km_2km_5km"
    ),
    fileFormat: "CSV"
  });

  print("Prepared road density export for:", label);
}

// =============================================
// 5. RUN ALL 6 FILES
// =============================================
for (var i = 0; i < assetConfigs.length; i++) {
  processAsset(
    assetConfigs[i].assetName,
    assetConfigs[i].label,
    assetConfigs[i].yearField
  );
}
