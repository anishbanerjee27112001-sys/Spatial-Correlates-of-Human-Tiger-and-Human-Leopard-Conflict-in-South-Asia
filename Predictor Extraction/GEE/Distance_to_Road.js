// ======================================================
// DISTANCE TO ROAD ONLY — MEMORY-SAFE VECTOR VERSION
// GRIP4 vector roads
// ======================================================

// =============================================
// 0. SETTINGS
// =============================================

var lonField = "Long";
var latField = "Lat";
var yearField = "Year";

var roadAsset = "projects/hwc-esa/assets/GRIP4";

var searchRadiusM = 61440; // 61.44 km
var maxErrorM = 30;

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

function prepFeature(f) {
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

// =============================================
// 3. LOAD ROADS
// =============================================
var roadsRaw = ee.FeatureCollection(roadAsset);

// =============================================
// 4. PROCESS EACH ASSET
// =============================================
function processAsset(assetName, label) {

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

  var studyRegion = pts.geometry().buffer(searchRadiusM);
  var roads = roadsRaw.filterBounds(studyRegion);

  var out = pts.map(function(f) {
    var geom = f.geometry();
    var searchArea = geom.buffer(searchRadiusM);

    var nearbyRoads = roads.filterBounds(searchArea);

    var roadDistances = nearbyRoads.map(function(r) {
      var d = geom.distance(r.geometry(), maxErrorM);
      return r.set("road_dist_m", d);
    });

    var minDistM = ee.Number(
      roadDistances.aggregate_min("road_dist_m")
    );

    var distRoadKm = ee.Algorithms.If(
      minDistM,
      minDistM.divide(1000),
      null
    );

    return f
      .set("road_search_radius_km_FINAL", searchRadiusM / 1000)
      .set("dist_to_road_km_FINAL", distRoadKm);
  });

  Export.table.toDrive({

    folder: exportFolder,
collection: out,
    description: sanitizeDescription(
      label + "_Distance_to_Road_GRIP4_Vector_FINAL"
    ),
    fileNamePrefix: sanitizeDescription(
      label + "_Distance_to_Road_GRIP4_Vector_FINAL"
    ),
    fileFormat: "CSV"
  });

  print("Prepared vector distance-to-road export for:", label);
}

// =============================================
// 5. RUN ALL 6 FILES
// =============================================
for (var i = 0; i < assetConfigs.length; i++) {
  processAsset(
    assetConfigs[i].assetName,
    assetConfigs[i].label
  );
}
