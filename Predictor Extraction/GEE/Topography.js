// ==========================================================
// TOPOGRAPHY EXTRACTION
// COPERNICUS DEM GLO-30
// Outputs:
// 1) Point elevation
// 2) Point slope
// 3) Point aspect
// 4) Mean elevation within 1, 2 and 5 km
// 5) Mean slope within 1, 2 and 5 km
// 6) Circular mean aspect within 1, 2 and 5 km
//
// IMPORTANT:
// Aspect is NOT averaged directly.
// Circular mean aspect is calculated from mean sine and cosine
// components and converted back to compass degrees.
// ==========================================================

// ===============================
// 0. SETTINGS
// ===============================

var lonField = "Long";
var latField = "Lat";

var scaleUse = 30;

var buf1 = 1000;
var buf2 = 2000;
var buf5 = 5000;

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
// 2. COPERNICUS DEM
// ===============================

// GLO-30 is supplied as an ImageCollection.
var copDemCol =
  ee.ImageCollection("COPERNICUS/DEM/GLO30")
    .select("DEM");

// Preserve native projection before terrain derivation.
var nativeProj =
  ee.Image(copDemCol.first())
    .projection();

// Mosaic DEM tiles and restore native projection.
var dem =
  copDemCol
    .mosaic()
    .rename("elevation")
    .setDefaultProjection(nativeProj);

// ===============================
// 3. TERRAIN DERIVATIVES
// ===============================

// Derive slope in degrees.
var slope =
  ee.Terrain
    .slope(dem)
    .rename("slope");

// Derive aspect in compass degrees:
// 0 = north
// 90 = east
// 180 = south
// 270 = west
var aspect =
  ee.Terrain
    .aspect(dem)
    .rename("aspect");

// ===============================
// 4. CIRCULAR-ASPECT COMPONENTS
// ===============================

// Convert aspect from degrees to radians.
var aspectRad =
  aspect
    .multiply(Math.PI / 180);

// Pixel-level sine and cosine components.
var aspectSin =
  aspectRad
    .sin()
    .rename("aspect_sin");

var aspectCos =
  aspectRad
    .cos()
    .rename("aspect_cos");

// Combined stack used for extraction.
var topoStack =
  dem
    .addBands(slope)
    .addBands(aspect)
    .addBands(aspectSin)
    .addBands(aspectCos);

// ===============================
// 5. HELPERS
// ===============================

function safeNum(val, defaultVal) {

  return ee.Number(
    ee.Algorithms.If(
      ee.Algorithms.IsEqual(
        val,
        null
      ),
      defaultVal,
      val
    )
  );
}

function sanitizeDescription(txt) {

  return String(txt)
    .replace(
      /[^A-Za-z0-9_]/g,
      "_"
    );
}

// Calculate mean values within a buffer.
function getBufferStats(
  img,
  geom,
  radiusM
) {

  return img.reduceRegion({

    reducer:
      ee.Reducer.mean(),

    geometry:
      geom.buffer(radiusM),

    scale:
      scaleUse,

    maxPixels:
      1e9,

    tileScale:
      4
  });
}

// ===============================
// 6. CIRCULAR MEAN ASPECT
// ===============================

function circularMeanAspect(
  meanSin,
  meanCos
) {

  meanSin =
    ee.Number(meanSin);

  meanCos =
    ee.Number(meanCos);

  // atan2 returns radians between -pi and +pi.
  var angleRad =
    meanSin.atan2(meanCos);

  // Convert radians to degrees.
  var angleDeg =
    angleRad
      .multiply(180 / Math.PI);

  // Convert negative angles to compass range 0–360°.
  var compassDeg =
    angleDeg
      .add(360)
      .mod(360);

  return compassDeg;
}

// ===============================
// 7. CORE FUNCTION
// ===============================

function processAsset(
  assetName,
  label
) {

  var raw =
    ee.FeatureCollection(
      assetName
    );

  // -------------------------------
  // Clean coordinates
  // -------------------------------

  var pts =
    raw

      .filter(
        ee.Filter.notNull(
          [lonField, latField]
        )
      )

      .filter(
        ee.Filter.neq(
          lonField,
          "NA"
        )
      )

      .filter(
        ee.Filter.neq(
          latField,
          "NA"
        )
      )

      .filter(
        ee.Filter.neq(
          lonField,
          ""
        )
      )

      .filter(
        ee.Filter.neq(
          latField,
          ""
        )
      )

      .map(function(f) {

        var lon =
          ee.Number.parse(
            ee.String(
              f.get(lonField)
            )
          );

        var lat =
          ee.Number.parse(
            ee.String(
              f.get(latField)
            )
          );

        return ee.Feature(
          ee.Geometry.Point(
            [lon, lat]
          ),
          f.toDictionary()
        );
      });

  // -------------------------------
  // Extract topographic variables
  // -------------------------------

  var out =
    pts.map(function(f) {

      var geom =
        f.geometry();

      // ===========================
      // POINT VALUES
      // ===========================

      var pointStats =
        topoStack.reduceRegion({

          reducer:
            ee.Reducer.first(),

          geometry:
            geom,

          scale:
            scaleUse,

          maxPixels:
            1e8
        });

      // ===========================
      // BUFFER VALUES
      // ===========================

      var b1 =
        getBufferStats(
          topoStack,
          geom,
          buf1
        );

      var b2 =
        getBufferStats(
          topoStack,
          geom,
          buf2
        );

      var b5 =
        getBufferStats(
          topoStack,
          geom,
          buf5
        );

      // ===========================
      // CIRCULAR MEAN ASPECT
      // ===========================

      var aspect1 =
        circularMeanAspect(
          safeNum(
            b1.get("aspect_sin"),
            0
          ),
          safeNum(
            b1.get("aspect_cos"),
            0
          )
        );

      var aspect2 =
        circularMeanAspect(
          safeNum(
            b2.get("aspect_sin"),
            0
          ),
          safeNum(
            b2.get("aspect_cos"),
            0
          )
        );

      var aspect5 =
        circularMeanAspect(
          safeNum(
            b5.get("aspect_sin"),
            0
          ),
          safeNum(
            b5.get("aspect_cos"),
            0
          )
        );

      // ===========================
      // SET OUTPUT ATTRIBUTES
      // ===========================

      return f

        // -----------------------
        // Point elevation
        // -----------------------

        .set(
          "point_elevation_m",

          safeNum(
            pointStats.get(
              "elevation"
            ),
            -999
          )
        )

        // -----------------------
        // Point slope
        // -----------------------

        .set(
          "point_slope_deg",

          safeNum(
            pointStats.get(
              "slope"
            ),
            -999
          )
        )

        // -----------------------
        // Point aspect
        // -----------------------

        .set(
          "point_aspect_deg",

          safeNum(
            pointStats.get(
              "aspect"
            ),
            -999
          )
        )

        // -----------------------
        // Mean elevation
        // -----------------------

        .set(
          "elev_mean_1km_buf",

          safeNum(
            b1.get(
              "elevation"
            ),
            -999
          )
        )

        .set(
          "elev_mean_2km_buf",

          safeNum(
            b2.get(
              "elevation"
            ),
            -999
          )
        )

        .set(
          "elev_mean_5km_buf",

          safeNum(
            b5.get(
              "elevation"
            ),
            -999
          )
        )

        // -----------------------
        // Mean slope
        // -----------------------

        .set(
          "slope_mean_1km_buf",

          safeNum(
            b1.get(
              "slope"
            ),
            -999
          )
        )

        .set(
          "slope_mean_2km_buf",

          safeNum(
            b2.get(
              "slope"
            ),
            -999
          )
        )

        .set(
          "slope_mean_5km_buf",

          safeNum(
            b5.get(
              "slope"
            ),
            -999
          )
        )

        // -----------------------
        // Circular mean aspect
        // -----------------------

        .set(
          "aspect_circmean_1km_buf",
          aspect1
        )

        .set(
          "aspect_circmean_2km_buf",
          aspect2
        )

        .set(
          "aspect_circmean_5km_buf",
          aspect5
        );
    });

  // ===============================
  // EXPORT
  // ===============================

  Export.table.toDrive({

    folder: exportFolder,
collection:
      out,

    description: sanitizeDescription(
        label +
        "_Topography_COPDEM_1_2_5km"
      ),
    fileNamePrefix: sanitizeDescription(
        label +
        "_Topography_COPDEM_1_2_5km"
      ),
    fileFormat:
      "CSV"
  });

  print(
    "Prepared topography export for:",
    label
  );
}

// ===============================
// 8. RUN ALL FILES
// ===============================

for (
  var i = 0;
  i < assetConfigs.length;
  i++
) {

  processAsset(
    assetConfigs[i].assetName,
    assetConfigs[i].label
  );
}
