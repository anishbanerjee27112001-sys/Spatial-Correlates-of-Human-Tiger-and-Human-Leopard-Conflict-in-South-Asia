// ==========================================================
// GLAD LAND COVER + FOREST DISTURBED/DEGRADED + SHDI
// ==========================================================
//
// Outputs:
//   1. Percentage cover of consolidated GLAD land-cover classes
//      at 1 km, 2 km and 5 km
//
//   2. Shannon's Diversity Index (SHDI)
//      at 1 km, 2 km and 5 km
//
//   3. Percentage disturbed/degraded forest
//      at 1 km, 2 km and 5 km
//
// Notes:
//   - SHDI is calculated across 8 consolidated land-cover classes.
//   - Natural logarithm is used:
//         SHDI = -Σ(p_i * ln(p_i))
// ==========================================================

// ==========================================================
// 0. BASIC SETTINGS
// ==========================================================

var lonField = "Long";
var latField = "Lat";
var yearField = "Year";

// Native / approximate GLAD analysis scale
var scaleUse = 30;

// ==========================================================
// 1. INPUT FILES
// ==========================================================

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

// ==========================================================
// 2. GLAD DATA
// ==========================================================

// GLAD ocean mask
var landmask =
  ee.Image("projects/glad/OceanMask")
    .lte(1);

// ------------------------------
// GLAD land-cover rasters
// ------------------------------

var gladImgs = {

  2000:
    ee.Image(
      "projects/glad/GLCLU2020/v2/LCLUC_2000"
    )
    .updateMask(landmask)
    .rename("GLAD_LC"),

  2005:
    ee.Image(
      "projects/glad/GLCLU2020/v2/LCLUC_2005"
    )
    .updateMask(landmask)
    .rename("GLAD_LC"),

  2010:
    ee.Image(
      "projects/glad/GLCLU2020/v2/LCLUC_2010"
    )
    .updateMask(landmask)
    .rename("GLAD_LC"),

  2015:
    ee.Image(
      "projects/glad/GLCLU2020/v2/LCLUC_2015"
    )
    .updateMask(landmask)
    .rename("GLAD_LC"),

  2020:
    ee.Image(
      "projects/glad/GLCLU2020/v2/LCLUC_2020"
    )
    .updateMask(landmask)
    .rename("GLAD_LC")

};

// ------------------------------
// GLAD forest type layer
// ------------------------------
//
// Class 4 = disturbed/degraded forest
//

var forestType =
  ee.Image(
    "projects/glad/GLCLU2020/Forest_type"
  )
  .updateMask(landmask)
  .rename("Forest_type");

// ==========================================================
// 3. CONSOLIDATED LAND-COVER GROUPS
// ==========================================================

// Cropland
var cropCodes =
  [1, 2, 3];

// Forest
var forestCodes =
  [5, 6, 7, 8, 9, 10, 11, 12, 13];

// Shrubland
var shrubCodes =
  [14, 15, 16, 21];

// Grassland
var grassCodes =
  [17];

// Built-up
var builtCodes =
  [26];

// Bare ground
var bareCodes =
  [27, 28, 29];

// Water
var waterCodes =
  [30];

// ------------------------------
// Codes already assigned to a
// consolidated category
// ------------------------------

var groupedCodes =
  cropCodes
    .concat(forestCodes)
    .concat(shrubCodes)
    .concat(grassCodes)
    .concat(builtCodes)
    .concat(bareCodes)
    .concat(waterCodes)
    .concat([0]);

// ==========================================================
// 4. HELPER FUNCTIONS
// ==========================================================

// ----------------------------------------------------------
// Sanitize export description
// ----------------------------------------------------------

function sanitizeDescription(txt) {

  return String(txt)
    .replace(
      /[^A-Za-z0-9_]/g,
      "_"
    );
}

// ----------------------------------------------------------
// Create binary mask from GLAD codes
//
// focal category = 1
// all remaining categories = 0
// ----------------------------------------------------------

function maskFromCodes(
  img,
  codes,
  name
) {

  return img
    .remap(
      codes,
      ee.List.repeat(
        1,
        codes.length
      ),
      0
    )
    .rename(name);

}

// ----------------------------------------------------------
// "Other" land-cover category
//
// Pixels not belonging to one of the
// explicitly grouped classes are assigned 1.
//
// Class 0 is not included in "other".
// ----------------------------------------------------------

function othersMask(img) {

  var known =
    maskFromCodes(
      img,
      groupedCodes,
      "known"
    );

  return img
    .neq(0)
    .and(
      known.eq(0)
    )
    .rename("others");

}

// ----------------------------------------------------------
// Temporal matching
//
// < 2000     -> -999 / null outputs
// 2000-2004  -> 2000
// 2005-2009  -> 2005
// 2010-2014  -> 2010
// 2015-2019  -> 2015
// >= 2020    -> 2020
// ----------------------------------------------------------

function assignGladYear(y) {

  y =
    ee.Number(y);

  return ee.Number(

    ee.Algorithms.If(
      y.lt(2000),
      -999,

      ee.Algorithms.If(
        y.lte(2004),
        2000,

        ee.Algorithms.If(
          y.lte(2009),
          2005,

          ee.Algorithms.If(
            y.lte(2014),
            2010,

            ee.Algorithms.If(
              y.lte(2019),
              2015,
              2020
            )
          )
        )
      )
    )
  );
}

// ----------------------------------------------------------
// Retrieve correct GLAD raster
// ----------------------------------------------------------

function getGladImage(y) {

  y =
    ee.Number(y);

  return ee.Image(

    ee.Algorithms.If(
      y.eq(2000),
      gladImgs[2000],

      ee.Algorithms.If(
        y.eq(2005),
        gladImgs[2005],

        ee.Algorithms.If(
          y.eq(2010),
          gladImgs[2010],

          ee.Algorithms.If(
            y.eq(2015),
            gladImgs[2015],
            gladImgs[2020]
          )
        )
      )
    )
  );
}

// ==========================================================
// 5. PERCENTAGE LAND-COVER FUNCTION
// ==========================================================
//
// Binary raster:
//     focal category = 1
//     other categories = 0
//
// Mean × 100 therefore gives percentage cover.
// ==========================================================

function pctInBuffer(
  binaryImg,
  geom,
  radiusM,
  bandName
) {

  var v =
    binaryImg
      .unmask(0)
      .reduceRegion({

        reducer:
          ee.Reducer.mean(),

        geometry:
          geom.buffer(radiusM),

        scale:
          scaleUse,

        maxPixels:
          1e8,

        tileScale:
          8

      })
      .get(bandName);

  return ee.Algorithms.If(

    ee.Algorithms.IsEqual(
      v,
      null
    ),

    null,

    ee.Number(v)
      .multiply(100)

  );
}

// ==========================================================
// 6. SHANNON'S DIVERSITY INDEX
// ==========================================================
//
// SHDI = -Σ(p_i × ln(p_i))
//
// p_i = proportional abundance of land-cover category i.
//
// Eight consolidated categories:
//
//   1. Forest
//   2. Cropland
//   3. Shrubland
//   4. Grassland
//   5. Built-up
//   6. Bare ground
//   7. Water
//   8. Other
//
// Categories with p_i = 0 contribute zero.
//
// Maximum possible SHDI with eight equally abundant
// categories:
//
//   ln(8) ≈ 2.079
//
// SHDI = 0 where only one represented land-cover
// category occurs.
//
// Proportions are normalized across the eight represented
// categories before SHDI is calculated.
// ==========================================================

function shdiInBuffer(
  classStack,
  geom,
  radiusM
) {

  // --------------------------------------------------------
  // Mean value of each binary band =
  // proportional cover of that category
  // --------------------------------------------------------

  var proportions =
    classStack.reduceRegion({

      reducer:
        ee.Reducer.mean(),

      geometry:
        geom.buffer(radiusM),

      scale:
        scaleUse,

      maxPixels:
        1e8,

      tileScale:
        8

    });

  // --------------------------------------------------------
  // Names of the eight consolidated classes
  // --------------------------------------------------------

  var bandNames =
    ee.List([

      "forest",
      "cropland",
      "shrubland",
      "grassland",
      "built",
      "bare",
      "water",
      "others"

    ]);

  // --------------------------------------------------------
  // Extract proportional cover values.
  //
  // Null values are converted to zero.
  // --------------------------------------------------------

  var pList =
    bandNames.map(
      function(name) {

        name =
          ee.String(name);

        var value =
          proportions.get(name);

        return ee.Number(

          ee.Algorithms.If(

            ee.Algorithms.IsEqual(
              value,
              null
            ),

            0,

            value
          )
        );
      }
    );

  // --------------------------------------------------------
  // Sum represented proportions
  // --------------------------------------------------------

  var total =
    ee.Number(
      pList.reduce(
        ee.Reducer.sum()
      )
    );

  // --------------------------------------------------------
  // Calculate:
  //
  //     -p_i × ln(p_i)
  //
  // after normalising p_i so represented categories
  // sum to 1.
  //
  // p = 0 contributes zero.
  // --------------------------------------------------------

  var shannonContributions =
    pList.map(
      function(p) {

        p =
          ee.Number(p);

        var pNorm =
          ee.Number(

            ee.Algorithms.If(

              total.gt(0),

              p.divide(total),

              0
            )
          );

        return ee.Number(

          ee.Algorithms.If(

            pNorm.gt(0),

            pNorm
              .multiply(
                pNorm.log()
              )
              .multiply(-1),

            0
          )
        );
      }
    );

  // --------------------------------------------------------
  // Sum contributions
  // --------------------------------------------------------

  var shdi =
    ee.Number(

      shannonContributions.reduce(
        ee.Reducer.sum()
      )

    );

  // --------------------------------------------------------
  // Return null where no valid represented
  // land-cover categories are available.
  // --------------------------------------------------------

  return ee.Algorithms.If(

    total.gt(0),

    shdi,

    null

  );
}

// ==========================================================
// 7. NULL OUTPUTS FOR PRE-2000 ROWS
// ==========================================================

function nullOutputs(f) {

  return f

    // ------------------------------------------------------
    // GLAD reference year
    // ------------------------------------------------------

    .set(
      "GLAD_year_used_Final",
      null
    )

    // ------------------------------------------------------
    // Forest
    // ------------------------------------------------------

    .set(
      "GLAD_pct_forest_1km_Final",
      null
    )

    .set(
      "GLAD_pct_forest_2km_Final",
      null
    )

    .set(
      "GLAD_pct_forest_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Cropland
    // ------------------------------------------------------

    .set(
      "GLAD_pct_cropland_1km_Final",
      null
    )

    .set(
      "GLAD_pct_cropland_2km_Final",
      null
    )

    .set(
      "GLAD_pct_cropland_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Shrubland
    // ------------------------------------------------------

    .set(
      "GLAD_pct_shrubland_1km_Final",
      null
    )

    .set(
      "GLAD_pct_shrubland_2km_Final",
      null
    )

    .set(
      "GLAD_pct_shrubland_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Grassland
    // ------------------------------------------------------

    .set(
      "GLAD_pct_grassland_1km_Final",
      null
    )

    .set(
      "GLAD_pct_grassland_2km_Final",
      null
    )

    .set(
      "GLAD_pct_grassland_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Built-up
    // ------------------------------------------------------

    .set(
      "GLAD_pct_built_1km_Final",
      null
    )

    .set(
      "GLAD_pct_built_2km_Final",
      null
    )

    .set(
      "GLAD_pct_built_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Bare ground
    // ------------------------------------------------------

    .set(
      "GLAD_pct_bare_1km_Final",
      null
    )

    .set(
      "GLAD_pct_bare_2km_Final",
      null
    )

    .set(
      "GLAD_pct_bare_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Water
    // ------------------------------------------------------

    .set(
      "GLAD_pct_water_1km_Final",
      null
    )

    .set(
      "GLAD_pct_water_2km_Final",
      null
    )

    .set(
      "GLAD_pct_water_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Other
    // ------------------------------------------------------

    .set(
      "GLAD_pct_others_1km_Final",
      null
    )

    .set(
      "GLAD_pct_others_2km_Final",
      null
    )

    .set(
      "GLAD_pct_others_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Shannon's Diversity Index
    // ------------------------------------------------------

    .set(
      "GLAD_SHDI_1km_Final",
      null
    )

    .set(
      "GLAD_SHDI_2km_Final",
      null
    )

    .set(
      "GLAD_SHDI_5km_Final",
      null
    )

    // ------------------------------------------------------
    // Disturbed/degraded forest
    // ------------------------------------------------------

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_1km_Final",
      null
    )

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_2km_Final",
      null
    )

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_5km_Final",
      null
    );
}

// ==========================================================
// 8. PROCESS EACH POINT
// ==========================================================

function processFeature(f) {

  // Point geometry
  var geom =
    f.geometry();

  // Temporally matched GLAD year
  var y =
    ee.Number(
      f.get(
        "GLAD_year_used_internal"
      )
    );

  // Retrieve corresponding GLAD image
  var lc =
    getGladImage(y);

  // ========================================================
  // Binary consolidated land-cover masks
  // ========================================================

  var forest =
    maskFromCodes(
      lc,
      forestCodes,
      "forest"
    );

  var crop =
    maskFromCodes(
      lc,
      cropCodes,
      "cropland"
    );

  var shrub =
    maskFromCodes(
      lc,
      shrubCodes,
      "shrubland"
    );

  var grass =
    maskFromCodes(
      lc,
      grassCodes,
      "grassland"
    );

  var built =
    maskFromCodes(
      lc,
      builtCodes,
      "built"
    );

  var bare =
    maskFromCodes(
      lc,
      bareCodes,
      "bare"
    );

  var water =
    maskFromCodes(
      lc,
      waterCodes,
      "water"
    );

  var others =
    othersMask(lc);

  // ========================================================
  // Stack consolidated categories for SHDI
  // ========================================================

  var landCoverStack =
    ee.Image.cat([

      forest,
      crop,
      shrub,
      grass,
      built,
      bare,
      water,
      others

    ]);

  // ========================================================
  // Disturbed / degraded forest
  // ========================================================

  var forestDisturbedDegraded =
    forestType
      .eq(4)
      .rename(
        "forest_disturbed_degraded"
      );

  // ========================================================
  // Attach outputs
  // ========================================================

  return f

    // ------------------------------------------------------
    // GLAD year
    // ------------------------------------------------------

    .set(
      "GLAD_year_used_Final",
      y
    )

    // ------------------------------------------------------
    // Forest percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_forest_1km_Final",
      pctInBuffer(
        forest,
        geom,
        1000,
        "forest"
      )
    )

    .set(
      "GLAD_pct_forest_2km_Final",
      pctInBuffer(
        forest,
        geom,
        2000,
        "forest"
      )
    )

    .set(
      "GLAD_pct_forest_5km_Final",
      pctInBuffer(
        forest,
        geom,
        5000,
        "forest"
      )
    )

    // ------------------------------------------------------
    // Cropland percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_cropland_1km_Final",
      pctInBuffer(
        crop,
        geom,
        1000,
        "cropland"
      )
    )

    .set(
      "GLAD_pct_cropland_2km_Final",
      pctInBuffer(
        crop,
        geom,
        2000,
        "cropland"
      )
    )

    .set(
      "GLAD_pct_cropland_5km_Final",
      pctInBuffer(
        crop,
        geom,
        5000,
        "cropland"
      )
    )

    // ------------------------------------------------------
    // Shrubland percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_shrubland_1km_Final",
      pctInBuffer(
        shrub,
        geom,
        1000,
        "shrubland"
      )
    )

    .set(
      "GLAD_pct_shrubland_2km_Final",
      pctInBuffer(
        shrub,
        geom,
        2000,
        "shrubland"
      )
    )

    .set(
      "GLAD_pct_shrubland_5km_Final",
      pctInBuffer(
        shrub,
        geom,
        5000,
        "shrubland"
      )
    )

    // ------------------------------------------------------
    // Grassland percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_grassland_1km_Final",
      pctInBuffer(
        grass,
        geom,
        1000,
        "grassland"
      )
    )

    .set(
      "GLAD_pct_grassland_2km_Final",
      pctInBuffer(
        grass,
        geom,
        2000,
        "grassland"
      )
    )

    .set(
      "GLAD_pct_grassland_5km_Final",
      pctInBuffer(
        grass,
        geom,
        5000,
        "grassland"
      )
    )

    // ------------------------------------------------------
    // Built-up percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_built_1km_Final",
      pctInBuffer(
        built,
        geom,
        1000,
        "built"
      )
    )

    .set(
      "GLAD_pct_built_2km_Final",
      pctInBuffer(
        built,
        geom,
        2000,
        "built"
      )
    )

    .set(
      "GLAD_pct_built_5km_Final",
      pctInBuffer(
        built,
        geom,
        5000,
        "built"
      )
    )

    // ------------------------------------------------------
    // Bare ground percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_bare_1km_Final",
      pctInBuffer(
        bare,
        geom,
        1000,
        "bare"
      )
    )

    .set(
      "GLAD_pct_bare_2km_Final",
      pctInBuffer(
        bare,
        geom,
        2000,
        "bare"
      )
    )

    .set(
      "GLAD_pct_bare_5km_Final",
      pctInBuffer(
        bare,
        geom,
        5000,
        "bare"
      )
    )

    // ------------------------------------------------------
    // Water percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_water_1km_Final",
      pctInBuffer(
        water,
        geom,
        1000,
        "water"
      )
    )

    .set(
      "GLAD_pct_water_2km_Final",
      pctInBuffer(
        water,
        geom,
        2000,
        "water"
      )
    )

    .set(
      "GLAD_pct_water_5km_Final",
      pctInBuffer(
        water,
        geom,
        5000,
        "water"
      )
    )

    // ------------------------------------------------------
    // Other land cover percentage
    // ------------------------------------------------------

    .set(
      "GLAD_pct_others_1km_Final",
      pctInBuffer(
        others,
        geom,
        1000,
        "others"
      )
    )

    .set(
      "GLAD_pct_others_2km_Final",
      pctInBuffer(
        others,
        geom,
        2000,
        "others"
      )
    )

    .set(
      "GLAD_pct_others_5km_Final",
      pctInBuffer(
        others,
        geom,
        5000,
        "others"
      )
    )

    // ======================================================
    // SHANNON'S DIVERSITY INDEX
    // ======================================================

    .set(
      "GLAD_SHDI_1km_Final",
      shdiInBuffer(
        landCoverStack,
        geom,
        1000
      )
    )

    .set(
      "GLAD_SHDI_2km_Final",
      shdiInBuffer(
        landCoverStack,
        geom,
        2000
      )
    )

    .set(
      "GLAD_SHDI_5km_Final",
      shdiInBuffer(
        landCoverStack,
        geom,
        5000
      )
    )

    // ======================================================
    // DISTURBED / DEGRADED FOREST
    // ======================================================

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_1km_Final",
      pctInBuffer(
        forestDisturbedDegraded,
        geom,
        1000,
        "forest_disturbed_degraded"
      )
    )

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_2km_Final",
      pctInBuffer(
        forestDisturbedDegraded,
        geom,
        2000,
        "forest_disturbed_degraded"
      )
    )

    .set(
      "GLAD_2000_2020_pct_forest_disturbed_degraded_5km_Final",
      pctInBuffer(
        forestDisturbedDegraded,
        geom,
        5000,
        "forest_disturbed_degraded"
      )
    );
}

// ==========================================================
// 9. PROCESS EACH ASSET
// ==========================================================

function processAsset(
  assetName,
  label
) {

  // --------------------------------------------------------
  // Load source FeatureCollection
  // --------------------------------------------------------

  var raw =
    ee.FeatureCollection(
      assetName
    );

  // --------------------------------------------------------
  // Keep rows with usable coordinates and year
  // --------------------------------------------------------

  var pts =
    raw

      .filter(
        ee.Filter.notNull([
          lonField,
          latField,
          yearField
        ])
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

      .filter(
        ee.Filter.neq(
          yearField,
          ""
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
          yearField,
          "NA"
        )
      )

      .map(
        function(f) {

          // Parse longitude
          var lon =
            ee.Number.parse(
              ee.String(
                f.get(
                  lonField
                )
              )
            );

          // Parse latitude
          var lat =
            ee.Number.parse(
              ee.String(
                f.get(
                  latField
                )
              )
            );

          // Parse first four characters of year
          var yr =
            ee.Number.parse(

              ee.String(
                f.get(
                  yearField
                )
              )
              .slice(
                0,
                4
              )

            );

          // Match observation year to GLAD raster year
          var gladYear =
            assignGladYear(yr);

          // Rebuild feature using parsed coordinates
          return ee.Feature(

            ee.Geometry.Point([
              lon,
              lat
            ]),

            f.toDictionary()

          )

          .set(
            "parsed_year",
            yr
          )

          .set(
            "GLAD_year_used_internal",
            gladYear
          );

        }
      );

  // ========================================================
  // QA: GLAD class histogram
  // ========================================================

  print(

    label +
      " GLAD class histogram QA",

    getGladImage(2020)
      .sampleRegions({

        collection:
          pts.limit(5000),

        scale:
          30,

        geometries:
          false

      })
      .aggregate_histogram(
        "GLAD_LC"
      )

  );

  // ========================================================
  // Process rows
  // ========================================================

  var out =
    pts.map(
      function(f) {

        return ee.Feature(

          ee.Algorithms.If(

            ee.Number(
              f.get(
                "GLAD_year_used_internal"
              )
            )
            .eq(-999),

            // Pre-2000
            nullOutputs(f),

            // 2000 onward
            processFeature(f)

          )
        );
      }
    );

  // ========================================================
  // Export
  // ========================================================

  Export.table.toDrive({

    folder: exportFolder,
collection:
      out,

    description: sanitizeDescription(
        label +
        "_GLAD_LANDCOVER_SHDI_FOREST_DEGRADED"
      ),
    fileNamePrefix: sanitizeDescription(
        label +
        "_GLAD_LANDCOVER_SHDI_FOREST_DEGRADED"
      ),
    fileFormat:
      "CSV"

  });

  // ========================================================
  // QA messages
  // ========================================================

  print(
    "Prepared GLAD landcover + SHDI + forest disturbed/degraded:",
    label
  );

  print(
    "Input rows:",
    raw.size()
  );

  print(
    "Output rows with valid coordinates:",
    out.size()
  );

}

// ==========================================================
// 10. RUN ALL ASSETS
// ==========================================================

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
