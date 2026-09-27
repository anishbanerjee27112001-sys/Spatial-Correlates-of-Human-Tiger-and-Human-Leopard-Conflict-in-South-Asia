# Spatial correlates of human–tiger and human–leopard conflict across South Asia

Code for background sampling, predictor extraction and full binomial GLM analyses. HLC denotes human–leopard conflict; HTC denotes human–tiger conflict. The reported statistical analyses used R 4.4.1. Google Earth Engine (GEE) scripts use JavaScript.

## Files

| Folder/file | Purpose |
|---|---|
| `background_sampling/` | Accessible Area, Bias Corrected and Uniform Random sampling for both species |
| `predictor_extraction/GEE/` | Topography, water, NDVI, tree-cover loss, land cover, road density and distances |
| `predictor_extraction/R/` | Bioclimate, population, livestock and forest fragmentation |
| `statistical_analysis/` | `HLC_GLM.R` and `HTC_GLM.R`; full models only |
| `config/paths.R` | Local data and raster paths; shared input helpers |
| `Setup_Packages.R` | Installs missing R packages |
| `Prepare_GEE_Inputs.R` | Converts the conflict workbooks and six extraction datasets to GEE input CSVs |

## Set paths once

Run R from this repository's root folder. Edit `config/paths.R` for your computer. The supplied data root is:

`C:/Users/Anish Banerjee/Desktop/HCC/Spatial Determinants of HCC/Final Datasets`

Conflict workbooks are `<species> Data/<species> Final Conflict Data.xlsx`. Each species requires these three CSVs in its `Background` folder:

| Species | Input filenames |
|---|---|
| HLC | `HLC_AccessibleArea_Merged_Final.csv`, `HLC_BiasCorrected_Merged_Final.csv`, `HLC_UniformRandom_Merged_Final.csv` |
| HTC | `HTC_AccessibleArea_Merged_Final.csv`, `HTC_BiasCorrected_Merged_Final.csv`, `HTC_UniformRandom_Merged_Final.csv` |

Extraction inputs require `Long`, `Lat`, `Year` and `Conflict`. GLM inputs additionally require the predictor columns named in the scripts. Check the raster locations in `config/paths.R`, especially `worldclim_file`; its default is a configurable location, not a verified file on your computer.

## Run the statistical analyses

If your six final merged datasets already contain the predictors:

```r
source("Setup_Packages.R")
source("statistical_analysis/HLC_GLM.R")
source("statistical_analysis/HTC_GLM.R")
```

Results are saved to `HLC Data/HLC_GLM_Final_Outputs` and `HTC Data/HTC_GLM_Final_Outputs`. Matching existing output files are overwritten. Only full models are produced; old reduced-model files are not deleted. Output numbering is retained for compatibility with figure/table scripts. Each run saves `16_SessionInfo.txt`.

## Run sampling and extraction

1. Run `Setup_Packages.R`. Run `Prepare_GEE_Inputs.R` to create conflict CSVs in each species' `GEE_Input` folder. Its second stage requires all six merged datasets; if they are absent, it stops after writing the conflict CSVs.
2. Upload the conflict CSVs as Earth Engine table assets named `HLC_Final_Conflict_Data` and `HTC_Final_Conflict_Data`. Edit `assetRoot` in each GEE script to your asset folder. The GRIP4 road asset referenced in the scripts must also be accessible.
3. Run each background script in the Earth Engine Code Editor and start its export tasks. Exports go to Google Drive.
4. Combine conflict and background records, assign background years using the study's original procedure, and prepare the six named datasets. Raw sampling exports use lowercase `conflict`; rename it to `Conflict` when preparing these inputs. The supplied scripts do not implement background-year assignment or the original HTC background-ratio adjustment.
5. Rerun `Prepare_GEE_Inputs.R` and upload the six merged input CSVs to Earth Engine, using each filename without `.csv` as its asset name. Run the relevant scripts in `predictor_extraction/GEE/`. Every extraction script is configured for all six datasets.
6. Run each R extraction script from the repository root, for example `source("predictor_extraction/R/Bioclimatic_Variables.R")`. Outputs go to each species' `Extracted` folder. Bioclimate extraction includes BIO3, BIO4, BIO5, BIO6, BIO12, BIO16 and BIO17.
7. Merge extraction outputs into the six final datasets using `Record_ID`, then run the GLMs. Preserve existing identifiers and row order when preparing inputs; generated identifiers depend on input row order. Merging component outputs is not automated here.

## Preparation and validation

Broken expressions introduced in the Word document were restored against the earlier intact scripts, including GLM summary references, PR-AUC arithmetic, livestock-year assignment and extraction filenames. Missing shared setup and the settlement-distance batch wrapper were restored. The duplicate placeholder WorldClim example was removed. Analytical calculations were preserved.

All 12 JavaScript files passed syntax checks. The nine R files passed delimiter checks; R execution and authenticated GEE execution were not available during preparation. This is a packaging cleanup, not a numerical revalidation.

Remaining analytical checks from the code review:
- The custom PR-AUC function does not explicitly add the initial zero-recall endpoint. Performance ROC-AUC uses automatic direction, while permutation importance fixes direction. These calculations are preserved and require review before finalizing results.
- VIF filtering can exit on errors; model convergence and five successful folds per configuration should be checked. `iteration = 100` searches fold assignments; it does not run 100 repeated cross-validations.
- The 2020 livestock filename denotes an ALL-LU product. Verify compatibility with the 2010/2015 animal-density rasters before describing all years as animals per square kilometre.

Source data, raster files and results are not included. No licence is imposed by this package; choose the intended code licence before public release.
