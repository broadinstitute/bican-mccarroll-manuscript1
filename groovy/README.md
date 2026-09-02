# QuPath Groovy scripts

Groovy scripts used for image analysis in [QuPath](https://qupath.github.io/) as part of the OPC
(oligodendrocyte precursor cell) manuscript analysis: cell detection, ROI generation, object
classification, and area/count quantification.

## Requirements

- [QuPath](https://qupath.github.io/) (tested with the current stable release — see
  [qupath.github.io](https://qupath.github.io/) for download/install instructions for your OS).
  No separate Groovy install is needed; QuPath bundles its own Groovy runtime and script editor.
- A QuPath project with images already imported.

## Running a script

1. Open your QuPath project.
2. Go to **Automate > Script editor** (or **Automate > Show script editor**).
3. **File > Open...** and select the desired `.groovy` file from this directory.
4. With an image open, run the script via **Run > Run** (single image) or **Run > Run for project**
   (batch over selected images in the project), depending on the script.
5. Output prints to the script editor's console; some scripts also return
   tab-separated tables meant to be pasted into Excel.

## Scripts

| Script | Purpose |
| --- | --- |
| `Full Image ROI Generator.groovy` | Creates a full-image annotation to use as the region of interest for downstream detection. |
| `Cell Detection Script (Modifiable).groovy` | Runs watershed cell detection (DAPI channel) over the selected annotation(s) and prints cell count / area / density for the current image. Detection parameters are hardcoded at the top of the script — adjust as needed per staining/imaging conditions. |
| `Run Classifier.groovy` | Runs a previously trained object classifier on the current image's detections. |
| `Parameter List.groovy` | Prints all measurement names available on the first detection object — useful for finding the exact class/measurement name to reference in other scripts. |
| `Print Count Proportion.groovy` | Prints per-image area, total/positive/unclassified detection counts, and percent positive for a given class (set `POS_CLASS` at the top of the script). |
| `ROI Area Calculator.groovy` | Prints the number of area annotations and total annotated area (mm²) for the current image. |
| `Image Area and Cell Detection Summation.groovy` | Run once via **Run for project**: walks every image in the project, sums annotated area and detection counts, and prints a formatted table plus manuscript-ready summary stats (median density, totals for Methods text). Supports an `INCLUDE`/`EXCLUDE` image-name filter at the top of the script. |

Most single-image scripts are meant to be run per-image (or via **Run for project** if they don't
aggregate state); `Image Area and Cell Detection Summation.groovy` is meant to be run once across
the whole project to produce manuscript summary numbers.
