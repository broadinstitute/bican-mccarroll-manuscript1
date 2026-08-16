# Gabitto 2024 Supplementary Table 5 conversion specification

## Goal

Convert all differential expression results in Gabitto 2024 Supplementary
Table 5 into a canonical "voom-like" file format that can be loaded
interchangeably with existing limma/voom/dream differential expression
outputs.

The conversion should use an adapter pattern. Gabitto-specific sheet names and
column names are interpreted only during conversion. Downstream code should see
the same file naming convention and DE column structure used for other
converted or native limma/voom/dream results.

The converted files are intended for downstream workflows including TRADE,
joint DE/mash-style analyses, and direct comparison of effect sizes.

## Input workbook

The input workbook is:

```text
/downloads/Gabitto_2024/original/gabitto_2024_supplementary_table_5.xlsx
```

The top-level conversion function must accept the workbook path as an argument;
the path above is the expected input for this analysis.

Supplementary Table 5 contains 555 DE-result sheets representing four analysis
types.

| Source sheet pattern | Canonical contrast | n |
|---|---|---:|
| `MTG <subclass>_<supertype>_across_Continuous_Pseudo-progression_Score_DE.csv CPS` | `ad_cps` | 139 |
| `early <subclass>_<supertype>_across_Continuous_Pseudo-progression_Score_DE.csv CPS` | `early_ad_cps` | 139 |
| `late <subclass>_<supertype>_across_Continuous_Pseudo-progression_Score_DE.csv CPS` | `late_ad_cps` | 139 |
| `MTG <supertype>_versus_all_DE.csv vs_all` | `versus_all` | 138 |

All sheets must be converted.

## Excel access

Use an R Excel-reading package rather than converting the workbook to another
format first. `readxl` is a suitable dependency.

Use explicit namespace qualification rather than attaching the package. For
example:

```r
readxl::excel_sheets(in_file)
readxl::read_excel(in_file, sheet = sheet_name)
```

A failure to list workbook sheets or read a sheet is a hard error and must stop
the conversion. An Excel read failure is treated as an I/O or workbook-level
problem that may affect other sheets.

## Public function interface

The conversion should expose two primary functions.

### Top-level conversion

```r
convert_gabitto_2024 <- function(in_file, out_dir) {
    ...
}
```

This is the normal entry point. It must:

1. Create `out_dir` recursively if it does not already exist.
2. Enumerate every sheet in the workbook.
3. Strictly parse every sheet name into canonical metadata.
4. Determine the expected output path for every sheet.
5. Detect output filename collisions before conversion begins.
6. Stop immediately if any sheet name cannot be parsed or if two sheets map to
   the same output filename.
7. Convert each sheet by calling `convert_gabitto_2024_sheet()`.
8. Overwrite existing output DE files.
9. Collect one manifest row for every source sheet.
10. Write the manifest to `<out_dir>/conversion_manifest.txt`.
11. Invisibly return the same manifest as a data frame.

Preflighting all sheet names and expected output paths before writing DE files
is preferred so parsing problems and output collisions are detected before the
conversion creates a partial result set.

### Single-sheet conversion

```r
convert_gabitto_2024_sheet <- function(in_file, sheet_name, out_dir) {
    ...
}
```

This function must convert exactly one Excel sheet into at most one canonical
DE output file.

It must:

1. Parse `sheet_name`.
2. Read that sheet from `in_file`.
3. Select the source columns appropriate for its analysis type.
4. Clean and validate the adapter fields.
5. Build the canonical voom-like table.
6. Construct the canonical output filename.
7. Write the file if usable rows remain.
8. Return a single manifest row describing the result.

The function should remain independently callable for debugging a single
sheet.

Additional helper functions should be added as needed. Keep source-specific
parsing, validation, duplicate handling, and output writing separated where
that improves clarity.

## Sheet-name parsing

Sheet-name parsing is strict.

Every sheet must match one of the known Gabitto Supplementary Table 5 patterns
and must produce:

- a SEA-AD supertype used as `celltype`
- the constant region `MTG`
- one of the four canonical contrasts

An unrecognized or ambiguously parsed sheet name is a hard error. Do not skip
such sheets. A parsing failure indicates that the specification or parser needs
to be refined before conversion continues.

Use the complete supertype name. Do not truncate suffixes or split names on
single underscores.

For filesystem normalization:

- replace spaces with `_`
- remove `/`
- otherwise preserve the parsed supertype text

Examples:

```text
L2/3 IT_1       -> L23_IT_1
L5/6 NP_2       -> L56_NP_2
Sst Chodl_1     -> Sst_Chodl_1
OPC_2_1-SEAAD   -> OPC_2_1-SEAAD
```

Single underscores are valid characters within cell type names. The canonical
filename convention uses double underscores to separate metadata fields, so
names such as `OPC_2_1-SEAAD` remain unambiguous.

## Source table structure

The first source column contains the gene name and has no header label. Treat
this column as the gene identifier.

The source `gene_id` column is a separate identifier and must not replace the
first unnamed gene column.

A representative CPS sheet contains columns including:

```text
<unnamed gene column>
logFC_(Intercept)
logFC_Age_at_Death_binned_codes
logFC_SexM
logFC_Genes_detected
logFC_Race_choice_WhiteUnchecked
logFC_method10xMulti
logFC_Continuous_Pseudo-progression_Score
se_(Intercept)
se_Age_at_Death_binned_codes
se_SexM
se_Genes_detected
se_Race_choice_WhiteUnchecked
se_method10xMulti
se_Continuous_Pseudo-progression_Score
p_(Intercept)
p_Age_at_Death_binned_codes
p_SexM
p_Genes_detected
p_Race_choice_WhiteUnchecked
p_method10xMulti
p_Continuous_Pseudo-progression_Score
gene_id
```

Do not extract coefficients for the intercept, age, sex, genes detected, race,
chemistry/method, or other model covariates.

## Source-to-adapter mapping

For the three CPS analysis types, adapt the following columns:

| Adapter field | CPS source |
|---|---|
| `gene` | first unnamed column |
| `beta` | `logFC_Continuous_Pseudo-progression_Score` |
| `se` | `se_Continuous_Pseudo-progression_Score` |
| `p_value` | `p_Continuous_Pseudo-progression_Score` |

For `versus_all`, adapt:

| Adapter field | `versus_all` source |
|---|---|
| `gene` | first unnamed column |
| `beta` | `logFC_comparison1` |
| `se` | `se_comparison1` |
| `p_value` | `p_comparison1` |

The Gabitto field called `logFC` is used as an effect-size/regression
coefficient. Preserve it directly as canonical `logFC`; do not apply an
additional fold-change transformation.

## Missing required columns

If a recognized sheet is missing any source column required for its adapter:

1. Emit a `logger::log_warn()` message identifying the sheet and missing
   columns.
2. Do not write a DE output file for that sheet.
3. Record the sheet as `skipped` in the manifest.
4. Continue processing the remaining sheets.

Missing required DE columns are considered a source-data problem rather than a
workbook I/O problem.

## Gene cleanup

Gene order should follow the source sheet after filtering. Do not sort genes.

### Blank gene names

Rows with missing, empty, or whitespace-only gene names must be removed.

Emit a `logger::log_warn()` message if any such rows are removed and record the
count in the manifest.

### Duplicate genes

Duplicate handling must use only the adapter fields:

```text
gene
beta
se
p_value
```

Do not use unused source columns to decide whether duplicate gene rows are
equivalent.

For each gene appearing more than once:

1. Compare `beta`, `se`, and `p_value` using floating-point tolerant numerical
   comparison rather than direct `==`.
2. An `all.equal()`-style comparison is appropriate. Numeric comparison should
   ignore irrelevant attributes and use an R floating-point tolerance suitable
   for values read from the same workbook.
3. If all duplicate rows are numerically equivalent in the adapter fields,
   retain only the first occurrence.
4. Preserve the position of that first occurrence in the original sheet
   ordering.
5. Count the removed equivalent duplicate rows in
   `n_exact_duplicates_removed`.
6. If duplicate rows for the same gene are not numerically equivalent, remove
   all rows for that gene.
7. Emit a `logger::log_warn()` message identifying that ambiguous duplicate
   genes were removed.
8. Record the number of rows removed for ambiguous duplicate genes in
   `n_ambiguous_gene_rows_removed`.

Do not arbitrarily choose among conflicting DE results for the same gene.

## Numeric validation

The adapter fields `beta`, `se`, and `p_value` must be numeric.

If a value is not numeric or cannot be cleanly converted to numeric:

1. Emit a `logger::log_warn()` message.
2. Remove the affected gene row.
3. Record the number of rows removed in
   `n_invalid_numeric_rows_removed`.

After numeric conversion, also require finite values where needed.

### Standard error validation

Remove rows where `se` is:

- `NA`
- `NaN`
- infinite
- zero
- negative

Warn with `logger::log_warn()` and record the number of rows removed in
`n_invalid_se_rows_removed`.

### P-value validation

Remove rows where `p_value` is:

- `NA`
- `NaN`
- infinite
- less than 0
- greater than 1

Warn with `logger::log_warn()` and record the number of rows removed in
`n_invalid_p_rows_removed`.

The beta must also be finite. Non-finite beta values should be handled as
invalid numeric rows.

## Empty sheets after filtering

If no usable gene rows remain after validation and filtering:

1. Emit a `logger::log_warn()` message.
2. Do not write an empty DE output file.
3. Record the sheet as `skipped` in the manifest.
4. Continue processing the remaining sheets.

## Canonical voom-like output schema

Each output DE file must contain exactly these seven columns, in this order:

```text
logFC
AveExpr
t
P.Value
adj.P.Val
B
z.std
```

The gene name must be written as the row name and must not be included as a
named `gene` column.

Construct the columns as:

| Canonical column | Value |
|---|---|
| `logFC` | adapter `beta` |
| `AveExpr` | `NA` |
| `t` | `beta / se` |
| `P.Value` | adapter `p_value` |
| `adj.P.Val` | BH-adjusted `P.Value`, independently within the sheet |
| `B` | `NA` |
| `z.std` | `NA` |

Use:

```r
stats::p.adjust(p_value, method = "BH")
```

for `adj.P.Val`.

Although the output field is named `t`, the Gabitto value is a Wald z-like
statistic derived from the reported effect and standard error. The canonical
name is retained for compatibility with existing DE-loading code.

Do not add a separate standard-error column. Downstream code reconstructs the
standard error using:

```r
dt[, log_fc_se := log_fc / t]
```

Because:

```text
t = beta / se
```

this reconstruction should return the original source standard error, subject
to normal floating-point precision.

## Output filenames

All output DE files must use the existing three-field naming convention:

```text
<celltype>__MTG__<contrast>_DE_results.txt
```

Examples:

```text
L23_IT_1__MTG__ad_cps_DE_results.txt
Sst_23__MTG__early_ad_cps_DE_results.txt
Vip_9__MTG__late_ad_cps_DE_results.txt
VLMC_2__MTG__versus_all_DE_results.txt
OPC_2_1-SEAAD__MTG__ad_cps_DE_results.txt
```

The generated filename must be parsable by the existing
`parse_de_result_filenames()` function and must yield:

```text
celltype    = <normalized complete supertype>
interaction = MTG
contrast    = <canonical contrast>
```

The output filename should be validated with the existing parser before the
file is written.

## Output collision handling

If two distinct workbook sheets resolve to the same output filename, stop the
entire conversion with an error.

Do not allow the later sheet to overwrite the earlier sheet in this case.

A collision indicates that the sheet-name parser has discarded information or
that the naming specification is incomplete and must be reviewed.

The top-level function should detect collisions during preflight, before
writing DE output files.

## File writing

Write one tab-delimited text file per successfully converted sheet.

Use behavior equivalent to:

```r
utils::write.table(
    df,
    out_file,
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
)
```

Requirements:

- overwrite an existing file with the same expected path
- retain all valid genes regardless of statistical significance
- preserve source gene order after filtering and duplicate cleanup
- write gene names as row names
- write exactly the seven canonical columns
- do not quote fields

## Manifest

The top-level converter must produce a conversion manifest with one row for
every source workbook sheet, including sheets that are skipped.

Write the manifest as a tab-delimited file:

```text
<out_dir>/conversion_manifest.txt
```

The same data frame should be invisibly returned by
`convert_gabitto_2024()`.

The manifest must include at least:

```text
sheet_name
cell_type
region
contrast
output_file
n_source_rows
n_output_rows
n_blank_gene_rows_removed
n_exact_duplicates_removed
n_ambiguous_gene_rows_removed
n_invalid_numeric_rows_removed
n_invalid_se_rows_removed
n_invalid_p_rows_removed
status
message
```

Use:

```text
region = MTG
```

for all successfully parsed sheets.

`status` should at minimum distinguish:

```text
written
skipped
```

`message` should contain a concise human-readable explanation for skipped
sheets and may be empty for successful sheets.

For a skipped sheet, `output_file` should still contain the output path that
would have been used if conversion had succeeded. This makes the manifest
easier to audit.

The removal-count fields should be zero when no rows were removed for that
reason.

## Logging behavior

Use `logger` for recoverable data-quality problems.

Warnings should identify the affected sheet and summarize what was removed or
why the sheet was skipped.

Recoverable warnings include:

- missing required DE columns
- blank gene names
- adapter-equivalent duplicate rows
- conflicting duplicate gene rows
- non-numeric adapter values
- invalid standard errors
- invalid p-values
- no usable rows remaining after filtering

Hard errors must stop the conversion and should not be converted into warnings.

Hard errors include:

- inability to enumerate workbook sheets
- inability to read an Excel sheet
- an unrecognized or ambiguous sheet-name pattern
- an output filename collision

## Downstream compatibility

The generated DE files must work with the existing DE infrastructure without a
Gabitto-specific branch after conversion.

In particular:

1. `parse_de_result_filenames()` must recover `celltype`, `MTG`, and the
   canonical contrast from each basename.
2. `parse_de_inputs()` must recover `gene` from the output row names.
3. The loaded table must contain `logFC`, `t`, `P.Value`, and `adj.P.Val`.
4. Joint DE/mash-style code must be able to reconstruct standard errors from
   `abs(logFC / t)`.
5. The TRADE formatting path must be able to map:
   - `interaction` to `region`
   - `contrast` to `test`
   - `logFC` to `log_fc`
   - `P.Value` to `p_value`
6. The TRADE path must be able to reconstruct:
   `log_fc_se = log_fc / t`.
7. No Gabitto-specific source column names should be needed downstream.

This is the adapter boundary: Gabitto-specific interpretation ends when the
canonical files are written.

## Validation and acceptance checks

For every successfully written file, verify at minimum:

- the output filename parses successfully
- parsed `celltype` matches the normalized complete source supertype
- parsed `interaction` is `MTG`
- parsed `contrast` matches the source analysis type
- gene row names are non-empty and unique
- gene order matches retained source order
- output columns are exactly:
  `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `z.std`
- `logFC` equals the selected source coefficient
- `P.Value` equals the selected source p-value
- `adj.P.Val` equals BH adjustment within that sheet
- for retained rows, `logFC / t` reproduces the source standard error within
  floating-point tolerance
- no significance filtering has been applied
- the output row count equals the number of valid retained genes
- the manifest counts agree with the filtering performed

The conversion is complete when all 555 sheets have either:

- produced a canonical DE file with `status = written`, or
- been intentionally skipped for a recoverable data-quality problem with
  `status = skipped`

and no hard parsing, collision, or I/O error has occurred.
