# Chromatin accessibility QTL mapping with SAIGE-QTL: usage

The same scripts used for single-cell eQTLs map caQTLs from single-cell or single-nucleus ATAC-seq and 10x Multiome data. Each peak is one phenotype, and each peak's cis window is the region tested. [caQTL_methods.md](caQTL_methods.md) gives the model, the reasons for the settings below and the calibration results. Everything here was checked at commit `d89c238a` with the `--tauFloor` option (PR #TBD) added.

## Inputs

**Phenotype file.** A tab-delimited file with one row per cell and these columns:
- donor ID;
- `log_offset`;
- donor-level covariates;
- cell-level covariates;
- one column of counts per peak.

Rules for the columns:
- **Counts.** Use paired-insertion counts, the default `counting_strategy` of SnapATAC2's `make_peak_matrix`.
- **`log_offset`.** The log of the cell's total in-peak count. It is used as given, so it must be on the log scale.
- **Peak column names.** They must be valid R names, for example `chr1_100600_101100`. Names such as `chr1:100600-101100` or `chr1-100600-101100` make step 1 fail.
- **Covariates.** All must be numeric; encode categorical covariates as 0/1 columns. Suggested donor covariates are sex, age, batch indicators, genotype PCs and PCs of donor-level accessibility. Suggested cell covariates are TSS enrichment and the leading spectral (LSI) components.

**Genotypes.**
- **Step 2:** PLINK `.bed/.bim/.fam` files, one set per chromosome. Step 2 reads the whole `.bim` at start-up.
- **Step 1:** a separate PLINK file of about 1,000 random autosomal variants with minor allele count ≥ 20, for the variance ratio. Its `.bim` must use numeric chromosome codes (`1`, not `chr1`).

**Peak windows.** One file per peak with three tab-separated columns and no header: chromosome, start, end. Coordinates are inclusive. The chromosome must be written exactly as in the step 2 `.bim`, or step 2 stops with an error. We use the peak ± 10 kb.

**Peak filter.** Peaks with few donors having nonzero counts carry little information. Dropping them first saves compute; calibration does not require it.

## Step 1: fit the null model (one run per peak)

```bash
Rscript extdata/step1_fitNULLGLMM_qtl.R \
  --phenoFile=pheno.tsv --phenoCol=chr1_100600_101100 \
  --sampleIDColinphenoFile=donor \
  --traitType=count \
  --offsetCol=log_offset \
  --covarColList=sex,age,gPC1,gPC2,aPC1,aPC2,tss,lsi2,lsi3 \
  --sampleCovarColList=sex,age,gPC1,gPC2,aPC1,aPC2 \
  --useGRMtoFitNULL=FALSE --useSparseGRMtoFitNULL=FALSE \
  --isCovariateOffset=FALSE --isCovariateTransform=TRUE --isRemoveZerosinPheno=FALSE \
  --tol=0.001 --tauFloor=0.000001 \
  --plinkFile=vr_markers \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix=out/chr1_100600_101100
```

Notes on these flags:
- `--traitType=count` is required. The default is `binary`, which fails on counts.
- List donor covariates first in `--covarColList` and repeat them in `--sampleCovarColList`. Step 1 fails without `--sampleCovarColList`. With `--isCovariateTransform=TRUE`, a cell covariate listed before the donor covariates corrupts the donor-level design that step 2 uses.
- Do not use `--qCovarColList` for donor covariates. Covariates declared there are silently dropped from the step 2 donor-level design; use 0/1 columns instead.
- `--tol` is the convergence threshold, and `--tauFloor` is the level below which the donor variance is set to 0.
  - Without `--tauFloor`, the floor equals `--tol`. The default 0.02, and even 1e-3, then sets small donor variances to zero, which inflates the tests for well-covered peaks.
  - With a version that lacks `--tauFloor`, use `--tol=0.00001` instead. It was validated for donors with up to about 800 expected counts per peak, but more sparse peaks then fall back or fail (methods).
- `--useGRMtoFitNULL` defaults to `TRUE`, so set it to `FALSE` explicitly. For related donors, add a sparse GRM (`--useSparseGRMtoFitNULL=TRUE` with `--sparseGRMFile` and `--sparseGRMSampleIDFile`).

**Check the outcome.** The script exits with status 0 even when the fit failed. The outcome is the line `Convergence Status: SUCCESS` or `FAILED` in `<outputPrefix>.status.txt`. If the full model fails or does not converge, the script refits with the covariates moved into the offset. Its log then contains `Retrying with all covariates as offset` or `Trying with all covariates as offset`. Record this flag and report calibration and results separately for such peaks.

**Outputs.** The model `<outputPrefix>.rda` stores several vectors of length equal to the number of cells. Delete it after step 2 when running many peaks.

## Step 2: test the variants in the window

```bash
Rscript extdata/step2_tests_qtl.R \
  --bedFile=geno/chr1.bed --bimFile=geno/chr1.bim --famFile=geno/chr1.fam \
  --chrom=1 --rangestoIncludeFile=windows/chr1_100600_101100.txt \
  --GMMATmodelFile=out/chr1_100600_101100.rda \
  --varianceRatioFile=out/chr1_100600_101100.varianceRatio.txt \
  --LOCO=FALSE --minMAC=20 --SPAcutoff=2 \
  --SAIGEOutputFile=out/chr1_100600_101100.txt
```

- Set `--LOCO=FALSE` to match step 1; the step 2 script's default is `TRUE`.
- `--minMAC=20` matches the variants used for the variance ratio.
- An optional `--minMAF` can be added.

## Step 3: peak-level p-value

```bash
Rscript extdata/step3_gene_pvalue_qtl.R \
  --assocFile=out/chr1_100600_101100.txt \
  --geneName=chr1_100600_101100 \
  --genePval_outputFile=out/chr1_100600_101100.acat.txt
```

This gives the ACAT p-value over the window, which is the same as `SAIGEQTL::get_CCT_pvalue(p.value)` on the step 2 output. Adjust across peaks with the Benjamini–Hochberg procedure or Storey's q-value.

## Running many peaks

- Starting R and loading the package dominate the run time of small fits, so run many peaks per job.
- When jobs run in parallel on one machine, set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.
- Run step 2 once per peak, or as a multi-trait run (`--GMMATmodel_varianceRatio_multiTraits_File`) over nearby peaks followed by filtering each peak's window. The multi-trait run gave results identical to per-peak runs in our tests without a GRM. With a sparse GRM, use one peak per run.
- Do not call the step 2 function twice in one R session. Its output streams stay open until the process exits, so later results are written to the first file.

## Pitfalls reproduced

| Setting | What happens |
|---|---|
| `--traitType` omitted | Step 1 fails: `y values must be 0 <= y <= 1` |
| No `--sampleCovarColList` | Step 1 fails |
| Character covariate column | Step 1 fails: `object '<name>' not found` |
| Categorical covariate in `--qCovarColList` | Step 1 runs, but the covariate is missing from the step 2 donor-level design |
| Cell covariate listed before donor covariates (`--isCovariateTransform=TRUE`) | Step 1 runs, but the step 2 donor-level design is wrong |
| Peak name with `:` or `-` | Step 1 fails: `object 'chr1' not found` |
| All-zero peak | Step 1 reports `SUCCESS`, with the donor variance at its floor (1e-6 with `--tauFloor=1e-6`) |
| `--tol=0.02` (default) without `--tauFloor` | Donor variances below 0.02 are set to 0 |
| `chr1`-style codes in the step 1 `--plinkFile` | Step 1 fails: `chromosome column in plink bim file is no numeric` |
| Step 2 `--chrom` differs from the `.bim` code | Step 2 stops with an error and writes only the header |
| Sparse peak whose counts all fall in some levels of a 0/1 donor covariate, such as batch | The full model diverges, is flagged as not converged and is refitted with the covariates as offset. Versions before this change report `SUCCESS` for the diverged fit (in a 12-count simulated peak the fitted values summed to 1.6e17), and step 2 then gives p = 0 for every variant; with those versions, check that the fitted values sum to the observed counts. |

In every step 1 failure above, the script still exits with status 0.
