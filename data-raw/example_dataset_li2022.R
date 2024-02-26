here::i_am("data-raw/new_example_dataset.R")

library(here)
library(tidyverse)
library(readxl)
library(janitor)
library(UpSetR)
library(XML)
library(xml2)
library(tabulizer)

options(timeout = 600)

## Helper functions -----------------------------------------------------------

#' Extract from a list of dataset-specific IDs the IDs that are unique to a
#' subset of these datasets.
unique_to <- function(id_list, x) {
  res <- list(
    id_list[x] |>
      reduce(intersect),
    id_list[-which(names(id_list) %in% x)] |>
      reduce(union)
  ) |>
    reduce(setdiff) |>
    sort()
  message(length(res))

  res
}

#' Finds the unique rows between two data-frames
diff_rows <- function(df1, df2) {
  bind_rows(
    setdiff(df1, df2) |>
      mutate(only_in = "df1"),
    setdiff(df2, df1) |>
      mutate(only_in = "df2")
  )
}

## Downloading files -----------------------------------------------------------

## Genotype, metabolites and samples information available in the borealis
## database (DOI doi:10.5683/SP3/ZETWNY)
## at https://doi.org/10.5683/SP3/ZETWNY

## RNAseq read counts table available in the GEO database (accession GSE217317)
## at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217317

## The genome annotation file is available through the Ensembl database
## (https://asia.ensembl.org/info/data/index.html) on their FTP site
## at https://ftp.ensembl.org/pub/release-110/gff3/bos_taurus/

dir_temp <- tempdir()

files_url <- c(
  genotype = "https://borealisdata.ca/api/access/datafile/411517",
  metabolite = "https://borealisdata.ca/api/access/datafile/411515",
  phenotype = "https://borealisdata.ca/api/access/datafile/411516",
  snp_info = "https://borealisdata.ca/api/access/datafile/411514",
  rnaseq = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217317&format=file&file=GSE217317%5F143%5FBRD%5FCON%5Freadcounts%2Ecsv%2Egz",
  gff = "https://ftp.ensembl.org/pub/release-110/gff3/bos_taurus/Bos_taurus.ARS-UCD1.2.110.gff3.gz",
  rnaseq_samples1 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217317/matrix/GSE217317-GPL23295_series_matrix.txt.gz",
  rnaseq_samples2 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217317/matrix/GSE217317-GPL26012_series_matrix.txt.gz",
  suppl_mat = "https://figshare.com/ndownloader/files/38514059"
)

files_ext <- c(
  genotype = ".csv",
  metabolite = ".tab",
  phenotype = ".tab",
  snp_info = ".txt",
  rnaseq = ".csv.gz",
  gff = ".gff3.gz",
  rnaseq_samples1 = ".txt.gz",
  rnaseq_samples2 = ".txt.gz",
  suppl_mat = ".xlsx"
)

files <- map2( ## Download the files
  files_url,
  files_ext,
  \(x, y) {
    out_file <- tempfile(tmpdir = dir_temp, fileext = y)
    download.file(x, out_file)
    return(out_file)
  }
) |>
  map_chr( ## Unzip the gz files
    \(x) {
      if (str_detect(x, "\\.gz$")) {
        cmd <- paste("gzip -d", x)
        system(cmd)
        x <- str_remove(x, "\\.gz")
      }
      return(x)
    }
  )


## Manually downloaded files

manuscript_file <- "/workspace/hrpoab/transfer_files/Li et al. - 2022 - Applying multi-omics data to study the genetic bac.pdf"

## Metabolites information available in the HMDB database (https://hmdb.ca/)
## Manually downloaded the files from database version 5.0
## (under Downloads > Metabolite and Protein Data (in XML format) > All Metabolites)
hmdb_file <- "/workspace/hrpoab/transfer_files/hmdb_metabolites.xml"

## Reading in datasets ---------------------------------------------------------

geno_raw_df <- read_csv(files[["genotype"]])

geno_fmeta_raw_df <- read_tsv(files[["snp_info"]], name_repair = make_clean_names)

metabo_raw_df <- read_tsv(files[["metabolite"]])

rnaseq_raw_df <- read_csv(files[["rnaseq"]])

samples_raw_df <- read_tsv(files[["phenotype"]])

rnaseq_smeta_raw_df <- files[paste0("rnaseq_samples", 1:2)] |>
  unname() |>
  map_dfr(
    ~ read_tsv(.x, col_names = FALSE) |>
      mutate(X1 = str_remove(X1, "^!")) |>
      filter(X1 %in% c("Sample_title", "Sample_geo_accession", "Sample_characteristics_ch1", "Sample_description")) |>
      separate_wider_delim(
        cols = X2,
        delim = "\t",
        names_sep = "_",
        too_few = "align_start"
      ) |>
      column_to_rownames("X1") |>
      t() |>
      as_tibble() |>
      rename_with(.fn = ~ str_remove(.x, "Sample_")) |>
      rename_with(.fn = ~ str_remove(.x, "_ch1")) |>
      mutate(
        title = str_remove(title, " \\[Re.+"),
        characteristics = str_remove(characteristics, "disease status: ")
      ),
    .id = "batch"
  )

paper_tables <- extract_tables(manuscript_file)

suppl_mat_raw <- c("cis eQTL", "trans eQTL") |>
  set_names() |>
  map(
    \(.x) {
      read_excel(files[["suppl_mat"]], sheet = .x, skip = 1)
    }
  ) |>
  list_rbind(names_to = "qtl_type") |>
  clean_names() |>
  ## Row added from Table 1 in the manuscript
  add_row(
    qtl_type = "QTL",
    snp = "BovineHD1800016801",
    gene_id = NA,
    b = 0.674,
    p_value = 7.65e-6,
    fdr = NA
  ) |>
  group_by(snp) |>
  slice_min(order_by = fdr, n = 1, with_ties = FALSE)


## Wrangling datasets ----------------------------------------------------------

## ---------------------- ##
##      Genotype data     ##
## ---------------------- ##

## Genotype data
dosage_values <- c("AA" = 0, "AB" = 1, "BB" = 2)
geno_df <- geno_raw_df |>
  remove_empty("cols") |>
  select(-Neogen.ID, -ID_RNA) |>
  rename(sample_id = Submitted.ID) |>
  pivot_longer(
    cols = -sample_id,
    names_to = "marker",
    values_to = "dosage"
  ) |>
  mutate(
    dosage = dosage_values[dosage]
  ) |>
  pivot_wider(
    names_from = sample_id,
    values_from = dosage
  )

## Genotypes features metadata
geno_fmeta_df <- geno_fmeta_raw_df |>
  select(-index) |>
  rename(marker = name) |>
  mutate(
    snp = str_remove_all(snp, "\\[|\\]"),
    chromosome = as.character(chromosome)
  ) |>
  separate_wider_delim(
    cols = snp,
    delim = "/",
    names = c("ref", "alt")
  ) |>
  left_join(
    suppl_mat_raw |>
      select(
        marker = snp,
        qtl_type,
        qtl_effect = b,
        p_value,
        fdr
      ),
    by = "marker"
  ) |>
  replace_na(list("qtl_type" = "non signif"))

## Filtering according to the paper; ie removing markers with:
## - > 10% more missing values
## - minor allele frequency < 5%
## - located on sex chromosomes
## Didn't do the filtering according to Hardy-Weinberg equilibrium test
markers_filtering <- geno_df |>
  left_join(
    select(geno_fmeta_df, marker, chromosome),
    by = "marker"
  ) |>
  pivot_longer(
    cols = -c(marker, chromosome),
    names_to = "sample_id",
    values_to = "dosage"
  ) |>
  group_by(marker, chromosome) |>
  summarise(
    n_tot = n(),
    missing = sum(is.na(dosage)),
    mas = sum(dosage, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    frac_missing = missing / n_tot,
    maf = mas / (2 * (n_tot - missing)),
    maf = case_when(
      maf > 0.5 ~ 1 - maf,
      TRUE ~ maf
    ),
    sex_chr = chromosome %in% c("MT", "X", "Y"),
    retained = (frac_missing <= 0.1) & (maf >= 0.05) & !sex_chr
  ) |>
  filter(retained) |>
  pull(marker)

## Random sub-sampling of markers to reduce the dataset size
## Keeping all markers detected as QTL or eQTL

markers_qtls <- geno_fmeta_df |>
  filter(!is.na(qtl_type)) |>
  pull(marker)

set.seed(36)
markers_subsampling <- geno_fmeta_df |>
  filter(qtl_type == "non signif", marker %in% markers_filtering) |>
  slice_sample(n = 23000) |>
  pull(marker)

markers_retained <- c(markers_qtls, markers_subsampling)

geno_df <- geno_df |>
  filter(marker %in% markers_retained)
geno_fmeta_df <- geno_fmeta_df |>
  filter(marker %in% markers_retained)


## ---------------------- ##
##      RNAseq data       ##
## ---------------------- ##

## RNAseq dataset
rnaseq_new_names <- rnaseq_smeta_raw_df |>
  select(title) |>
  mutate(
    new_name = str_remove(title, "_((BRD)|(NB))"),
    new_name = case_when(
      new_name == "P4647" ~ "P4687",
      str_detect(new_name, "^\\d{4}$") ~ paste0("R", new_name),
      TRUE ~ new_name
    )
  ) |>
  select(new_name, title)

rnaseq_df <- rnaseq_raw_df |>
  select(-(Chr:Length)) |>
  rename(gene_id = Geneid, all_of(deframe(rnaseq_new_names)))

## Removing genes with only zero values
genes_to_remove <- rnaseq_df |>
  pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "read_count"
  ) |>
  group_by(gene_id) |>
  summarise(frac_zero = sum(read_count == 0) / n()) |>
  filter(frac_zero >= 0.9) |>
  pull(gene_id)

rnaseq_df <- rnaseq_df |>
  filter(!(gene_id %in% genes_to_remove))

## Genome annotation file
## To make the GFF3 file smaller, we'll remove all information about the
## exons from it.
new_gff_file <- tempfile(tmpdir = dir_temp, fileext = ".gff3")
cmd <- paste(
  "grep -v -E \"ensembl[[:space:]](exon|CDS)\"",
  files[["gff"]],
  ">",
  new_gff_file
)
system(cmd)


## ---------------------------- ##
##      Metabolomics data       ##
## ---------------------------- ##

## Metabolomics dataset
metabo_df <- metabo_raw_df |>
  select(-(Genomicbreedingcomposition1:Dayonfeed)) |>
  rename(sample_id = AnimalID)

## Confirming that there are no metabolites with only 0s or NAs
metabo_df |>
  pivot_longer(
    cols = -sample_id,
    names_to = "feature_id",
    values_to = "value"
  ) |>
  replace_na(list(value = 0)) |>
  group_by(feature_id) |>
  summarise(total = sum(value)) |>
  filter(total == 0)


## ---------------------- ##
##  Samples information   ##
## ---------------------- ##

## Samples metadata
smeta_df <- samples_raw_df |>
  bind_rows(
    metabo_raw_df |>
      select(AnimalID:Dayonfeed) |>
      filter(!(AnimalID %in% samples_raw_df$AnimalID)) |>
      mutate(
        Gender = case_when(
          Gender == 1 ~ 1,
          Gender == 2 ~ 0
        )
      )
  ) |>
  clean_names() |>
  rename(
    rnaseq_batch = batchof_rn_asequencing,
    day_on_feed = dayonfeed
  ) |>
  rename_with(
    .fn = \(x) str_replace(x, "genomicbreedingcomposition", "geno_comp_"),
    .cols = contains("genomicbreedingcomposition")
  ) |>
  relocate(status, .after = "gender") |>
  relocate(rnaseq_batch, .after = "day_on_feed") |>
  ## Adding info about RNAseq samples. We know the ones missing from here
  ## are all males from feedlot 2, and from RNAseq batch 1
  bind_rows(
    rnaseq_smeta_raw_df |>
      full_join(rnaseq_new_names, by = "title") |>
      filter(!(new_name %in% union(samples_raw_df$AnimalID, metabo_raw_df$AnimalID))) |>
      mutate(
        batch = as.numeric(batch),
        characteristics = case_when(
          characteristics == "BRD" ~ "BRD",
          TRUE ~ "Control"
        ),
        feedlot = 2,
        gender = 0
      ) |>
      select(
        animal_id = new_name,
        rnaseq_batch = batch,
        status = characteristics,
        feedlot,
        gender
      )
  ) |>
  ## Recoding factors
  mutate(
    gender = case_when(
      gender == 0 ~ "male",
      gender == 1 ~ "female"
    ),
    feedlot = paste0("F", feedlot),
    rnaseq_batch = paste0("B", rnaseq_batch)
  )


## Checking match between sample IDs
list(
  geno = colnames(geno_df)[-1],
  metabo = metabo_df$sample_id,
  rnaseq = colnames(rnaseq_df)[-1],
  smeta = smeta_df$animal_id
) |>
  reduce(intersect) |>
  length()

## Extracting metabolites information ------------------------------------------

metabo_properties <- c("accession", "name", "chemical_formula",
                       "monisotopic_molecular_weight", "cas_registry_number",
                       "smiles","inchikey", "kegg_id")

extract_metabo_info <- function(x) {
  c(
    x[c("accession", "name", "chemical_formula",
        "monisotopic_molecular_weight", "cas_registry_number",
        "smiles","inchikey", "kegg_id")],
    x[["taxonomy"]][c("direct_parent", "super_class")]
  ) |>
    list_flatten() |>
    discard(is.null)
}

metabo_data <- read_xml(hmdb_file)
metabo_children <- xml_children(metabo_data)

## Function to extract all HMDB IDs for each metabolite
get_metab_ids <- function(x) {

  x <- str_split(x, "\\<\\/secondary_accessions\\>")[[1]][[1]]

  main_id <- str_extract(x, "HMDB\\d+")
  all_ids <- str_extract_all(x, "HMDB\\d+")[[1]]

  tibble(
    main_id = main_id,
    id = all_ids
  )
}

## Lookup table linking each metabolite to its HMDB IDs
metabo_accs <- xml_children(metabo_data) |>
  map(
    ~.x |>
      as.character() |>
      get_metab_ids()
  ) |>
  list_rbind(names_to = "index")

## Checking that all IDs from the dataset are in the database
all(colnames(metabo_df)[-1] %in% metabo_accs$id)

## Extract index of metabolites that are in the dataset
metab_indx <- metabo_accs |>
  filter(id %in% colnames(metabo_df)[-1])

## Checking that there's a unique match
nrow(metab_indx) == (ncol(metabo_df) - 1)

metabo_fmeta_df <- metabo_children[metab_indx$index] |>
  set_names(metab_indx$id) |>
  map(
    ~ .x |>
      as_list() |>
      extract_metabo_info() |>
      as_tibble()
  ) |>
  list_rbind(names_to = "feature_id") |>
  rename(hmdb_id = accession) |>
  mutate(
    monisotopic_molecular_weight = as.numeric(monisotopic_molecular_weight)
  )


## Extracting gene GO annotation -----------------------------------------------
library(biomaRt)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "btaurus_gene_ensembl",
  version = 110
)

go_annot_df <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
  filters = "ensembl_gene_id",
  values = rnaseq_df$gene_id,
  mart = ensembl
) |>
  as_tibble() |>
  filter(go_id != "", namespace_1003 != "") |>
  rename(
    gene_id = ensembl_gene_id,
    go_name = name_1006,
    go_domain = namespace_1003
  ) |>
  mutate(
    go_domain = str_replace(go_domain, "_", " "),
    go_domain = str_to_sentence(go_domain)
  )

## Running the RNAseq DE analysis ----------------------------------------------

library(edgeR)

select <- dplyr::select

rnaseq_mat <- rnaseq_df |>
  column_to_rownames(var = "gene_id") |>
  as.matrix()

## Average genomics composition to be used to replace
## missing values
geno_comp_means <- smeta_df |>
  select(matches("geno_comp_\\d")) |>
  pivot_longer(
    cols = everything(),
    names_to = "var",
    values_to = "val"
  ) |>
  group_by(var) |>
  summarise(
    mean = mean(val, na.rm = TRUE)
  ) |>
  deframe() |>
  as.list()

rnaseq_samples_df <- smeta_df |>
  mutate(
    status = factor(status, levels = c("Control", "BRD"))
  ) |>
  replace_na(geno_comp_means) |>
  column_to_rownames("animal_id") |>
  as.data.frame()

rnaseq_samples_df <- rnaseq_samples_df[colnames(rnaseq_mat), ]

y <- DGEList(
  counts = rnaseq_mat,
  group = rnaseq_samples_df$status,
  samples = rnaseq_samples_df
)

## Step 1: filter out genes with CPM < 0.5 in at least 63 samples (count per million)
##         done in edgeR with filterByExpr
genes_to_keep <- filterByExpr(y, min.count = 0.5)
y <- y[genes_to_keep, , keep.lib.sizes = FALSE]

## Step 2: TMM normalisation
y <- calcNormFactors(y) ## should be the equivalent of normLibSizes function

## Step 3: DE analysis with a negative binomial GLM with the following fixed effects:
##         feedlot, genomic breed composition, sequencing batch
design <- with(
  rnaseq_samples_df,
  model.matrix(~ status + feedlot + geno_comp_1 + geno_comp_2 + geno_comp_3 + rnaseq_batch)
)
y2 <- estimateDisp(y, design)
fit <- glmQLFit(y2, design)
comp <- glmQLFTest(fit, coef = 2)
de_res <- topTags(comp, n = nrow(rnaseq_mat))$table |>
  as_tibble(rownames = "gene_id") |>
  clean_names() |>
  mutate(
    de_signif = case_when(
      fdr < 0.01 & abs(log_fc) > 2 & log_cpm > 2 ~ "DE",
      TRUE ~ "Not DE"
    ),
    de_status = case_when(
      de_signif == "Not DE" ~ "Not DE",
      de_signif == "DE" & log_fc < 0 ~ "downregulated",
      de_signif == "DE" & log_fc > 0 ~ "upregulated"
    )
  )

## Adding genes missing from analysis
de_res <- de_res |>
  bind_rows(
    tibble(
      gene_id = setdiff(rnaseq_df$gene_id, de_res$gene_id),
      de_signif = "Not DE",
      de_status = "Not DE"
    )
  )


## Running the metabolomics DE analysis ----------------------------------------------

metabo_smeta <- smeta_df |>
  filter(animal_id %in% metabo_df$sample_id)

metabo_models <- colnames(metabo_df)[-1] |>
  set_names() |>
  map(
    ~ metabo_df |>
      select(sample_id, value = all_of(.x)) |>
      left_join(
        metabo_smeta,
        by = c("sample_id" = "animal_id")
      ) |>
      mutate(
        value = case_when(
          abs(value) == 0 ~ NA,
          TRUE ~ value
        ),
        value = log(value)) |>
      filter(!is.na(value)) |>
      column_to_rownames("sample_id") |>
      as.data.frame()
  ) |>
  map(
    ~ lm(
      value ~ feedlot + gender + geno_comp_1 + geno_comp_2 + geno_comp_3,
      data = .x
    )
  )

metabo_corr_df <- metabo_models |>
  imap(
    ~ residuals(.x) |>
      enframe(name = "sample_id", value = .y)
  ) |>
  reduce(full_join, by = "sample_id")

metabo_de_res <- colnames(metabo_corr_df)[-1] |>
  set_names() |>
  map(
    ~ metabo_corr_df |>
      select(sample_id, value = all_of(.x)) |>
      filter(!is.na(value)) |>
      left_join(
        select(metabo_smeta, sample_id = animal_id, status),
        by = "sample_id"
      ) |>
      mutate(
        status = factor(status, levels = c("Control", "BRD"))
      ) |>
      group_by(status) |>
      nest() |>
      mutate(
        data = map(data, deframe)
      ) |>
      deframe()
  ) |>
  map(
    ~ t.test(.x[["BRD"]], .x[["Control"]],
             alternative = "two.sided",
             var.equal = TRUE)
  ) |>
  map(
    ~ tibble(
      t_value = .x[["statistic"]],
      p_value = .x[["p.value"]]
    )
  ) |>
  list_rbind(names_to = "feature_id") |>
  mutate(
    padj = p.adjust(p_value, method = "fdr"),
    de_signif = case_when(
      padj < 0.05 ~ "DE",
      TRUE ~ "Not DE"
    ),
    de_status = case_when(
      de_signif == "Not DE" ~ "Not DE",
      de_signif == "DE" & t_value < 0 ~ "downregulated",
      de_signif == "DE" & t_value > 0 ~ "upregulated"
    )
  )

metabo_fmeta_df <- metabo_fmeta_df |>
  left_join(metabo_de_res, by = "feature_id")


## Extracting DE results from manuscript ---------------------------------------

de_genes <- paper_tables[4:6] |>
  map(
    ~ .x |>
      as_tibble() |>
      row_to_names(row_number = 1) |>
      clean_names()
  ) |>
  list_rbind() |>
  mutate(
    across(
      .cols = c(p_value, fdr),
      .fns = ~ str_replace(., " × 10", "e")
    ),
    across(
      log2_fold_change:fdr,
      .fns = ~ str_replace(., "−", "-")
    ),
    across(
      log2_fold_change:fdr,
      .fns = as.numeric
    )
  )

de_metabo <- paper_tables[[8]] |>
  row_to_names(row_number = 1) |>
  as_tibble() |>
  select(-matches("v\\d")) |>
  clean_names() |>
  filter(metabolite_t_value != "") |>
  mutate(
    wrong_order = case_when(
      p_value == "" ~ TRUE,
      TRUE ~ FALSE
    ),
    extract_val = case_when(
      wrong_order ~ str_extract(metabolite_t_value, "^[\\d\\.\\s×−]+"),
      TRUE ~ NA_character_
    ),
    p_value = case_when(
      wrong_order ~ extract_val,
      TRUE ~ p_value
    ),
    metabolite_t_value = case_when(
      wrong_order ~ str_remove(metabolite_t_value, "^[\\d\\.\\s×−]+"),
      TRUE ~ metabolite_t_value
    )
  ) |>
  select(-wrong_order, -extract_val) |>
  separate_wider_regex(
    cols = metabolite_t_value,
    patterns = c(
      "metabolite" = "[\\w-\\s]+",
      "\\s",
      "t_value" = "−?[\\d\\.]+"
    )
  ) |>
  mutate(
    across(
      .cols = c(p_value, fdr),
      .fns = ~ str_replace(., " × 10", "e")
    ),
    across(
      t_value:fdr,
      .fns = ~ str_replace(., "−", "-")
    ),
    across(
      t_value:fdr,
      .fns = as.numeric
    ),
    metabolite = str_replace(metabolite, " -", "-"),
    metabolite = case_when(
      metabolite == "2-Hydroxyisovalerate" ~ "2-Hydroxy-3-methylbutyric acid",
      metabolite == "Isopropanol" ~ "Isopropyl alcohol",
      metabolite == "L-Aspartate" ~ "L-Aspartic acid",
      metabolite == "Malonate" ~ "Malonic acid",
      metabolite == "Methionine" ~ "L-Methionine",
      metabolite == "Tyrosine" ~ "L-Tyrosine",
      metabolite == "Valine" ~ "L-Valine",
      TRUE ~ metabolite
    )
  )

## There must have been a formatting issue for the table
## of metabolomics DE results because the t-values don't
## match with the reported p-values:
## (note that the df migth be slightly different for some
## metabolites due to missing values)
de_metabo |>
  mutate(exp_tval = qt(p_value / 2, 137)) |>
  View()


## Clustering samples based on genomics composition results --------------------

geno_comp_df <- smeta_df |>
  filter(!is.na(geno_comp_1)) |>
  select(animal_id, matches("geno_comp_\\d")) |>
  column_to_rownames("animal_id") |>
  as.data.frame()

set.seed(685)
kmeans_res <- kmeans(geno_comp_df, centers = 3, nstart = 10)
smeta_df <- smeta_df |>
  mutate(
    geno_comp_cluster = kmeans_res$cluster[animal_id],
    geno_comp_cluster = str_c("K", geno_comp_cluster)
  )

## Writing datasets ------------------------------------------------------------

data_path <- here("inst/extdata/")

write_extdata <- function(df, file) {
  write_csv(df, file = paste0(data_path, "/", file))
}

write_extdata(geno_df, "genomics_dataset.csv")
write_extdata(geno_fmeta_df, "genomics_features_info.csv")

write_extdata(rnaseq_df, "transcriptomics_dataset.csv")
file.copy(new_gff_file, paste0(data_path, "/bos_taurus_gene_model.gff3"))
write_extdata(de_res, "transcriptomics_de_results.csv")
write_extdata(go_annot_df, "transcriptomics_go_annotation.csv")

write_extdata(metabo_df, "metabolomics_dataset.csv")
write_extdata(metabo_fmeta_df, "metabolomics_features_info.csv")

write_extdata(smeta_df, "samples_info.csv")
