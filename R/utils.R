# Assumes l = list(a, b, ...) where a, b, ... themselves are lists, and are
# *rows/observations*. Handles cases where, if a = list(x, y, ...), x is a
# character vector. Useful for combining gene filter objects into data frames
# with rows as gene filters.
.df_from_list_of_lists = function(l) {
  init_df = data.frame(t(sapply(l,c)))
  # Unlist df types that are nested
  for (column in names(init_df)) {
    if (all(init_df[,column] %>% purrr::map(length) == 1)) {
      # Get type of first elt in column; CURRENTLY LOSES S3 CLASS!
      base_type = typeof(init_df[,column][[1]])
      init_df[,column] = as.vector(init_df[,column], mode = base_type)
    }
  }
  init_df
}

# Eliminate columns with at most n unique values from a df, and return the a
# list of the simplified df and the eliminated column names and unique values.
# A further improvement here would be to store the IDs corresponding to the
# unique values (group by, essentially) in list_of_elims.
.elim_unif_cols = function(df, cols, n = 1) {
  elims = list()
  for (column in names(df)) {
    if (column %in% cols) {
      unique_vals = unique(df[,column])
      if(length(unique_vals) <= n) {
        df[,column] = NULL # remove column
        elims[[column]] = unique_vals
      }
    }
  }
  list(df = df, elims = elims)
}

.param_sentence = function(param_list) {
  names(df_and_elims$elims) %>% purrr::map_chr(function(name) {
    glue::glue("{.sentence_case(name, cap = F, abbreviations = c('FDR'), repl_list = list(thresh = 'threshold',
               pval = 'p-value'))} (", glue::glue(df_and_elims$elims[[name]], .sep = ", "), ")")
  }) %>% glue::collapse(sep = ", ", last = " and ")
}

.sentence_case = function(str_vec, abbreviations = NULL, repl_list = NULL, cap = TRUE) {
  lc_words = snakecase::to_any_case(str_vec, case = 'snake', sep_out = ' ', abbreviations = abbreviations)
  # Map shorthand to longhand
  for (repl in names(repl_list)) {
    lc_words = gsub(repl, repl_list[[repl]], lc_words)
  }
  # Recapitalize abbreviations (doesn't deal with abbreviations within words)
  for (abbrev in abbreviations) {
    lc_words = gsub(tolower(abbrev), toupper(abbrev), lc_words)
  }
  ifelse(rep(cap, length(lc_words)), paste0(toupper(substr(lc_words, 1, 1)), substr(lc_words, 2, nchar(lc_words))), lc_words)
}

.geom_text = function(...) {
  geom_text(..., family = 'Palatino')
}

.cache_ = function(name, directory, fn) {
  full_filename = file.path(directory, paste0(name, '.rds'))
  if (file.exists(full_filename)) {
    readRDS(full_filename)
  } else {
    res = match.fun(fn)()
    if (!dir.exists(directory)) {
      dir.create(directory)
    }
    saveRDS(res, full_filename)
    res
  }
}

.filterGSE = function(gse, classesOfInterest) {
  includeIndices = which(gse$pheno$group %in% classesOfInterest)
  gse$pheno = gse$pheno[includeIndices,]
  gse$expr = gse$expr[,includeIndices]
  rownames(gse$pheno) = colnames(gse$expr)
  if ("class" %in% colnames(gse) && !is.null(gse$class)) {
    gse$class = gse$class[includeIndices]
    names(gse$class) = rownames(gse$pheno)
  }
  gse
}

.addClassVec = function(gse, caseClasses) {
  gse$class = ifelse(gse$pheno$group %in% caseClasses, 1, 0)
  names(gse$class) = rownames(gse$pheno)
  gse
}

.removeLOCMIRFromFilter = function(geneFilter) {
  geneFilter$posGeneNames = discard(geneFilter$posGeneNames, ~grepl("(^(LOC|MIR)[[:digit:]]+|orf)", .))
  geneFilter$negGeneNames = discard(geneFilter$negGeneNames, ~grepl("(^(LOC|MIR)[[:digit:]]+|orf)", .))
  geneFilter
}

.vec_to_csl = function(vec) {
  ret = collapse(vec, sep = ", ", last = " and ")
  if (length(ret) == 0) {
    ""
  } else {
    ret
  }
}

.shortDescribeGeneSig = function(geneFilter) {
  glue("{.num_genes_in_sig(geneFilter)}-gene signature")
}

.describeGeneSig = function(geneFilter, case_classes = NULL, control_classes = NULL) {
  desc = glue("{.num_genes_in_sig(geneFilter)}-gene signature")
  if (length(case_classes) > 0) {
    desc = glue("{desc} for {.vec_to_csl(case_classes)}")
  }
  if (length(control_classes) > 0) {
    desc = glue("{desc} vs. {.vec_to_csl(control_classes)}")
  }
}

.describeGeneSigVerbose = function(geneFilter, case_classes = NULL, control_classes = NULL,
                                   datasets = NULL,
                                   include_params = c("effectSizeThresh", "FDRThresh", "numberStudiesThresh", "isLeaveOneOut")) {
  desc = .describeGeneSig(geneFilter, case_classes, control_classes)
  if (length(datasets) > 0) {
    desc = glue("{desc} discovered in {.vec_to_csl(datasets)}")
  }
  if (length(include_params) > 0) {
    params_list = geneFilter[include_params]
    desc = glue("{desc} (params: {map2_chr(names(params_list), params_list, ~ glue('{.x}={.y}')) %>% .vec_to_csl})")
  }
  desc = glue("{desc} (pos.: {.vec_to_csl(geneFilter$posGeneNames)};",
              " neg.: {.vec_to_csl(geneFilter$negGeneNames)})")
  desc
}

.num_genes_in_sig = function(geneFilter) {
  return(length(geneFilter$posGeneNames) + length(geneFilter$negGeneNames))
}

.first_num_in_string = function(str_vec) {
  as.double(stringr::str_match(pattern = "([0-9]+)", string = str_vec)[,1])
}

.avg_expr_per_gene = function(gse) {
  gse %>%
    { dplyr::inner_join(as.data.frame(.$expr, stringsAsFactors = FALSE) %>% tibble::rownames_to_column("probe"),
                        data.frame(symbol=.$keys, stringsAsFactors = FALSE) %>% tibble::rownames_to_column("probe") %>% dplyr::filter(!is.na(symbol))) %>%
        dplyr::select(-probe) } %>%
    splitstackshape::cSplit("symbol", direction = "long") %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarize_all(dplyr::funs(mean)) %>%
    tibble::column_to_rownames(var = "symbol") %>%
    t %>%
    { dplyr::inner_join(as.data.frame(., stringsAsFactors = FALSE) %>% tibble::rownames_to_column("sample"),
                        data.frame(gse$pheno[, "group", drop=FALSE], stringsAsFactors = FALSE) %>% tibble::rownames_to_column("sample"))}
}

.excl_named_list = function(a_list, names_vec) {
  names_to_keep = names(a_list) %>% purrr::discard(function(name) {
    name %in% names_vec
  })
  a_list[names_to_keep]
}

.gen_metrics_df = function(metaObject, filtersList) {
  plyr::ldply(filtersList, function(gene_filter) { .gen_metrics_df_per_filter(metaObject, gene_filter) }, .id = "gene_filter")
}

.gen_metrics_df_per_filter = function(metaObject, filterObject) {
  groupName = "group"
  predName = "score"
  roc_obj_list = MetaIntegrator:::.rocObjList(metaObject, filterObject)
  roc_plot_data = plyr::llply(roc_obj_list, function(x) with(x, MetaIntegrator:::.rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  roc_stats = plyr::ldply(roc_plot_data, function(ds) { ds$stats }, .id = "dataset") %>%
    dplyr::rename(metric_val = auc, metric_upper = ci.upper, metric_lower = ci.lower) %>%
    dplyr::select(-p.value) %>%
    dplyr::mutate(metric_name = "auroc")
  prc_plot_data = plyr::llply(metaObject$originalData, function(dataset) {MetaIntegrator:::.prcData(dataset, filterObject)})
  prc_stats = plyr::ldply(prc_plot_data, function(ds) {
    stat_obj = ds$auprcInfo
    data.frame(metric_val = stat_obj$auprc, metric_upper = stat_obj$auprc.CI[2], metric_lower = stat_obj$auprc.CI[1])
  }, .id = "dataset") %>%
    dplyr::mutate(metric_name = "auprc")
  dplyr::bind_rows(roc_stats, prc_stats)
}

.gen_cohort_labels_df = function(originalData) {
  plyr::ldply(originalData, function(ds) { ds$pheno[, "group",drop=F] }, .id = "dataset")
}

.get_gene_filter = function(datetime = "latest") {
  file_pattern = glue::glue("geneFilters_{if (datetime == 'latest') '' else datetime}")
  files = list.files(path = here::here("data"), pattern = file_pattern, recursive = TRUE)
  filename = sort(files, decreasing = TRUE)[1]
  readRDS(here::here("data", filename))
}

.abbreviate_if = function(str_vec, condition = function(x) {TRUE}, ...) {
  cond_fn = match.fun(condition)
  str_vec[cond_fn(str_vec)] = abbreviate(str_vec[cond_fn(str_vec)], ...)
  str_vec
}
