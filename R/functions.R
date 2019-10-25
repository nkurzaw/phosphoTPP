#' @export
getProb <- function(x, p_thres = 0.75){
    tmp <- strsplit(gsub("[^0-9,\\.]", ",", x), split = ",")
    sapply(tmp, function(y){
        any(as.numeric(y[which(y != "")]) >  p_thres)
    })
}

#' @export
countGoodProb <- function(x, p_thres = 0.75){
    tmp <- strsplit(gsub("[^0-9,\\.]", ",", x), split = ",")
    sapply(tmp, function(y){
        length(which(as.numeric(y[which(y != "")]) >  p_thres))
    })
}

#' @export
getTm <- function(h1_model, x_in = seq(35, 75, by = 0.001), tolerance = 0.01){
    predicted_vals <- predict(h1_model, newdata = data.frame(
        temperature = x_in))
    tm_id <- ifelse(any(abs(predicted_vals-0.5) < tolerance),
                    which(abs(predicted_vals-0.5) ==
                              min(abs(predicted_vals-0.5), na.rm = T))[1],
                    NA)
    tm <-  ifelse(!is.na(tm_id), x_in[tm_id], NA)
    return(tm)
}

fitMeltcurve <- function(df, curve_func = drc::LL.4(fixed = c(NA, NA, 1, NA))){
    h1_model <- try(drc::drm(df$rel_value ~ df$temperature,
                             data = df, fct = curve_func),
                    silent = TRUE)
    return(h1_model)
}

#' @export
fitMeltcurveModelAndEval <- function(df, curve_func = drc::LL.4(fixed = c(NA, NA, 1, NA))){
    h0_model <- lm(df$rel_value ~ 1)
    h1_model <- fitMeltcurve(df = df, curve_func = curve_func)
    if((class(h1_model) %in% c("drc")) & !(class(h1_model) %in% c("try-error"))){
        rssH0 <- sum(residuals(h0_model)^2)
        rssH1 <- sum(residuals(h1_model)^2)
        r_sq <- 1 - (rssH1/rssH0)
        predicted_vals <- predict(h1_model, newdata = data.frame(
            temperature = seq(40, 70, by = 0.001)))
        tm <-  getTm(h1_model = h1_model)
        data.frame(rSq = r_sq,
                   meltPoint = tm,
                   slopeFactor = coef(h1_model)[1],
                   plateau = coef(h1_model)[2],
                   #maximum = coef(h1_model)[3],
                   inflectionPoint = coef(h1_model)[3])
    }else{
        data.frame(rSq = NA,
                   meltPoint = NA,
                   slopeFactor = NA,
                   plateau = NA,
                   maximum = NA,
                   inflectionPoint = NA)
    }
}

#' @export
getPhosphoSiteId <- function(gene_name, string, seq){
    seq_match <- str_locate(seq, gsub("[p,_]", "", string))[1]
    ids_mapped <- str_locate_all(pattern = "p", string = string)[[1]][,1]
    corrected_ids <- seq_match + ids_mapped - seq_len(length(ids_mapped)) -1
    phospho_aa <- paste0("p", substring(string, ids_mapped + 1, ids_mapped + 1))
    return(paste(gene_name, paste0(phospho_aa, corrected_ids, collapse = "" ), sep = "_"))
}

getHqSubset <- function(df, n_rep = 3,
                        filter_criteria = list("rel_fc_127H_low" = 0,
                                               "rel_fc_129H_low" = 0.4,
                                               "rel_fc_129H_high" = 0.6,
                                               "rel_fc_130H_low" = 0,
                                               "rel_fc_130H_high" = 0.3,
                                               "rel_fc_131L_low" = 0,
                                               "rel_fc_131L_high" = 0.2)){
    jointP <- df %>% 
        group_by(protein_id, sequence, modifications) %>% 
        filter(length(unique(dataset_id)) >= n_rep) %>% 
        ungroup() %>% 
        filter(rel_fc_127H > filter_criteria[["rel_fc_127H_low"]],
               rel_fc_129H > filter_criteria[["rel_fc_129H_low"]],
               rel_fc_129H < filter_criteria[["rel_fc_129H_high"]],
               rel_fc_130H > filter_criteria[["rel_fc_130H_low"]],
               rel_fc_130H < filter_criteria[["rel_fc_130H_high"]],
               rel_fc_131L > filter_criteria[["rel_fc_131L_low"]],
               rel_fc_131L < filter_criteria[["rel_fc_131L_high"]],
               !grepl("##", protein_id)) %>% 
        rowwise() %>% 
        mutate(id = paste(protein_id, sequence, modifications, sep = "_")) %>% 
        ungroup         
    
    rep_tab <- table(jointP$dataset_id)
    mostP_rep <- names(rep_tab[rep_tab == max(rep_tab)])
    
    normP <- df %>% 
        rowwise() %>% 
        mutate(id = paste(protein_id, sequence, modifications, sep = "_")) %>% 
        ungroup %>% 
        filter(id %in% filter(jointP, dataset_id == mostP_rep)$id) 
    
    return(normP)
    
}

fitHqSubset <- function(normP, temperature_anno){
    normP_median <- normP %>% 
        group_by(dataset_id) %>% 
        summarize(rel_fc_126 = median(rel_fc_126, na.rm = TRUE),
                  rel_fc_127L = median(rel_fc_127L, na.rm = TRUE),
                  rel_fc_127H = median(rel_fc_127H, na.rm = TRUE),
                  rel_fc_128L = median(rel_fc_128L, na.rm = TRUE),
                  rel_fc_128H = median(rel_fc_128H, na.rm = TRUE),
                  rel_fc_129L = median(rel_fc_129L, na.rm = TRUE),
                  rel_fc_129H = median(rel_fc_129H, na.rm = TRUE),
                  rel_fc_130L = median(rel_fc_130L, na.rm = TRUE),
                  rel_fc_130H = median(rel_fc_130H, na.rm = TRUE),
                  rel_fc_131L = median(rel_fc_131L, na.rm = TRUE)) %>% 
        ungroup
    
    fit_df <- gather(normP_median, fc_channel, rel_value, matches("rel_fc")) %>% 
        left_join(temperature_anno %>% 
                      dplyr::select(fc_channel, temperature),
                  by = "fc_channel")
    
    meltcurve_fits <-  fit_df %>% 
        group_by(dataset_id) %>% 
        do({
            suppressWarnings(fitMeltcurveModelAndEval(df = .))
        }) %>% 
        ungroup 
    
    best_fit_instance <- meltcurve_fits %>% 
        filter(rSq == max(meltcurve_fits$rSq, na.rm = TRUE))
    
    best_fit_df <- fit_df %>% 
        filter(dataset_id == best_fit_instance$dataset_id[1])
    
    best_fit <- fitMeltcurve(df = best_fit_df)
    
    return(list("best_fit" = best_fit,
                "median_df" = fit_df,
                "R_sq" = best_fit_instance$rSq))
}

getNormFactors <- function(meltcurve_fit_list){
    dataset_ids <- unique(meltcurve_fit_list$median_df$dataset_id)
    out_list <- lapply(dataset_ids, function(rep_i){
        out_df <- meltcurve_fit_list[["median_df"]] %>% 
            filter(dataset_id == rep_i) %>% 
            mutate(predicted = predict(meltcurve_fit_list[["best_fit"]])) %>% 
            mutate(norm_factor = predicted / rel_value) %>% 
            dplyr::select(temperature, median_rel_value = rel_value,
                          norm_factor, predicted)
        return(out_df)
    })
    return(out_list)
}

#' @export
getNormFactors4Data <- function(pep_tab, temperature_anno, n_rep,
                                filter_criteria = list("rel_fc_127H_low" = 0,
                                                       "rel_fc_129H_low" = 0.4,
                                                       "rel_fc_129H_high" = 0.6,
                                                       "rel_fc_130H_low" = 0,
                                                       "rel_fc_130H_high" = 0.3,
                                                       "rel_fc_131L_low" = 0,
                                                       "rel_fc_131L_high" = 0.2)){
    hq_df <- getHqSubset(df = pep_tab, n_rep = n_rep, 
                         filter_criteria = filter_criteria)
    meltcurve_fit_list <- fitHqSubset(hq_df, temperature_anno)
    norm_factor_list <- getNormFactors(meltcurve_fit_list)
    return(norm_factor_list)
}

decideOnGeneName <- function(in_vec, ref_vec){
    vapply(in_vec, FUN = function(x){
        alts <- strsplit(x, split = "\\|")[[1]]
        if(length(which(alts %in% ref_vec)) == 1){
            alts[which(alts %in% ref_vec)]
        }else{
            alts[1]
        }
    } , FUN.VALUE = "tst", USE.NAMES = FALSE)
}

#' @export
chooseCoherentGeneNames <- function(in_df, their_data){
    their_gene_names <- unique(gsub("_.+", "", their_data$Gene_pSite))
    out_df <- in_df %>%
        within(gene_name[grep("\\|", gene_name)] <-
                   decideOnGeneName(in_vec = gene_name[grep("\\|", gene_name)],
                                    ref_vec = their_gene_names))
    return(out_df)
}

#' @export
lowerFn <- function(data, mapping, dot.size) {
    p <- ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.15) +
        geom_line(aes(x,y), data = tibble(x = 40:70, y = 40:70),
                  color = "darkgray") +
        ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
        theme_bw()
    p
}

#' @export
countAtLeast <- function(x, val_with_sign, min_n = 2, thres = 0.1){
    ids <- which(x <= thres)
    non_na_ids <- which(!is.na(x))
    (length(ids) >= min_n) & (length(unique(sign(val_with_sign[non_na_ids]))) < 2)
}

#' @export
plotDiffPhosphoMeltcurve <- function(g_name, p_seq, 
                                     predict_vec = seq(35, 70, by = 0.1),
                                     nbf_df, phospho_df, Gene_pSite,
                                     in_ylim = c(0, 1.25)){
    prot_df <- bind_rows(
        filter(phospho_df, gene_name == g_name,
               grepl(p_seq, mod_sequence)) %>% 
            mutate(group = "phosphopeptide"),
        filter(nbf_df, gene_name == g_name)%>% 
            mutate(group = "all unmodified")) 
    
    prot_summarized_df <- prot_df %>% 
        group_by(group, temperature) %>% 
        summarize(mean_rel_value = mean(rel_value, na.rm = TRUE),
                  #sd_rel_value  = sd(rel_value, na.rm = TRUE),
                  se_rel_value  = sd(rel_value, na.rm = TRUE)/
                      sqrt(n())) %>% 
        ungroup()
    
    prot_nbf_fit <- fitMeltcurve(df = filter(prot_df, group == "all unmodified"))
    prot_phospho_fit <- fitMeltcurve(df = filter(prot_df, group == "phosphopeptide"))
    
    prot_fit_df <- tibble(
        temperature = rep(predict_vec, 2),
        mean_rel_value = c(
            predict(prot_nbf_fit, newdata = data.frame(
                temperature = predict_vec)),
            predict(prot_phospho_fit, newdata = data.frame(
                temperature = predict_vec))
        ),
        group = c(
            rep("all unmodified", length(predict_vec)),
            rep("phosphopeptide", length(predict_vec)))
    )
    
    ggplot(prot_summarized_df, aes(temperature, mean_rel_value, color = group)) +
        geom_point(shape = 16, size = 0.5) +
        geom_errorbar(aes(ymin = mean_rel_value - se_rel_value, 
                          ymax = mean_rel_value + se_rel_value), 
                      width=.1, size = 0.25) +
        # geom_smooth(method = "lm", 
        #             formula = y ~ splines::ns(x, df = 4),
        #             se = FALSE) +
        geom_line(data = prot_fit_df, size = 0.25) +
        scale_color_manual("", values = c("all unmodified" = "chartreuse3", "phosphopeptide" = "purple")) +
        theme_bw() +
        labs(x = expression('Temperature ('*~degree*C*')'),
             y = "fraction non-denatured") +
        theme(legend.position = "bottom") +
        ggtitle(gsub("_", " ", Gene_pSite)) +
        coord_cartesian(ylim = in_ylim) +
        theme_paper
}

#' @export
plotDiffPhosphoSpline <- function(g_name, p_seq, spline_df = 4,
                                  predict_vec = seq(35, 70, by = 0.1),
                                  nbf_df, phospho_df, Gene_pSite){
    prot_df <- bind_rows(
        filter(phospho_df, gene_name == g_name,
               grepl(p_seq, mod_sequence)) %>% 
            mutate(group = "phosphopeptide"),
        filter(nbf_df, gene_name == g_name)%>% 
            mutate(group = "all unmodified")) 
    
    prot_summarized_df <- prot_df %>% 
        group_by(group, temperature) %>% 
        summarize(mean_rel_value = mean(rel_value, na.rm = TRUE),
                  #sd_rel_value  = sd(rel_value, na.rm = TRUE),
                  se_rel_value  = sd(rel_value, na.rm = TRUE)/
                      sqrt(n())) %>% 
        ungroup()
    
    ggplot(prot_summarized_df, aes(temperature, mean_rel_value, color = group)) +
        geom_point() +
        geom_errorbar(aes(ymin = mean_rel_value - se_rel_value, 
                          ymax = mean_rel_value + se_rel_value), 
                      width=.1) +
        geom_smooth(method = "lm",
                    formula = y ~ splines::ns(x, df = spline_df),
                    se = FALSE, size = 0.5) +
        #geom_line(data = prot_fit_df) +
        scale_color_manual("", values = c("all unmodified" = "chartreuse3", "phosphopeptide" = "purple")) +
        theme_bw() +
        labs(x = expression('Temperature ('*~degree*C*')'),
             y = "fraction non-denatured") +
        theme(legend.position = "bottom") +
        ggtitle(gsub("_", " ", Gene_pSite))
}