
get_summ <- function(results, rm_bounds = c(-Inf, Inf), quant_bound = c(0, 1)){
  
  
  if(sum(rm_bounds != c(-Inf, Inf)) != 0){
    results <- results %>%
      rowwise() %>%
      mutate(sdr = ifelse( sdr < rm_bounds[2] & sdr > rm_bounds[1], sdr, NA),
             tml = ifelse( tml < rm_bounds[2] & tml > rm_bounds[1], tml, NA)) %>%
      # filter(!(sdr > rm_bounds[2] | sdr < rm_bounds[1]),
      #        !(tml > rm_bounds[2] | tml < rm_bounds[1])) %>%
      ungroup()
  } else if(sum(quant_bound != c(0, 1)) != 0) {
    results <- results %>%
      group_by(n, scen) %>%
      mutate(sdr = ifelse( sdr < quantile(sdr, quant_bound[2]) & sdr > quantile(sdr, quant_bound[1]), sdr, NA),
             tml = ifelse( tml < quantile(tml, quant_bound[2]) & tml > quantile(tml, quant_bound[1]), tml, NA)) %>%
      ungroup()
  }
  
  
  ress_update <- results %>%
    mutate(
      ## Bias ##
      # U1 bias #
      sdrU1_bias_k = sdrU1 - trueU1,
      tmlU1_bias_k = tmlU1 - trueU1,
      
      # biassub = sub - true,
      biassdr_k = sdr - true,
      biastml_k = tml - true,
      biasipw_k = IPW - true,
      
      ## abs bias ##
      # biassub = sub - true,
      abs_biassdr_k = abs(sdr - true),
      abs_biastml_k = abs(tml - true),
      abs_biasipw_k = abs(IPW - true),
      
      ## relative bias ##
      # rel_biassub = (sub - true)/true,
      rel_biassdr_k = (sdr - true)/true,
      rel_biastml_k = (tml - true)/true,
      rel_biasipw_k = (IPW - true)/true,
      
      ## abs relative bias ##
      # rel_biassub = (sub - true)/true,
      abs_rel_biassdr_k = abs((sdr - true)/true),
      abs_rel_biastml_k = abs((tml - true)/true),
      abs_rel_biasipw_k = abs((IPW - true)/true),
      
      
      coversdr_k = true >= sdrlo & true <= sdrhi,
      covertml_k = true >= tmllo & true <= tmlhi
      
    ) %>%
    ungroup() %>%
    group_by(scen, n, type) %>%
    summarise(
      
      ## mean bias ##
      U1biassdr = mean(sdrU1_bias_k, na.rm = TRUE),
      U1biastml = mean(tmlU1_bias_k, na.rm = TRUE),
      
      ## mean bias ##
      biassdr = mean(biassdr_k, na.rm = TRUE),
      biastml = mean(biastml_k, na.rm = TRUE),
      biasipw = mean(biasipw_k, na.rm = TRUE),
      
      ## median bias ##
      med_biassdr = median(biassdr_k, na.rm = TRUE),
      med_biastml = median(biastml_k, na.rm = TRUE),
      med_biasipw = median(biasipw_k, na.rm = TRUE),
      
      ## mean abs bias ##
      abs_biassdr = mean(abs_biassdr_k, na.rm = TRUE),
      abs_biastml = mean(abs_biastml_k, na.rm = TRUE),
      abs_biasipw = mean(abs_biasipw_k, na.rm = TRUE),
      
      ## mean rel bias ##
      rel_biassdr = mean(rel_biassdr_k, na.rm = TRUE),
      rel_biastml = mean(rel_biastml_k, na.rm = TRUE),
      rel_biasipw = mean(rel_biasipw_k, na.rm = TRUE),
      
      ## mean abs rel bias ##
      abs_rel_biassdr = mean(abs_rel_biassdr_k, na.rm = TRUE),
      abs_rel_biastml = mean(abs_rel_biastml_k, na.rm = TRUE),
      abs_rel_biasipw = mean(abs_rel_biasipw_k, na.rm = TRUE),
      
      ## Mean squared error ##
      MSEsdr = mean(biassdr_k^2, na.rm = TRUE),
      MSEtml = mean(biastml_k^2, na.rm = TRUE),
      MSEipw = mean(biasipw_k^2, na.rm = TRUE),
      
      ## variance of estimates ##
      varsdr_est = var(sdr, na.rm = TRUE),
      vartml_est = var(tml, na.rm = TRUE),
      varipw_est = var(IPW, na.rm = TRUE),
      
      ## variance of bias ##
      vb_sdr = var(biassdr_k, na.rm = TRUE),
      vb_tml = var(biastml_k, na.rm = TRUE),
      vb_ipw = var(biasipw_k, na.rm = TRUE),
      
      ## Coverage ##
      covsdr = mean(coversdr_k, na.rm = TRUE),
      covtml = mean(covertml_k, na.rm = TRUE)
    ) %>% ungroup()
  
  return(ress_update)
}



get_plots <- function(summ, results, rm_bounds, quant_bound){
  
  ### U1 ###
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # subset #
  if(sum(rm_bounds != c(-Inf, Inf)) != 0){
    results <- results %>%
      rowwise() %>%
      mutate(sdr = ifelse( sdr < rm_bounds[2] & sdr > rm_bounds[1], sdr, NA),
             tml = ifelse( tml < rm_bounds[2] & tml > rm_bounds[1], tml, NA)) %>%
      # filter(!(sdr > rm_bounds[2] | sdr < rm_bounds[1]),
      #        !(tml > rm_bounds[2] | tml < rm_bounds[1])) %>%
      ungroup()
  } 
  if(sum(quant_bound != c(0, 1)) != 0) {
    results <- results %>%
      group_by(n, scen) %>%
      mutate(sdr = ifelse( sdr < quantile(sdr, quant_bound[2]) & sdr > quantile(sdr, quant_bound[1]), sdr, NA),
             tml = ifelse( tml < quantile(tml, quant_bound[2]) & tml > quantile(tml, quant_bound[1]), tml, NA)) %>%
      ungroup()
  } 
  if(sum(quant_bound != c(0, 1)) != 0 || sum(rm_bounds != c(-Inf, Inf)) != 0) {
    rm_num <- results %>% 
      group_by(n, scen) %>% 
      summarise(n_it = n(), 
                rm_sdr = sum(is.na(sdr)),
                rm_tml = sum(is.na(tml))) %>% 
      ungroup()
  } else{
    rm_num <- NULL
  }
  
  # process data for U1 #
  processed_data <- results %>%
    select(scen, n, sdrU1, tmlU1) %>% 
    group_by(scen, n) %>% 
    mutate(mean_sdr = mean(sdrU1, na.rm = TRUE),
           mean_tml = mean(tmlU1, na.rm = TRUE)) %>% 
    ungroup()
  
  # Melt the data
  melted_data <- tidyr::pivot_longer(
    processed_data,
    cols = c(sdrU1, tmlU1),
    names_to = "variable",
    values_to = "value"
  )
  
  
  trueU1 <- unique(results$trueU1)
  pU1 <- ggplot(melted_data, aes(x = value, fill = variable)) +
    geom_density(aes(y = ..density..), alpha = 0.3, adjust = 1.5, size = 0.5) +
    geom_vline(data = subset(melted_data, variable == "sdrU1"), aes(xintercept = mean_sdr), linetype = "dashed", color = "blue") +
    geom_vline(data = subset(melted_data, variable == "tmlU1"), aes(xintercept = mean_tml), linetype = "dashed", color = "red") +
    facet_wrap(~ scen + n, nrow = 5) +
    geom_vline(xintercept = trueU1, linetype = "dashed") +
    theme_bw() +
    scale_fill_manual(values = c("blue", "red"), labels = c("SDR", "TMLE")) + xlab("U1") + ylab("Density") + xlim(c(trueU1-1,trueU1+1))
  
  
  
  
  ### U2 ###
  #### raw values ####
  # process data for U1 #
  processed_data_U2 <- results %>%
    select(scen, n, sdr, tml) %>% 
    group_by(scen, n) %>% 
    mutate(mean_sdr = mean(sdr, na.rm = TRUE),
           mean_tml = mean(tml, na.rm = TRUE)) %>% 
    ungroup()
  
  # Melt the data
  melted_data_U2 <- tidyr::pivot_longer(
    processed_data_U2,
    cols = c(sdr, tml),
    names_to = "variable",
    values_to = "value"
  )
  
  trueU2 <- unique(results$true)
  pU2 <- ggplot(melted_data_U2, aes(x = value, fill = variable)) +
    geom_density(aes(y = ..density..), alpha = 0.3, adjust = 1.5, size = 0.5) +
    geom_vline(data = subset(melted_data_U2, variable == "sdr"), aes(xintercept = mean_sdr), linetype = "dashed", color = "blue") +
    geom_vline(data = subset(melted_data_U2, variable == "tml"), aes(xintercept = mean_tml), linetype = "dashed", color = "red") +
    facet_wrap(~ scen + n, nrow = 5) +
    geom_vline(xintercept = trueU2, linetype = "dashed") +
    theme_bw() +
    scale_fill_manual(values = c("blue", "red"), labels = c("SDR", "TMLE")) + xlab("U1") + ylab("Density") + xlim(c(trueU2-0.2,trueU2+0.2))
  
  
  ## U1 plot ##
  pU1bias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("U1bias", variable), !grepl("abs|rel|med",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  pU1bias_slope_n <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("U1bias", variable), !grepl("abs|rel|med",variable)) %>%
    mutate(value = value * sqrt(n)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  
  #### summary plots ####
  pmedbias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("med_bias", variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  pbias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("bias", variable), !grepl("abs|rel|med|U1",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  pbias_slope_n <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("bias", variable), !grepl("abs|rel|med|U1",variable)) %>%
    mutate(value = value * sqrt(n)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias * sqrt(n)") +
    labs(variable = "Method")
  p_absbias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("abs_bias", variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  p_absbias_slope_n <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("abs_bias", variable)) %>%
    mutate(value = value * sqrt(n)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias * sqrt(n)") +
    labs(variable = "Method")
  p_relbias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("rel_bias", variable), !grepl("abs",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias") +
    labs(variable = "Method")
  p_relbias_slope_n <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("rel_bias", variable), !grepl("abs",variable)) %>%
    mutate(value = value * sqrt(n)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Bias * sqrt(n)") +
    labs(variable = "Method")
  pMSE_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("MSE", variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) +
    geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("MSE") +
    labs(variable = "Method")
  pMSE_slope_n <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("MSE", variable)) %>%
    mutate(value = value * n) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("MSE * n") +
    labs(variable = "Method")
  p_var_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("var", variable), #| grepl("vm", variable),
           !grepl("cover",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Variance of estimated values") +
    labs(variable = "Method")
  p_varbias_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("vb_", variable), #| grepl("vm", variable),
           !grepl("cover",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black", "red"), labels = c("SDR","TMLE","IPW")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1) + ylab("Variance of estimated values") +
    labs(variable = "Method")
  p_cov_slope <- melt(summ, id = c("scen", "n", "type")) %>%
    filter(type == "slope", grepl("cov", variable),!grepl("cover",variable)) %>%
    ggplot(aes(x = n, y = value, group = variable, color = variable)) + geom_point(aes(group = variable)) +
    scale_color_manual(values=c("blue", "black"), labels = c("SDR","TMLE")) +
    geom_line(aes(group = variable)) + theme_bw() + facet_wrap(~scen,nrow = 1)  + ylab("Coverage") + coord_cartesian(ylim = c(0.6, 1)) +
    labs(variable = "Method") + geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 0.5)
  
  
  return(list(
    pU1 = pU1,
    pU2 = pU2,
    pmedbias_slope = pmedbias_slope,
    pbias_slope = pbias_slope,
    pbias_slope_n = pbias_slope_n,
    p_absbias_slope = p_absbias_slope,
    p_absbias_slope_n = p_absbias_slope_n,
    p_relbias_slope = p_relbias_slope,
    p_relbias_slope_n = p_relbias_slope_n,
    pMSE_slope = pMSE_slope,
    pMSE_slope_n = pMSE_slope_n,
    p_var_slope = p_var_slope,
    p_varbias_slope = p_varbias_slope,
    p_cov_slope = p_cov_slope,
    rm_num = rm_num,
    pU1bias_slope = pU1bias_slope,
    pU1bias_slope_n = pU1bias_slope_n
  ))
}
