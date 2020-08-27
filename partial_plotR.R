#####Function to make partial dependence plots of all variables wanted
##
##uses the packages patchwork, tidyverse, gtable, magrittr
##arguments required: 
#data - data.frame of full output from xgboost
#hist.data - data.frame of frequencies for plotting histograms (column names currently are variable.name, varaible.value for the x axis and value for the y axis frequency)
#vars - a vector of the variables of interest. It may be simplest to simply do unique(data$variable.name) to grab all variables
#type - whether you want the plots to show the mean and 95% confidence interval (as derived from the t-distribution) or mean and all bootstrap results
#histogram - a logical argument with TRUE meaning that histogram plots are made and overlayed with the partial dependency plots. 

##example: partial_plot(pd_out, out, vars = levels(out$variable.name), type = "mean")
partial_plot <- function(data, hist.data, vars, type = c("mean", "all"), histogram = T, ...) {
  PLT <- lapply(vars, function(vars) {
    BOOT <- data %>% dplyr::filter(variable.name == vars) %$% bootstrap_run %>% unique
    if(length(BOOT) > 1){
      ROWN <- sapply(BOOT, function(j) data %>% dplyr::filter(variable.name == vars & bootstrap_run == j) %$% x %>% length)
      ROWN <- do.call(c, lapply(ROWN, function(i) 1:i))
      UNIQ <- data %>% dplyr::filter(variable.name == vars) %$% x %>% unique %>% as.character %>% unique %>% as.numeric
      TAB <- data %>% dplyr::filter(variable.name == vars) %$% x %>% base::table()
      NUMS <- base::sort(UNIQ[which(UNIQ %in% names(sort(TAB, decreasing = T)[1:max(ROWN)]))])
      if(length(NUMS) != max(ROWN)) {
        NUMA <- data %>% dplyr::filter(variable.name == vars) %$% x %>% table/max(BOOT)
        if(length(NUMA) == 1) NUMS <- rep(names(NUMA), NUMA) else {
        NUMS <- do.call(c, sapply(1:length(NUMA), function(x) rep(names(NUMA)[x], NUMA[x])))
        }
      }
      if(type == "all") {
        data %>% dplyr::filter(variable.name == vars) %>% dplyr::mutate(rown = ROWN) %>% dplyr::group_by(rown) %>% dplyr::summarize(mean = mean(yhat), conf25 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[1]), conf95 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[2])) %>% dplyr::mutate(var = NUMS) %>% ggplot() +
          geom_line(data = data %>% filter(variable.name == vars), aes(x = x, y = yhat, group = bootstrap_run), color = "grey50") +
          #geom_ribbon(aes(x = var, ymin = conf25, ymax = conf95), alpha = 0.5) +
          geom_line(aes(x = var, y = mean), color = "black", size = 1.5) +
          labs(x = vars) +
          theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8), panel.grid = element_line(color = "transparent"), axis.title = element_text(face = "bold"))} else {
            data %>% dplyr::filter(variable.name == vars) %>% dplyr::mutate(rown = ROWN) %>% dplyr::group_by(rown) %>% dplyr::summarise(mean = mean(yhat), conf25 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[1]), conf95 = ifelse(n() == 1, NA, t.test(yhat)$conf.int[2])) %>% dplyr::mutate(var = NUMS) %>% ggplot() +
              geom_ribbon(aes(x = var, ymin = conf25, ymax = conf95), alpha = 0.5) +
              geom_line(aes(x = var, y = mean), color = "black", size = 1.5) +
              labs(x = vars) +
              theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8), panel.grid = element_line(color = "transparent"), axis.title = element_text(face = "bold"))
          }} else {
            data %>% dplyr::filter(variable.name == vars) %>% ggplot() +
              geom_line(data = data %>% dplyr::filter(variable.name == vars), aes(x = x, y = yhat), color = "black", size = 1.5) +
              labs(x = vars) +
              theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8), panel.grid = element_line(color = "transparent"), axis.title = element_text(face = "bold"))  
          }})
  if(isTRUE(histogram)) {
    PLT2 <- lapply(vars, function(vars) {
      hist.data %>% dplyr::filter(variable.name == vars) %>% ggplot(aes(x = variable.value, y = value)) +
        geom_bar(stat = "identity", color = "black", fill = "grey85") +
        scale_y_continuous(position = "right") +
        labs(x = vars, y = "Frequency") +
        theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 0.8), panel.grid.major = element_line(color = "grey90"), axis.title = element_text(face = "bold"))
    })
    PLTS <- lapply(1:length(PLT2), function(j) {
      base_g <- ggplot_gtable(ggplot_build(PLT2[[j]]))
      overlay_g <- ggplot_gtable(ggplot_build(PLT[[j]]))
      plt_panel = c(subset(base_g$layout, name == "panel", se = t:r))
      pnl_ind = which(overlay_g$layout$name == "panel")
      leg_ind = which(overlay_g$layout$name == "axis-l")
      lab_ind <- which(overlay_g$layout$name == "ylab-l")
      final_grob = gtable_add_grob(base_g,
                                   overlay_g$grobs[[pnl_ind]],
                                   plt_panel$t,
                                   plt_panel$l,
                                   plt_panel$b,
                                   plt_panel$r, name = "a")
      final_grob <- gtable_add_cols(final_grob, widths = unit(20, "points"), pos = 0)
      final_grob <- gtable_add_cols(final_grob, widths = unit(20, "points"), pos = 0)
      final_grob = gtable_add_grob(final_grob,
                                   overlay_g$grobs[[leg_ind]],
                                   7, #plt_axis$t,
                                   3, #plt_axis$l,
                                   7, #plt_axis$b,
                                   3, #plt_axis$r,
                                   name = "b")
      final_grob <- gtable_add_grob(final_grob,
                                    overlay_g$grob[[lab_ind]],
                                    7,
                                    1,
                                    7,
                                    1,
                                    name = "c")
      final_grob
    })
    wrap_plots(PLTS)} else wrap_plots(PLT)}
