library(ggplot2)
load("pval_des.RData")
load("pval_lasso.RData")
load("pval_MCP.RData")
load("pval_oracle.RData")
load("pval_SCAD.RData")
load("pval_subset.RData")
load("pval_group.RData")

sigma_all = seq(from = 0.5, to = 1, by = 0.05)
nsigma =length(sigma_all)
signal = seq(from = 0, to = 1, by = 0.2)

alpha = 0.1
N = 1000

power_oracle = apply(Pval_oracle, c(2), function(x) sum(as.numeric(x<alpha))/N)
power_des = apply(Pval_des, c(2), function(x) sum(as.numeric(x<alpha))/N)
power_lasso = apply(Pval_lasso, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_subset = apply(Pval_subset, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_group = apply(Pval_group, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_SCAD = apply(Pval_SCAD, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_MCP = apply(Pval_MCP, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)

# print("power_oracle")
# print(power_oracle)
# print("power_des")
# print(power_des)
# print("power_lasso")
# print(power_lasso)
# print("power_SCAD")
# print(power_SCAD)
# print("power_MCP")
# print(power_MCP)
# print("power_subset")
# print(power_subset)
# print("power_group")
# print(power_group))

# --- Custom legend labels ---
label_names <- c(
  "lasso"          = "aCSS(lasso)",
  "SCAD"           = "aCSS(SCAD)",
  "MCP"            = "aCSS(MCP)",
  "subset"         = "aCSS(subset)",
  "group-SCAD"     = "aCSS(group-SCAD)",
  "oracle"         = "oCRT",
  "debiased-lasso" = "debiased-lasso"
)

# --- Distinct colors for all methods ---
color_set <- c(
  "oCRT"              = "#D4AC0D",  # Gold
  "debiased-lasso"    = "#E67E22",  # Orange
  "aCSS(lasso)"       = "#4E79A7",  # Blue
  "aCSS(SCAD)"        = "#F28E2B",  # Light Orange
  "aCSS(MCP)"         = "#59A14F",  # Green
  "aCSS(subset)"      = "#AF7AC5",  # Purple
  "aCSS(group-SCAD)"  = "#17BECF"   # Teal
)

# --- Same linetype within group: solid = oracle-based, dashed = aCSS ---
linetype_set <- c(
  "oCRT"              = "solid",
  "debiased-lasso"    = "solid",
  "aCSS(lasso)"       = "dashed",
  "aCSS(SCAD)"        = "dashed",
  "aCSS(MCP)"         = "dashed",
  "aCSS(subset)"      = "dashed",
  "aCSS(group-SCAD)"  = "dashed"
)

# --- Power Plot: Loop Over sigma_all ---
for (i in 1:length(sigma_all)) {
  
  plot_data <- data.frame(
    signal = rep(signal, times = 7),
    power = c(
      unlist(power_oracle),
      unlist(power_des),
      unlist(power_lasso[, i]),
      unlist(power_SCAD[, i]),
      unlist(power_MCP[, i]),
      unlist(power_group[, i]),
      unlist(power_subset[, i])
    ),
    method = rep(c("oracle", "debiased-lasso", "lasso", "SCAD", "MCP", "group-SCAD", "subset"), each = length(signal))
  )
  
  # Assign legend labels
  plot_data$label <- label_names[as.character(plot_data$method)]
  plot_data$label <- factor(plot_data$label, levels = unique(label_names))
  
  # Confidence ribbon
  plot_data$delta <- sqrt(plot_data$power * (1 - plot_data$power) / N)
  
  # Plot
  p <- ggplot(plot_data, aes(x = signal, y = power, color = label, linetype = label, fill = label)) +
    geom_line(size = 0.8) +
    geom_ribbon(aes(ymin = power - delta, ymax = power + delta), alpha = 0.3, color = NA) +
    scale_color_manual(values = color_set) +
    scale_fill_manual(values = color_set) +
    scale_linetype_manual(values = linetype_set) +
    labs(x = expression(beta[0]), y = "Power", color = "Method", linetype = "Method", fill = "Method") +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "red", size = 0.5) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 11),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.text = element_text(size = 8)
    )+guides(col = guide_legend(nrow = 3, theme = theme(legend.byrow = FALSE)))
  
  print(p)
  ggsave(filename = paste0("myplot", sigma_all[i], ".pdf"), plot = p, height = 5, width = 6, units = "in")
}

# --- Type I Error Plot ---
data <- data.frame(
  sigma = rep(sigma_all, 5),
  type1_error = c(
    power_lasso[1, ],
    power_SCAD[1, ],
    power_MCP[1, ],
    power_group[1, ],
    power_subset[1, ]
  ),
  method = rep(c( "lasso", "SCAD", "MCP", "group-SCAD", "subset"), each = length(sigma_all))
)

# Assign custom labels
data$label <- label_names[as.character(data$method)]
data$label <- factor(data$label, levels = unique(label_names))

# Confidence ribbon
data$delta <- sqrt(data$type1_error * (1 - data$type1_error) / N)

# Plot
plot <- ggplot(data, aes(x = sigma, y = type1_error, color = label, linetype = label, fill = label)) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = type1_error - delta, ymax = type1_error + delta), alpha = 0.3, color = NA) +
  scale_color_manual(values = color_set) +
  scale_fill_manual(values = color_set) +
  scale_linetype_manual(values = linetype_set) +
  labs(x = expression(sigma), y = "Type I Error", color = "Method", linetype = "Method", fill = "Method") +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "red", size = 0.5) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )+ 
  guides(col = guide_legend(nrow = 2, theme = theme(legend.byrow = TRUE)))

print(plot)
ggsave('type1_sparse.pdf', plot, width = 9, height = 4, units = "in")
ggsave('type1_sparse.png', plot, width = 7, height = 3, units = "in")



