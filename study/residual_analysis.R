

# run this script after running study.R

library("ggplot2")
library("cowplot")
library("tikzDevice")

residuals <- resid(fit); residuals <- residuals / sd(residuals)
g2 <- ggplot() + theme_minimal() + coord_equal(1) + 
  geom_qq(aes(sample = residuals)) +
  theme(panel.grid.minor = element_line(color = NA)) + 
  geom_abline(slope = 1, intercept = 0, colour = "red") + 
  ylab("residuals quantiles")

ranef <- as.numeric(getME(fit, "b")); ranef <- ranef / sd(ranef)
g3 <- ggplot() + theme_minimal() + coord_fixed(1) + 
  geom_qq(aes(sample = ranef)) +
  theme(panel.grid.minor = element_line(color = NA)) + 
  geom_abline(slope = 1, intercept = 0, colour = "red") + 
  ylab("estimated random effect quantiles")

plot_grid(g2, g3, align = "h", nrow = 1, 
          rel_heights = c(1/2, 1/2))


g4 <- ggplot(prisons, aes(x = county_log_mortality, 
                    y = resid(fit) )) + theme_minimal() + #coord_fixed(1) + 
  geom_point() + labs(x = "sdf")

df <- tibble(x = 1:45, y = ranef(fit)[[1]][[1]])
g5 <- ggplot(df, aes(x, y)) + theme_minimal() +geom_point()


plot_grid(g4, g5, align = "h", nrow = 1, 
          rel_heights = c(1/2, 1/2))
