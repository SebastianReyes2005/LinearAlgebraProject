original <- c(-0.041002, -0.028534, 0.000848, 0.009175, -0.006356,
              0.007206, 0.049740, 0.057683, 0.024539, 0.041740)

fitted <- c(-0.040999, -0.028522, 0.000815, 0.009258, -0.006479,
            0.007330, 0.049659, 0.057718, 0.024531, 0.041742)

days <- 1:10

df <- data.frame(
  Day = days,
  Original = original,
  Fitted = fitted
)

library(ggplot2)

df$RelativeError <- abs(df$Original - df$Fitted) / abs(df$Original) * 100

scale_factor <- max(abs(c(df$Original, df$Fitted))) / max(df$RelativeError)
df$RelativeErrorScaled <- df$RelativeError * scale_factor

ggplot(df, aes(x = Day)) +
  geom_line(aes(y = Original, color = "Original"), size = 1.2) +
  geom_point(aes(y = Original, color = "Original"), size = 2) +
  geom_line(aes(y = Fitted, color = "Fitted"), linetype = "dashed", size = 1.2) +
  geom_point(aes(y = Fitted, color = "Fitted"), shape = 1, size = 2.5) +
  geom_line(aes(y = RelativeErrorScaled, color = "Relative Error (%)"), linetype = "dotdash", size = 1.2) +
  scale_color_manual(values = c("Original" = "steelblue", "Fitted" = "firebrick", "Relative Error (%)" = "darkgreen")) +
  scale_x_continuous(breaks = days) +
  scale_y_continuous(
    name = "Overnight Movement",
    sec.axis = sec_axis(~ . / scale_factor, name = "Relative Error (%)")
  ) +
  labs(
    title = "Original vs Fitted vs Relative Error",
    x = "Day",
    color = "Legend"
  ) +
  theme_minimal()