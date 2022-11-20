x <- fread('~/Documents/data/CBGdata/abstract_exports/export_admissionDuration_4_days_hypothresh_NAs_included_3_ratio_2_WITH_LOO_PRED.csv')

x_no_prior_hypo <- x[training_min > 3]
x_hr <- x_no_prior_hypo[prediction > 0.4]
boxplot(x_hr[, -c('id')], ylim = c(0, 28), las=3, varwidth=T)

export=x

plot(export$cV, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))

plot(export$min_by_day_20, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))
plot(export$min_by_day_3, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))
plot(export$min_by_day_1, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))

boxplot(export$prediction ~ cut(export$min_by_day_20, 100), varwidth=T)
boxplot(export$prediction ~ cut(export$min_by_day_13, 100), varwidth=T)
boxplot(export$prediction ~ cut(export$min_by_day_1, 100), varwidth=T)

boxplot(export$prediction ~ cut(export$max_by_day_13, 100), varwidth=T)
boxplot(export$prediction ~ cut(export$max_by_day_7, 100), varwidth=T)
boxplot(export$prediction ~ cut(export$max_by_day_1, 100), varwidth=T)

boxplot(export$prediction ~ export$day_N_13, varwidth=T)
boxplot(export$prediction ~ export$day_N_7, varwidth=T)
boxplot(export$prediction ~ export$day_N_1, varwidth=T)

boxplot(export$prediction ~ cut(export$gradient, 100), varwidth=T)
boxplot(export$prediction ~ cut(export$cV, 100), varwidth=T)
