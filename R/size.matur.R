#' A nonparametric method to estimate size at first maturity of marine exploited species.
#'
#' This function fits a Generalized Additive Model (GAM) to size vs maturity (in proportion) data. The function returns a table with the estimated size at first maturity from a percentage of interest of the user (usually 50%; L50, but also suited for 75%, 95%, etc) with the associated 95% confidence intervals. It also plots the data and the fitted curve with 95% confidence bands.
#' @param size a vector of size classes (total length, shell length, etc.).
#' @param prop.mat Proportion of mature individuals for each size class (must be the same lenght as size).
#' @param xlab Character string for the x-axis label of the generated scatter plot (for example total length (cm))
#' @param ylab same as xlab but for the y-axis label.
#' @param simsize.incr Increase intervals for the simulated vector of sizes, which is used to predict the curve from the fitted GAM. A small interval will create a smoother curve but may take longer time to compute.
#' @param p.car A pch value for the marker for the observed maturity data in the generated plot. For example p.car = 19 uses a solid black circle.
#' @param perc.size The percentage of the size to estimate the first maturity at. For example, perc.size = 0.50 estimates the size at first maturity of the 50% of the curve (L50).
#' @keywords Size at first maturity; Nonparametric statistics; Marine exploited species
#' @export
#' @examples attach(generosa) #geoduck clam data from (Aragon-Noriega 2015)
#' @examples size.matur(size = generosa$size.mm, prop.mat = generosa$mat.prop, xlab = "Shell size (mm)", ylab = "Proportion of matures", simsize.incr = 0.1, p.car = 19, perc.size = 0.5)
size.matur <- function(size, prop.mat, xlab, ylab,simsize.incr,p.car,perc.size){


  # main plot #
  plot(size, prop.mat, xlab=xlab, ylab=ylab, pch=p.car, cex.axis=1.3, cex.lab=1.3, las=1)

  legend("topleft", c("Fit", "95% C.I."), lty=c(1,2), lwd=c(2,1), cex=1.5, bty = "n")

#fit the GAM #
mod.fit <- mgcv::gam(prop.mat ~ s(size))


#create vector of simulated sizes #
size.sim <- seq(min(size), max(size),simsize.incr)
mod.pred <- mgcv::predict.gam(mod.fit, newdata = data.frame(size=size.sim), type="response", se.fit = T)
#Constrain the prediction and the C.I. to be [0:1]
mod.pred$fit <- ifelse(mod.pred$fit < 1, mod.pred$fit, 1)
mod.pred$fit <- ifelse(mod.pred$fit > 0, mod.pred$fit, 0)
low.ic <- mod.pred$fit-mod.pred$se.fit
hi.ic <- mod.pred$fit+mod.pred$se.fit

low.ic <- ifelse(low.ic > 0, low.ic, 0)
hi.ic <- ifelse(hi.ic < 1, hi.ic, 1)

#add lines of the fitted model and associated Bayesian credibility bands #
lines(size.sim, mod.pred$fit, lwd=2)
lines(size.sim, hi.ic, lty=2)
lines(size.sim, low.ic, lty=2)

#find the size at first maturity #
close.min <-  which.min(abs(mod.pred$fit - perc.size))
close.low <- which.min(abs((hi.ic) - perc.size))
close.hi <- which.min(abs((low.ic) - perc.size))

size.mat <-  size.sim[close.min]
low.size.mat <- size.sim[close.low]
hi.size.mat <- size.sim[close.hi]

#add the lines of the estimate size at first maturity and associated 95% confidence intervals#
segments(x0=size.mat,y0=-1,x1=size.mat,y1=mod.pred$fit[close.min], lty=1, col="gray31")
segments(x0=0,y0=perc.size,x1=size.sim[close.min],y1=perc.size, lty=1, col="gray31")
segments(x0=low.size.mat,y0=-1,x1=low.size.mat,y1=(hi.ic)[close.low], lty=2, col="gray31")

segments(x0=hi.size.mat,y0=-1,x1=hi.size.mat,y1=(low.ic)[close.hi], lty=2, col="gray31")


#create the a table with the results of the analysis #
df.1 <- data.frame(low.CI=low.size.mat,x=size.mat, high.CI=hi.size.mat, dev.exp = summary(mod.fit)$dev.expl)
names(df.1) <-c( "0.25%",paste("L",perc.size*100,sep = ""),"0.975%", "Exp. Dev.")
print(df.1)


}
