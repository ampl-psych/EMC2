rm(list = ls())
load("~/Documents/UVA/2022/EMC2/test_new_obj.RData")
devtools::load_all()

get_drift_angle <- function(x0, x1, y0, y1){
  width <- x1 - x0
  height <- y1 - y0
  w_plot <- par("pin")[1] / diff(par("usr")[1:2])
  h_plot <- par("pin")[2] / diff(par("usr")[3:4])
  return(atan((h_plot/w_plot) * (height/width)) * 180/pi)
}

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun <- function(d)d$S==d$lR
design <- make_design(data = forstmann,model=LBA,matchfun=matchfun,
                          formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                          contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
p_vector <- c(v=1.5,v_lMd=.75, sv_lMTRUE = log(.5), B = log(1.2), B_Eneutral = .2,
           B_Eaccuracy = .25, B_lRright = 0.3, A = log(.4), t0=log(.2))

pars <- mapped_par(p_vector, design)


data <- make_data(p_vector, design, n_trials = 100)
dadm <- design_model(data,design , verbose = FALSE)
dadm <- dadm[dadm$winner,]

get_cells <- function(df, fnames){
  cells <- df[,fnames,drop=FALSE]
  for (i in fnames) cells[,i] <- paste(i,cells[,i],sep="=")
  return(apply(cells,1,paste,collapse=" "))
}

get_cols <- function(cols, dups, i, name){
  if(sum(dups) == 1) return('black')
  else(
    return(cols[name])
  )
}




get_unq_names <- function(pars){
  dups <- duplicated(pars)
  if(sum(!dups) == 1) return(pars)
  nams <- names(pars)
  unq_nams <- nams[!dups]
  unq_split <- strsplit(unq_nams, split = " ")
  split <- strsplit(nams, split = " ")
  out_names <- rep("", length(split))
  for(i in 1:length(split[[1]])){
    is_unq <- length(unique(sapply(unq_split, function(x) x[i]))) != 1
    if(is_unq){
      if(i == 1){
        out_names <- paste0(out_names, sapply(split, function(x) x[i]))
      } else{
        out_names <- paste0(out_names, " ", sapply(split, function(x) x[i]))
      }
    }
  }
  names(pars) <- out_names
  return(pars)
}

plot_factors <- c("E")
vary_factor <- c("lM")
cells <- get_cells(dadm, plot_factors)
parcells <- get_cells(pars, plot_factors)
ucells <- unique(cells)

# To do, just choose which parameter to vary
make_lba_plot_helper <- function(data, pars, factors, main, fList = NULL){
  squash_dens <- .75
  alpha <- .45
  pR <- as.data.frame(table(data[,factors])/nrow(data))
  if(length(factors) == 1) colnames(pR)[1] <- factors
  pR <- pR$Freq[order(get_cells(pR, factors))]
  cells <- get_cells(data, factors)
  parcells <- get_cells(pars, factors)
  ucells <- unique(cells)
  vs <- list()
  bs <- list()
  t0s <- list()
  all_data <- list()
  for (i in sort(ucells)) {
    tmp_data <- data[cells == i,]
    tmp_pars <- pars[parcells == i,]
    vs[[i]] <- mean(tmp_pars$v)
    bs[[i]] <- mean(tmp_pars$b)
    t0s[[i]] <- mean(tmp_pars$t0)
    all_data[[i]] <- tmp_data
  }
  vs <- get_unq_names(unlist(vs))
  bs <- get_unq_names(unlist(bs))
  t0s <- get_unq_names(unlist(t0s))
  free_v <- !duplicated(vs)
  free_b <- !duplicated(bs)
  all_id_names <- unique(c(names(vs[free_v]), names(bs[free_b])))
  cols <- c("#8C2E89","#116214", "#D68E07","#200BA5",  "#02C6B2", "#C60202", "#788714", "#D33A7C") #colorspace::qualitative_hcl(length(all_id_names), "Pastel 1", c = 120)
  names(cols) <- all_id_names
  denses <- lapply(all_data, function(x) density(x$rt))
  t0 <- mean(t0s)
  x <- c(0, t0)
  y <- c(0, 0)
  tops <- numeric(length(denses))
  for(i in 1:length(denses)){
    denses[[i]]$y <- denses[[i]]$y*pR[i]*squash_dens + bs[i] + max(bs)/75
  }
  ylim <- c(-.1, max(sapply(denses, function(x) max(x$y))))
  max_rt <- quantile(data$rt, .975)
  plot(x, y, type = "l", xlim = c(0, max_rt), xlab = "", ylab = "", axes = F, ylim = ylim, lwd = 3, main = main)
  arrows(x0 = 0, x1 = max_rt, y0 = -.1, y1 = -.1, lwd = 3, length = .1)
  if(sum(free_v) != 1){
    dens_free <- free_v
    names_dens <- names(vs)
  } else{
    dens_free <- free_b
    names_dens <- names(bs)
  }

  for(i in 1:length(all_data)){
    if(free_b[i]) lines(c(0, max_rt), c(bs[i], bs[i]), lwd = 3,
                                 lty = 2, col = get_cols(cols, free_b, i, names(bs)[i]))
    rt <- min(bs)/vs[i] + t0
    if(free_v[i]) arrows(x0 = t0, x1 = rt, y0 = 0, y1 = min(bs)-.05, lwd = 3,
                                  length = .1, col = get_cols(cols, free_v, i, names(vs)[i]))
    lines(denses[[i]], lwd = 3, col = adjustcolor(get_cols(cols, dens_free, i, names_dens[i])))
    polygon(denses[[i]], col = adjustcolor(get_cols(cols, dens_free, i, names_dens[i]), alpha.f = alpha))
  }
  if(sum(free_v) != 1) legend("left", legend = names(vs[free_v]), lty = 1, lwd = 3, col = cols[names(vs[free_v])],
         cex = 1, bty = "n")
  if(sum(free_b) != 1 & !identical(names(bs)[free_b], names(vs)[free_v])){
    legend("bottomright", legend = names(bs[free_b]), lty = 1, lwd = 3, col = cols[names(bs[free_b])],
                                               cex = 1, bty = "n", inset = c(.05, .05))
  }
  title(xlab="RT", line=0, cex.lab = 1.5)
}


# debug(get_unq_names)

# debug(make_lba_plot_helper)
# par(mfrow = c(2,2))
# make_lba_plot_helper(dadm, pars, c("lM", "lR"), main = "test")

par(mfrow = c(1,2))
for (i in sort(ucells)) {
  tmp_data <- dadm[cells == i,]
  tmp_pars <- pars[parcells == i,]
  make_lba_plot_helper(tmp_data, tmp_pars, vary_factor, main = i)
}

# undebug(make_lba_plot_helper)


#
# LBA_plothelper <- function(props, pars){
#   rtsT <- data$rt[data$lM]
#   rtsF <- data$rt[!data$lM]
#   rtsT <- rtsT[rtsT < quantile(rtsT, .975)]
#   rtsF <- rtsF[rtsF < quantile(rtsF, .975)]
#   dT <- density(rtsT)
#   dF <- density(rtsF)
#   dT$y <- dT$y * mean(data$lM)
#   dF$y <- dF$y * (1-mean(data$lM))
#   rtT <- B/vT + t0
#   rtF <- B/vF + t0
#
#   x <- c(0, t0)
#   y <- c(0, 0)
#   ylim <- c(-.1, B + max(c(dT$y, dF$y)))
#   dT$y <- dT$y*squash + B
#   dF$y <- dF$y*squash + B
#
#
#   max_rt <- max(c(rtsT, rtsF))
#
#
#   plot(x, y, type = "l", xlim = c(0, max_rt), xlab = "", ylab = "", axes = F, ylim = ylim, lwd = 3)
#   arrows(x0 = 0, x1 = max_rt, y0 = -.1, y1 = -.1, lwd = 3, length = .1)
#   arrows(x0 = t0, x1 = rtT, y0 = 0, y1 = B-.05, lwd = 3, length = .1, col = "purple")
#   arrows(x0 = t0, x1 = rtF, y0 = 0, y1 = B-.05, lwd = 3, length = .1, col = "darkgreen")
#   lines(c(0, max_rt), c(B, B), lwd = 3, lty = 2)
#   title(xlab="RT", line=0, cex.lab = 1.5)
#   text(t0/2, .15, "t0", cex = 1.5)
#   text(t0/5, B + .15, "B", cex = 1.5)
#   legend(x = t0/5, y = B/1.05, legend = c("v_TRUE", "v_FALSE"), lty = 1, lwd = 3, col = c("purple", "darkgreen"),
#          cex = 1.5, bty = "n")
#
#   lines(dT, lwd = 3)
#   polygon(dT, col = adjustcolor("purple", alpha.f = alpha))
#
#   lines(dF, lwd = 3)
#   polygon(dF, col = adjustcolor("darkgreen", alpha.f = alpha))
#
# }
#
#
#
#
#
#
#
#
#
