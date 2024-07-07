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

design <- make_design(data = forstmann,model=DDM,
                          formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                          constants=c(s=log(1)))

p_vector <- c(v_Sleft = -2, v_Sright = 2, a = log(2), a_Eneutral = .3, a_Eaccuracy = .3, t0 = log(.2),
           Z = qnorm(.5), sv = 1, SZ = qnorm(.3))

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
vary_factor <- c("S")
cells <- get_cells(dadm, plot_factors)
parcells <- get_cells(pars, plot_factors)
ucells <- unique(cells)

# To do, just choose which parameter to vary
make_DDM_plot_helper <- function(data, pars, factors, main, fList = NULL){
  squash_dens <- .75
  alpha <- .45
  pR <- as.data.frame(table(data[,factors])/nrow(data))
  if(length(factors) == 1) colnames(pR)[1] <- factors
  pR <- pR$Freq[order(get_cells(pR, factors))]
  cells <- get_cells(data, factors)
  parcells <- get_cells(pars, factors)
  ucells <- unique(cells)
  vs <- list()
  as <- list()
  zs <- list()
  t0s <- list()
  all_data <- list()
  is_bot <- c()
  for (i in sort(ucells)) {
    tmp_data <- data[cells == i,]
    tmp_pars <- pars[parcells == i,]
    if(tmp_pars$S[1] == levels(data$S)[1]){
      is_bot[i] <- T
    } else{
      is_bot[i] <- F
    }
    vs[[i]] <- mean(tmp_pars$v)
    as[[i]] <- mean(tmp_pars$a)
    t0s[[i]] <- mean(tmp_pars$t0)
    zs[[i]] <- mean(tmp_pars$z)
    all_data[[i]] <- tmp_data
  }
  vs <- get_unq_names(unlist(vs))
  as <- get_unq_names(unlist(as))
  t0s <- get_unq_names(unlist(t0s))
  zs <- get_unq_names(unlist(zs))
  free_v <- !duplicated(vs)
  free_a <- !duplicated(as)
  all_id_names <- unique(c(names(vs[free_v]), names(as[free_a])))
  cols <- c("#8C2E89","#116214", "#D68E07","#200BA5",  "#02C6B2", "#C60202", "#788714", "#D33A7C") #colorspace::qualitative_hcl(length(all_id_names), "Pastel 1", c = 120)
  names(cols) <- all_id_names
  denses <- lapply(all_data, function(x) density(x$rt))
  t0 <- mean(t0s)
  Z <- mean(zs)
  x <- c(0, t0)
  y <- c(Z, Z)
  for(i in 1:length(denses)){
    if(is_bot[i]){
      denses[[i]]$y <- -1 * denses[[i]]$y*pR[i]*squash_dens - max(as)/75
    } else{
      denses[[i]]$y <- denses[[i]]$y*pR[i]*squash_dens + as[i] + max(as)/75
    }
  }
  ylim <- c(min(sapply(denses, function(x) min(x$y))) - max(as)/50, max(sapply(denses, function(x) max(x$y))))
  max_rt <- quantile(data$rt, .975)
  plot(x, y, type = "l", xlim = c(0, max_rt), xlab = "", ylab = "", axes = F, ylim = ylim, lwd = 3, main = main)
  arrows(x0 = 0, x1 = max_rt, y0 = ylim[1] - max(as)/50, y1 = ylim[1]  - max(as)/50, lwd = 3, length = .1)
  if(sum(free_v) != 1){
    dens_free <- free_v
    names_dens <- names(vs)
  } else{
    dens_free <- free_a
    names_dens <- names(as)
  }

  for(i in 1:length(all_data)){
    rt <- min(as)/abs(vs[i]) + t0
    if(is_bot[i]){
      if(free_a[i]) lines(c(0, max_rt), c(0, 0), lwd = 3,
                          lty = 2, col = get_cols(cols, free_a, i, names(as)[i]))
      if(free_v[i]) arrows(x0 = t0, x1 = rt, y0 = Z, y1 = 0 + .05, lwd = 3,
                           length = .1, col = get_cols(cols, free_v, i, names(vs)[i]))
    } else{
      if(free_a[i]) lines(c(0, max_rt), c(as[i], as[i]), lwd = 3,
                          lty = 2, col = get_cols(cols, free_a, i, names(as)[i]))
      rt <- min(as)/abs(vs[i]) + t0
      if(free_v[i]) arrows(x0 = t0, x1 = rt, y0 = Z, y1 = min(as)-.05, lwd = 3,
                           length = .1, col = get_cols(cols, free_v, i, names(vs)[i]))
    }

    lines(denses[[i]], lwd = 3, col = adjustcolor(get_cols(cols, dens_free, i, names_dens[i])))
    polygon(denses[[i]], col = adjustcolor(get_cols(cols, dens_free, i, names_dens[i]), alpha.f = alpha))
  }
  if(sum(free_v) != 1){
    legend("topright", legend = names(vs[free_v]), title = "v", lty = 1, lwd = 3, col = cols[names(vs[free_v])],
           cex = 1, bty = "n", inset = c(0, .15), title.cex = 1.25, title.adj = .2)
  } else{
    text(rt*.6 - min(as)/50, (rt-t0)*.6*vs[1] + min(as)/50, "v", cex = 1.5, srt = get_drift_angle(t0, rt, 0, min(as) - .05))
  }
  text(t0/2, min(as)/10, "t0", cex = 1.5)
  if(sum(free_a) != 1 & !identical(names(as)[free_a], names(vs)[free_v])){
    legend("bottomright", legend = names(as[free_a]), title = "b", lty = 1, lwd = 3, col = cols[names(as[free_a])],
           cex = 1, bty = "n", inset = c(.05, .05), title.adj = .2, title.cex = 1.5)
  }
  title(xlab="RT", line=0, cex.lab = 1.5)
}




# debug(get_unq_names)

# debug(make_lba_plot_helper)
# par(mfrow = c(2,2))
# make_lba_plot_helper(dadm, pars, c("lM", "lR"), main = "test")

undebug(make_DDM_plot_helper)
par(mfrow = c(1,2))
for (i in sort(ucells)) {
  tmp_data <- dadm[cells == i,]
  tmp_pars <- pars[parcells == i,]
  make_DDM_plot_helper(tmp_data, tmp_pars, vary_factor, main = i)
}

