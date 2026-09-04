compare(list(samples_LNR), BayesFactor = FALSE, type = "conditional")
compare(list(samples_LNR), BayesFactor = FALSE, type = "conditional", use_best_fit = FALSE)


ModelCompare <- compare(

  list(samples_LNR), BayesFactor = FALSE, type = "joint"

)



ModelCompare_Group <- compare(

  list(samples_LNR), BayesFactor = FALSE, type = "joint", group_only=TRUE

)



ModelCompare$DIC <- ModelCompare$DIC - ModelCompare_Group$DIC

ModelCompare$BPIC <- ModelCompare$BPIC - ModelCompare_Group$BPIC



getp <- function(IC) {

  IC <- -(IC - min(IC))/2

  exp(IC)/sum(exp(IC))

}



ModelCompare$wDIC <- getp(ModelCompare$DIC)

ModelCompare$wBPIC <- getp(ModelCompare$BPIC)
