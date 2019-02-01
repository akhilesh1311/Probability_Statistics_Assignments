library(package=ISLR, lib.loc = "~/Documents/probStats/finals/")
library(UpSetR)

FDR <- function(pvalues) {
  alpha <- seq(0, 0.2, 0.001)
  alpha <- sort(alpha, decreasing = TRUE)
  R <- 101
  x <- 0
  while(R>100) {
    x <- x+1
    a <- alpha[x]
    ps <- sort(pvalues, decreasing = FALSE)
    m <- length(pvalues)
    l <- c((1:m)*a/m)
    com <- ps - l
    R <- which(com<=0)[length(which(com<=0))]
    if(length(R) == 0) {
      R <- 0
    }
  }
  Tvalue <- ps[R]
  rej <- which(pvalues<Tvalue)
  # print(rej)
  return(rej)
  # print(R)
  # print(x)
}

PValue_Generator <- function(lis) {
  mu_list <- apply(lis[2:nrow(lis),], 2, mean)
  sd_list <- apply(lis[2:nrow(lis),], 2, sd)
  tlist <- mu_list/(sd_list/sqrt(nrow(lis) - 1))
  tlist <- na.omit(tlist)
  plist.suppressed <- pt(tlist, nrow(lis)-2)
  plist.active <- 1 - pt(tlist, nrow(lis)-2)
  list(plist.active = plist.active, plist.suppressed = plist.suppressed)
}

Main_function <- function() {
  breast.mat <- matrix(ncol = 6830)
  breast.rows <- which(NCI60$labs %in% "BREAST")
  for (i in breast.rows) {
    breast.mat <- rbind(breast.mat, NCI60$data[i,])
  }
  
  cns.mat <- matrix(ncol = 6830)
  cns.rows <- which(NCI60$labs %in% "CNS")
  for (i in cns.rows) {
    cns.mat <- rbind(cns.mat, NCI60$data[i,])
  }
  
  colon.mat <- matrix(ncol = 6830)
  colon.rows <- which(NCI60$labs %in% "COLON")
  for (i in colon.rows) {
    colon.mat <- rbind(colon.mat, NCI60$data[i,])
  }
  
  ka.repro.mat <- matrix(ncol = 6830)
  ka.repro.rows <- which(NCI60$labs %in% "K562A-repro")
  for (i in ka.repro.rows) {
    ka.repro.mat <- rbind(ka.repro.mat, NCI60$data[i,])
  }
  
  kb.repro.mat <- matrix(ncol = 6830)
  kb.repro.rows <- which(NCI60$labs %in% "K562B-repro")
  for (i in kb.repro.rows) {
    kb.repro.mat <- rbind(kb.repro.mat, NCI60$data[i,])
  }
  
  leukemia.mat <- matrix(ncol = 6830)
  leukemia.rows <- which(NCI60$labs %in% "LEUKEMIA")
  for (i in leukemia.rows) {
    leukemia.mat <- rbind(leukemia.mat, NCI60$data[i,])
  }
  
  ma.repro.mat <- matrix(ncol = 6830)
  ma.repro.rows <- which(NCI60$labs %in% "MCF7A-repro")
  for (i in ma.repro.rows) {
    ma.repro.mat <- rbind(ma.repro.mat, NCI60$data[i,])
  }
  
  md.repro.mat <- matrix(ncol = 6830)
  md.repro.rows <- which(NCI60$labs %in% "MCF7D-repro")
  for (i in md.repro.rows) {
    md.repro.mat <- rbind(md.repro.mat, NCI60$data[i,])
  }
  
  melanoma.mat <- matrix(ncol = 6830)
  melanoma.rows <- which(NCI60$labs %in% "MELANOMA")
  for (i in melanoma.rows) {
    melanoma.mat <- rbind(melanoma.mat, NCI60$data[i,])
  }
  
  nsclc.mat <- matrix(ncol = 6830)
  nsclc.rows <- which(NCI60$labs %in% "NSCLC")
  for (i in nsclc.rows) {
    nsclc.mat <- rbind(nsclc.mat, NCI60$data[i,])
  }
  
  ovarian.mat <- matrix(ncol = 6830)
  ovarian.rows <- which(NCI60$labs %in% "OVARIAN")
  for (i in ovarian.rows) {
    ovarian.mat <- rbind(ovarian.mat, NCI60$data[i,])
  }
  
  prostate.mat <- matrix(ncol = 6830)
  prostate.rows <- which(NCI60$labs %in% "PROSTATE")
  for (i in prostate.rows) {
    prostate.mat <- rbind(prostate.mat, NCI60$data[i,])
  }
  
  renal.mat <- matrix(ncol = 6830)
  renal.rows <- which(NCI60$labs %in% "RENAL")
  for (i in renal.rows) {
    renal.mat <- rbind(renal.mat, NCI60$data[i,])
  }
  
  unknown.mat <- matrix(ncol = 6830)
  unknown.rows <- which(NCI60$labs %in% "UNKNOWN")
  for (i in unknown.rows) {
    unknown.mat <- rbind(unknown.mat, NCI60$data[i,])
  }
  
  breast.supp <- PValue_Generator(breast.mat)$plist.suppressed
  breast.active <- PValue_Generator(breast.mat)$plist.active
  print("breast.suppressed")
  breast.supp.ans <- FDR(breast.supp)
  print("breast.active")
  breast.active.ans <- FDR(breast.active)
  
  cns.supp <- PValue_Generator(cns.mat)$plist.suppressed
  cns.active <- PValue_Generator(cns.mat)$plist.active
  print("cns.suppressed")
  cns.supp.ans <- FDR(cns.supp)
  print("cns.active")
  cns.active.ans <- FDR(cns.active)
  
  colon.supp <- PValue_Generator(colon.mat)$plist.suppressed
  colon.active <- PValue_Generator(colon.mat)$plist.active
  print("colon.suppressed")
  colon.supp.ans <- FDR(colon.supp)
  print("colon.active")
  colon.active.ans <- FDR(colon.active)
  
  leukemia.supp <- PValue_Generator(leukemia.mat)$plist.suppressed
  leukemia.active <- PValue_Generator(leukemia.mat)$plist.active
  print("leukemia.suppressed")
  leukemia.supp.ans <- FDR(leukemia.supp)
  print("leukemia.active")
  leukemia.active.ans <- FDR(leukemia.active)
  
  melanoma.supp <- PValue_Generator(melanoma.mat)$plist.suppressed
  melanoma.active <- PValue_Generator(melanoma.mat)$plist.active
  print("melanoma.suppressed")
  melanoma.supp.ans <- FDR(melanoma.supp)
  print("melanoma.active")
  melanoma.active.ans <- FDR(melanoma.active)
  
  nsclc.supp <- PValue_Generator(nsclc.mat)$plist.suppressed
  nsclc.active <- PValue_Generator(nsclc.mat)$plist.active
  print("nsclc.suppressed")
  nsclc.supp.ans <- FDR(nsclc.supp)
  print("nsclc.active")
  nsclc.active.ans <- FDR(nsclc.active)
  
  ovarian.supp <- PValue_Generator(ovarian.mat)$plist.suppressed
  ovarian.active <- PValue_Generator(ovarian.mat)$plist.active
  print("ovarian.suppressed")
  ovarian.supp.ans <- FDR(ovarian.supp)
  print("ovarian.active")
  ovarian.active.ans <- FDR(ovarian.active)
  
  prostate.supp <- PValue_Generator(prostate.mat)$plist.suppressed
  prostate.active <- PValue_Generator(prostate.mat)$plist.active
  print("prostate.suppressed")
  prostate.supp.ans <- FDR(prostate.supp)
  print("prostate.active")
  prostate.active.ans <- FDR(prostate.active)
  
  renal.supp <- PValue_Generator(renal.mat)$plist.suppressed
  renal.active <- PValue_Generator(renal.mat)$plist.active
  print("renal.suppressed")
  renal.supp.ans <- FDR(renal.supp)
  print("renal.active")
  renal.active.ans <- FDR(renal.active)
  
  bar.names <- c("Breast", "CNS", "Colon", "Leukemia", "Melanoma", "Nsclc", "Ovarian", "Prostate", "Renal")
  barplot(c(length(breast.active.ans), length(cns.active.ans), length(colon.active.ans), length(leukemia.active.ans), length(melanoma.active.ans), length(nsclc.active.ans),
            length(ovarian.active.ans), length(prostate.active.ans), length(renal.active.ans)), names.arg = bar.names, 
          main = "Bar Plot for number of active genes causing cancer of given type",xlab = "Cancer type", ylab = "Number of genes")
  
  barplot(c(length(breast.supp.ans), length(cns.supp.ans), length(colon.supp.ans), length(leukemia.supp.ans), length(melanoma.supp.ans), length(nsclc.supp.ans),
            length(ovarian.supp.ans), length(prostate.supp.ans), length(renal.supp.ans)), names.arg = bar.names,
          main = "Bar plot for number of suppressed genes for cancer of given type", xlab = "Cancer type", ylab = "Number of genes")
  
  # n1 <- max(length(breast.active.ans), length(cns.active.ans), length(colon.active.ans), length(leukemia.active.ans), length(melanoma.active.ans), length(nsclc.active.ans),
  #           length(ovarian.active.ans), length(prostate.active.ans), length(renal.active.ans))
  # 
  # 
  # length(breast.active.ans) <- n1
  # length(cns.active.ans) <- n1
  # length(colon.active.ans) <- n1
  # length(leukemia.active.ans) <- n1
  # length(melanoma.active.ans) <- n1
  # length(nsclc.active.ans) <- n1
  # length(ovarian.active.ans) <- n1
  # length(prostate.active.ans) <- n1
  # length(renal.active.ans) <- n1
  # 
  # n2 <- max(length(breast.supp.ans), length(cns.supp.ans), length(colon.supp.ans), length(leukemia.supp.ans), length(melanoma.supp.ans), length(nsclc.supp.ans), 
  #           length(ovarian.supp.ans), length(prostate.supp.ans), length(renal.supp.ans))
  # 
  # length(breast.supp.ans) <- n2
  # length(cns.supp.ans) <- n2
  # length(colon.supp.ans) <- n2
  # length(leukemia.supp.ans) <- n2
  # length(melanoma.supp.ans) <- n2
  # length(nsclc.supp.ans) <- n2
  # length(ovarian.supp.ans) <- n2
  # length(prostate.supp.ans) <- n2
  # length(renal.supp.ans) <- n2
  # 
  # active.df <<- data.frame(breast.active.ans, cns.active.ans, colon.active.ans, leukemia.active.ans, melanoma.active.ans, nsclc.active.ans,
  #                         ovarian.active.ans, prostate.active.ans, renal.active.ans)
  # supp.df <<- data.frame(breast.supp.ans, cns.supp.ans, colon.supp.ans, leukemia.supp.ans, melanoma.supp.ans, nsclc.supp.ans,
  #                       ovarian.supp.ans, prostate.supp.ans, renal.supp.ans)


  # write.csv(active.df, "~/Documents/probStats/finals/active.csv")
  # write.csv(supp.df, "~/Documents/probStats/finals/supp.csv")
  #finding the common elements between two cancer types
  # common.two.active <<- sapply(seq_len(length(active.df)), function(x)
  #   sapply(seq_len(length(active.df)), function(y) intersect(unlist(active.df[x]), unlist(active.df[y]))))
  # 
  # common.two.supp <<- sapply(seq_len(length(supp.df)), function(x)
  #   sapply(seq_len(length(supp.df)), function(y) intersect(unlist(supp.df[x]), unlist(supp.df[y]))))
  
  active.data <- read.csv("~/Documents/probStats/finals/active12lol.csv", header=T, sep="," )
  upset(active.data, nsets = 9, sets = c("breast.active.ans", "renal.active.ans", "prostate.active.ans", "cns.active.ans", "colon.active.ans", "melanoma.active.ans", "leukemia.active.ans", "nsclc.active.ans"), sets.bar.color = "#56B4E9")
  

  supp.data <- read.csv("~/Documents/probStats/finals/supp1lol.csv", header=T, sep="," )
  upset(supp.data,  sets.bar.color = "#56B4E9")
}