################################
## Graphs results simulations ##
################################

library(ggplot2)
library(cowplot)

###################

tLMs <- c("s = 2","s = 4")
tHors <- c("w = 1","w = 2")

############
## Small1 ##
############

load("Scenarios/Small1/Small1_resu.RData")

BS.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  BS.list <- lapply(res, FUN = function(x) return(x$BS[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(BS.list[[1]]), each = nrow(BS.list[[1]])), length(BS.list)),
                   "BS" = unlist(BS.list),
                   "tLM" = rep(tLMs[i], ncol(BS.list[[1]])*length(BS.list)),
                   "tHor" = rep(tHors, ncol(BS.list[[1]])*length(BS.list)))

  BS.mat.ggplot <- rbind(BS.mat.ggplot, df)

}

BS.mat.ggplot <- BS.mat.ggplot[-which(BS.mat.ggplot$Submethod=="pred_KM"),]

BS.mat.ggplot$Method <- NA
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_dyn1","pred_dyn2","pred_dyn3","pred_dyn4"))] <- "DynForest"
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("surv_JM_intslope"))] <- "JMbayes"

###################

p.BS <- ggplot(aes(x = Submethod, y = BS), data = BS.mat.ggplot[which(BS.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c("DynForest (mtry=1)","DynForest (mtry=2)","DynForest (mtry=3)","DynForest (mtry=4)","JMbayes")) +
  scale_y_continuous(breaks=c(0.02, 0.05, 0.08, 0.11)) +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  xlab("") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold"))


AUC.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  AUC.list <- lapply(res, FUN = function(x) return(x$AUC[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(AUC.list[[1]]), each = nrow(AUC.list[[1]])), length(AUC.list)),
                   "AUC" = unlist(AUC.list),
                   "tLM" = rep(tLMs[i], ncol(AUC.list[[1]])*length(AUC.list)),
                   "tHor" = rep(tHors, ncol(AUC.list[[1]])*length(AUC.list)))

  AUC.mat.ggplot <- rbind(AUC.mat.ggplot, df)

}

AUC.mat.ggplot <- AUC.mat.ggplot[-which(AUC.mat.ggplot$Submethod=="pred_KM"),]

AUC.mat.ggplot$Method <- NA
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_dyn1","pred_dyn2","pred_dyn3","pred_dyn4"))] <- "DynForest"
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("surv_JM_intslope"))] <- "JMbayes"


p.AUC <- ggplot(aes(x = Submethod, y = AUC), data = AUC.mat.ggplot) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c("DynForest (mtry=1)","DynForest (mtry=2)","DynForest (mtry=3)","DynForest (mtry=4)","JMbayes")) +
  scale_y_continuous(breaks=c(0.5, 0.6, 0.7, 0.8)) +
  facet_grid(tLM ~ tHor, scales = "free") +
  xlab("") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold")) +
  coord_flip()

pB2 <- plot_grid(p.BS, p.AUC, nrow = 1, rel_widths = c(1.3,1))

ggsave(filename = "Graphs/Small1_resu.pdf",
       #width = 170, height = 130, units = "mm",
       width = 7, height = 4,
       dpi = 300,
       plot = pB2,
       device = "pdf")

############
## Small2 ##
############

load("Scenarios/Small2/Small2_resu.RData")

BS.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  BS.list <- lapply(res, FUN = function(x) return(x$BS[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(BS.list[[1]]), each = nrow(BS.list[[1]])), length(BS.list)),
                   "BS" = unlist(BS.list),
                   "tLM" = rep(tLMs[i], ncol(BS.list[[1]])*length(BS.list)),
                   "tHor" = rep(tHors, ncol(BS.list[[1]])*length(BS.list)))

  BS.mat.ggplot <- rbind(BS.mat.ggplot, df)

}

BS.mat.ggplot <- BS.mat.ggplot[-which(BS.mat.ggplot$Submethod=="pred_KM"),]

BS.mat.ggplot$Method <- NA
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_dyn1","pred_dyn2","pred_dyn3","pred_dyn4"))] <- "DynForest"
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("surv_JM_intslope"))] <- "JMbayes"

###################

p.BS <- ggplot(aes(x = Submethod, y = BS), data = BS.mat.ggplot[which(BS.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c("DynForest (mtry=1)","DynForest (mtry=2)","DynForest (mtry=3)","DynForest (mtry=4)","JMbayes")) +
  scale_y_continuous(breaks=c(0.06, 0.09, 0.12, 0.15)) +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  xlab("") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold"))


AUC.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  AUC.list <- lapply(res, FUN = function(x) return(x$AUC[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(AUC.list[[1]]), each = nrow(AUC.list[[1]])), length(AUC.list)),
                   "AUC" = unlist(AUC.list),
                   "tLM" = rep(tLMs[i], ncol(AUC.list[[1]])*length(AUC.list)),
                   "tHor" = rep(tHors, ncol(AUC.list[[1]])*length(AUC.list)))

  AUC.mat.ggplot <- rbind(AUC.mat.ggplot, df)

}

AUC.mat.ggplot <- AUC.mat.ggplot[-which(AUC.mat.ggplot$Submethod=="pred_KM"),]

AUC.mat.ggplot$Method <- NA
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_dyn1","pred_dyn2","pred_dyn3","pred_dyn4"))] <- "DynForest"
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("surv_JM_intslope"))] <- "JMbayes"


p.AUC <- ggplot(aes(x = Submethod, y = AUC), data = AUC.mat.ggplot) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c("DynForest (mtry=1)","DynForest (mtry=2)","DynForest (mtry=3)","DynForest (mtry=4)","JMbayes")) +
  #scale_y_continuous(breaks=c(0.5, 0.6, 0.7, 0.8)) +
  facet_grid(tLM ~ tHor, scales = "free") +
  xlab("") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold")) +
  coord_flip()

pC <- plot_grid(p.BS, p.AUC, nrow = 1, rel_widths = c(1.3,1))

ggsave(filename = "Graphs/Small2_resu.pdf",
       #width = 170, height = 130, units = "mm",
       width = 7, height = 4,
       dpi = 300,
       plot = pC,
       device = "pdf")

###

pBC <- plot_grid(pB2, pC, ncol = 1, rel_widths = c(1,1), labels="AUTO")

ggsave(filename = "Graphs/Small_resu.pdf",
       #width = 170, height = 130, units = "mm",
       width = 7, height = 7,
       dpi = 300,
       plot = pBC,
       device = "pdf")


############
## Large1 ##
############

load("Scenarios/Large1/Large1_resu.RData")

BS.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  BS.list <- lapply(res, FUN = function(x) return(x$BS[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(BS.list[[1]]), each = nrow(BS.list[[1]])), length(BS.list)),
                   "BS" = unlist(BS.list),
                   "tLM" = rep(tLMs[i], ncol(BS.list[[1]])*length(BS.list)),
                   "tHor" = rep(tHors, ncol(BS.list[[1]])*length(BS.list)))

  BS.mat.ggplot <- rbind(BS.mat.ggplot, df)

}

BS.mat.ggplot <- BS.mat.ggplot[-which(BS.mat.ggplot$Submethod%in%c("pred_KM","pred_2step_PACE_dyn")),]

BS.mat.ggplot$Method <- NA
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_dyn"))] <- "DynForest"
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_2step_dyn"))] <- "Regression calibration"


p.BS <- ggplot(aes(x = Submethod, y = BS), data = BS.mat.ggplot[which(BS.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c(paste0("Regression", "\n" ,"calibration"),"DynForest")) +
  scale_y_continuous(breaks=c(0.04, 0.06, 0.08, 0.10, 0.12)) +
  xlab("") +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold"))


AUC.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  AUC.list <- lapply(res, FUN = function(x) return(x$AUC[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(AUC.list[[1]]), each = nrow(AUC.list[[1]])), length(AUC.list)),
                   "AUC" = unlist(AUC.list),
                   "tLM" = rep(tLMs[i], ncol(AUC.list[[1]])*length(AUC.list)),
                   "tHor" = rep(tHors, ncol(AUC.list[[1]])*length(AUC.list)))

  AUC.mat.ggplot <- rbind(AUC.mat.ggplot, df)

}

AUC.mat.ggplot <- AUC.mat.ggplot[-which(AUC.mat.ggplot$Submethod%in%c("pred_KM","pred_2step_PACE_dyn")),]

AUC.mat.ggplot$Method <- NA
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_dyn"))] <- "DynForest"
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_2step_dyn"))] <- "Regression calibration"


p.AUC <- ggplot(aes(x = Submethod, y = AUC), data = AUC.mat.ggplot[which(AUC.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  xlab("") +
  ylim(c(0.35, 0.95)) +
  scale_x_discrete(labels = c(paste0("Regression", "\n" ,"calibration"),"DynForest")) +
  #scale_y_continuous(breaks=c(0.4, 0.6, 0.8)) +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"))

pD2 <- plot_grid(p.BS, p.AUC, nrow = 1, rel_widths = c(1.2,1))

ggsave(filename = "Graphs/Large1_resu.pdf",
       #width = 170, height = 110, units = "mm",
       width = 7, height = 4,
       dpi = 300,
       plot = pD2,
       device = "pdf")

############
## Large2 ##
############

load("Scenarios/Large2/Large2_resu.RData")

BS.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  BS.list <- lapply(res, FUN = function(x) return(x$BS[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(BS.list[[1]]), each = nrow(BS.list[[1]])), length(BS.list)),
                   "BS" = unlist(BS.list),
                   "tLM" = rep(tLMs[i], ncol(BS.list[[1]])*length(BS.list)),
                   "tHor" = rep(tHors, ncol(BS.list[[1]])*length(BS.list)))

  BS.mat.ggplot <- rbind(BS.mat.ggplot, df)

}

BS.mat.ggplot <- BS.mat.ggplot[-which(BS.mat.ggplot$Submethod%in%c("pred_KM","pred_2step_PACE_dyn")),]

BS.mat.ggplot$Method <- NA
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_dyn"))] <- "DynForest"
BS.mat.ggplot$Method[which(BS.mat.ggplot$Submethod%in%c("pred_2step_dyn"))] <- "Regression calibration"


p.BS <- ggplot(aes(x = Submethod, y = BS), data = BS.mat.ggplot[which(BS.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  scale_x_discrete(labels = c(paste0("Regression", "\n" ,"calibration"),"DynForest")) +
  scale_y_continuous(breaks=c(0.07, 0.10, 0.13, 0.16)) +
  xlab("") +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold"))


AUC.mat.ggplot <- NULL

for (i in seq(length(tLMs))){

  AUC.list <- lapply(res, FUN = function(x) return(x$AUC[[i]]))

  df <- data.frame("Submethod" = rep(rep(colnames(AUC.list[[1]]), each = nrow(AUC.list[[1]])), length(AUC.list)),
                   "AUC" = unlist(AUC.list),
                   "tLM" = rep(tLMs[i], ncol(AUC.list[[1]])*length(AUC.list)),
                   "tHor" = rep(tHors, ncol(AUC.list[[1]])*length(AUC.list)))

  AUC.mat.ggplot <- rbind(AUC.mat.ggplot, df)

}

AUC.mat.ggplot <- AUC.mat.ggplot[-which(AUC.mat.ggplot$Submethod%in%c("pred_KM","pred_2step_PACE_dyn")),]

AUC.mat.ggplot$Method <- NA
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_dyn"))] <- "DynForest"
AUC.mat.ggplot$Method[which(AUC.mat.ggplot$Submethod%in%c("pred_2step_dyn"))] <- "Regression calibration"


p.AUC <- ggplot(aes(x = Submethod, y = AUC), data = AUC.mat.ggplot[which(AUC.mat.ggplot$Submethod!="pred_2step_PACE_dyn"),]) +
  geom_boxplot(aes(fill = Method)) +
  xlab("") +
  scale_x_discrete(labels = c(paste0("Regression", "\n" ,"calibration"),"DynForest")) +
  scale_y_continuous(breaks=c(0.4, 0.6, 0.8)) +
  guides(fill = FALSE) +
  facet_grid(tLM ~ tHor, scales = "free") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"))

pE <- plot_grid(p.BS, p.AUC, nrow = 1, rel_widths = c(1.2,1))

ggsave(filename = "Graphs/Large2_resu.pdf",
       #width = 170, height = 110, units = "mm",
       width = 7, height = 4,
       dpi = 300,
       plot = pE,
       device = "pdf")

###

pDE <- plot_grid(pD2, pE, ncol = 1, rel_widths = c(1,1), labels="AUTO")

ggsave(filename = "Graphs/Large_resu.pdf",
       #width = 170, height = 130, units = "mm",
       width = 7, height = 7,
       dpi = 300,
       plot = pDE,
       device = "pdf")

###