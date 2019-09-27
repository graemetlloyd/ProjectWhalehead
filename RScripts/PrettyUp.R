p <- ggplot(fi, aes(factor(Bayes), MIG))
b <- p + geom_boxplot() + theme_bw() +scale_x_discrete(breaks = c(0,1), labels = c("Parsimony", "Bayes")) + xlab(NULL)
 q <- p + geom_violin(scale = "count") + geom_jitter(height = 0, width = 0.2, alpha = .4) +theme_bw() +scale_x_discrete(breaks = c(0,1), labels = c("Parsimony", "Bayes")) + ylab(NULL) + xlab(NULL)

 yp <- custom_plots(CombinedTrees = "Data/CombinedTrees/Yates_2003a.nex", StratFile = "Data/StratCongruenceResults/InputTreeResults/Yates_2003a.txt")
 
 grid.arrange(arrangeGrob(b,q, ncol=2, nrow=1), arrangeGrob(yp, ncol=1, nrow=1), heights = c(4,4))
