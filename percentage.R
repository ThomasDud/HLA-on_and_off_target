library(ggplot2)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
genes <- c("HLA-A","HLA-B","HLA-C","HLA-DQA1","HLA-DQB1","HLA-DRB1")

df_target <- lapply(genes, function(gene) {
    df <- read.table(paste0("on_targets/", gene, "_target.tsv"), sep="\t")
    df[order(df[[1]]),]
})

df_off <- lapply(genes, function(gene) {
    df <- read.table(paste0("off_targets/",gene, "_off_sum.tsv"), sep="\t")
    df[order(df[[1]]),]
})

df_perc_list <- mapply(function(gene, on, off) {
    df_perc <- data.frame(
    samples = on[,1],
    on = on[,2],
    off= off[,2],
    percentage = on[,2] / (on[,2] + off[,2]),
    Gene=gene
    )
    write.table(
        df_perc,
        file = paste0("tables/", gene,"_on_target_percentage.tsv"),
        sep = "\t",
        row.names=FALSE
    )
    return(df_perc)
}, genes, df_target, df_off,
SIMPLIFY=FALSE
)

# Distribution of the percentages for each of the genes

mapply(function(df, gene) {
    plot <- ggplot(df, aes(x=samples,y=percentage)) + 
    geom_bar(stat="identity", fill="skyblue", color="black") +
    labs(
        title = paste0("distribution of the percentage of ", gene, " on target"),
        x="Samples",
        y="Percentage on target"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    ggsave(paste0("plots/", gene, "_distribution.png"), plot = plot, width = 8, height = 6, dpi = 300)
}, df_perc_list, genes)

# Boxplot to summarize the distributions

df_all <- do.call(rbind, df_perc_list)

boxplot <- ggplot(df_all, aes(x = Gene, y = percentage)) +
geom_boxplot(fill = "lightgreen") +
labs(
    title = "On-target % per gene",
    x = "Gene",
    y = "Percentage"
  ) +
  theme_minimal()

ggsave("plots/boxplot.png", plot = boxplot, width = 8, height = 6, dpi = 300)