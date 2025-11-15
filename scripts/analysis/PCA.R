rm(list=ls())      # Se vacía el entorno de trabajo

# Se cargan las librerías necesarias
library(DESeq2)   # Librería para la transformación VST
library(stats)    # Librería para el PCA
library(ggplot2)  # Librería para hacer el gráfico

# Se abren los datos a analizar
datos <- read.csv("data/processed_data/GSE_unificado_ordenado.csv", header = TRUE, row.names = 1, check.names = FALSE)
datos <- datos[, !startsWith(colnames(datos), "anno_")] # Se eliminanan las columnas no necesarias

#Agrupamos las muestras por tipo de tratamiento y género
muestras <- c(rep("F_control", 3), rep("F_40μg/L_DBAN", 3), rep("F_200μg/L_DBAN", 3), rep("M_control", 3), rep("M_40μg/L_DBAN", 3)) 

colData <- data.frame(row.names = colnames(datos), Grupo = muestras)
dds <- DESeqDataSetFromMatrix(countData = datos,
                              colData = colData,
                              design = ~ Grupo)

# Se filtran los genes con muy pocos counts
dds <- dds[rowSums(counts(dds)) > 10, ]

vsd <- vst(dds, blind = TRUE)  # blind=TRUE para PCA exploratoria

matriz <- assay(vsd)       # Se extraee la matriz transformada
datos_tras <- t(matriz)    # filas = muestras, columnas = genes


# Se escalan los datos y se hace la PCA
dat_scaled <- scale(datos_tras)
pca.results <- prcomp(dat_scaled, center = TRUE, scale. = TRUE)

# Resultado de las componentes principales
pca.df <- data.frame(pca.results$x)

# Varianza (cuadrado de la desviacion tipica)
varianzas <- pca.results$sdev^2

# Total de la varianza de los datos
total.varianza <- sum(varianzas)

# Varianza explicada por cada componente principal
varianza.explicada <- varianzas/total.varianza

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Se guarda el resultado de la PCA en un .CSV
write.csv(pca.df, "results/tables/PCA_results_vst.csv", row.names = TRUE)

# Representación gráfica de las primeras dos componentes principales respecto a los datos
p <- ggplot(pca.df, aes(x = PC1, y = PC2, color = muestras)) +
  geom_point(size=4) +
  scale_color_manual(values=c('#F4A6C6', '#E75480', '#C2185B', '#89CFF0', '#0056B3')) +
  labs(title = "PCA de expresión génica", x=x_label, y=y_label, color='Muestras') +
  theme_classic() +
  theme(
        legend.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(color="gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5)
      )
# Se guarda el grafico como png en la carpeta figures
ggsave(filename = "results/figures/PCA_plot.png",
       plot = p,
       width = 8,       # ancho en pulgadas
       height = 6,      # alto en pulgadas
       dpi = 300)       # resolución