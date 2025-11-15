rm(list=ls())      # Se vacía el entorno de trabajo
library(pheatmap) # Librería para el heatmap

# Se abren los datos a analizar
datos <- read.csv("data/processed_data/GSE_unificado_ordenado.csv", header = TRUE, row.names = 1, check.names = FALSE)
datos <- datos[, !startsWith(colnames(datos), "anno_")] # Se eliminanan las columnas no necesarias

# Transformación logaritmica de los counts
# Se usa log2(datos + 1) para:
#  - Reducir la asimetría de los datos (genes con counts muy altos vs muy bajos)
#  - Evitar valores infinitos al aplicar log a cero (se suma 1)
datos_log <- log2(datos + 1)
# Filtrado de genes de varianza 0
datos_log <- datos_log[apply(datos_log, 1, var) != 0, ]

# Escalado de los datos por cada gen (filas)
# Se centra y normaliza cada gen restando la media y dividiendo por la desviación estándar
# La doble transposición (t()) se usa porque scale() opera por columnas.
datos_scaled <- t(scale(t(datos_log)))

# Verificar valores infinitos o NA
datos_scaled[is.na(datos_scaled)] <- 0
datos_scaled[is.infinite(datos_scaled)] <- 0

# Cambia los nombres a los tipos de tratamiento
muestras <- c(rep("F control", 3), rep("F 40 μg/L DBAN", 3), rep("F 200 μg/L DBAN", 3), rep("M control", 3), rep("M 40 μg/L DBAN", 3))
names(muestras) <- colnames(datos_scaled)

# Se crea la anotación para las columnas
annotation_col <- data.frame(Muestras = muestras)
rownames(annotation_col) <- colnames(datos_scaled)

# Se define la paleta de colores de las muestras
colores_grupo <- c(
  "F control" = "#F4A6C6",
  "F 40 μg/L DBAN"   = "#E75480",
  "F 200 μg/L DBAN"   = "#C2185B",
  "M control" = "#89CFF0",
  "M 40 μg/L DBAN"    = "#0056B3"
)

# Se guarda el heatmap en un PDF a un tamaño en el que se puedan leer los genes
pdf("results/figures/heatmap_genes.pdf", width = 10, height = 700)  # Se define la altura y anchura

# Se genera el heatmap
set.seed(2001)
pheatmap(
  datos_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_colnames = FALSE,
  fontsize_row = 5,  # Tamaño de los nombres de las filas
  annotation_col = annotation_col, # Se asigna las anotaciones de las columnas al heatmap
  annotation_colors = list(Muestras = colores_grupo),
  clustering_distance_rows = "euclidean", # Se calcula la distancia euclidiana entre los perfiles de expresión
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", # Se define el método de enlace para el cluster
  main = "Heatmap de expresión génica",
  treeheight_row = 10, # Tamaños de los dendrogramas
  treeheight_col = 0
)

dev.off() # Se cierra el pdf