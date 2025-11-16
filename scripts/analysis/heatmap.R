rm(list=ls())      # Se vacía el entorno de trabajo
library(pheatmap) # Librería para el heatmap

source("scripts/data_processing/DESeq2.R") # Llama y ejecuta el archivo DESeq2

# Crear vector de grupos
muestras <- c(rep("F_control", 3), rep("F_40μg/L_DBAN", 3), rep("F_200μg/L_DBAN", 3),
              rep("M_control", 3), rep("M_40μg/L_DBAN", 3))

# Preparar DESeqDataSet desde CSV preprocesado
dds <- preparar_dds("data/processed_data/GSE_unificado_ordenado.csv")

# Para heatmap o PCA
vsd <- obtener_vst(dds)
datos_scaled <- assay(vsd)

# Escalamos por gen, como scale() opera sobre las columnas transponemos para tener genes en columnas escalamos y transponemos.
datos_scaled <- t(scale(t(datos_scaled)))

# Reemplazamos posibles NA o Inf
datos_scaled[is.na(datos_scaled)] <- 0
datos_scaled[is.infinite(datos_scaled)] <- 0

# Se crea la anotación para las columnas
annotation_col <- data.frame(Grupo = muestras)
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
