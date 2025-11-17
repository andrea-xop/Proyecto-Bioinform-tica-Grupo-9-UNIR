# DESeq2.R
library(DESeq2)

# Función para preparar DESeq2
preparar_dds <- function(file_csv, col_names_anno = "anno_", min_counts = 10, diseño = "~ Grupo") {
  # Leer datos
  datos <- read.csv(file_csv, header = TRUE, row.names = 1, check.names = FALSE)
  datos <- datos[, !startsWith(colnames(datos), col_names_anno)]
  
  # Crear vector de grupos
  muestras <- c(rep("F_control", 3), rep("F_40μg/L_DBAN", 3), rep("F_200μg/L_DBAN", 3),
                rep("M_control", 3), rep("M_40μg/L_DBAN", 3))
  
  colData <- data.frame(row.names = colnames(datos), Grupo = factor(muestras))
  
  # Crear DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = datos,
                                colData = colData,
                                design = as.formula(diseño))
  
  # Filtrar genes con pocos counts
  dds <- dds[rowSums(counts(dds)) > min_counts, ]
  
  return(dds)
}

# Función para obtener la transformación VST
obtener_vst <- function(dds, blind = TRUE) {
vsd <- vst(dds, blind = blind)
return(vsd)
}