rm(list=ls())      # Limpia el entorno de trabajo

# Se cargan las librerías necesarias
library(ggplot2)   # Librería para la creación de gráficos
library(ggrepel)   # Librería para añadir etiquetas al gráfico evitando solapamientos

source("scripts/data_processing/DESeq2.R") # Llama y ejecuta el archivo DESeq2

# Crear vector de grupos
muestras <- c(rep("F_control", 3), rep("F_40μg/L_DBAN", 3), rep("F_200μg/L_DBAN", 3),
              rep("M_control", 3), rep("M_40μg/L_DBAN", 3))

# Preparar DESeqDataSet desde CSV preprocesado
dds <- preparar_dds("data/processed_data/GSE_unificado_ordenado.csv")

# Ejecutamos el análisis de expresión diferencial con DESeq2
dds <- DESeq(dds)

# Se eligen las parejas de datos que se van a representar
comparaciones <- list(
  "F_40μg/L_DBAN_vs_F_control" = c("Grupo", "F_40μg/L_DBAN", "F_control"),
  "F_200μg/L_DBAN_vs_F_control" = c("Grupo", "F_200μg/L_DBAN", "F_control"),
  "M_40μg/L_DBAN_vs_M_control" = c("Grupo", "M_40μg/L_DBAN", "M_control")
)

# Creamos una función para realizar los tres volcano_plots
crear_volcano <- function(res, titulo, carpeta_salida = NULL,  ruta_csv = NULL) {
  # Se añade el nombre de gen como columna
  res$gene <- rownames(res)
  
  # Se calcula el -log10 del p-valor
  res$negLogP <- -log10(res$pvalue)
  
  # Umbrales de significancia
  logFC_threshold <- 1
  pvalue_threshold <- 0.05
  
  # Tamaño del punto proporcional al p-valor (limitado a 6)
  res$PointSize <- pmin(res$negLogP / 2, 6)
  
  # Clasificamos los genes según su significancia
  res$Significance <- "No significativo"
  res$Significance[res$log2FoldChange > logFC_threshold & res$pvalue < pvalue_threshold] <- "Positivo"
  res$Significance[res$log2FoldChange < -logFC_threshold & res$pvalue < pvalue_threshold] <- "Negativo"
  
  # Creamos el gráfico de tipo volcano
  p <- ggplot(res, aes(x = log2FoldChange, y = negLogP)) +
    geom_point(aes(color = Significance, size = PointSize), alpha = 0.8) +  # puntos coloreados según significancia
    scale_color_manual(values = c("Negativo" = "#377eb8",
                                  "No significativo" = "grey70",
                                  "Positivo" = "#e41a1c")) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold),
               linetype = "dashed", color = "white") +                      # líneas verticales de umbral
    geom_hline(yintercept = -log10(pvalue_threshold),
               linetype = "dashed", color = "white") +                      # línea horizontal de p-valor
    geom_text_repel(data = subset(res, Significance != "No significativo" & negLogP > 10),
                    aes(label = gene),
                    size = 3, box.padding = 0.5, point.padding = 0.5,
                    max.overlaps = 15, color = "white") +                   # etiquetas de genes significativos
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "#222222"),                   # fondo oscuro
      plot.background = element_rect(fill = "#222222"),
      panel.grid = element_line(color = "#444444"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(hjust = 0.5, face = "bold", color = "white"),
      legend.title = element_blank(),
      legend.text = element_text(color = "white")
    ) +
    labs(title = titulo, x = "log2 Fold Change", y = "-log10(p-value)") +
    guides(size = "none")  # oculta la leyenda del tamaño de punto
  
  # Se guardan los gráficos si se indica una carpeta de salida
  if (!is.null(carpeta_salida)) {
    if (!dir.exists(carpeta_salida)) dir.create(carpeta_salida, recursive = TRUE)
    
    # Limpia el nombre del archivo
    nombre_archivo <- titulo
    nombre_archivo <- gsub(" ", "_", nombre_archivo)       # reemplaza espacios por "_"
    nombre_archivo <- gsub("[/\\\\]", "-", nombre_archivo) # reemplaza "/" y "\" por "-"
    nombre_archivo <- gsub("μ", "u", nombre_archivo)       # reemplaza "μ" por "u"
    nombre_archivo <- paste0(nombre_archivo, ".png")
    
    # Guarda los gráficos como PNGs
    ggsave(filename = file.path(carpeta_salida, nombre_archivo),
           plot = p, width = 8, height = 6, dpi = 300, device = "png")
  }
  # Guardar CSV en otra ruta
  if (!is.null(ruta_csv)) {
    if (!dir.exists(ruta_csv)) dir.create(ruta_csv, recursive = TRUE)
    
    nombre_csv <- paste0("Resultados_", gsub("[/\\\\]", "_", titulo), ".csv")
    write.csv(res, file = file.path(ruta_csv, nombre_csv), row.names = TRUE)
  }
  return(p)
}


# Ejecutamos todas las comparaciones y guardamos los gráficos
for (nombre in names(comparaciones)) {
  contraste <- comparaciones[[nombre]]                        # obtiene el contraste
  res <- results(dds, contrast = contraste)                   # ejecuta DESeq2 para esa comparación
  print(crear_volcano(res,                                 # Llama a la función para crear el volcano plot
                      titulo = nombre, 
                      carpeta_salida = "results/figures", 
                      ruta_csv = "results/tables"))
}
