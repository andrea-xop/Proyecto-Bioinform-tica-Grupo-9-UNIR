rm(list=ls())      # Se vacía el entorno de trabajo
setwd("C:\\Users") # Se introduce el directorio en el que se va a trabajar

# Se carga la librería del tidyverse
library(dplyr)


# Carga de los archivos a procesar
archivo1 <- read.csv("GSE306907_gene_count_matrix.csv", header = TRUE, stringsAsFactors = FALSE)
archivo2 <- read.csv("GSE306907_gene_count_2_matrix.csv", header = TRUE, stringsAsFactors = FALSE)

# Como tienen las mismas filas se unen por columnas
datos_unidos <- cbind(archivo1, archivo2)

# Algunas columnas se repiten por lo que se deben eliminar.

duplicadas <- duplicated(as.list(datos_unidos)) # Se convierte el dataframe en una lista de columnas y se evalua logicamente cuales estan duplicadas guardando sus valores en un vector 
datos_unicos <- datos_unidos[, !duplicadas]     # Se crea un nuevo dataframe en el que solo se guardan las columnas no duplicadas

# Las columnas quedan desordenadas por lo que se deben ordenar.

# Se extraen los nombres de las columnas
nombres <- names(datos_unicos)

# Se crean vectores de columnas según el orden en el que las queremos poner
gene_col    <- nombres[grepl("^gene_id$", nombres, ignore.case = TRUE)] # Empieza por gene_id y termina
F1_cols     <- nombres[grepl("^F1", nombres)]                           # Empieza por F1
F4_cols     <- nombres[grepl("^F4", nombres)]                           # Empieza por F4
F5_cols     <- nombres[grepl("^F5", nombres)]                           # Empieza por F5
M_cols      <- nombres[grepl("^M", nombres)]                            # Empieza por M
anno_cols   <- nombres[grepl("^anno_(chr|start|end)$", nombres)]        # Empieza por anno_ seguido de chr o start o end y termina

# Se combinan en el orden deseado
orden_final <- c(gene_col, F1_cols, F4_cols, F5_cols, M_cols, anno_cols)

# Se reordena el dataframe
datos_ordenados <- datos_unicos[, orden_final]

# Se guarda el resultado en un nuevo archivo .csv
write.csv(datos_ordenados, "GSE_unificado_ordenado.csv", row.names = FALSE)