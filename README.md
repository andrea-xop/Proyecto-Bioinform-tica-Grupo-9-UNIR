# Análisis de expresión génica diferencial 

**Asignatura:** Introducción a la Programación Científica  

**Universidad Internacional de La Rioja (UNIR)**  

**Integrantes del grupo:** 
- **Daniel Carrera Cabezuelo** (Dancar96)
- **David Carrera Natera** (b62canad)
- **Eloy Gomez Ayerbe** (Eloy27)
- **Andrea Moriñigo Velaz** (andre-mori)
- **Cristina Ochoa Varela** (Crisov31)
- **Andrea Unzu Redín** (andrea_xop)

 ---
## *Descripción general*
 El presente proyecto tiene como finalidad realizar un **análisis completo de expresión génica diferencial** a partir de datos públicos del estudio **[GSE306907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306907)** disponible en el repositorio Gene Expression Omnibus (GEO) de NCBI. 
 
 Este estudio emplea la tecnología de secuenciación masiva de RNA (RNA-seq) para caracterizar los cambios transcripcionales en muestras biológicas bajo diferentes condiciones experimentales. Analiza los efectos neurotóxicos del **dibromoacetonitrilo (DBAN)**, en el pez cebra (*Danio rerio*) mediante datos de RNA-seq. 
 
 Los peces fueron expuestos durante ocho semanas a distintas concentraciones de DBAN y se evaluaron cambios de comportamiento y expresión génica. El análisis transcriptómico mostró alteraciones en vías neuronales clave, con efectos más marcados en machos. 
 
 A partir de los datos de **RNA-seq**, se evaluarán las diferencias de expresión entre grupos experimentales mediante el paquete **DESeq2**, con el fin de indentificar genes diferencialmente expresados y generar visualizaciones que resuman los resultados biológicos. 
 
 ---
## *Objetivos*

El objetivo principal de este proyecto es desarrollar un **pipeline reproducible de análisis de expresión génica diferencial** a partir del conjunto de datos público **GSE306907**, empleando R y una estructura organizada de scripts, datos y resultados.

### Objetivos específicos

1. **Adquirir y organizar los datos del estudio**
   - Descargar la matriz de conteos y los metadatos utilizando scripts (por ejemplo, `scripts/data_processing/Descarga_rawData.R`).
   - Almacenar los datos originales en `data/raw_data/` y separarlos de los datos procesados en `data/processed_data/`.

2. **Realizar el preprocesamiento y limpieza de la matriz de expresión**
   - Filtrar genes de baja expresión, revisar la integridad de las muestras y efectuar cualquier transformación necesaria para garantizar la calidad del análisis.
   - Implementar este procesamiento en los scripts correspondientes del directorio `scripts/data_processing/`.

3. **Ejecutar el análisis de expresión diferencial**
   - Construir el diseño experimental a partir de los metadatos almacenados en `metadata/info_estudio.md`.
   - Crear el objeto DESeqDataSet y ejecutar el pipeline de DESeq2 para identificar genes diferencialmente expresados.
   - Exportar los resultados estructurados en tablas dentro de `results/tables/`.

4. **Generar visualizaciones interpretables**
   - Producir gráficos clave como PCA, heatmap y volcano plot mediante los scripts alojados en `scripts/analysis/` (`PCA.R`, `heatmap.R`, `Volcano_plot.R`).
   - Guardar las figuras generadas en `results/figures/`.

5. **Documentar el proyecto y garantizar su reproducibilidad**
   - Integrar todo el flujo de análisis, resultados visuales y conclusiones en el documento `docs/report.Rmd`.
   - Mantener una trazabilidad clara mediante un sistema de control de versiones, commits descriptivos y trabajo colaborativo mediante ramas en GitHub.

 ---
 ## *Instrucciones de ejecución del proyecto*
Esta sección explica **cómo ejecutar el proyecto paso a paso en tu propio ordenador**.
Incluye el orden de ejecución y las acciones prácticas que cualquier usuario debe realizar para reproducir el análisis completo.
Piensa en esto como un manual de ejecución. 

**Requisitos del sistema** (TODOS LOS PAQUETES QUE USO PARA COMENTAR Y QUITAR DESPUÉS DE AQUI!!!)
Para reproducir el análisis se recomienda un entorno Linux o WSL con las siguientes herramientas:
   - **R ≥ 4.2.0** (para el análisis)
   - **SRA Toolkit** (para descargar archivos FASTQ)
   - **FASTQC** (para el análisis de control de calidad)
   - **fastp** (para la limpieza de secuencias)
   - **Salmon ≥ 1.10.0** (para el análisis de la expresión génica)
   - **Paquetes de R: "tidyverse", "readr", "ggplot2", "pheatmap", "BioCManager", "DESeq2", "GEOquery", "tixmport"** (para la comparación y visualización de los resultados)

### 1º Descarga de datos desde GEO
   1. Accede a la página Gene Expression Omnibus (GEO): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306907
   2. En el apartado "Selector de ejecución SRA" identificar qué runs SRR corresponden a cada condición
   3. Mediante SRA Toolkit descarga los archivos FASTQ:
```
      #Crear carpeta para almacenar datos en bruto
      mkdir -p data/raw

      #Descargar los archivos .sra de todas las muestras
        prefetch SRR35217944 SRR35217945 SRR35217946 SRR35217947 SRR35217948 \
        SRR35217949 SRR35217950 SRR35217951 SRR35217952 SRR35217953 \
        SRR35217954 SRR35217955 SRR35217956 SRR35217957 SRR35217958

      #Convertir .sra -> FASTQ
      #Sustituye "SRR35217944 SRR35217945 ..." por la lista real si cambia
           for srr in SRR35217944 SRR35217945 SRR35217946 SRR35217947 SRR35217948 \
           SRR35217949 SRR35217950 SRR35217951 SRR35217952 SRR35217953 \
           SRR35217954 SRR35217955 SRR35217956 SRR35217957 SRR35217958; do
   
       fasterq-dump $srr -O data/raw --split-files
       done
```

### 2º Preprocesamiento
Consiste en:
   - Análisis de calidad de las muestras: interfaz gráfica FASTQC

```
     mkdir -p results/fastqc
     fastqc -o results/fastqc data/raw/fastq/*.fastq
```

   - Limpieza de datos: fastp

```
     mkdir -p data/processed/fastq_clean
for fq in data/raw/fastq/*_1.fastq; do
  base=$(basename $fq _1.fastq)
  fastp -i data/raw/fastq/${base}_1.fastq -I data/raw/fastq/${base}_2.fastq \
        -o data/processed/fastq_clean/${base}_1.clean.fastq \
        -O data/processed/fastq_clean/${base}_2.clean.fastq \
        -h results/fastp/${base}_fastp.html -j results/fastp/${base}_fastp.json
done
```

### 3º Cuantificación de la expresión génica
   1. Crear índice de transcriptoma de zebrafish
 ```
      salmon index -t transcripts.fa -i index_salmon --type quasi -k 31
```  
   2. Cuantificación de la expresión: Salmon
```
      mkdir -p results/salmon
for sample in $(ls data/processed/fastq_clean/*_1.clean.fastq | sed 's/_1.clean.fastq//' ); do
  base=$(basename $sample)
  salmon quant -i index_salmon -l A \
    -1 ${sample}_1.clean.fastq -2 ${sample}_2.clean.fastq \
    -p 8 --validateMappings -o results/salmon/${base}_quant
done
```
   3. Preparar tabla de conteos de genes en R:: "readr", "DESeq2" y "tixmport"
      - Crear tx2gene (tabla que relaciona transcritos con genes)
      - Importar cuantificación   

```
#Preparar tabla de conteos de genes desde cuantificación Salmon
#src/03_import_tximport.R (fragmento)
      library(tximport)
      library(readr)

#CARGAR TABLA DE MUESTRAS
#tabla_samples.csv debe tener columnas: sample, condition, path_to_quant
     samples <- read_csv("data/raw/samples_table.csv")

#Crear bvector con rutas a quant.sf
     files <- file.path(samples$path_to_quant, "quant.sf")
     names(files) <- samples$sample

#Comprobar que todos los archivos existen
     missinf <- files[!file.exists(files)]
     if(length(missing) > 0){
  stop("ERROR: Las siguientes rutas no existen:\n",
       paste(missing, collapse = "\n"))
}

#Cargar tx2gene (transcrito --> gen)
     tx2gene <- read_csv("data/raw/tx2gene.csv")
# crear a partir de GTF o descargado

#mportar cuantificación con tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
     counts <- txi$counts
     write.csv(counts, "data/processed/counts_matrix.csv")

#Guardar matrices de salida
#Matriz de conteos por gen
write_csv(counts, "data/processed/counts_matrix.csv")

#(Opcional) Matriz TPM — útil para exploración
write_csv(as.data.frame(txi$abundance), "data/processed/TPM_matrix.csv")

#(Opcional) Longitudes efectivas por gen
write_csv(as.data.frame(txi$length), "data/processed/gene_lengths.csv")

message("✔ Importación completada: counts_matrix.csv generado")
```

### 4º Análisis diferencial
   1. Mirar los genes más expresados: "dplyr"
   3. Filtrado de genes de baja expresión: "DESeq2"
   4. Transformación logarítmica: "DESeq2"

```
#REVISAR LOS GENES MÁS EXPRESADOS
message("Most expressed genes:")

top_genes <- counts %>%
  mutate(gene = rownames(counts)) %>%
  rowwise() %>%
  mutate(mean_expression = mean(c_across(where(is.numeric)))) %>%
  arrange(desc(mean_expression)) %>%
  slice(1:20)

print(top_genes %>% select(gene, mean_expression))

write_csv(top_genes, "results/tables/top_expressed_genes.csv")

#FILTRADO DE GENES DE BAJA EXPRESIÓN (DESeq2)
meta <- read_csv("data/raw/metadata_for_DE.csv")

#metadata_for_DE.csv debe tener:

#sample, condition

#Asegurar que el orden coincide:
all(colnames(counts) == meta$sample)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ condition)

#Filtrar genes con muy pocos conteos:
#Mantener genes con 10 o más conteos totales
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(paste("Genes tras filtrado:", nrow(dds)))

#TRANSFORMACIÓN LOGARÍTMICA (normalización)
#VST (más rápido)
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "data/processed/vst_normalized.rds")

#rlog (opcional, más lento)

#rld <- rlog(dds, blind = FALSE)

#saveRDS(rld, "data/processed/rlog_normalized.rds")
message("Transformación logarítmica completada. Archivos guardados.")
```

### 5º Interpretación y visualización 
   Interpretación:
   - Criterios: padj < 0.05 y |log2FC| > 1
   - PCA: comprobación de si las muestras se separan por condición; si no, revisar la calidad
   - Mapping rate baja: puede indicar mala extracción
   - Genes de interés: si no aparecen significativos, revisar expresión; puede que no se expresen en ese tejido.

   Visualización:
   - Tabla resumen para cada gen de interés: mean counts por condición + log2FC + padj
   - Boxplot: para observar la distribución por muestras
   - Heatmap: muestra varios genes de interes (rows) vs muestras (cols)

### Resultados esperados
   - Identificación de genes y rutas neuronales alteradas por la exposición a dibromoacetonitrilo 
   - Evidencia de disrrupción de vías sinápticas en el grupo expuesto
   - Gráficos y tablas para informe o publicación 
 
 ---
## *Requisitos técnicos y del entorno*
Para la recreación de este proyecto será necesario tener en cuenta la versión de R utilizada, así como los paquetes que necesitarán ser instalados, estando enlistado a continuación:

### 1º Versión de R necesaria:
Debe ser una versión de R igual o superior a 4.3

(En este caso se ha utilizado el entorno gráfico proporcionado por RStudio, pero no es necesario).

### 2º Paquetes de R necesarios:
```r
“tidyverse”
“readr”
“janitor”
“ggplot2”
“pheatmap”
“BiocManager”
”DESeq2”
“GEOquery”
```

### 3º Elaboración del entorno:
```r
install.packages(c("tidyverse", "readr", "janitor", "ggplot2", "pheatmap"))

if (!requireNamespace("BiocManager", quietly=TRUE)) 
install.packages("BiocManager") BiocManager::install(c("DESeq2", "GEOquery"))
```

 ---
## *Visión conceptual del flujo de análisis*
Mientras que la sección anterior (*Instrucciones de ejecución del proyecto*) detalla los pasos técnicos necesarios para ejecutar cada herramienta del proyecto, esta sección ofrece una **visión conceptual del proceso completo**, explicando **qué ocurre en cada fase del análisis y con qué finalidad**, sin entrar en comandos ni parámetros específicos. Su propósito es ayudar a entender el sentido global del pipeline.

 ### 1º Obtención de los datos
El análisis comienza con la recuperación del dataset correspondiente al estudio GSE306907, incluyendo tanto la información de expresión como los metadatos asociados a cada muestra. Estos metadatos permiten definir las condiciones experimentales y son esenciales para cualquier análisis comparativo.

### 2º Control y preparación de la información
Los datos brutos se someten a una revisión de calidad para asegurar que cumplen estándares mínimos. A continuación, se organizan y procesan con el fin de eliminar ruido técnico y estructurar la información en un formato adecuado para los análisis posteriores.

### 3º Cuantificación de la expresión génica
Una vez depurados, los datos se emplean para estimar la abundancia relativa de cada transcrito o gen en las muestras del estudio. Esta cuantificación da lugar a matrices de expresión que servirán como base para identificar diferencias entre condiciones.

### 4º Análisis de expresión diferencial
A través de modelos estadísticos se comparan las distintas condiciones experimentales. Este proceso permite detectar genes que presentan una expresión significativamente mayor o menor en respuesta al tratamiento o estímulo, revelando posibles mecanismos biológicos implicados.

### 5º Visualización e interpretación de resultados
Los resultados se representan mediante técnicas gráficas —como PCA, volcano plots o heatmaps— que facilitan la interpretación del comportamiento global de las muestras y la identificación de genes relevantes. Estas visualizaciones ayudan a detectar patrones y validar la coherencia del análisis.

### 6. Informe reproducible
Finalmente, todo el análisis se integra en un documento reproducible (por ejemplo, un archivo RMarkdown) que reúne los métodos, resultados y visualizaciones. Este informe permite revisar de forma ordenada y transparente todo el proceso analítico.

 ---
## *Referencias*
-	GEO Accession: GSE306907 – NCBI GEO
-	Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.Genome Biology (2014).
-	R Core Team (2024). R: A language and environment for statistical computing.

 ---
## *Licencia*
Este proyecto se distribuye bajo la licencia **MIT**, permitiendo su uso, modificación y distribución libre con atribución al autor.
