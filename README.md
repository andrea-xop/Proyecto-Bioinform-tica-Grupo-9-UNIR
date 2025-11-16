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
 ## *Instrucciones de uso*
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

### Descarga de datos desde GEO
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

### Preprocesamiento
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

### Cuantificación de la expresión génica
   1. Crear índice de transcriptoma de zebrafish
```
      salmon index -t transcripts.fa -i index_salmon --type quasi -k 31
```

   3. Cuantificación de la expresión: Salmon

```
      mkdir -p results/salmon
for sample in $(ls data/processed/fastq_clean/*_1.clean.fastq | sed 's/_1.clean.fastq//' ); do
  base=$(basename $sample)
  salmon quant -i index_salmon -l A \
    -1 ${sample}_1.clean.fastq -2 ${sample}_2.clean.fastq \
    -p 8 --validateMappings -o results/salmon/${base}_quant
done


      4. Importar cuantificaciones y preparar matriz de conteo (R)
      
```
      #src/03_import_tximport.R (fragmento)
      library(tximport)
      library(readr)

     #tabla_samples.csv debe tener columnas: sample, condition, path_to_quant
     samples <- read_csv("data/raw/samples_table.csv")

     files <- file.path(samples$path_to_quant, "quant.sf")
     names(files) <- samples$sample
 
     tx2gene <- read_csv("data/raw/tx2gene.csv") # crear a partir de GTF o descargado
     txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
     counts <- txi$counts
     write.csv(counts, "data/processed/counts_matrix.csv")
```

### Análisis diferencial
   1. Preparar tabla de conteos de genes en R: "readr", "DESeq2" y "tixmport"
      - Crear tx2gene (tabla que relaciona transcritos con genes)
      - Importar cuantificación
   2. Mirar los genes más expresados: "dplyr"
   3. Filtrado de genes de baja expresión: "DESeq2"
   4. Transformación logarítmica: "DESeq2"

### Interpretación y visualización 
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
 ## *Flujo de trabajo*
Esta sección describe el **proceso lógico del análisis**, paso a paso, desde la descarga de datos hasta la obtención de resultados finales..
Piensa en esto como una **explicación conceptual del pipeline**.

 ### 1º Descargar los datos
El script `01_download_data.R` usa el paquete **GEOquery** para descargar la matriz de expresión y los metadatos del estudio GSE306907.
 ```r
   library(GEOquery)
   gse <- getGEO("GSE306907", GSEMatrix = TRUE)
   expr <- exprs(gse[[1]])
   meta <- pData(gse[[1]])
   write.csv(expr, "data/raw/expression_matrix.csv")
   write.csv(meta, "data/raw/metadata.csv")
```
### 2º Preprocesamiento
- Limpieza de datos.
- Filtrado de genes de baja expresión.
- Transformaciones necesarias para trabajar con RNA-Seq.

### 3º Análisis de expresión diferencial
Con DESeq2 se modelan los datos y se identifican genes sobreexpresados o infraexpresados entre condiciones experimentales.
```r
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta,
                              design    = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```
### 4º Visualización
Generación de gráficos clave:
- PCA: variabilidad global entre muestras.
- Volcano plot: genes significativos (p-adj y log2FC).
- Heatmap: principales genes diferencialmente expresados.

### 5º Informe reproducible
Todo el análisis se integra en un archivo RMarkdown para documentar métodos, resultados y visualizaciones de forma clara y formal.
El informe `docs/report.Rmd` integrará estos resultados para su presentación final.

 ---
## *Referencias*
-	GEO Accession: GSE306907 – NCBI GEO
-	Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.Genome Biology (2014).
-	R Core Team (2024). R: A language and environment for statistical computing.

 ---
## *Licencia*
Este proyecto se distribuye bajo la licencia **MIT**, permitiendo su uso, modificación y distribución libre con atribución al autor.
