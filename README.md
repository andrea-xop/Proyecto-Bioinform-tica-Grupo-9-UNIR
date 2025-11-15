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
 ## *Instrucciones de uso*
Esta sección explica **cómo ejecutar el proyecto paso a paso en tu propio ordenador**.
Incluye los comandos necesarios, el orden de ejecución de los scripts y las acciones prácticas que cualquier usuario debe realizar para reproducir el análisis completo.
Piensa en esto como un manual de ejecución.
 A rellenar
 
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
