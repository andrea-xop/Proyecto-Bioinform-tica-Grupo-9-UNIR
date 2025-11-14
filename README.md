# Proyecto-Bioinformática-Grupo-9-UNIR
Implementación de un repositorio base en el marco del proyecto de bioinformática del Grupo 9 (UNIR).

## Miembros del grupo 9
- **Daniel Carrera Cabezuelo** ()
- **David Carrera Natera** ()
- **Eloy Gomez Ayerbe** ()
- **Andrea Moríñigo Velaz** (andre-mori)
- **Cristina Ochoa Varela** ()
- **Andrea Unzu Redín** (andrea_xop)


**Instrucciones de uso**
Este proyecto analiza las diferencias en la expresión génica de *Danio rerio* (zebrafish) adultos expuestos a dibromoacetonitrilo 

**Requisitos del sistema**
Para reproducir el análisis se recomienda un entorno Linux o WSL con las siguientes herramientas:
   - **R ≥ 4.2.0** (para el análisis)
   - **SRA Toolkit** (para descargar archivos FASTQ)
   - **FASTQC** (para el análisis de control de calidad)
   - **fastp** (para la limpieza de secuencias)
   - **Salmon ≥ 1.10.0** (para el análisis de la expresión génica)
   - **Paquetes de R: "dplyr", "readr", "ggplot2", "pheatmap", "DESeq2", "tixmport"** (para la comparación y visualización de los resultados)

  
**Descarga de datos desde GEO**
   1. Accede a la página Gene Expression Omnibus (GEO): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306907
   2. En el apartado "Selector de ejecución SRA" identificar qué runs SRR corresponden a cada condición
   3. Mediante SRA Toolkit descarga los archivos FASTQ:
      ```bash
      #Crear carpeta
      mkdir -p data/faw

      #Descargar los datos del proyecto
      prefetch SRR35217944 SRR35217945 SRR35217946 SRR35217947 SRR35217948 SRR35217949 SRR35217950 SRR35217951 SRR35217952 SRR35217953 SRR35217954 SRR35217955 SRR35217956 SRR35217957 SRR35217958
      fasterq-dump SRRXXXXXXXX* -0 data/raw

      Sustituye SRRXXXXXXXX* por el nombre del acceso real

**Preprocesamiento**
Consiste en:
   - Análisis de calidad de las muestras: interfaz gráfica FASTQC 
   - Limpieza de datos: fastp

**Cuantificación de la expresión génica**
   1. Crear índice de transcriptoma de zebrafish
   2. Cuantificación de la expresión con Salmon

**Análisis diferencial**
   1. Preparar tabla de conteos de genes en R: "readr", "DESeq2" y "tixmport"
      - Crear tx2gene (tabla que relaciona transcritos con genes)
      - Importar cuantificación
   2. Mirar los genes más expresados: "dplyr"
   3. Filtrado de genes de baja expresión: "DESeq2"
   4. Transformación logarítmica: "DESeq2"

**Interpretación y visualización**


**Resultados esperados**
   - Identificación de genes y rutas neuronales alteradas por la exposición a dibromoacetonitrilo 
   - Evidencia de disrrupción de vías sinápticas en el grupo expuesto
   - Gráficos y tablas para informe o publicación 

**Referencias**
   1. GEO Accession:
   2. Publicación asociada:
   3. Organismo:
