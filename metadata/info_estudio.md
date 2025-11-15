# Información del estudio GSE306907

Este archivo contiene información complementaria del estudio utilizado en este proyecto de análisis de expresión génica diferencial.

## 1. Organización del experimento
- Organismo: *Danio rerio* (pez cebra)
- Sustancia: Dibromoacetonitrilo (DBAN), subproducto de la desinfección del agua
- Duración de la exposición: 8 semanas
- Tipo de muestras: peces adultos (machos y hembras)
- Tecnología: RNA-seq

## 2. Condiciones experimentales
Se establecieron los siguientes grupos de exposición:

| Grupo | Concentración DBAN (μg/L) | Sexo | Observaciones |
|-------|-----------------------------|------|---------------|
| Control | 0 | male/female | Sin exposición |
| DBAN_1.6 | 1.6 | male/female | Exposición baja |
| DBAN_8 | 8 | male/female | Exposición moderada |
| DBAN_40 | 40 | male/female | Nivel donde se observan efectos transcriptómicos claros |
| DBAN_200 | 200 | male/female | Afecta comportamiento |
| DBAN_1000 | 1000 | male/female | Afecta supervivencia |

## 3. Archivos recomendados para incluir aquí
- Metadatos descargados desde GEO  
- Diccionario de variables de las tablas  
- Notas del estudio  
- Archivos .txt con descripción del diseño experimental  

## 4. Información útil para el análisis
- Factores de diseño: `sex`, `dose`
- Tipo de análisis: DESeq2
- Comparaciones principales: exposición vs control
