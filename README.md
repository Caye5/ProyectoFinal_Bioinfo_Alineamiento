# 🔬 Proyecto Final - Bioinformática

**Autores:** Cayetana Reátegui & Daniela Guillén  
**Curso:** Principios de Programación  
**Tema:** Automatización del alineamiento de secuencias y post-procesamiento básico

---

## 🧠 Problema Biológico

El alineamiento de secuencias es un paso fundamental en los análisis genómicos realizados con tecnologías de secuenciación de alto rendimiento (NGS). Este proceso, cuando se realiza manualmente, puede resultar muy lento y propenso a errores, especialmente al trabajar con múltiples muestras.

Este proyecto busca resolver ese problema mediante un **script automatizado en Bash** que ejecute el alineamiento de lecturas contra un genoma de referencia y realice el post-procesamiento básico del resultado.

---

## 🎯 Objetivos

- Automatizar el alineamiento de lecturas FASTQ con un genoma de referencia.  
- Realizar el post-procesamiento: conversión, ordenamiento y eliminación de duplicados.  
- Generar un resumen con estadísticas clave del alineamiento final.

---

## 🛠️ Herramientas Utilizadas

- **BWA**: Herramienta de alineamiento para lecturas cortas.  
- **Samtools**: Suite para manipulación de archivos SAM/BAM.  
- **SRA Toolkit** *(opcional)*: Descarga de lecturas desde NCBI.

---

## 🗂️ Estructura del Proyecto

```text
ProyectoFinal_CayeDani/
├── script_proyecto.sh              # Script principal en Bash
├── ecoli_ref.fa                    # Genoma de referencia (FASTA)
├── ecoli_reads.fastq.gz           # Lecturas (FASTQ comprimido)
├── resultados/                     # Carpeta con archivos generados
└── resumen_alineamiento.txt       # Estadísticas del alineamiento
```

---

## ⚙️ Pasos Ejecutados

### 1. Preparar entorno de trabajo

```bash
mkdir ProyectoFinal_CayeDani
cd ProyectoFinal_CayeDani
```

---

### 2. Descargar el genoma de referencia (Escherichia coli)

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O ecoli_ref.fa.gz
gunzip ecoli_ref.fa.gz
```

---

### 3. Descargar lecturas FASTQ desde NCBI

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2584861/SRR2584861_1.fastq.gz -O ecoli_reads.fastq.gz
```

---

### 4. Indexar el genoma de referencia

```bash
bwa index ecoli_ref.fa
```

---

### 5. Ejecutar el script de alineamiento

```bash
bash script_proyecto.sh
```

---

## 📜 Descripción del Script (`script_proyecto.sh`)

Este script automatiza todo el flujo de alineamiento y post-procesamiento:

1. Alinea las lecturas con `bwa mem`  
2. Convierte `.sam` a `.bam` con `samtools view`  
3. Ordena el `.bam` con `samtools sort`  
4. Indexa el archivo ordenado con `samtools index`  
5. Elimina duplicados con `samtools markdup`  
6. Genera un resumen con `samtools flagstat`

---

### 🔧 Contenido del script

```bash
#!/bin/bash

# Script de alineamiento y post-procesamiento
# Autoras: Cayetana Reátegui y Daniela Guillén

# ================================================
# BIENVENIDA PERSONALIZADA
# ================================================
echo "Bienvenido/a al proyecto de alineamiento automatizado."
read -p "Por favor, ingresa tu nombre: " nombre
echo "Hola, $nombre. Iniciando el análisis..."

# ================================================
# CREACIÓN DE DIRECTORIOS
# ================================================
mkdir -p ProyectoFinal_CayeDani
cd ProyectoFinal_CayeDani
mkdir resultados
mkdir datos

# ================================================
# DESCARGA DEL GENOMA DE REFERENCIA
# ================================================
echo "Descargando genoma de referencia (Escherichia coli)..."
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O datos/ecoli_ref.fa.gz
gunzip datos/ecoli_ref.fa.gz

# ================================================
# DESCARGA DE LECTURAS
# ================================================
echo "Descargando lecturas FASTQ (SRR2584861)..."
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2584861/SRR2584861_1.fastq.gz -O datos/ecoli_reads.fastq.gz

# ================================================
# INDEXACIÓN DEL GENOMA
# ================================================
echo "Indexando el genoma..."
bwa index datos/ecoli_ref.fa

# ================================================
# ALINEAMIENTO
# ================================================
echo "Alineando lecturas contra el genoma..."
bwa mem datos/ecoli_ref.fa datos/ecoli_reads.fastq.gz > resultados/aligned_reads.sam

# ================================================
# POST-PROCESAMIENTO
# ================================================
echo "Convirtiendo SAM a BAM..."
samtools view -S -b resultados/aligned_reads.sam > resultados/aligned_reads.bam

echo "Ordenando archivo BAM..."
samtools sort resultados/aligned_reads.bam -o resultados/aligned_reads_sorted.bam

echo "Indexando archivo BAM ordenado..."
samtools index resultados/aligned_reads_sorted.bam

echo "Eliminando duplicados..."
samtools markdup -r resultados/aligned_reads_sorted.bam resultados/aligned_reads_nodup.bam

# ================================================
# RESUMEN DEL ALINEAMIENTO
# ================================================
echo "Generando resumen del alineamiento..."
samtools flagstat resultados/aligned_reads_nodup.bam > resumen_alineamiento.txt

# ================================================
# MENSAJE FINAL
# ================================================
echo "Proceso completado, $nombre. Los resultados se encuentran en la carpeta 'ProyectoFinal_CayeDani/resultados'."

```

---

## 📊 Salida Esperada

- `aligned_reads.sam`: alineamiento original  
- `aligned_reads_sorted.bam`: archivo ordenado  
- `aligned_reads_nodup.bam`: archivo sin duplicados  
- `resumen_alineamiento.txt`: estadísticas de mapeo  

---

## ✅ Observaciones Finales

- El flujo puede adaptarse fácilmente a múltiples muestras.  
- Mejora la reproducibilidad del análisis bioinformático.  
- Se puede extender para incluir análisis de variantes, expresión génica u otros estudios.

---
