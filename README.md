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

# Colores para mensajes
GREEN="\e[32m"
RESET="\e[0m"

echo -e "${GREEN}==> Iniciando alineamiento y post-procesamiento...${RESET}"

# Crear carpeta de resultados
mkdir -p resultados

# Paso 1: Alineamiento
echo -e "${GREEN}==> Ejecutando BWA MEM...${RESET}"
bwa mem ecoli_ref.fa ecoli_reads.fastq.gz > resultados/aligned_reads.sam

# Paso 2: Conversión SAM a BAM
echo -e "${GREEN}==> Convirtiendo SAM a BAM...${RESET}"
samtools view -S -b resultados/aligned_reads.sam > resultados/aligned_reads.bam

# Paso 3: Ordenar BAM
echo -e "${GREEN}==> Ordenando archivo BAM...${RESET}"
samtools sort resultados/aligned_reads.bam -o resultados/aligned_reads_sorted.bam

# Paso 4: Indexar BAM ordenado
echo -e "${GREEN}==> Indexando archivo BAM ordenado...${RESET}"
samtools index resultados/aligned_reads_sorted.bam

# Paso 5: Eliminar duplicados
echo -e "${GREEN}==> Eliminando duplicados con samtools markdup...${RESET}"
samtools markdup -r resultados/aligned_reads_sorted.bam resultados/aligned_reads_nodup.bam

# Paso 6: Generar resumen del alineamiento
echo -e "${GREEN}==> Generando resumen con samtools flagstat...${RESET}"
samtools flagstat resultados/aligned_reads_nodup.bam > resumen_alineamiento.txt

echo -e "${GREEN}==> Proceso completado. Resultados guardados en la carpeta 'resultados'.${RESET}"
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
