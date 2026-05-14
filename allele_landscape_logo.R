#!/usr/bin/env Rscript

# ==========================================
# BacGWAS-Visualizer (Sprint 3: Módulo de Renderizado en R)
# This script is automatically called by bacwas_visualizer.py.
# It generates the Genetic Landscape and Differential logos graphics.
# ==========================================

# 1. Cargar librerías silenciosamente
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Biostrings)
  # library(ggseqlogo) # Descomentar si usas ggseqlogo para tus matrices
})

# 2. Capturar argumentos desde la terminal (Enviados por Python)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Faltan argumentos. Uso: Rscript generar_paisaje_logos.R <matriz_integrada.csv> <rois.csv> <fasta_referencia> <outdir>", call.=FALSE)
}

archivo_matriz <- args[1]
archivo_rois   <- args[2]
archivo_ref    <- args[3]
outdir         <- args[4]

cat("\n[R-Script] Iniciando renderizado del Paisaje Genético y Logos...\n")

# 3. Cargar Datos
df_matriz <- read.csv(archivo_matriz)
df_rois   <- read.csv(archivo_rois)
referencia <- readDNAStringSet(archivo_ref)
secuencia_ref <- as.character(referencia[[1]])

# Nota de diseño: Para que el paisaje (línea azul) se dibuje correctamente a lo largo de 
# todo el gen, df_matriz necesita tener una columna de 'coordenada_x'. 
# (Asumiremos que tu script de Python mapea todos los unitigs, o que tienes una lógica 
# para asignar el eje X en el paisaje).

# ==========================================
# 📊 GRÁFICO A: EL PAISAJE COMPLETO
# ==========================================
cat("[R-Script] Dibujando Paisaje Diferencial...\n")

# Recreando la estética de tu imagen original
plot_paisaje <- ggplot() +
  # Aquí iría tu geom_line o geom_segment con los datos de df_matriz
  # geom_line(data = df_matriz, aes(x = posicion, y = delta_freq), color = "steelblue") +
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  
  # Dibujar dinámicamente las bandas amarillas basadas en las ROIs
  geom_rect(data = df_rois, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "gold", alpha = 0.4) +
  
  theme_classic() +
  labs(title = "A) Paisaje Diferencial", x = "", y = "GC (+) . Frecuencia . NAG (-)") +
  theme(plot.title = element_text(face = "bold", size = 14))

# ==========================================
# 🧬 GRÁFICO B: GENERACIÓN DINÁMICA DE LOGOS
# ==========================================
cat(sprintf("[R-Script] Procesando %d Regiones de Interés (ROIs) para logos...\n", nrow(df_rois)))

lista_logos <- list()

for (i in 1:nrow(df_rois)) {
  roi_start <- df_rois$start[i]
  roi_end   <- df_rois$end[i]
  
  # Extraer la secuencia específica de la región desde el FASTA
  secuencia_roi <- substr(secuencia_ref, roi_start, roi_end)
  
  # -----------------------------------------------------
  # AQUÍ INSERTAS TU CÓDIGO ORIGINAL DEL LOGO DIFFERENCIAL
  # -----------------------------------------------------
  # Ejemplo conceptual de cómo envolverlo en un plot de ggplot:
  
  plot_logo <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = paste("Logo ROI", i, "\n", roi_start, "-", roi_end), size = 6) +
    theme_void() +
    labs(title = sprintf("Detalle ROI %d (%d - %d)", i, roi_start, roi_end)) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  
  lista_logos[[i]] <- plot_logo
}

# ==========================================
# 🧩 ENSAMBLAJE FINAL (PATCHWORK) Y GUARDADO
# ==========================================
cat("[R-Script] Ensamblando panel final...\n")

# Patchwork es magia: 'plot_paisaje / plot_logos' pone uno arriba y los otros abajo
if (length(lista_logos) == 1) {
  panel_final <- plot_paisaje / lista_logos[[1]] + plot_layout(heights = c(1, 1))
} else {
  # Si hay múltiples logos, los pone en fila en la parte inferior usando 'wrap_plots'
  panel_logos <- wrap_plots(lista_logos, nrow = 1)
  panel_final <- plot_paisaje / panel_logos + plot_layout(heights = c(1, 1))
}

ruta_salida <- file.path(outdir, "Paisaje_y_Logos_Diferenciales.png")
ggsave(ruta_salida, plot = panel_final, width = 12, height = 8, dpi = 300)

cat(sprintf("[R-Script] ¡Éxito! Gráfico guardado en: %s\n\n", ruta_salida))