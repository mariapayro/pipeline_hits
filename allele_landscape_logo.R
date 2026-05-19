#!/usr/bin/env Rscript
# =============================================================================
# Motor de Renderizado R (BacGWAS-Visualizer)
# Construye dinámicamente matrices de ggseqlogo a partir de hits de GWAS
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggseqlogo)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Faltan argumentos. Uso: Rscript allele_landscape_logo.R <matriz.csv> <rois.csv> <fasta> <outdir>", call.=FALSE)
}

archivo_matriz <- args[1]
archivo_rois   <- args[2]
archivo_ref    <- args[3]
outdir         <- args[4]

cat("\n[R-Script] Iniciando renderizado del Paisaje Genético y Logos...\n")

df_matriz <- read.csv(archivo_matriz, check.names = FALSE)
df_rois   <- read.csv(archivo_rois, check.names = FALSE)

# Asegurarnos de que las columnas están en mayúsculas para evitar errores
names(df_matriz) <- toupper(names(df_matriz))
names(df_rois) <- tolower(names(df_rois))

# ==========================================
# 📈 PANEL A: PAISAJE GENÓMICO GLOBAL
# ==========================================
cat("[R-Script] Dibujando Paisaje Diferencial...\n")

p_global <- ggplot(df_matriz, aes(x = START, y = DELTA_FREQ)) +
  # Sombreado amarillo dinámico para todas las ROIs detectadas
  geom_rect(data = df_rois, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "#FFF200", alpha = 0.5, inherit.aes = FALSE) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  # Usamos geom_segment (tipo Manhattan) porque son mutaciones puntuales/kmers
  geom_segment(aes(xend = START, yend = 0), color = "steelblue", linewidth = 0.6) +
  
  theme_classic() +
  labs(
    title = "A) Paisaje Diferencial Global de Variantes",
    x = "", 
    y = "Δ Frecuencia\nNAG (-)                      GC (+)"
  ) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

# ==========================================
# 🧬 PANEL B: CONSTRUCCIÓN DINÁMICA DE LOGOS
# ==========================================
cat(sprintf("[R-Script] Procesando %d Regiones de Interés (ROIs) para logos...\n", nrow(df_rois)))

paleta_adn <- make_col_scheme(
  chars  = c('A', 'C', 'G', 'T', '-'),
  groups = c('A', 'C', 'G', 'T', '-'),
  cols   = c('#109648', '#255C99', '#F7B32B', '#D62839', 'grey50')
)

lista_logos <- list()

for (i in 1:nrow(df_rois)) {
  # Extraemos coordenadas y le damos un "colchón" de 2 bases para que se vea más bonito
  roi_start <- df_rois$start[i] - 2
  roi_end   <- df_rois$end[i] + 2
  longitud  <- roi_end - roi_start + 1
  
  # 1. Crear matriz vacía para ggseqlogo
  mat_zoom <- matrix(0, nrow = 5, ncol = longitud)
  rownames(mat_zoom) <- c('A', 'C', 'G', 'T', '-')
  colnames(mat_zoom) <- as.character(roi_start:roi_end)
  
  # 2. Filtrar mutaciones que caen en esta ROI
  df_roi <- df_matriz %>% filter(START >= roi_start & START <= roi_end)
  
  # 3. Rellenar la matriz interpretando la mutación (Ej. A>T)
  if(nrow(df_roi) > 0 && "SNP" %in% names(df_roi)) {
    for(j in 1:nrow(df_roi)) {
      mutacion <- as.character(df_roi$SNP[j])
      posicion <- as.numeric(df_roi$START[j])
      delta    <- as.numeric(df_roi$DELTA_FREQ[j])
      
      # Calcular en qué columna de la matriz cae esta posición
      col_idx <- posicion - roi_start + 1
      
      if(grepl(">", mutacion) && col_idx > 0 && col_idx <= longitud) {
        partes <- strsplit(mutacion, ">")[[1]]
        ref <- partes[1]
        alt <- partes[2]
        
        # El alelo alterno (mutación) toma el valor del delta
        if(alt %in% rownames(mat_zoom)) mat_zoom[alt, col_idx] <- delta
        # El alelo de referencia baja proporcionalmente
        if(ref %in% rownames(mat_zoom)) mat_zoom[ref, col_idx] <- -delta
      }
    }
  }
  
  # 4. Dibujar el logo usando tu configuración original
  # Ajustamos los saltos del eje X dinámicamente según el tamaño de la ROI
  salto_x <- max(1, floor(longitud / 4)) 
  
  p_zoom <- ggseqlogo(mat_zoom, method = 'custom', seq_type = "dna", col_scheme = paleta_adn) +
    theme_logo() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 10),
      axis.title.x = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = sprintf("ROI %d (%d - %d)", i, roi_start, roi_end),
      x = "Posición",
      y = if(i==1) "Δ Frecuencia" else "" # Solo poner la etiqueta Y en el primer logo
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    scale_x_continuous(breaks = seq(1, ncol(mat_zoom), by = salto_x), 
                       labels = seq(roi_start, roi_end, by = salto_x))
  
  lista_logos[[i]] <- p_zoom
}

# ==========================================
# 🛠️ ENSAMBLAJE VECTORIAL Y GUARDADO
# ==========================================
cat("[R-Script] Ensamblando panel final...\n")

# Si hay múltiples logos, los pone uno al lado del otro
panel_logos <- wrap_plots(lista_logos, nrow = 1)
figura_final <- p_global / panel_logos + plot_layout(heights = c(1, 1.2))

ruta_salida_png <- file.path(outdir, "Figura_Paisaje_Genetico_Diferencial.png")
ruta_salida_pdf <- file.path(outdir, "Figura_Paisaje_Genetico_Diferencial.pdf")

ggsave(ruta_salida_png, plot = figura_final, width = 16, height = 10, dpi = 300)
ggsave(ruta_salida_pdf, plot = figura_final, width = 16, height = 10, dpi = 300, device = pdf)

cat(sprintf("[R-Script] ¡Éxito total! Figuras guardadas en: %s\n", outdir))