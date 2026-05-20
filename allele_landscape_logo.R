#!/usr/bin/env Rscript
# =============================================================================
# Motor de Renderizado R (BacGWAS-Visualizer)
# Modo Dual: Dibuja SNPs puntuales o Paisajes de K-mers (Unitigs) superpuestos
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggseqlogo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Faltan argumentos. Uso: Rscript allele_landscape_logo.R <matriz.csv> <rois.csv> <fasta> <outdir> <tipo>", call.=FALSE)
}

archivo_matriz <- args[1]
archivo_rois   <- args[2]
archivo_ref    <- args[3] # (Se podría usar en el futuro si se necesita extraer más contexto)
outdir         <- args[4]
tipo_variante  <- args[5] # "snp" o "unitig"

cat(sprintf("\n[R-Script] Iniciando Motor Gráfico en Modo: %s...\n", toupper(tipo_variante)))

df_matriz <- read.csv(archivo_matriz, check.names = FALSE)
df_rois   <- read.csv(archivo_rois, check.names = FALSE)

names(df_matriz) <- toupper(names(df_matriz))
names(df_rois) <- tolower(names(df_rois))

# ==========================================
# 📈 PANEL A: PAISAJE GENÓMICO GLOBAL
# ==========================================
p_global <- ggplot(df_matriz, aes(x = START, y = DELTA_FREQ)) +
  geom_rect(data = df_rois, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "#FFF200", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_segment(aes(xend = START, yend = 0), color = "steelblue", linewidth = 0.6) +
  theme_classic() +
  labs(
    title = sprintf("A) Paisaje Diferencial Global (%s)", toupper(tipo_variante)),
    x = "", 
    y = "Δ Frecuencia\nNAG (-)                      GC (+)"
  ) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

# ==========================================
# 🧬 PANEL B: CONSTRUCCIÓN DINÁMICA DE LOGOS
# ==========================================
paleta_adn <- make_col_scheme(
  chars  = c('A', 'C', 'G', 'T', '-'),
  groups = c('A', 'C', 'G', 'T', '-'),
  cols   = c('#109648', '#255C99', '#F7B32B', '#D62839', 'grey50')
)

lista_logos <- list()

for (i in 1:nrow(df_rois)) {
  # Colchón dinámico dependiendo si es SNP o Unitig
  colchon <- ifelse(tipo_variante == "snp", 2, 0) 
  roi_start <- df_rois$start[i] - colchon
  roi_end   <- df_rois$end[i] + colchon
  longitud  <- roi_end - roi_start + 1
  
  mat_zoom <- matrix(0, nrow = 5, ncol = longitud)
  rownames(mat_zoom) <- c('A', 'C', 'G', 'T', '-')
  colnames(mat_zoom) <- as.character(roi_start:roi_end)
  
  df_roi <- df_matriz %>% filter(START >= roi_start & START <= roi_end)
  
  if(nrow(df_roi) > 0) {
    for(j in 1:nrow(df_roi)) {
      delta <- as.numeric(df_roi$DELTA_FREQ[j])
      posicion_inicio <- as.numeric(df_roi$START[j])
      
      # ------------------------------------------------
      # LÓGICA MODO SNP
      # ------------------------------------------------
      if (tipo_variante == "snp" && "SNP" %in% names(df_roi)) {
        mutacion <- as.character(df_roi$SNP[j])
        col_idx <- posicion_inicio - roi_start + 1
        if(grepl(">", mutacion) && col_idx > 0 && col_idx <= longitud) {
          partes <- strsplit(mutacion, ">")[[1]]
          alt <- partes[2]
          if(alt %in% rownames(mat_zoom)) mat_zoom[alt, col_idx] <- delta
        }
      } 
      # ------------------------------------------------
      # LÓGICA MODO UNITIG (K-MERS LARGOS)
      # ------------------------------------------------
      else if (tipo_variante == "unitig" && "SECUENCIA_REAL" %in% names(df_roi)) {
        secuencia <- as.character(df_roi$SECUENCIA_REAL[j])
        bases <- strsplit(secuencia, "")[[1]]
        
        for (k in seq_along(bases)) {
          base <- bases[k]
          # La posición real de esta letra en el genoma
          pos_actual <- posicion_inicio + k - 1 
          col_idx <- pos_actual - roi_start + 1
          
          if (col_idx > 0 && col_idx <= longitud && base %in% rownames(mat_zoom)) {
            # Nos quedamos con la señal más fuerte si varios unitigs se enciman
            if (abs(delta) > abs(mat_zoom[base, col_idx])) {
              mat_zoom[base, col_idx] <- delta
            }
          }
        }
      }
    }
  }
  
  salto_x <- max(1, floor(longitud / 5)) 
  
  p_zoom <- ggseqlogo(mat_zoom, method = 'custom', seq_type = "dna", col_scheme = paleta_adn) +
    theme_logo() +
    coord_cartesian(ylim = c(-0.6, 0.6)) + # Evita que se estiren las letras feo
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = sprintf("ROI %d (%d - %d)", i, roi_start, roi_end),
      x = "Posición",
      y = if(i==1) "Δ Frecuencia" else ""
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    scale_x_continuous(breaks = seq(1, ncol(mat_zoom), by = salto_x), 
                       labels = seq(roi_start, roi_end, by = salto_x))
  
  lista_logos[[i]] <- p_zoom
}

# ==========================================
# 🛠️ ENSAMBLAJE VECTORIAL Y GUARDADO
# ==========================================
panel_logos <- wrap_plots(lista_logos, nrow = 1)
figura_final <- p_global / panel_logos + plot_layout(heights = c(1, 1.2))

ruta_salida_png <- file.path(outdir, "Figura_Paisaje_Genetico_Diferencial.png")
ruta_salida_pdf <- file.path(outdir, "Figura_Paisaje_Genetico_Diferencial.pdf")

ggsave(ruta_salida_png, plot = figura_final, width = 16, height = 10, dpi = 300)
ggsave(ruta_salida_pdf, plot = figura_final, width = 16, height = 10, dpi = 300, device = pdf)

cat(sprintf("[R-Script] ¡Éxito total! Figuras guardadas en: %s\n", outdir))