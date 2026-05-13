#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BacGWAS-Visualizer (Fase 1: Motor de Magnitud)
Pipeline modular para el análisis y visualización de GWAS bacteriano.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gzip
from matplotlib.lines import Line2D

def parse_arguments():
    parser = argparse.ArgumentParser(description="BacGWAS-Visualizer: Motor de cálculo de magnitud y gráficos.")
    
    # Inputs OBLIGATORIOS
    parser.add_argument('-v', '--vcf', required=True, help="Archivo VCF de entrada (.vcf o .vcf.gz).")
    parser.add_argument('-p', '--pheno', required=True, help="Archivo de fenotipos (.txt). Formato: ID  Fenotipo (1/0).")
    parser.add_argument('-g', '--gwas', required=True, help="Output del GWAS de Pyseer (.tsv).")
    
    # Parámetros CUSTOMIZABLES
    parser.add_argument('-t', '--type', choices=['snp', 'unitig', 'kmer'], default='unitig', help="Tipo de variante (afecta los títulos).")
    parser.add_argument('-c', '--case-name', default="Casos (1)", help="Etiqueta para fenotipo 1 (Ej. GC)")
    parser.add_argument('-n', '--control-name', default="Controles (0)", help="Etiqueta para fenotipo 0 (Ej. NAG)")
    
    # Outputs y Filtros
    parser.add_argument('-o', '--outdir', default="BacGWAS_Output", help="Carpeta de salida.")
    parser.add_argument('--pval', type=float, default=1e-5, help="Umbral de significancia (defecto: 1e-5).")
    
    return parser.parse_args()

def load_phenotypes(pheno_file):
    print(">>> [1/5] Cargando y limpiando metadatos fenotípicos...")
    try:
        df_pheno = pd.read_csv(pheno_file, sep=r'\s+', header=None)
        if isinstance(df_pheno.iloc[0, 1], str) and not df_pheno.iloc[0, 1].replace('.', '', 1).isdigit():
            df_pheno = df_pheno.iloc[1:].copy()
            
        df_pheno[0] = df_pheno[0].astype(str).str.strip()
        df_pheno[1] = pd.to_numeric(df_pheno[1])
        
        casos_ids = set(df_pheno[df_pheno[1] == 1][0].tolist())
        controles_ids = set(df_pheno[df_pheno[1] == 0][0].tolist())
        
        print(f"    ✔️ Casos ({len(casos_ids)}) | Controles ({len(controles_ids)})")
        return casos_ids, controles_ids
    except Exception as e:
        print(f"❌ ERROR leyendo fenotipos: {e}")
        exit(1)

def process_vcf(vcf_file, casos_ids, controles_ids):
    print(">>> [2/5] Extrayendo matemáticas desde la matriz VCF...")
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    casos_indices, controles_indices = [], []
    resultados = []
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'): continue
            
            # Mapear columnas a cepas
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                for idx, sample in enumerate(samples):
                    if sample in casos_ids: casos_indices.append(idx + 9)
                    elif sample in controles_ids: controles_indices.append(idx + 9)
                print(f"    ✔️ Columnas VCF emparejadas con éxito.")
                continue
            
            # Conteo súper rápido por fila
            parts = line.strip().split('\t')
            var_id = parts[2]
            if var_id == '.': var_id = f"{parts[0]}_{parts[1]}" 
            
            conteo_casos = sum(1 for i in casos_indices if '1' in parts[i].split(':')[0])
            conteo_controles = sum(1 for i in controles_indices if '1' in parts[i].split(':')[0])
            
            freq_case = conteo_casos / len(casos_indices) if casos_indices else 0
            freq_control = conteo_controles / len(controles_indices) if controles_indices else 0
            
            resultados.append({
                'variant': var_id,
                'freq_casos': freq_case,
                'freq_controles': freq_control,
                'delta_freq': freq_case - freq_control,
                'abs_delta_freq': abs(freq_case - freq_control)
            })
            
    return pd.DataFrame(resultados)

def generar_graficos(df_plot, args):
    print(">>> [4/5] Renderizando gráficos de alta resolución...")
    paleta = {args.case_name: '#d63031', args.control_name: '#0984e3'}
    titulo_eje = fr"Cambio de la variante (|$\Delta$ Frecuencia|)"
    
    plt.figure(figsize=(13, 8))
    sns.regplot(data=df_plot, x='abs_delta_freq', y='log10_p', scatter=False, color='gray', line_kws={'linestyle':'-.', 'linewidth': 1.5, 'alpha': 0.6})
    
    scatter = sns.scatterplot(data=df_plot, x='abs_delta_freq', y='log10_p', size='af', sizes=(20, 300), hue='Grupo Dominante', palette=paleta, alpha=0.7, edgecolor='black', linewidth=0.4)
    
    plt.axvline(0, color='black', linewidth=1.5)
    plt.axhline(-np.log10(args.pval), color='black', linestyle='--', linewidth=1.5)
    
    plt.title(f"Magnitud Absoluta vs Significancia ({args.type.upper()})", fontsize=16, fontweight='bold', pad=15)
    plt.xlabel(titulo_eje, fontsize=13)
    plt.ylabel(r"$-\log_{10}$(p-value)", fontsize=13)
    
    handles, labels = scatter.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='black', linestyle='--', linewidth=1.5))
    labels.append(f'Umbral (p={args.pval})')
    plt.legend(handles=handles, labels=labels, loc='upper left', bbox_to_anchor=(1.02, 1), title="Simbología", frameon=True, shadow=True)
    
    plt.grid(True, which='major', linestyle=':', linewidth=0.5, alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "Plot_Magnitud_Absoluta.png"), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    args = parse_arguments()
    if not os.path.exists(args.outdir): os.makedirs(args.outdir)
        
    casos_ids, controles_ids = load_phenotypes(args.pheno)
    df_vcf_freqs = process_vcf(args.vcf, casos_ids, controles_ids)
    
    print(">>> [3/5] Integrando datos con estadísticas GWAS...")
    try:
        df_gwas = pd.read_csv(args.gwas, sep="\t")
        df_gwas['variant'] = df_gwas['variant'].astype(str).str.strip()
        df_vcf_freqs['variant'] = df_vcf_freqs['variant'].astype(str).str.strip()
        
        # --- RAYOS X: ¿QUÉ ESTÁ LEYENDO PYTHON? ---
        print(f"    🔍 Total filas en GWAS: {len(df_gwas)}")
        print(f"    🔍 Total filas en VCF: {len(df_vcf_freqs)}")
        
        if len(df_gwas) > 0: 
            print(f"    🔍 Ejemplo ID GWAS (Pyseer): '{df_gwas['variant'].iloc[0]}'")
        if len(df_vcf_freqs) > 0: 
            print(f"    🔍 Ejemplo ID VCF: '{df_vcf_freqs['variant'].iloc[0]}'")
            
        df_plot = pd.merge(df_gwas, df_vcf_freqs, on='variant', how='inner')
        print(f"    ✔️ Variantes que SÍ coincidieron: {len(df_plot)}")
        
        if df_plot.empty:
            print("\n    ❌ ERROR: El cruce dio 0 resultados. Los IDs no coinciden.")
            exit(1)
            
    except Exception as e:
        print(f"❌ ERROR cruzando archivos: {e}")
        exit(1)

    df_plot['log10_p'] = -np.log10(df_plot['lrt-pvalue'])
    df_plot['Grupo Dominante'] = np.where(df_plot['delta_freq'] > 0, args.case_name, args.control_name)
    df_plot = df_plot.sample(frac=1, random_state=42).reset_index(drop=True) # Mezclar puntos
    
    df_plot.to_csv(os.path.join(args.outdir, "GWAS_Matriz_Integrada.csv"), index=False)
    
    generar_graficos(df_plot, args)
    
    # ==========================================
    # 🚀 EL PUENTE HACIA EL SPRINT 2
    # ==========================================
    print(">>> [5/5] Preparando hand-off para el Módulo de Mapeo...")
    df_sig = df_plot[df_plot['lrt-pvalue'] < args.pval].copy()
    
    if not df_sig.empty:
        ruta_hits = os.path.join(args.outdir, "HITS_SIGNIFICATIVOS.csv")
        df_sig.to_csv(ruta_hits, index=False)
        print(f"    ✔️ {len(df_sig)} hits superaron el umbral (p < {args.pval}).")
        print(f"    ✔️ Guardados en: {ruta_hits}")
    else:
        print(f"    ⚠️ Ningún variante superó el umbral (p < {args.pval}).")

    print("\n✨ ¡Fase 1 completada con éxito! Revisa tu carpeta de salida.")

if __name__ == "__main__":
    main()
