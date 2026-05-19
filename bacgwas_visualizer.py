#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BacGWAS-Visualizer (Sprint 1 & 2: Magnitud, Mapeo y Clustering)
Pipeline for Allele Frequency Analysis for Bacterial GWAS.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gzip
from matplotlib.lines import Line2D
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="BacGWAS-Visualizer: Motor de magnitud, mapeo y clustering.")
    
    # Inputs OBLIGATORIOS
    parser.add_argument('-v', '--vcf', required=True, help="Archivo VCF de entrada (.vcf o .vcf.gz).")
    parser.add_argument('-p', '--pheno', required=True, help="Archivo de fenotipos (.txt). Formato: ID  Fenotipo (1/0).")
    parser.add_argument('-g', '--gwas', required=True, help="Output del GWAS de Pyseer (.tsv completo).")
    
    # Inputs OPCIONALES (Para el Sprint 2)
    parser.add_argument('-r', '--ref', required=False, help="Genoma/Gen de referencia (.fasta) para mapear y agrupar unitigs.")
    
    # Parámetros CUSTOMIZABLES
    parser.add_argument('-t', '--type', choices=['snp', 'unitig', 'kmer'], default='unitig', help="Tipo de variante.")
    parser.add_argument('-c', '--case-name', default="Casos (1)", help="Etiqueta para fenotipo 1")
    parser.add_argument('-n', '--control-name', default="Controles (0)", help="Etiqueta para fenotipo 0")
    
    # Outputs y Filtros
    parser.add_argument('-o', '--outdir', default="BacGWAS_Output", help="Carpeta de salida.")
    parser.add_argument('--pval', type=float, default=1e-5, help="Umbral de significancia (defecto: 1e-5).")
    
    return parser.parse_args()

def load_phenotypes(pheno_file):
    print("\n>>> [1/5] Cargando y limpiando metadatos fenotípicos...")
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
            if line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                for idx, sample in enumerate(samples):
                    if sample in casos_ids: casos_indices.append(idx + 9)
                    elif sample in controles_ids: controles_indices.append(idx + 9)
                print(f"    ✔️ Columnas VCF emparejadas con éxito.")
                continue
            
            parts = line.strip().split('\t')
            var_id = f"{parts[0]}_{parts[1]}"  # Forzamos a usar Cromosoma_Posición en lugar de la secuencia para hacer match con Pyseer
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

# ==========================================
# 🗺️ MÓDULO SPRINT 2: MAPEO Y CLUSTERING
# ==========================================
def mapear_y_agrupar_unitigs(df_hits, ref_file, outdir):
    print("\n>>> [5/5] Iniciando Mapeo Espacial y Clustering...")
    try:
        ref_record = next(SeqIO.parse(ref_file, "fasta"))
        ref_seq = str(ref_record.seq).upper()
        print(f"    ✔️ Referencia cargada: {ref_record.id} ({len(ref_seq)} pb)")
    except Exception as e:
        print(f"    ❌ ERROR leyendo el FASTA de referencia: {e}")
        return

    # 1. Mapeo Bidireccional
    print("    ⏳ Mapeando hits en la secuencia...")
    mapeados = []
    for _, row in df_hits.iterrows():
        unitig = row['variant'].upper()
        
        # Intento Forward
        pos = ref_seq.find(unitig)
        strand = '+'
        
        # Intento Reverse Complement
        if pos == -1:
            unitig_rc = str(Seq(unitig).reverse_complement())
            pos = ref_seq.find(unitig_rc)
            strand = '-'
            
        if pos != -1:
            mapeados.append({
                'variant': unitig,
                'start': pos + 1,  # Biología usa base 1 (no base 0 como Python)
                'end': pos + len(unitig),
                'strand': strand,
                'pvalue': row['lrt-pvalue'],
                'grupo': row['Grupo Dominante']
            })
            
    df_map = pd.DataFrame(mapeados)
    if df_map.empty:
        print("    ⚠️ Ningún hit mapeó en esta referencia. Revisa si es la correcta.")
        return
        
    print(f"    ✔️ {len(df_map)}/{len(df_hits)} hits encontrados en la referencia.")
    df_map.to_csv(os.path.join(outdir, "Hits_Mapeados_Detalle.csv"), index=False)

    # 2. Clustering (Fusión de Intervalos superpuestos)
    print("    ⏳ Agrupando hits encimados en Regiones de Interés (ROIs)...")
    df_map = df_map.sort_values('start').reset_index(drop=True)
    
    rois = []
    current_roi = None
    
    for _, row in df_map.iterrows():
        if current_roi is None:
            current_roi = {'start': row['start'], 'end': row['end'], 'hits_contenidos': 1}
        else:
            # Si el inicio de este unitig cae "dentro" o "pegado" al anterior (margen de 5 bases)
            if row['start'] <= current_roi['end'] + 5:
                current_roi['end'] = max(current_roi['end'], row['end']) # Expandir la burbuja
                current_roi['hits_contenidos'] += 1
            else:
                rois.append(current_roi) # Guardar ROI completada
                current_roi = {'start': row['start'], 'end': row['end'], 'hits_contenidos': 1} # Empezar nueva
                
    if current_roi: rois.append(current_roi) # Guardar la última
        
    df_rois = pd.DataFrame(rois)
    df_rois.index = [f"ROI_{i+1}" for i in range(len(df_rois))]
    df_rois.index.name = "ID_Region"
    
    ruta_rois = os.path.join(outdir, "ROIs_para_R.csv")
    df_rois.to_csv(ruta_rois)
    print(f"    ✨ Clustering exitoso: Se crearon {len(df_rois)} Regiones de Interés.")
    print(f"    📄 Archivo puente para R listo: {ruta_rois}")

def main():
    args = parse_arguments()
    if not os.path.exists(args.outdir): os.makedirs(args.outdir)
        
    casos_ids, controles_ids = load_phenotypes(args.pheno)
    df_vcf_freqs = process_vcf(args.vcf, casos_ids, controles_ids)
    
    print(">>> [3/5] Integrando datos con estadísticas GWAS...")
    try:
        df_gwas = pd.read_csv(args.gwas, sep="\t")
        
        columnas_gwas = df_gwas.columns.tolist()
        col_buscada = 'variant'
        
        if col_buscada not in columnas_gwas:
            encontrada = False
            for col in columnas_gwas:
                if col.strip().lower() == col_buscada:
                    df_gwas.rename(columns={col: col_buscada}, inplace=True)
                    print(f"    ✔️ Columna '{col}' estandarizada a '{col_buscada}'.")
                    encontrada = True
                    break
            
            if not encontrada:
                primera_col = columnas_gwas[0]
                df_gwas.rename(columns={primera_col: col_buscada}, inplace=True)
                print(f"    ✔️ Usando la primera columna '{primera_col}' como ID de variante.")
        
        df_gwas = df_gwas.dropna(subset=['variant'])
        
        # --- NUEVO: Limpieza Quirúrgica de IDs ---
        # El VCF tiene '26695_1319583', pero Pyseer tiene '26695_1319583_A_T'.
        # Esta función corta todo lo que esté después del segundo guion bajo '_'
        def limpiar_id_pyseer(id_sucio):
            partes = str(id_sucio).strip().split('_')
            # Si tiene más de 2 partes (ej. Chrom_Pos_Ref_Alt), nos quedamos solo con las 2 primeras
            if len(partes) >= 2:
                return f"{partes[0]}_{partes[1]}"
            return str(id_sucio).strip()

        df_gwas['variant'] = df_gwas['variant'].apply(limpiar_id_pyseer)
        df_vcf_freqs['variant'] = df_vcf_freqs['variant'].astype(str).str.strip()
        
        df_plot = pd.merge(df_gwas, df_vcf_freqs, on='variant', how='inner')
        if df_plot.empty:
            print("\n    ❌ ERROR: El cruce dio 0 resultados.")
            print(f"    🔍 Ejemplo GWAS: '{df_gwas['variant'].iloc[0]}' | Ejemplo VCF: '{df_vcf_freqs['variant'].iloc[0]}'")
            exit(1)
            
        print(f"    ✔️ ¡Cruce exitoso! {len(df_plot)} variantes emparejadas.")
        
    except Exception as e:
        print(f"❌ ERROR cruzando archivos: {e}")
        exit(1)

    df_plot['log10_p'] = -np.log10(df_plot['P'])
#    df_plot['log10_p'] = -np.log10(df_plot['lrt-pvalue']) # si se llama lrt-pvalue
    df_plot['Grupo Dominante'] = np.where(df_plot['delta_freq'] > 0, args.case_name, args.control_name)
    df_plot = df_plot.sample(frac=1, random_state=42).reset_index(drop=True)
    df_plot.to_csv(os.path.join(args.outdir, "GWAS_Matriz_Integrada.csv"), index=False)
    
    generar_graficos(df_plot, args)
    
    df_sig = df_plot[df_plot['P'] < args.pval].copy()
    #df_sig = df_plot[df_plot['lrt-pvalue'] < args.pval].copy() # si se llama lrt-pvalue
    if not df_sig.empty:
        ruta_hits = os.path.join(args.outdir, "HITS_SIGNIFICATIVOS.csv")
        df_sig.to_csv(ruta_hits, index=False)
        print(f"    ✔️ {len(df_sig)} hits superaron el umbral (p < {args.pval}).")
        
        # --- ACTIVADOR DEL SPRINT 2 ---
        if args.ref:
            mapear_y_agrupar_unitigs(df_sig, args.ref, args.outdir)
        else:
            print("\n>>> [5/5] Mapeo Omitido (No se proporcionó --ref).")
    else:
        print(f"    ⚠️ Ningún variante superó el umbral (p < {args.pval}).")

    print("\n🚀 ¡Pipeline ejecutado exitosamente!")

if __name__ == "__main__":
    main()