# Bacterial GWAS hits visualizer

Pipeline for GWAS hits study.

Input:

```python bacgwas_visualizer.py \
  --vcf ruta/a/tus/unitigs.vcf \
  --pheno ruta/al/phenotype.txt \
  --gwas ruta/al/pyseer_output.tsv \
  --type unitig \
  --case-name "Cancer Gastrico (GC)" \
  --control-name "Gastritis (NAG)" \
  --outdir Resultados_BabA \
  --pval 1e-5```
