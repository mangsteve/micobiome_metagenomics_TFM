process mergeCentrifuge {
  label 'mg18_merge_centrifuge'
  conda params.mergeCentrifuge.conda
  cpus 1
  memory '4 GB'
  publishDir "$results_dir/mg18_centrifuge", mode: 'copy'

  input:
    path centrifuge_reports, each: "$params.mergeCentrifuge.output_dir/*.txt"

  output:
    path 'merged_centrifuge_results.txt'

  script:
  """
  python3 << EOF
import pandas as pd
import glob

# Obtener todos los reportes de centrifuge
files = glob.glob('$params.mergeCentrifuge.output_dir/*.txt')

# Inicializar listas de DataFrames para lecturas únicas y abundancia
dfs_reads = []
dfs_abundance = []

for file in files:
    sample_name = file.split('/')[-1].replace('_centrifuge_report.txt', '')
    
    # Leer el archivo de salida de Centrifuge y extraer los parámetros clave
    df = pd.read_csv(file, sep='\\t', usecols=['name', 'numUniqueReads', 'abundance'])
    
    # Crear DataFrames separados para numUniqueReads y abundance
    df_reads = df[['name', 'numUniqueReads']].copy()
    df_reads.columns = ['Species', f'{sample_name}_reads']
    dfs_reads.append(df_reads)
    
    df_abundance = df[['name', 'abundance']].copy()
    df_abundance.columns = ['Species', f'{sample_name}_abundance']
    dfs_abundance.append(df_abundance)

# Merge de DataFrames para lecturas únicas
merged_reads = dfs_reads[0]
for df in dfs_reads[1:]:
    merged_reads = pd.merge(merged_reads, df, on='Species', how='outer')

# Merge de DataFrames para abundancia
merged_abundance = dfs_abundance[0]
for df in dfs_abundance[1:]:
    merged_abundance = pd.merge(merged_abundance, df, on='Species', how='outer')

# Combinar ambas tablas (lecturas y abundancia) en una sola
merged_final = pd.merge(merged_reads, merged_abundance, on='Species', how='outer')

# Rellenar NaNs con 0s y guardar el archivo combinado
merged_final.fillna(0, inplace=True)
merged_final.to_csv('merged_centrifuge_results.txt', sep='\\t', index=False)

EOF
  """
}
