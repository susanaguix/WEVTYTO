import os
import re
import subprocess
import csv
import pandas as pd

SCRIPT_DIR = os.getcwd()
SAMPLE_DIR = os.getcwd()
REFERENCE_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(os.getcwd(), "results")

os.makedirs(OUTPUT_DIR, exist_ok=True)

def filter_fastq(f):
    sample = re.sub(r'\.fastq\.gz$', '', os.path.basename(f))  
    filtered_fastq = os.path.join(OUTPUT_DIR, f"{sample}_filtered.fastq") 
    if os.path.isfile(filtered_fastq):
        print(f"Skipping filtering for sample {sample} as the output file already exists")
        return 0
    print(f"Filtering sample {sample}")
    command = [
        "vsearch",
        "--fastx_filter",
        os.path.join(SAMPLE_DIR, f),
        "--fastq_stripleft",
        "20",
        "--fastq_stripright",
        "20",
        "--fastq_qmax",
        "90",
        "--fastqout",
        filtered_fastq, 
    ]
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        print(f"Error: filtering failed for sample {sample}")
        return 1
    return 0

def convert_to_fasta(f):
    sample = os.path.splitext(f)[0]
    filtered_fasta = os.path.join(OUTPUT_DIR, f"{sample}.fasta")
    if os.path.isfile(filtered_fasta):
        print(f"Skipping conversion to fasta for {sample} as the output file already exists")
        return 0
    print(f"Converting {sample} to fasta")
    command = ["seqtk", "seq", "-a", os.path.join(OUTPUT_DIR, f)]
    try:
        with open(filtered_fasta, "w") as output_file:
            subprocess.check_call(command, stdout=output_file)
    except subprocess.CalledProcessError:
        print(f"Error: converting to fasta failed for {sample}")
        return 1
    return 0

def cluster_fasta(f):
    sample = os.path.splitext(f)[0]
    cluster_fasta = os.path.join(OUTPUT_DIR, f"{sample}_clusters_consensus.fasta")
    if os.path.isfile(cluster_fasta):
        print(f"Skipping clustering for {sample} as the output file already exists")
        return 0
    print(f"Clustering {sample}")
    command = [
        "vsearch",
        "--cluster_fast",
        os.path.join(OUTPUT_DIR, f),
        "--id",
        "0.95",
        "--consout",
        cluster_fasta,
        "--minseqlength",
        "800",
        "--maxseqlength",
        "1200",
        "--sizeout",
        "--sizein",
        "--clusterout_sort",
        "--strand",
        "both",
        "--threads",
        "4",
    ]
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        print(f"Error: clustering failed for {sample}")
        return 1
    return 0

def blast_search(f, output_matched_fasta=True):
    sample = os.path.splitext(f)[0]
    blast_results = os.path.join(OUTPUT_DIR, f"{sample}_results.txt")
    if os.path.isfile(blast_results):
        print(f"Skipping BLAST search for {sample} as the output file already exists")
        return 0

    print(f"Performing BLAST search for {sample}")

    alignment_output = os.path.join(OUTPUT_DIR, f"{sample}_alignments.txt")
    matched_fasta = os.path.join(OUTPUT_DIR, f"{sample}_matched.fasta")
    updated_fasta = os.path.join(OUTPUT_DIR, f"{sample}_matched_updated.fasta")
    notmatched_fasta = os.path.join(OUTPUT_DIR, f"{sample}_notmatched.fasta")

    command = [
        "vsearch",
        "--usearch_global",
        os.path.join(OUTPUT_DIR, f),
        "--db",
        os.path.join(REFERENCE_DIR, "ev_reference_sequences.fasta.gz"),
        "--id",
        "0.80",
        "--strand",
        "both",
        "--alnout",
        alignment_output,
        "--blast6out",
        blast_results,
        "--notmatched",
        notmatched_fasta,
        "--threads",
        "4",
        "--top_hits_only",
        "--mincols",
        "200", 
        "--maxhits",
        "1"
    ]

    if output_matched_fasta:
        command += [
            "--matched",
            matched_fasta,
            "--sizein"
        ]

    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        print(f"Error: BLAST search failed for {sample}")
        return 1

    if os.stat(matched_fasta).st_size == 0 or os.stat(blast_results).st_size == 0:
        print(f"Warning: The matched.fasta or results.txt file for {sample} is empty. Skipping updating process.")
        return 0

    return 0


def calculate_abundance(f, max_rows_per_sheet=100000):
    sample = os.path.splitext(f)[0]
    sample_name = sample.split('_filtered')[0]
    blast_results = os.path.join(OUTPUT_DIR, f"{sample}_results.txt")
    excel_results_dir = os.path.join(OUTPUT_DIR, "Excel_results")
    os.makedirs(excel_results_dir, exist_ok=True)

    if os.path.isfile(os.path.join(excel_results_dir, f"{sample_name}_enterovirus_abundance_aggregated.xlsx")):
        print(f"Skipping enterovirus abundance calculation for {sample_name} as the output file already exists")
        return 0

    print(f"Calculating enterovirus abundance for {sample_name}")

    species_list = []
    seq_count_list = []
    blast_id_list = []
    size_80_94_list = []
    size_95_100_list = []

    with open(blast_results, 'r') as blast_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                sequence_info = fields[0].split(';')
                seqs = int(sequence_info[1].split('=')[1])
                blast_id = float(fields[2])  
                species = fields[1]

                if 80 <= blast_id < 95:
                    size_80_94_list.append(seqs)
                    size_95_100_list.append(0)
                elif 95 <= blast_id <= 100:
                    size_80_94_list.append(0)
                    size_95_100_list.append(seqs)
                else:
                    size_80_94_list.append(0)
                    size_95_100_list.append(0)

                species_list.append(species)
                seq_count_list.append(seqs)
                blast_id_list.append(blast_id)

    df = pd.DataFrame({
        'Reference': species_list,
        'Sequence Count': seq_count_list,
        'BLAST id': blast_id_list,
        'Seq BLAST id (80-94.9%)': size_80_94_list,
        'Seq BLAST id (95-100%)': size_95_100_list
    })

    df['BLAST id'] = pd.to_numeric(df['BLAST id'], errors='coerce')
    df['Reference'] = df['Reference'].astype(str)
    df['EV Type'] = df['Reference'].str.split('_').str[-1]

    df['BLAST id * Seq Count'] = df['BLAST id'] * df['Sequence Count']

    df_aggregated = df.groupby('EV Type').agg({
        'Sequence Count': 'sum',
        'BLAST id': ['min', 'max'],
        'Seq BLAST id (80-94.9%)': 'sum',
        'Seq BLAST id (95-100%)': 'sum',
        'BLAST id * Seq Count': 'sum'  
    }).reset_index()

    total_sequence_count_sum = df_aggregated['Sequence Count'].sum()

    df_aggregated.insert(2, 'Percentage', (df_aggregated['Sequence Count'] / total_sequence_count_sum) * 100)

    df_aggregated.columns = ['EV Type', 'Total Sequence Count', 'Percentage', 'Min BLAST id', 'Max BLAST id',
                             'Seq BLAST id (80-94.9%)', 'Seq BLAST id (95-100%)', 'BLAST id * Seq Count']

    df_aggregated['Mean Weighted BLAST id'] = df_aggregated['BLAST id * Seq Count'] / df_aggregated['Total Sequence Count']

    aggregated_file = os.path.join(excel_results_dir, f"{sample_name}_enterovirus_abundance_aggregated.xlsx")
    with pd.ExcelWriter(aggregated_file, engine='xlsxwriter') as writer:
        df_aggregated.to_excel(writer, index=False, sheet_name='Aggregated Data')

    num_rows = len(df)
    num_sheets = -(-num_rows // max_rows_per_sheet)  
    for i in range(num_sheets):
        start_idx = i * max_rows_per_sheet
        end_idx = min((i + 1) * max_rows_per_sheet, num_rows)
        chunk_df = df.iloc[start_idx:end_idx]
        sheet_name = f'Sheet_{i + 1}'  
        excel_file = os.path.join(excel_results_dir, f"{sample_name}_enterovirus_abundance_{sheet_name}_original.xlsx")
        chunk_df.to_excel(excel_file, index=False)

    return 0
    
def rename_matched_fasta_headers(results_dir="."):

    renamed_dir = os.path.join(results_dir, "renamed_fasta")
    os.makedirs(renamed_dir, exist_ok=True)

    fasta_files = [file for file in os.listdir(results_dir) if file.endswith('_filtered_clusters_consensus_matched.fasta')]

    centroids_mapping = {}
    for results_file in os.listdir(results_dir):
        if results_file.endswith('_filtered_clusters_consensus_results.txt'):

            sample = results_file.split('_filtered_clusters_consensus_results.txt')[0]
            with open(os.path.join(results_dir, results_file), 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:  
                        centroid_info = parts[0]
                        identifier = parts[1]
                        centroids_mapping[centroid_info] = identifier

    for fasta_file in fasta_files:

        sample = fasta_file.split('_filtered_clusters_consensus_matched.fasta')[0]

        if os.path.getsize(os.path.join(results_dir, fasta_file)) == 0:
            print(f"Skipping empty file: {fasta_file}")
            continue

        output_file = os.path.join(renamed_dir, f'{sample}_renamed.fasta')

        with open(os.path.join(results_dir, fasta_file), 'r') as input_file, open(output_file, 'w') as output_file:
            for line in input_file:
                if line.startswith('>'):

                    current_header = line.strip()[1:]

                    if current_header in centroids_mapping:

                        identifier = centroids_mapping[current_header]
                        new_header = f'>{current_header}_{identifier}\n'
                        output_file.write(new_header)
                    else:
                        output_file.write(line)
                else:
                    output_file.write(line)

        print(f"Renaming successful for {sample}, output: {output_file}")

    print("Header renaming complete.")



fastq_files = [f for f in os.listdir(SAMPLE_DIR) if f.endswith(".fastq.gz")]  
filter_results = [filter_fastq(f) for f in fastq_files]

filtered_fastq_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith("_filtered.fastq")]
convert_results = [convert_to_fasta(f) for f in filtered_fastq_files]

filtered_fasta_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith("_filtered.fasta")]
cluster_results = [cluster_fasta(f) for f in filtered_fasta_files]

consensus_fasta_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith("_clusters_consensus.fasta")]
blast_results = [blast_search(f) for f in consensus_fasta_files]
abundance_results = [calculate_abundance(f) for f in consensus_fasta_files]

rename_matched_fasta_headers(results_dir=OUTPUT_DIR)
