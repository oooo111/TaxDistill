import sys
import pandas as pd

def load_gtdb_mapping(metadata_path):
    print(f"正在加载映射表: {metadata_path} ...")
    # 修改此处列名
    df = pd.read_csv(metadata_path, sep='\t', low_memory=False, usecols=['ncbi_taxid', 'gtdb_taxonomy'])
    
    # 建立映射字典
    mapping = {}
    for tid, line in zip(df['ncbi_taxid'], df['gtdb_taxonomy']):
        if pd.notnull(tid):
            mapping[str(int(tid))] = line
    return mapping

def convert_raw_kraken_to_gtdb(kraken_file, mapping, output_file):
    results = []
    with open(kraken_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            status = parts[0]
            contig_id = parts[1]
            taxid = parts[2]
            
            if status == 'U' or taxid == '0':
                gtdb_label = "Unclassified"
            else:
                gtdb_label = mapping.get(taxid, "Unclassified")
            
            results.append({'anonymous_contig_id': contig_id, 'ground_truth_label': gtdb_label})
    
    out_df = pd.DataFrame(results)
    # 输出为 TSV 格式，包含您要求的列名表头
    out_df.to_csv(output_file, sep='\t', index=False)
    print(f"转换完成！输出文件: {output_file}")

if __name__ == "__main__":
    # 请确保路径正确
    metadata_file = "/disk1/yerongye/TaxDistill/Taxnometer/vamb/eval/bac120_metadata_r226.tsv"
    kraken_raw_input = "/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Airways/centrifuge_Airways_output.txt"
    output_tsv = "/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Airways/centrifuge_Airways.tsv"

    mapping_dict = load_gtdb_mapping(metadata_file)
    if mapping_dict:
        convert_raw_kraken_to_gtdb(kraken_raw_input, mapping_dict, output_tsv)