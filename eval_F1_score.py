import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


gt_path = './marine_ground_truth_v226_lineage_strict.csv'


pred_path = './results_taxometer_fixed.tsv'

CONFIDENCE_THRESHOLD = 0.80

def extract_species(tax_str):

    if pd.isna(tax_str) or tax_str == '':
        return None
    parts = str(tax_str).split(';')
    for part in parts:
        part = part.strip()
        if part.startswith('s__') or part.startswith('s_'):
            try:
                content = part.split('_', 1)[1].lstrip('_')
                return content.strip()
            except IndexError:
                return None
    return None

def extract_species_score(score_str):

    if pd.isna(score_str) or score_str == '':
        return 0.0
    
    try:
        scores = str(score_str).strip().split(';')
        species_score = float(scores[-1])
        return species_score
    except (ValueError, IndexError):
        return 0.0

def normalize_name(name):

    if name is None:
        return ""
    return name.lower().replace(' ', '_')

def main():
    try:
        df_gt = pd.read_csv(gt_path)
    except Exception as e:
        print(f"读取GT文件失败: {e}")
        return
    try:
        df_pred = pd.read_csv(pred_path, sep='\t')
    except Exception as e:
        print(f"读取预测文件失败: {e}")
        return

    df_merged = pd.merge(df_pred, df_gt, 
                         left_on='contigs', 
                         right_on='anonymous_contig_id', 
                         how='inner')
    
    print(f"共有 {len(df_merged)} 个 Contigs 匹配到 Ground Truth。")
    df_merged['gt_species'] = df_merged['ground_truth_label'].apply(extract_species)
    df_merged['pred_species'] = df_merged['predictions'].apply(extract_species)
    df_merged['species_score'] = df_merged['scores'].apply(extract_species_score)

    # --------------------------------

    def classify(row):
        pred = row['pred_species']
        gt = row['gt_species']
        score = row['species_score']
        
        if pred is None or pred == '':
            return 'No species label'
        
        if score <= CONFIDENCE_THRESHOLD:
            return 'No species label (Low Score)'
    
        if normalize_name(pred) == normalize_name(gt):
            return 'Correct'
        else:
            return 'Wrong'
    df_merged['category'] = df_merged.apply(classify, axis=1)
    counts = df_merged['category'].value_counts()
    
    correct_count = counts.get('Correct', 0)
    wrong_count = counts.get('Wrong', 0)
    no_label_base = counts.get('No species label', 0)
    no_label_low_score = counts.get('No species label (Low Score)', 0)
    no_label_total = no_label_base + no_label_low_score
    
    total = correct_count + wrong_count + no_label_total

    precision = correct_count / (correct_count + wrong_count) if (correct_count + wrong_count) > 0 else 0
    
    recall = correct_count / total if total > 0 else 0
    
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print(f"Total Contigs: {total}")
    print(f"Correct (Blue): {correct_count}")
    print(f"Wrong (Red):    {wrong_count}")
    print(f"No label (Orange): {no_label_total} (Empty: {no_label_base}, Low Score: {no_label_low_score})")
    print("-" * 30)
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-score:  {f1:.4f}")

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(6, 8))

    categories_plot = ['Correct', 'Wrong', 'No Label']
    values_plot = [correct_count, wrong_count, no_label_total]
    colors = ['#1f77b4', '#d62728', '#ff7f0e'] 
    
    bottom = 0
    bar_width = 0.6
    
    for val, color, label in zip(values_plot, colors, categories_plot):
        ax.bar('Taxnometer', val, bottom=bottom, color=color, label=label, width=bar_width, edgecolor='white')
        if val > 0:
            ax.text('Taxnometer', bottom + val/2, str(val), ha='center', va='center', color='white', fontweight='bold')
        bottom += val

    ax.set_ylabel('Number of Contigs')
    ax.set_title(f'Taxnometer Performance (Score > {CONFIDENCE_THRESHOLD})')
    ax.legend(loc='upper right', frameon=True)
    
    output_img = 'taxnometer_evaluation.png'
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)

if __name__ == '__main__':
    main()
