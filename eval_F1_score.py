import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ================= 配置路径 =================
# Ground Truth 路径 (保持不变)
gt_path = '/disk1/yerongye/Taxnometer/vamb/data/Marine/marine_ground_truth_v226_lineage_strict.csv'

# Taxnometer 预测结果路径
pred_path = '/disk1/yerongye/KD_Taxnometer_ntnf/vamb/Ablation_study/parameter_ablation/fix_temp_design_alpha/kraken2_marine_0.3_ntnf_KD_100M/results_taxometer_fixed.tsv'

CONFIDENCE_THRESHOLD = 0.80
# ===========================================

def extract_species(tax_str):
    """
    从分类字符串中提取种名 (兼容 s__ 和 s_ 前缀)
    例如: s_Bordetella pertussis -> Bordetella pertussis
    """
    if pd.isna(tax_str) or tax_str == '':
        return None
    
    # 按分号分割层级
    parts = str(tax_str).split(';')
    for part in parts:
        part = part.strip()
        # 匹配 s__ 或 s_ 前缀
        if part.startswith('s__') or part.startswith('s_'):
            # split('_', 1) 确保只分割第一个下划线
            try:
                content = part.split('_', 1)[1].lstrip('_')
                return content.strip()
            except IndexError:
                return None
    return None

def extract_species_score(score_str):
    """
    从分数各个层级中提取最后一个分数 (Species level score)
    输入示例: "1.0;1.0;1.0;1.0;1.0;1.0;0.99999"
    输出: 0.99999 (float)
    """
    if pd.isna(score_str) or score_str == '':
        return 0.0
    
    try:
        # 按分号分割，取最后一个值
        scores = str(score_str).strip().split(';')
        species_score = float(scores[-1])
        return species_score
    except (ValueError, IndexError):
        return 0.0

def normalize_name(name):
    """
    标准化名称以便比较：转小写，空格转下划线
    """
    if name is None:
        return ""
    return name.lower().replace(' ', '_')

def main():
    print("正在读取文件...")
    # 1. 读取 Ground Truth
    try:
        df_gt = pd.read_csv(gt_path)
    except Exception as e:
        print(f"读取GT文件失败: {e}")
        return

    # 2. 读取 Taxnometer 预测文件
    try:
        # Taxnometer 输出可能有 header，根据你提供的示例包含 'contigs', 'predictions', 'scores'
        df_pred = pd.read_csv(pred_path, sep='\t')
    except Exception as e:
        print(f"读取预测文件失败: {e}")
        return

    print("正在匹配 Contigs...")
    # 3. 合并数据
    df_merged = pd.merge(df_pred, df_gt, 
                         left_on='contigs', 
                         right_on='anonymous_contig_id', 
                         how='inner')
    
    print(f"共有 {len(df_merged)} 个 Contigs 匹配到 Ground Truth。")

    # 4. 提取数据列
    print("正在提取物种名和置信度分数...")
    df_merged['gt_species'] = df_merged['ground_truth_label'].apply(extract_species)
    df_merged['pred_species'] = df_merged['predictions'].apply(extract_species)
    df_merged['species_score'] = df_merged['scores'].apply(extract_species_score)

    # --- 调试信息：查看提取是否正确 ---
    print("\n[Debug] 数据提取示例 (前5行):")
    print(df_merged[['predictions', 'pred_species', 'scores', 'species_score']].head())
    print("-" * 30)
    # --------------------------------

    # 5. 定义分类逻辑 (加入阈值判断)
    def classify(row):
        pred = row['pred_species']
        gt = row['gt_species']
        score = row['species_score']
        
        # 情况 A: 预测结果为空
        if pred is None or pred == '':
            return 'No species label'
        
        # 情况 B: 分数过低 (过滤掉) -> 视为没有有效标签
        # 你的要求： > 0.95 才算数
        if score <= CONFIDENCE_THRESHOLD:
            return 'No species label (Low Score)'
        
        # 情况 C: 分数达标，进行比对
        if normalize_name(pred) == normalize_name(gt):
            return 'Correct'
        else:
            return 'Wrong'

    # 应用分类
    df_merged['category'] = df_merged.apply(classify, axis=1)

    # 6. 统计结果
    counts = df_merged['category'].value_counts()
    
    correct_count = counts.get('Correct', 0)
    wrong_count = counts.get('Wrong', 0)
    # 将 "No species label" 和 "No species label (Low Score)" 合并统计
    no_label_base = counts.get('No species label', 0)
    no_label_low_score = counts.get('No species label (Low Score)', 0)
    no_label_total = no_label_base + no_label_low_score
    
    total = correct_count + wrong_count + no_label_total

    # 7. 计算指标
    # Precision = Correct / (Correct + Wrong) -> 只看置信度够且做出了预测的部分
    precision = correct_count / (correct_count + wrong_count) if (correct_count + wrong_count) > 0 else 0
    
    # Recall = Correct / Total -> 低分数的也被视为“漏召回”
    recall = correct_count / total if total > 0 else 0
    
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\n=== 评估结果 (Taxnometer, Threshold > {CONFIDENCE_THRESHOLD}) ===")
    print(f"Total Contigs: {total}")
    print(f"Correct (Blue): {correct_count}")
    print(f"Wrong (Red):    {wrong_count}")
    print(f"No label (Orange): {no_label_total} (Empty: {no_label_base}, Low Score: {no_label_low_score})")
    print("-" * 30)
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-score:  {f1:.4f}")

    # 8. 绘图
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
    print(f"\n图表已保存为: {output_img}")

if __name__ == '__main__':
    main()