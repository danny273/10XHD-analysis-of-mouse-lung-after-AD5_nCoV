###总体情况
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.sparse import issparse
import scanpy as sc
import glob
import os
from tqdm import tqdm
import warnings

# 忽略特定警告
warnings.filterwarnings("ignore", category=FutureWarning)




def add_plaque_info(adata, plaque_info):
    """添加斑块信息到AnnData对象"""
    print(f"为样本添加斑块信息...")

    # 创建初始列
    adata.obs['pla_nich'] = 'neg'
    adata.obs['plaque'] = np.nan

    # 处理细胞ID格式
    plaque_info['cell_id'] = [int(s.split('_')[1]) for s in plaque_info['cell_id']]
    plaque_info['cell_id'] = plaque_info['cell_id'].astype(str)
    adata.obs['object_id'] = adata.obs.index.astype(str)

    # 标记斑块微环境细胞
    plaque_mask = adata.obs['object_id'].isin(plaque_info['cell_id'])
    adata.obs.loc[plaque_mask, 'pla_nich'] = 'pos'

    # 添加斑块名称
    plaque_map = plaque_info.set_index('cell_id')['pla_id'].to_dict()
    adata.obs.loc[plaque_mask, 'plaque'] = adata.obs.loc[plaque_mask, 'object_id'].map(plaque_map)

    # 统计信息
    n_plaque_cells = plaque_mask.sum()
    plaque_counts = adata.obs.loc[plaque_mask, 'plaque'].value_counts()

    print(f"斑块微环境细胞数量: {n_plaque_cells}/{len(adata)}")
    print(f"斑块数量: {len(plaque_counts)}")
    print(f"平均每个斑块细胞数: {n_plaque_cells / len(plaque_counts):.1f}")

    return adata


def calculate_distance_to_plaque(adata):
    """计算所有细胞到最近斑块微环境细胞的距离"""
    # 获取pos细胞坐标
    pos_mask = adata.obs["pla_nich"] == "pos"
    if pos_mask.sum() == 0:
        print("没有找到pos细胞，无法计算距离")
        return None

    pos_coords = adata.obs.loc[pos_mask, ["x", "y"]].values
    all_coords = adata.obs[["x", "y"]].values

    # 计算所有细胞到最近pos细胞的距离
    pos_kdtree = cKDTree(pos_coords)
    distances, _ = pos_kdtree.query(all_coords, k=1)  # 最近距离

    return distances


def preprocess_samples(file_list, plaque_info):
    """预处理所有样本，添加斑块信息并计算距离"""
    processed_samples = {}

    for file_path in tqdm(file_list, desc="预处理样本"):
        sample_name = os.path.basename(file_path).split('.')[0]
        print(f"\n处理样本: {sample_name}")

        # 加载数据
        adata = sc.read_h5ad(file_path)

        # 筛选当前样本的斑块信息
        sample_plaque_info = plaque_info[plaque_info["source"] == sample_name]

        # 如果样本没有斑块信息，跳过
        if sample_plaque_info.empty:
            print(f"样本 {sample_name} 没有斑块信息，跳过")
            continue

        # 添加斑块信息
        adata = add_plaque_info(adata, sample_plaque_info)

        # 计算距离
        distances = calculate_distance_to_plaque(adata)

        if distances is not None:
            # 存储处理后的样本数据
            processed_samples[sample_name] = {
                'adata': adata,
                'distances': distances
            }

    return processed_samples


def analyze_genes(processed_samples, genes_of_interest, max_dist=800, bin_size=50):
    """
    分析指定基因的空间表达梯度

    参数:
    processed_samples: 预处理后的样本字典
    genes_of_interest: 要分析的基因列表
    max_dist: 最大分析距离
    bin_size: 距离分组大小

    返回:
    dict: 包含分析结果的字典
    """
    # 准备距离分组
    bins = np.arange(0, max_dist + bin_size, bin_size)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bin_centers)

    # 初始化结果存储
    results = {
        gene: {
            "all_samples_mean": np.zeros(n_bins),
            "all_samples_sem": np.zeros(n_bins),
            "sample_counts": np.zeros(n_bins),
            "sample_results": []  # 存储每个样本的原始结果
        }
        for gene in genes_of_interest
    }

    # 添加细胞计数结果
    results["cell_counts"] = {
        "all_samples_mean": np.zeros(n_bins),
        "all_samples_sem": np.zeros(n_bins),
        "sample_counts": np.zeros(n_bins),
        "sample_results": []
    }

    # 处理每个样本
    for sample_name, sample_data in processed_samples.items():
        adata = sample_data['adata']
        distances = sample_data['distances']

        # 计算每个距离区间的细胞数量
        cell_counts = np.zeros(len(bins) - 1)
        for i in range(len(bins) - 1):
            bin_mask = (distances >= bins[i]) & (distances < bins[i + 1])
            cell_counts[i] = np.sum(bin_mask)

        # 存储细胞计数结果
        results["cell_counts"]["sample_results"].append(cell_counts)

        # 分析每个基因
        for gene in genes_of_interest:
            # 获取基因表达数据
            if gene in adata.var_names:
                expr_vector = adata[:, gene].X
                if issparse(expr_vector):
                    gene_expr = expr_vector.toarray().flatten()
                else:
                    gene_expr = expr_vector.flatten()
            else:
                print(f"基因 {gene} 在样本 {sample_name} 中不存在，跳过")
                # 添加NaN数组保持对齐
                results[gene]["sample_results"].append(np.full(len(bins) - 1, np.nan))
                continue

            # 初始化基因结果存储
            gene_mean_expr = np.zeros(len(bins) - 1) * np.nan

            # 按距离分组计算基因表达
            for i in range(len(bins) - 1):
                bin_mask = (distances >= bins[i]) & (distances < bins[i + 1])
                if np.sum(bin_mask) > 0:
                    bin_expr = gene_expr[bin_mask]
                    gene_mean_expr[i] = np.mean(bin_expr)

            # 存储基因结果
            results[gene]["sample_results"].append(gene_mean_expr)

    # 汇总所有样本的结果
    for gene in genes_of_interest + ["cell_counts"]:
        all_samples = np.array(results[gene]["sample_results"])

        # 计算每个距离bin的平均值
        for bin_idx in range(n_bins):
            bin_values = all_samples[:, bin_idx]
            valid_values = bin_values[~np.isnan(bin_values)]

            if len(valid_values) > 0:
                results[gene]["all_samples_mean"][bin_idx] = np.mean(valid_values)
                results[gene]["all_samples_sem"][bin_idx] = np.std(valid_values) / np.sqrt(len(valid_values))
                results[gene]["sample_counts"][bin_idx] = len(valid_values)
            else:
                results[gene]["all_samples_mean"][bin_idx] = np.nan
                results[gene]["all_samples_sem"][bin_idx] = np.nan
                results[gene]["sample_counts"][bin_idx] = 0

    # 添加bin信息
    results["bins"] = bins
    results["bin_centers"] = bin_centers

    return results


# def plot_results(results, genes_of_interest, output_dir):
#     """可视化分析结果"""
#     bin_centers = results["bin_centers"]
#
#     # 基因表达梯度图（高度标准化）
#     plt.figure(figsize=(14, 10))
#     colors = plt.cm.tab10(np.linspace(0, 1, len(genes_of_interest)))
#
#     for i, gene in enumerate(genes_of_interest):
#         # 获取原始均值和标准误
#         mean_expr = results[gene]["all_samples_mean"]
#         sem = results[gene]["all_samples_sem"]
#
#         # 高度标准化：归一化到0-1范围
#         valid_mask = ~np.isnan(mean_expr)
#         if np.any(valid_mask):
#             min_val = np.min(mean_expr[valid_mask])
#             max_val = np.max(mean_expr[valid_mask])
#             if max_val > min_val:  # 避免除以零
#                 norm_mean = (mean_expr - min_val) / (max_val - min_val)
#                 norm_sem = sem / (max_val - min_val)  # 标准化标准误
#             else:
#                 norm_mean = mean_expr
#                 norm_sem = sem
#         else:
#             norm_mean = mean_expr
#             norm_sem = sem
#
#         # 绘制高度标准化后的曲线
#         plt.plot(bin_centers, norm_mean, 'o-', color=colors[i], label=gene, alpha=0.8)
#         plt.fill_between(bin_centers,
#                          norm_mean - norm_sem,
#                          norm_mean + norm_sem,
#                          alpha=0.2, color=colors[i])
#
#     plt.axvline(0, linestyle='--', color='gray', alpha=0.5)
#     plt.xlabel('距最近斑块微环境细胞的距离 (μm)', fontsize=14)
#     plt.ylabel('标准化表达量 (0-1)', fontsize=14)
#     plt.title('基因空间表达梯度 (高度标准化)', fontsize=16)
#     plt.legend(fontsize=12, loc='best')
#     plt.grid(alpha=0.2)
#     plt.xlim(0, bin_centers[-1])
#     plt.ylim(-0.1, 1.1)  # 固定y轴范围便于比较
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, 'multi_gene_spatial_gradient_normalized.png'), dpi=300)
#     plt.show()
#
#     # 细胞数量分布图（保持原始尺度）
#     plt.figure(figsize=(14, 6))
#     mean_counts = results["cell_counts"]["all_samples_mean"]
#     sem_counts = results["cell_counts"]["all_samples_sem"]
#
#     plt.bar(bin_centers, mean_counts, width=results["bins"][1] - results["bins"][0],
#             color='steelblue', alpha=0.7,
#             yerr=sem_counts, error_kw=dict(elinewidth=1, ecolor='darkred'))
#
#     plt.xlabel('距最近斑块微环境细胞的距离 (μm)', fontsize=14)
#     plt.ylabel('平均细胞数量', fontsize=14)
#     plt.title('不同距离区间的细胞分布', fontsize=16)
#     plt.grid(axis='y', alpha=0.2)
#     plt.xlim(0, bin_centers[-1])
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, 'cell_distribution_distance_bins.png'), dpi=300)
#     plt.show()
#
#     # 创建热图展示基因表达模式
#     plt.figure(figsize=(12, 8))
#     heatmap_data = []
#
#     for gene in genes_of_interest:
#         # 高度标准化处理
#         values = results[gene]["all_samples_mean"]
#         if np.isnan(values).all():
#             continue
#
#         # 归一化处理 (0-1)
#         min_val = np.nanmin(values)
#         max_val = np.nanmax(values)
#         if max_val > min_val:  # 避免除以零
#             norm_values = (values - min_val) / (max_val - min_val)
#         else:
#             norm_values = values * 0  # 如果所有值相同，设为0
#
#         heatmap_data.append(norm_values)
#
#     if heatmap_data:
#         heatmap_data = np.array(heatmap_data)
#         sns.heatmap(heatmap_data,
#                     cmap='viridis',
#                     xticklabels=[f"{int(c)}" for c in bin_centers],
#                     yticklabels=genes_of_interest,
#                     cbar_kws={'label': '标准化表达量 (0-1)'},
#                     vmin=0, vmax=1)  # 固定颜色范围
#
#         plt.xlabel('距斑块的距离 (μm)', fontsize=12)
#         plt.ylabel('基因', fontsize=12)
#         plt.title('斑块周围的空间表达模式', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir, 'gene_expression_heatmap_normalized.png'), dpi=300)
#         plt.show()
#     else:
#         print("没有足够的数据生成热图")


def save_results(results, genes_of_interest, output_dir):
    """保存分析结果到CSV文件"""
    bin_centers = results["bin_centers"]
    bins = results["bins"]

    # 创建结果数据框
    results_df = pd.DataFrame({
        'distance_bin_center': bin_centers,
        'distance_range': [f"{bins[i]}-{bins[i + 1]}" for i in range(len(bins) - 1)]
    })

    # 添加基因结果（原始值）
    for gene in genes_of_interest:
        results_df[f'{gene}_mean'] = results[gene]["all_samples_mean"]
        results_df[f'{gene}_sem'] = results[gene]["all_samples_sem"]
        results_df[f'{gene}_n_samples'] = results[gene]["sample_counts"]

    # 添加细胞计数结果
    results_df['cell_count_mean'] = results["cell_counts"]["all_samples_mean"]
    results_df['cell_count_sem'] = results["cell_counts"]["all_samples_sem"]
    results_df['cell_count_n_samples'] = results["cell_counts"]["sample_counts"]

    # 保存到CSV
    csv_path = os.path.join(output_dir, 'spatial_gradient_results.csv')
    results_df.to_csv(csv_path, index=False)
    print(f"结果已保存至: {csv_path}")

    # 创建标准化结果数据框
    norm_results_df = pd.DataFrame({
        'distance_bin_center': bin_centers
    })

    for gene in genes_of_interest:
        values = results[gene]["all_samples_mean"]
        if np.any(~np.isnan(values)):
            min_val = np.nanmin(values)
            max_val = np.nanmax(values)
            if max_val > min_val:
                norm_values = (values - min_val) / (max_val - min_val)
            else:
                norm_values = np.zeros_like(values)
        else:
            norm_values = np.full_like(values, np.nan)

        norm_results_df[f'{gene}_normalized'] = norm_values

    # 保存标准化结果
    norm_csv_path = os.path.join(output_dir, 'normalized_spatial_gradient_results.csv')
    norm_results_df.to_csv(norm_csv_path, index=False)
    print(f"标准化结果已保存至: {norm_csv_path}")

    return norm_results_df

# 主程序

# 加载数据文件列表
file_list = sorted([f for f in glob.glob('D:/atlas/corrected_matrices/*.h5ad') if '_with_clusters' not in f])[:8]
print("待处理文件列表:")
for f in file_list:
    print(f)

# 加载斑块信息
plaque_info = pd.read_csv("D:/atlas/pla/pla_id.csv")

print(f"斑块信息包含 {len(plaque_info)} 行数据")


# 预处理样本（只执行一次）
print("\n预处理样本中...")
processed_samples = preprocess_samples(file_list, plaque_info)
print(f"成功预处理 {len(processed_samples)} 个样本")
# 设置要分析的基因
#genes_of_interest = ["Grn", "Tnfsf13", "Tnfsf13b", "Gas6","Mki67"]
# 设置输出目录
output_dir = "D:/atlas/pla/轨迹"
os.makedirs(output_dir, exist_ok=True)
genes_of_interest = ["Mki67","S1pr1",'C3', "Cxcl12",'Ccl8', 'Cxcl16', 'Ccl11', 'Ccl12', 'Cxcl10', 'Cxcl9', 'Ccl19', 'Ccl20',
                     'Ccl21a', 'Flt3l', 'Tgfb3', 'Tnf', 'Tnfsf11', 'Lta', 'Pdgfc', 'Igf1', 'Gas6', 'Grn',
                     'Tnfsf13b', 'Tnfsf13b', 'Tnfsf13', 'Il16']
genes_of_interest = ["S1pr1","Cxcl12","Mki67", "Tnfsf13","Tnfsf13b", "Tnfrsf13b", "Tnfrsf17","Sdc1","Ifng",'Cxcl9', "Cxcl10","Tnf","Ccl5"]
# 分析基因（可多次调用，无需重新计算斑块信息）
results = analyze_genes(
    processed_samples=processed_samples,
    genes_of_interest=genes_of_interest,
    max_dist=1200,
    bin_size=50
)
# 可视化结果
# plot_results(results, genes_of_interest, output_dir)

save_results(results, genes_of_interest, output_dir)
# 保存结果

print("分析完成!")

####################################################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.colors import LinearSegmentedColormap

# -------------------------------
# Step 1: 读取标准化后的基因空间梯度结果
# 说明：使用标准化结果画热图，便于比较不同基因沿 DNBS 距离的相对变化
# -------------------------------
output_dir = "D:/atlas/pla/轨迹"
norm_csv_path = os.path.join(output_dir, "normalized_spatial_gradient_results.csv")
data = pd.read_csv(norm_csv_path)

df = data.set_index("distance_bin_center")
df_transposed = df.T

# -------------------------------
# Step 2: 设置 Nature 风格绘图参数
# 说明：使用小字号、可编辑 PDF 字体和克制配色，适合放入论文拼图
# -------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

sns.set_theme(style="white", context="paper")

nature_cmap = LinearSegmentedColormap.from_list(
    "nature_blue",
    ["#F7F7F7", "#D6E5F3", "#6BAED6", "#08519C"]
)

# -------------------------------
# Step 3: 绘制热图
# 说明：去掉格内数字，减少视觉噪音；保留细边框帮助识别距离分箱
# -------------------------------
fig, ax = plt.subplots(figsize=(7.2, 3.2))

hm = sns.heatmap(
    df_transposed,
    cmap=nature_cmap,
    vmin=0,
    vmax=1,
    annot=False,
    linewidths=0.25,
    linecolor="white",
    cbar_kws={
        "label": "Normalized expression",
        "shrink": 0.75,
        "pad": 0.02
    },
    ax=ax
)

ax.set_xlabel("Distance Bin Center (μm)", fontsize=8)
ax.set_ylabel("")
ax.set_title("Spatial gene-expression gradient", fontsize=9, pad=8, weight="bold")

# 只显示部分 x 轴标签，避免距离分箱过密
x_labels = [int(float(x.get_text())) for x in ax.get_xticklabels()]
show_every = 2
ax.set_xticklabels(
    [str(x) if i % show_every == 0 else "" for i, x in enumerate(x_labels)],
    rotation=45,
    ha="right",
    fontsize=7
)

# 基因名用斜体，更接近论文图风格
ax.set_yticklabels(
    ax.get_yticklabels(),
    rotation=0,
    fontsize=7,
    fontstyle="italic"
)

# 优化 colorbar 样式
cbar = hm.collections[0].colorbar
cbar.ax.tick_params(labelsize=7, length=2, width=0.5)
cbar.outline.set_linewidth(0.5)
cbar.set_label("Normalized expression", fontsize=7)

sns.despine(left=True, bottom=True)
plt.tight_layout()

plt.savefig(os.path.join(output_dir, "gene_expression_heatmap_nature.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(output_dir, "gene_expression_heatmap_nature.png"), dpi=600, bbox_inches="tight")
plt.show()
######################################
output_dir = "D:/atlas/pla/轨迹"
csv_path = os.path.join(output_dir, 'spatial_gradient_results.csv')
data=pd.read_csv(csv_path)
data.head(5)
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 创建保存图像的目录
plot_dir = os.path.join(output_dir, 'gene_expression_plots')
os.makedirs(plot_dir, exist_ok=True)

# 获取所有基因名称（去除后缀）
genes = [col.split('_')[0] for col in data.columns if '_mean' in col and col != 'cell_count_mean']
print(f"找到 {len(genes)} 个基因: {genes}")

# 设置图形参数
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.figsize': (12, 6),
    'figure.dpi': 150,
    'font.family': 'DejaVu Sans',  # 支持中文的字体
})

# 为每个基因创建折线图
# 设置支持中文的字体（Windows系统使用Microsoft YaHei）
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'KaiTi', 'Arial Unicode MS', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
for gene in genes:
    # 创建新图形
    fig, ax = plt.subplots()

    # 提取相关数据
    x = data['distance_bin_center']
    y_mean = data[f'{gene}_mean']
    y_sem = data[f'{gene}_sem']

    # 绘制折线图
    ax.plot(x, y_mean, 'o-', linewidth=2.5, markersize=8, color='#1f77b4', label='mean')

    # 添加误差线（SEM）
    ax.fill_between(x,
                    y_mean - y_sem,
                    y_mean + y_sem,
                    alpha=0.25,
                    color='#1f77b4',
                    label='SEM')

    # 添加标题和标签
    ax.set_title(f'{gene} ', pad=20)
    ax.set_xlabel('distance (Pixel)', labelpad=15)
    ax.set_ylabel('expresion', labelpad=15)

    # 添加网格线
    ax.grid(True, linestyle='--', alpha=0.3)

    # 添加图例
    ax.legend(loc='best')

    # 优化x轴刻度
    ax.set_xticks(np.arange(0, max(x) + 1, 200))

    # 添加零线
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

    # 保存图像
    plot_path = os.path.join(plot_dir, f'{gene}_expression.pdf')
    plt.tight_layout()
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()

    print(f"已保存: {plot_path}")

print(f"\n所有基因表达图已保存至: {plot_dir}")
