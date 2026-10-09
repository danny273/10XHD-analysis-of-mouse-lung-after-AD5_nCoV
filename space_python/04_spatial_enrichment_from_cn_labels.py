#基于Delaunay三角网
import numpy as np
from scipy.spatial import Delaunay
import pandas as pd
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
import scanpy as sc
from joblib import Parallel, delayed
from tqdm import tqdm  # 可选进度条支持
def spatial_enrichment_analysis(df, cluster_col='clu', n_simulations=100, adjust_method='fdr_bh', seed=123,n_jobs=-1):
    # 重置索引确保安全访问
    df = df.reset_index(drop=True)

    # 生成Delaunay三角网（保持不变）
    points = df[['x', 'y']].values
    tri = Delaunay(points)

    # 提取唯一边（保持不变）
    edges = set()
    for simplex in tri.simplices:
        for i in range(3):
            a, b = sorted([simplex[i], simplex[(i + 1) % 3]])
            if a != b:
                edges.add((a, b))
    edges = list(edges)

    # 统计实际连接数（保持不变）
    actual_counts = defaultdict(int)
    for i, j in edges:
        clu_i = df.iloc[i][cluster_col]
        clu_j = df.iloc[j][cluster_col]
        clu_pair = tuple(sorted([clu_i, clu_j]))
        actual_counts[clu_pair] += 1

    # 并行化模拟函数
    def run_simulation(rng_seed):
        """
        单次模拟的独立函数
        - rng_seed: 独立随机种子
        返回: 本次模拟的聚类对计数字典
        """
        # 创建独立的随机状态
        rng = np.random.RandomState(rng_seed)
        shuffled_clu = rng.permutation(df[cluster_col])

        # 统计模拟连接
        sim = defaultdict(int)
        for i, j in edges:
            clu_pair = tuple(sorted([shuffled_clu[i], shuffled_clu[j]]))
            sim[clu_pair] += 1
        return sim

    # 生成随机种子序列（保证可重复性）
    seed_seq = np.random.SeedSequence(seed)
    child_seeds = seed_seq.generate_state(n_simulations)

    # 并行执行模拟
    sim_counts = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(run_simulation)(child_seeds[i])
        for i in tqdm(range(n_simulations), desc='Running simulations')
    )

    # 收集所有可能的聚类对
    all_pairs = set(actual_counts.keys())
    for sim in sim_counts:
        all_pairs.update(sim.keys())
    all_pairs = list(all_pairs)

    # 准备数据表格
    results = []
    for pair in all_pairs:
        # 实际观测值
        obs = actual_counts.get(pair, 0)
        # 模拟值
        sim_vals = [s.get(pair, 0) for s in sim_counts]
        # 计算p值
        n_higher = sum((obs + 1) > (s + 1) for s in sim_vals)
        p_high = 1 - (n_higher / n_simulations)

        n_lower = sum((obs + 1) < (s + 1) for s in sim_vals)
        p_low = 1 - (n_lower / n_simulations)

        # 计算富集分数
        mean_sim = np.mean(sim_vals)
        enrich = np.log2((obs + 1) / (mean_sim + 1))

        results.append({
            'cluster_pair': pair,
            'observed': obs,
            'sim_mean': mean_sim,
            'enrichment': enrich,
            'p_high': p_high,
            'p_low': p_low
        })

    # 转换为DataFrame
    results_df = pd.DataFrame(results)

    # p值校正
    pvals_high = results_df['p_high'].values
    pvals_low = results_df['p_low'].values

    _, p_adj_high, _, _ = multipletests(pvals_high, method=adjust_method)
    _, p_adj_low, _, _ = multipletests(pvals_low, method=adjust_method)

    results_df['p_adj_high'] = p_adj_high
    results_df['p_adj_low'] = p_adj_low

    # 计算PI值
    results_df['PI'] = 0.0
    for idx, row in results_df.iterrows():
        if row['p_adj_high'] <= row['p_adj_low']:
            pi = -np.log10(row['p_adj_high'] + 1 / n_simulations) * row['enrichment']
        else:
            pi = -np.log10(row['p_adj_low'] + 1 / n_simulations) * row['enrichment']
        results_df.at[idx, 'PI'] = pi

    # 按PI值排序
    results_df = results_df.sort_values('PI', ascending=False).reset_index(drop=True)
    results_df['rank'] = results_df.index + 1

    return results_df


# 使用示例

# 示例数据 (需替换为真实数据)
leiden = pd.read_csv("D:/space/d211/leiden_labels_with_index.csv", index_col=0)
adata = sc.read_h5ad("D:/atlas/d211.h5ad")

# 过滤特定数据源
leiden = leiden[leiden["source"] == 4]
adata.obs["leiden"] = leiden["merged_cluster"].values.astype(str)  # 确保分类为字符串
data = {
    'x':adata.obs["x"],
    'y':adata.obs["y"],
    'leiden':adata.obs["leiden"]
}
df = pd.DataFrame(data)

# 运行分析
results = spatial_enrichment_analysis(df,cluster_col='leiden')
results.to_csv("D:/atlas/d211_clu.csv", index=False)
# 查看结果
print(results.head())
##############################################################################################################################
import seaborn as sns
import matplotlib.pyplot as plt


def plot_enrichment_heatmap(results_df,
                            value_col='enrichment',
                            pval_threshold=0.05,
                            figsize=(20, 20),
                            cmap='coolwarm',
                            annot_fontsize=8):
    """
    修复版本的热图绘制函数，支持层次聚类
    """
    # 生成聚类列表
    clusters = sorted(set([c for pair in results_df['cluster_pair'] for c in pair]))

    # 构建矩阵
    value_matrix = pd.DataFrame(np.zeros((len(clusters), len(clusters))),
                                index=clusters, columns=clusters)
    sig_matrix = pd.DataFrame(np.full((len(clusters), len(clusters)), ''),
                              index=clusters, columns=clusters)

    # 填充数据
    for _, row in results_df.iterrows():
        c1, c2 = row['cluster_pair']
        val = row[value_col]
        p_high = row['p_adj_high']
        p_low = row['p_adj_low']

        value_matrix.loc[c1, c2] = val
        value_matrix.loc[c2, c1] = val

        # 显著性标记
        sig = []
        if p_high <= pval_threshold:
            sig.append('*' * (1 + int(p_high < 0.01)))
        if p_low <= pval_threshold:
            sig.append('*' * (1 + int(p_low < 0.01)))
        sig = '\n'.join(sig)

        sig_matrix.loc[c1, c2] = sig
        sig_matrix.loc[c2, c1] = sig

    # 组合标注
    annot_matrix = value_matrix.round(2).astype(str) + '\n' + sig_matrix

    # 绘制带聚类的热图
    g = sns.clustermap(
        value_matrix,
        cmap=cmap,
        annot=annot_matrix.values,
        fmt='',
        center=0,
        dendrogram_ratio=0.1,  # 调整树状图比例
        figsize=figsize,
        annot_kws={'size': annot_fontsize, 'color': 'black'},
        cbar_kws={'label': value_col},
        linewidths=0.5,
        tree_kws={'linewidths': 1}  # 树状图线宽
    )

    # 设置标题和标签（修正部分）
    plt.subplots_adjust(top=0.9)  # 为标题留出空间
    g.fig.suptitle(f'Spatial Enrichment Heatmap ({value_col})', y=0.95)

    # 设置轴标签（通过热图坐标轴）
    g.ax_heatmap.set_xlabel('Cluster', fontsize=12)
    g.ax_heatmap.set_ylabel('Cluster', fontsize=12)

    # 调整颜色条标签
    g.cax.set_ylabel(value_col, rotation=-90, va='bottom')

    return g
# 基础调用
plot_enrichment_heatmap(results)
plt.show()
plt.savefig('D:/atlas/d211_cluheat.png')
##imm#########################################################################################################
def plot_clustered_subset_heatmap(results_df,
                                  subset_clusters,
                                  value_col='enrichment',
                                  pval_threshold=0.05,
                                  figsize=(15, 13),
                                  cmap='coolwarm',
                                  annot_fontsize=9):
    """
    绘制指定niche子集的聚类热图

    参数:
    results_df: 原始结果DataFrame
    subset_clusters: 目标niche列表，如['niche16','niche12',...]
    value_col: 显示值列名（enrichment/PI）
    pval_threshold: 显著性阈值
    figsize: 图像尺寸
    cmap: 颜色映射
    annot_fontsize: 标注字体大小
    """
    # 数据过滤
    target_pairs = [(c1, c2) for c1 in subset_clusters
                    for c2 in subset_clusters if c1 <= c2]
    filtered_df = results_df[results_df['cluster_pair'].isin(target_pairs)]

    # 构建子矩阵
    clusters = sorted(subset_clusters)
    n = len(clusters)

    # 初始化矩阵
    value_matrix = pd.DataFrame(np.nan, index=clusters, columns=clusters)
    sig_matrix = pd.DataFrame('', index=clusters, columns=clusters)

    # 填充数据
    for _, row in filtered_df.iterrows():
        c1, c2 = row['cluster_pair']
        val = row[value_col]
        p_high = row['p_adj_high']
        p_low = row['p_adj_low']

        # 对称填充
        value_matrix.loc[c1, c2] = val
        value_matrix.loc[c2, c1] = val

        # 显著性标记
        sig = []
        if p_high <= pval_threshold:
            sig.append(f'↑{"*" * (1 + int(p_high < 0.01))}')
        if p_low <= pval_threshold:
            sig.append(f'↓{"*" * (1 + int(p_low < 0.01))}')
        sig = '\n'.join(sig)

        sig_matrix.loc[c1, c2] = sig
        sig_matrix.loc[c2, c1] = sig

    # 填充缺失值为0（无显著交互）
    value_matrix = value_matrix.fillna(0)

    # 组合标注
    annot_matrix = value_matrix.round(2).astype(str) + '\n' + sig_matrix

    # 绘制聚类热图
    g = sns.clustermap(
        value_matrix,
        cmap=cmap,
        annot=annot_matrix.values,
        fmt='',
        center=0 if value_col == 'enrichment' else 1,
        dendrogram_ratio=0.15,
        figsize=figsize,
        annot_kws={'size': annot_fontsize, 'color': 'black'},
        cbar_kws={'label': value_col},
        linewidths=0.5,
        tree_kws={'linewidths': 1.5}
    )

    # 设置标注
    g.ax_heatmap.set_xlabel('Cluster', fontsize=12)
    g.ax_heatmap.set_ylabel('Cluster', fontsize=12)
    g.fig.suptitle(f'Clustered Subset Heatmap ({value_col})\nNiches: {subset_clusters}',
                   y=0.98, fontsize=14)

    # 调整颜色条
    g.cax.set_ylabel(value_col, rotation=-90, va='bottom')

    # 优化坐标标签
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    return g
# 绘制子集热图
target_niches = ['niche16','niche12','niche24','niche27',
                'niche26','niche33','niche37','niche36','niche29']

g = plot_clustered_subset_heatmap(
    results,
    subset_clusters=target_niches,
    value_col='PI',
    cmap='RdBu_r'
)
plt.show()
# 保存高清图片
plt.savefig('D:/atlas/d211_imm-cluheat.png', dpi=300, bbox_inches='tight')
