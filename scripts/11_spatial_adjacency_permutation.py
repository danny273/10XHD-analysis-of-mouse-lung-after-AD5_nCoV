import os
os.environ["OMP_NUM_THREADS"] = "1"  # 必须在所有import之前设置
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np
from scipy.spatial import cKDTree
import pandas as pd
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
import psutil
from threadpoolctl import threadpool_limits
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

# 为热图定义一个红蓝发散的颜色映射 (红表示正相关扎堆，蓝表示负相关排斥)
cmap_custom = sns.diverging_palette(240, 10, as_cmap=True)
# 设置 PDF 导出时保留文字的字符属性 (Type 42 是 TrueType 字体)
matplotlib.rcParams['pdf.fonttype'] = 42
def spatial_enrichment_analysis(df, cluster_col='clu', n_simulations=100,
                                adjust_method='fdr_bh', seed=123, n_jobs=-1,
                                distance_threshold=120):
    """移除进度条的高效版本"""

    # 数据预处理
    df = df.copy()
    df['_encoded_clu'], clu_mapping = pd.factorize(df[cluster_col])
    cluster_labels = df['_encoded_clu'].values

    # 阶段1：构建空间邻域网络
    points = df[['x', 'y']].values
    tree = cKDTree(points)
    edges = list(tree.query_pairs(r=distance_threshold))

    # 统计实际连接数
    actual_counts = defaultdict(int)
    for i, j in edges:
        clu_pair = tuple(sorted([cluster_labels[i], cluster_labels[j]]))
        actual_counts[clu_pair] += 1

    # 阶段2：并行模拟
    def run_simulation(rng_seed):
        """优化后的模拟函数"""
        with threadpool_limits(limits=1):
            rng = np.random.RandomState(rng_seed)
            shuffled_clu = rng.permutation(cluster_labels)
            edges_array = np.array(edges, dtype=np.int64)
            pairs = np.sort(shuffled_clu[edges_array], axis=1)
            unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
            return dict(zip(map(tuple, unique_pairs), counts))

    seed_seq = np.random.SeedSequence(seed)
    child_seeds = seed_seq.generate_state(n_simulations)

    # 简化并行调用
    sim_counts = Parallel(n_jobs=n_jobs, verbose=0, batch_size=10)(  # 增加批处理大小
        delayed(run_simulation)(child_seeds[i])
        for i in range(n_simulations)
    )
    # 阶段3：结果整合
    # 收集所有可能的聚类对
    all_pairs = set(actual_counts.keys())
    for sim in sim_counts:
        all_pairs.update(sim.keys())
    all_pairs = list(all_pairs)

    # 准备结果数据
    results = []
    for pair in all_pairs:
        obs = actual_counts.get(pair, 0)
        sim_vals = [s.get(pair, 0) for s in sim_counts]

        # 计算统计指标
        n_higher = sum((obs + 1) > (s + 1) for s in sim_vals)
        p_high = 1 - (n_higher / n_simulations)
        n_lower = sum((obs + 1) < (s + 1) for s in sim_vals)
        p_low = 1 - (n_lower / n_simulations)
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

    # 多重检验校正
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

    # 最终排序
    results_df = results_df.sort_values('PI', ascending=False).reset_index(drop=True)
    results_df['rank'] = results_df.index + 1
    results_df['cluster_pair'] = results_df['cluster_pair'].apply(
        lambda x: tuple(sorted(clu_mapping.take(list(x))))
    )
    return results_df
##################################################################################################################
# 在空间转录组学和计算生物学中，你提供的这段代码实现了一种基于**置换检验（Permutation Test）和空间邻域网络（Spatial Neighborhood Network）**的统计方法。它通过对比“真实物理分布”与“随机打乱分布”来量化不同细胞类型（Clusters）之间的空间相关性。
# 具体计算过程可以分为以下四个核心步骤：
# 1. 构建空间邻域网络（定义“邻居”）
# 计算的第一步是确定哪些细胞在空间上是“相邻”的。
# 坐标检索：利用 cKDTree（一种高效的空间索引数据结构）对所有细胞的 (x,y)(x,y) 坐标进行索引13。
# 距离阈值：通过 distance_threshold=120 定义一个物理半径。如果细胞 A 和细胞 B 之间的欧几里得距离小于此阈值，算法就认为它们之间存在一条“边（Edge）”。
# 结果：最终生成一个庞大的连接矩阵，记录了组织切片中所有物理位置靠近的细胞对。
# 2. 统计实际连接频率（观测值 ObsObs）
# 在建立了邻域网络后，算法会根据每个细胞的标签（Cluster）统计每种“组合”出现的次数。
# 配对计数：遍历所有边，记录每对 Cluster（如 (T细胞, 浆细胞)）在真实组织中相邻出现的总次数。
# 物理意义：如果 (T细胞, 浆细胞) 的实际连接数非常多，直观上说明这两种细胞在空间上倾向于“扎堆”分布，即存在正相关性48。
# 3. 生成随机背景分布（模拟值 SimSim）
# 这是统计显著性检验的核心。为了判断观测到的“扎堆”是由于生物学逻辑导致的，还是仅仅因为细胞数量多而产生的随机巧合，算法进行了模拟：
# 标签置换（Permutation）：保持细胞的物理位置和邻域网络不动，但随机打乱所有细胞的类型标签（Cluster Labels）17。
# 重复模拟：在 n_simulations=100 次模拟中，每次打乱后都会重新计算一遍每种 Cluster 对的连接数。
# 期望值：这 100 次模拟结果的平均值（sim_mean）代表了在完全随机分布的情况下，两种细胞理论上应该相邻出现的“期望次数”。
# 4. 计算相关性指标（Enrichment & PI）
# 最后，通过对比“现实”与“随机”，得出最终的相关性得分：
# 富集倍数 (Enrichment Score)：
# Enrichment=log2(Observed+1Simulated_Mean+1)Enrichment=log2​(Simulated_Mean+1Observed+1​)
# 如果值 >0>0，说明实际观察到的相邻次数多于随机预期，存在空间正相关（共定位）。
# 如果值 <0<0，说明实际观察到的相邻次数少于随机预期，存在空间负相关（排斥）410。
# 显著性 (P-value)：看实际观测值在 100 次模拟分布中所处的分位数。如果只有极少数（或没有）模拟值能超过观测值，则 PP 值极小，说明结果高度显著8。
# 优先指数 (PI值)：
# PI=−log⁡10(Padj)imesEnrichmentPI=−log10​(Padj)imesEnrichment
# 这是一个综合指标，类似于 Volcano Plot（火山图）的逻辑，同时考虑了相关性的强度（Enrichment）和结果的可靠性（P-value）。PI 值越大，说明这种空间共定位关系在生物学上越重要。
#################################################################################################################################
def process_pairs(results_df):
    """解析聚类对 (cluster_pair) 并生成显著性标记"""
    records = []
    for _, row in results_df.iterrows():
        # 解包聚类对，比如 ('AT1', 'Mac')
        cluster1, cluster2 = row['cluster_pair']
        # 选取最小的 p-value（判断是正相关显著还是负相关显著）
        p_val = min(row['p_adj_high'], row['p_adj_low'])
    # 生成显著性星号
        sig = ""
        if p_val < 0.001:
            sig = "***"
        elif p_val < 0.01:
            sig = "**"
        elif p_val < 0.05:
            sig = "*"
        records.append({
            'source': cluster1,
            'target': cluster2,
            'enrichment': row['enrichment'],
            'sig': sig
        })
    return pd.DataFrame(records)
def generate_heatmap(results, source, save_dir="D:/atlas"):
    """生成带显著性标记的聚类热图"""
    # 预处理结果
    filtered_results = results.iloc[:-1, :]  # 移除最后一行汇总数据
    pair_df = process_pairs(filtered_results)

    # 获取所有细胞类型
    cell_types = sorted(set(pair_df['source']).union(pair_df['target']))

    # 构建富集矩阵
    enrich_matrix = pd.DataFrame(0, index=cell_types, columns=cell_types)
    sig_matrix = pd.DataFrame("", index=cell_types, columns=cell_types)

    for _, row in pair_df.iterrows():
        a, b = row['source'], row['target']
        enrich_matrix.loc[a, b] = row['enrichment']
        sig_matrix.loc[a, b] = row['sig']
        enrich_matrix.loc[b, a] = row['enrichment']
        sig_matrix.loc[b, a] = row['sig']

    # 生成注释文本
    annot_matrix = enrich_matrix.round(1).astype(str) + sig_matrix
    annot_matrix = annot_matrix.replace("0.0", "")

    # 动态调整图像尺寸
    n_types = len(cell_types)
    figsize = (max(12, n_types * 0.8), max(10, n_types * 0.6))

    # 绘制热图
    plt.figure(figsize=figsize)
    sns.set(font_scale=1.2)
    vmax = np.nanmax(np.abs(enrich_matrix.values))

    g = sns.clustermap(
        enrich_matrix,
        cmap=cmap_custom,
        annot=annot_matrix,
        fmt="",
        norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax),
        linewidths=0.5,
        linecolor="lightgray",
        dendrogram_ratio=(0.15, 0.1),
        cbar_pos=(0.95, 0.3, 0.03, 0.4),
        figsize=figsize,
        annot_kws={"fontsize":20, "color": "black"}
    )

    # 调整颜色条
    # cbar = g.ax_heatmap.collections[0].colorbar
    # cbar.set_ticks([-vmax, 0, vmax])
    # cbar.set_ticklabels(["Negative", "Zero", "Positive"])

    # 调整标签样式
    plt.setp(g.ax_heatmap.get_xticklabels(),
             rotation=45, ha='right', fontsize=20)
    plt.setp(g.ax_heatmap.get_yticklabels(),
             rotation=0, fontsize=20)

    # 添加标题并保存
    g.fig.suptitle(
        f"Spatial Interaction - {source}",
        y=1.02,  # 控制垂直位置
        x=0.5,    # 水平居中
        fontsize=30,
        ha='center'  # 水平对齐方式
    )
    save_path = os.path.join(save_dir, f"{source}_interaction_heatmap.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    # 读取完整数据
    df = pd.read_csv('D:/atlas/leiden_26.csv')

    # 获取所有唯一source
    sources = df["source"].unique()

    # 遍历处理每个source
    for source in sources:
        print(f"正在处理: {source}")
        df_source = df[df["source"] == source]

        # 运行空间富集分析
        results = spatial_enrichment_analysis(
            df_source,
            cluster_col='nich4',
            distance_threshold=120,
            n_simulations=400,
            n_jobs=18
        )

        # 保存结果
        result_path = os.path.join("D:/atlas/空间相关性", f"{source}_results.csv")
        results.to_csv(result_path, index=False)

        # 生成热图
        generate_heatmap(results, source, save_dir="D:/atlas/空间相关性")
        print(f"完成处理: {source}")

    print("全部处理完成！")
############################################################################################################################
import os
import pandas as pd

# 配置路径
base_dir = "D:/atlas/空间相关性"
save_path = "D:/atlas/combined_results.csv"

def auto_merge():
    all_dfs = []
    # 遍历目录自动检测结果文件
    for filename in os.listdir(base_dir):
        if filename.endswith("_results.csv"):
            source = filename.replace("_results.csv", "")
            filepath = os.path.join(base_dir, filename)
            try:
                # 读取时指定列类型防止类型不匹配
                df = pd.read_csv(filepath, dtype={
                    'cluster_pair': str,
                    'observed': int,
                    'sim_mean': float,
                    'p_adj_low': float,
                    'PI': float,
                    'rank': int
                })
                # 添加来源标识列
                df['source'] = source
                all_dfs.append(df)
                print(f"已加载: {filename}")
            except Exception as e:
                print(f"加载 {filename} 失败: {str(e)}")
                continue

    # 验证是否找到足够文件
    if len(all_dfs) < 8:
        print(f"警告: 只找到 {len(all_dfs)} 个文件，不足8个")

    # 合并数据
    combined_df = pd.concat(all_dfs, axis=0, ignore_index=True)

    # 数据校验
    print("\n合并后数据概况:")
    print(f"总行数: {len(combined_df)}")
    print(f"唯一source数量: {combined_df['source'].nunique()}")
    print(f"列类型:\n{combined_df.dtypes}")

    # 保存结果
    combined_df.to_csv(save_path, index=False)
    print(f"\n合并结果已保存至: {save_path}")
    return combined_df
# 执行自动合并
final_df = auto_merge()

# 展示样例数据
print("\n合并数据示例:")
print(final_df.head(3))



import numpy as np
from ast import literal_eval


def generate_balt_heatmap(combined_df, save_dir="D:/atlas", target_cluster="BALT"):
    # 数据预处理流程（修改所有BALT为target_cluster参数）
    processed = (
        combined_df
        .assign(cluster_pair=lambda x: x['cluster_pair'].apply(literal_eval))
        .loc[lambda x: x['cluster_pair'].apply(lambda p: p[0] != p[1])]
        .assign(
            sorted_pair=lambda x: x['cluster_pair'].apply(
                lambda p: tuple(sorted(p)) if p[0] < p[1] else tuple(reversed(sorted(p)))
            )
        )
        .drop_duplicates(subset=['source', 'sorted_pair'])
        .loc[lambda x: x['sorted_pair'].apply(lambda p: target_cluster in p)]  # 修改点1
        .assign(
            interaction=lambda x: x['sorted_pair'].apply(
                lambda p: f"{target_cluster}-{p[1]}" if p[0] == target_cluster else f"{p[0]}-{target_cluster}"  # 修改点2
            )
        )
    )

    # 构建矩阵体系（保持原样）
    cell_types = processed['interaction'].unique()
    time_points = processed['source'].unique()

    # 初始化空矩阵
    enrich_matrix = pd.DataFrame(0, index=cell_types, columns=time_points)
    sig_matrix = pd.DataFrame('', index=cell_types, columns=time_points)

    # 填充数据矩阵（保持原样）
    for _, row in processed.iterrows():
        interaction = row['interaction']
        time_point = row['source']
        enrich_matrix.loc[interaction, time_point] = row['enrichment']
        sig = "**" if row['p_adj_high'] < 0.01 else "*" if row['p_adj_high'] < 0.05 else ""
        sig_matrix.loc[interaction, time_point] = sig

    # 生成组合标注（保持原样）
    annot_matrix = enrich_matrix.round(1).astype(str) + sig_matrix
    annot_matrix = annot_matrix.replace("0.0", "")

    # 可视化参数设置（保持原样）
    plt.style.use('seaborn-v0_8-white')
    cmap_custom = LinearSegmentedColormap.from_list(
        'diverging', ['#2b62d9', 'white', '#d93b2b']
    )
    vmax = np.nanmax(np.abs(enrich_matrix.values))

    # 动态尺寸计算（保持原样）
    n_rows = len(cell_types)
    base_height = 10
    row_height = 0.6
    figsize = (12, max(base_height, n_rows * row_height))

    # 绘制高级热图（保持原样）
    ax = sns.heatmap(
        enrich_matrix,
        cmap=cmap_custom,
        annot=annot_matrix,
        fmt="",
        norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax),
        linewidths=0.5,
        linecolor="lightgray",
        cbar_kws={
            "label": "Enrichment Score",
            "ticks": [-vmax, 0, vmax]
        },
        annot_kws={
            "fontsize": 12,
            "color": "black",
            "fontweight": "bold"
        }
    )

    # 调整标签样式（保持原样）
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=12)

    # 添加标题（修改标题动态显示目标cluster）
    plt.title(
        f"{target_cluster} Spatial Interaction Dynamics",  # 修改点3
        y=1.02,
        fontsize=18,
        fontweight='bold'
    )

    # 保存输出（文件名动态化）
    plt.savefig(
        f"{save_dir}/{target_cluster}_interaction_clustermap.png",  # 修改点4
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()
generate_balt_heatmap(final_df, target_cluster="Plasma")

# 使用示例
generate_balt_heatmap(final_df)
combined_df=final_df