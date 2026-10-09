import glob
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import numpy as np
import scipy.stats as stats

def build_neighbor_df(adata,type="type",radius=60):
    """动态生成包含当前样本所有细胞类型的邻域统计"""
    # 获取坐标和细胞类型
    cell_coords = adata.obs[['x', 'y']].values
    adata.obs[type] = adata.obs[type].astype(str)
    cell_types = adata.obs[type].values

    # 构建KDTree进行半径近邻搜索
    tree = cKDTree(cell_coords)
    neighbor_indices = tree.query_ball_point(cell_coords, r=radius)

    # 生成当前样本的类型映射
    unique_types = np.unique(cell_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    cell_type_indices = np.array([type_to_idx[t] for t in cell_types])

    # 统计邻域类型数量
    n_cells = cell_coords.shape[0]
    n_types = len(unique_types)
    neighbor_counts = np.zeros((n_cells, n_types), dtype=int)

    for i in range(n_cells):
        neighbors = [n for n in neighbor_indices[i] if n != i]
        counts = np.bincount(cell_type_indices[neighbors], minlength=n_types)
        neighbor_counts[i] = counts

    # 创建带类型前缀的列
    return pd.DataFrame(
        neighbor_counts,
        columns=[f'n_{t}' for t in unique_types],
        index=adata.obs.index
    )


# 获取所有h5ad文件
file_paths = glob.glob('D:/atlas/*.h5ad')[:8]
all_dfs = []
leiden = pd.read_csv("D:/space/d211/leiden_labels_with_index.csv", index_col=0)
for file_idx, path in enumerate(file_paths):
    print(f"Processing {file_idx + 1}/{len(file_paths)}: {path}")
    adata = sc.read_h5ad(path)
    # 提取当前样本对应的 leiden_cluster 并确保转换为字符串类型
    current_leiden = leiden.loc[
        leiden['source'] == file_idx, 'merged_cluster'
    ].astype(str).fillna('unknown')

    adata.obs['leiden'] = current_leiden.values
    # 保存修正后的数据（注意修正文件名）
    adata.write_h5ad(f'leiden_{file_idx}.h5ad')  # 根据需求调整文件名
    # 生成当前样本的统计
    df = build_neighbor_df(adata,type="leiden",radius=60)
    # 添加元数据
    df["leiden"] = adata.obs["leiden"].astype(str)
    df['source'] = file_idx
    all_dfs.append(df)

# 智能合并：自动对齐列，缺失列填充0
combined_neighbor = pd.concat(all_dfs, axis=0).fillna(0)

# 类型列排序（将n_*列排前，元数据排后）
n_columns = [c for c in combined_neighbor if c.startswith('n_niche')]
combined_neighbor [n_columns] = combined_neighbor [n_columns].astype('int32')
grouped = combined_neighbor.groupby('leiden')[n_columns].mean()
# 首次运行 clustermap 以获取行聚类顺序
g = sns.clustermap(
    grouped,
    row_cluster=True,
    col_cluster=False,  # 关闭列聚类
    standard_scale=1,
    cmap='viridis'
)
plt.close()  # 关闭临时图形

# 获取行聚类后的顺序（按树状图重新排列的索引）
row_order = g.dendrogram_row.reordered_ind
leiden_order = grouped.index[row_order].tolist()  # 例如 ['1', '0', '2']
# 生成列顺序：根据行聚类后的 leiden 顺序匹配 nich_ 列
column_order = [f'n_{leiden}' for leiden in leiden_order]

# 检查列是否存在，避免 KeyError
valid_columns = [col for col in column_order if col in grouped.columns]
grouped_sorted = grouped[valid_columns]
plt.figure(figsize=(20, 20))
sns.clustermap(
    grouped_sorted,
    cmap='viridis',
    row_cluster=True,  # 禁用行重新聚类（已手动排序）
    col_cluster=False,  # 列不聚类
    standard_scale=1,   # 行标准化
    dendrogram_ratio=0.1,
    figsize=(20, 20)
)
plt.suptitle("Clustered Heatmap with Aligned Row-Column Order", y=1.02)
plt.show()
