import glob
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import numpy as np
import scipy.stats as stats
##
import os

# 获取当前工作目录
current_dir = os.getcwd()
print("当前工作目录：", current_dir)

# 设置新的工作目录
new_dir = "D:/atlas"
os.chdir(new_dir)

# 再次获取当前工作目录
updated_dir = os.getcwd()
print("更新后的工作目录：", updated_dir)
def build_neighbor_df(adata, radius=60):
    """动态生成包含当前样本所有细胞类型的邻域统计"""
    # 获取坐标和细胞类型
    cell_coords = adata.obs[['x', 'y']].values
    adata.obs['type'] = adata.obs['type'].astype(str)
    cell_types = adata.obs['type'].values

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

for file_idx, path in enumerate(file_paths):
    print(f"Processing {file_idx + 1}/{len(file_paths)}: {path}")
    adata = sc.read_h5ad(path)

    # 生成当前样本的统计
    df = build_neighbor_df(adata,radius=60)

    # 添加元数据
    df['type'] = adata.obs['type'].astype(str)
    df['source'] = file_idx

    all_dfs.append(df)

# 智能合并：自动对齐列，缺失列填充0
combined_neighbor = pd.concat(all_dfs, axis=0).fillna(0)

# 类型列排序（将n_*列排前，元数据排后）
n_columns = [c for c in combined_neighbor if c.startswith('n_')]
meta_columns = ['type', 'source']
combined_neighbor = combined_neighbor[n_columns + meta_columns]
combined_neighbor [n_columns] = combined_neighbor [n_columns].astype('int32')
output_dir="combined_plots60_30_2_0.9_cosin"
# 保存结果
from pathlib import Path
Path(output_dir).mkdir(exist_ok=True)
combined_neighbor.to_csv(f"{output_dir}/combined_neighbor.csv")
column_sums = combined_neighbor.iloc[:, :50].sum(axis=1)
print(column_sums.mean(),column_sums.median())
###################################################################################################
# 移除非特征列并标准化
import hnswlib
import numpy as np
import pandas as pd
import leidenalg
import igraph as ig
from scipy.sparse import csr_matrix
import os
import time
from sklearn.preprocessing import StandardScaler
import scanpy as sc
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import scipy.stats as stats
import json
###################
#combined_neighbor.to_csv(f"{output_dir}/combined_neighbor.csv")
X=pd.read_csv(f'{output_dir}/combined_neighbor.csv',index_col=0)
X=combined_neighbor
source=X["source"].values
type=X["type"].values
X=X.iloc[:,:-2]
k=30
ef_construction=1500
M=48
num_threads=30
# 阶段2：HNSW索引构建
print(f"[{time.ctime()}] 阶段2：构建HNSW索引 (数据量: {X.shape})")
dim = X.shape[1]
index = hnswlib.Index(space='cosine', dim=dim)  # 使用L2距离
index.init_index(max_elements=X.shape[0], ef_construction=ef_construction, M=M)
index.add_items(X, num_threads=num_threads)
index.set_ef(max(k * 4, 1200))  # 动态调整ef参数

# 阶段3：KNN查询
print(f"[{time.ctime()}] 阶段3：KNN查询 (k={k})")
labels, distances = index.knn_query(X, k=k, num_threads=num_threads)

# 阶段4：构建稀疏邻接矩阵
print(f"[{time.ctime()}] 阶段4：构建邻接矩阵")
n_samples = X.shape[0]
rows = np.repeat(np.arange(n_samples), k)
cols = labels.flatten()
data = np.ones(n_samples * k)

# 创建对称邻接矩阵
adj_csr = csr_matrix((data, (rows, cols)), shape=(n_samples, n_samples))
adj_csr = adj_csr.maximum(adj_csr.T)  # 使矩阵对称

# 阶段5：Leiden聚类（内存优化版）
print(f"[{time.ctime()}] 阶段5：Leiden聚类")

# 转换为igraph图对象
sources, targets = adj_csr.nonzero()
edgelist = list(zip(sources.tolist(), targets.tolist()))
# 分块构建图（内存优化）
graph = ig.Graph(n=n_samples, edges=edgelist)
# 配置优化参数
optimiser = leidenalg.Optimiser()
# partition = leidenalg.CPMVertexPartition(
#     graph,
#     resolution_parameter=0.1,  # 根据目标聚类数量调整
#     n_iterations=15
# )
partition = leidenalg.find_partition(
    graph,
    leidenalg.RBConfigurationVertexPartition,
    resolution_parameter=2,  # 可调整聚类粒度
    n_iterations=15,
    seed=123
)

leiden_labels=np.array(partition.membership)
# 转换为带索引的 DataFrame
leiden_df = pd.DataFrame(
    leiden_labels,
    index=X.index,  # 确保索引与原始数据一致
    columns=["leiden_cluster"]
)
del(graph,edgelist,adj_csr,index)
# 添加来源信息（可选）
#combined_neighbor=pd.read_csv(combined_path,index_col=0)
leiden_df["source"] = combined_neighbor["source"]
df=leiden_df
df["type"]=type
df

# 保存结果
leiden_df.to_csv(f'{output_dir}/leiden_labels_with_index.csv')
#df=pd.read_csv("leiden_labels_with_index.csv")


import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
def calculate_type_ratio(df, cluster_col='leiden_cluster', type_col='type'):
    """计算每个cluster的细胞类型比例"""
    # 创建交叉表
    ct_matrix = pd.crosstab(df[cluster_col], df[type_col])
    # 计算比例并添加伪计数
    ct_ratio = ct_matrix.div(ct_matrix.sum(axis=1), axis=0) + 1e-6
    return ct_ratio

def get_similarity_matrix(ratio_matrix, threshold=0.95):
    """生成相似性矩阵并筛选候选对"""
    sim_matrix = cosine_similarity(ratio_matrix)
    np.fill_diagonal(sim_matrix, 0)  # 忽略自相似
    candidates = np.argwhere(sim_matrix >= threshold)
    return sim_matrix, candidates
def merge_clusters(df, ratio_matrix, candidates, cluster_col='leiden_cluster'):
    """基于候选对进行迭代合并，并按大小重命名为niche1、niche2..."""
    cluster_sizes = df[cluster_col].value_counts()
    merge_map = {}
    # 按相似对排序（先合并大群）
    sorted_pairs = sorted(
        [(i, j) for i, j in candidates],
        key=lambda x: cluster_sizes.loc[x[0]] + cluster_sizes.loc[x[1]],
        reverse=True
    )
    # 合并处理
    for clus1, clus2 in sorted_pairs:
        if clus1 in merge_map or clus2 in merge_map:
            continue

        # 保留更大的cluster编号
        if cluster_sizes[clus1] >= cluster_sizes[clus2]:
            merge_map[clus2] = clus1
        else:
            merge_map[clus1] = clus2
    # 生成最终映射
    final_map = {c: c for c in ratio_matrix.index}
    final_map.update(merge_map)

    # 应用映射
    df['merged_cluster'] = df[cluster_col].map(final_map)

    # 按合并后的簇大小排序并重命名
    cluster_order = df['merged_cluster'].value_counts().sort_values(ascending=False).index
    niche_mapping = {cluster: f'niche{i + 1}' for i, cluster in enumerate(cluster_order)}
    df['merged_cluster'] = df['merged_cluster'].map(niche_mapping)
    return df

type_ratio = calculate_type_ratio(df)
sim_matrix, merge_candidates = get_similarity_matrix(type_ratio,threshold=0.9)
merged_df = merge_clusters(df, type_ratio, merge_candidates)
print(f'原始群数量: {len(np.unique(merged_df["leiden_cluster"]))}')
print(f'最终群数量: {len(np.unique(merged_df["merged_cluster"]))}')
merged_df.to_csv(f'{output_dir}/leiden_labels_with_index.csv')
#######################################################################################################
# 绘制特定样本的细胞类型分布

subset = merged_df[merged_df['source'] == 6]
print(f'最终群数量: {len(np.unique(subset["merged_cluster"]))}')
print("每个聚类的样本数:", np.bincount(subset["merged_cluster"]))
subset=combined_neighbor[combined_neighbor['source'] == 6]
subset["n_GCB"]
subset["n_IM1"]
