
########################################################################################################################
import glob
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import numpy as np
import scipy.stats as stats


def build_neighbor_df(adata, radius=60):
    #     """构建基于最短距离的邻域统计DataFrame"""
    #     # 获取坐标和细胞类型
    cell_coords = adata.obs[['x', 'y']].values
    adata.obs['type'] = adata.obs['type'].astype(str)
    cell_types = adata.obs['type'].values

    # 构建KDTree进行半径近邻搜索
    tree = cKDTree(cell_coords)
    neighbor_indices = tree.query_ball_point(cell_coords, r=radius)  # 获取每个点的半径内邻居索引（包括自身）

    # 将细胞类型转换为数值索引
    unique_types = np.unique(cell_types)
    type_to_idx = {t: i for i, t in enumerate(unique_types)}
    cell_type_indices = np.array([type_to_idx[t] for t in cell_types])

    # 统计邻域类型数量
    n_cells = cell_coords.shape[0]
    n_types = len(unique_types)
    neighbor_counts = np.zeros((n_cells, n_types), dtype=int)

    for i in range(n_cells):
        # 排除自身，仅统计周围邻居
        neighbors = [n for n in neighbor_indices[i] if n != i]
        counts = np.bincount(cell_type_indices[neighbors], minlength=n_types)
        neighbor_counts[i] = counts

    # 创建DataFrame
    neighbor_df = pd.DataFrame(
        neighbor_counts,
        columns=[f'n_{t}' for t in unique_types],
        index=adata.obs.index
    )
    return neighbor_df

# 获取所有h5ad文件（根据实际情况修改路径）
path = 'D:/atlas\\d451.h5ad'  # 假设处理前8个文件

adata = sc.read_h5ad(path)

# 构建邻域统计
neighbor_df = build_neighbor_df(adata)

# 添加元信息
neighbor_df['type'] = adata.obs['type']  # 细胞类型
neighbor_df['source'] = 4  # 来源文件标识（从0开始）

# 绘制并保存热图
adata.obs = pd.concat([adata.obs, neighbor_df], axis=1)
# 合并所有数据
combined_neighbor =neighbor_df
# 保存合并的邻域数据
combined_neighbor.to_csv('combined_neighbor.csv')
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

def merge_clusters_by_similarity(similarity_matrix, cluster_labels, threshold=0.95):
    assert similarity_matrix.shape[0] == similarity_matrix.shape[1], "需要方阵"
    n_clusters = similarity_matrix.shape[0]

    # 初始化并查集结构
    parent = list(range(n_clusters))

    def find(u):
        """带路径压缩的查找"""
        while parent[u] != u:
            parent[u] = parent[parent[u]]  # 路径压缩
            u = parent[u]
        return u

    def union(u, v):
        """按秩合并"""
        root_u = find(u)
        root_v = find(v)
        if root_u != root_v:
            # 总是将小的聚类合并到大的聚类（可选策略）
            parent[root_v] = root_u

    # 转换为成对比较格式
    np.fill_diagonal(similarity_matrix, 0)  # 忽略自相似性
    triu_indices = np.triu_indices_from(similarity_matrix, k=1)
    similarities = similarity_matrix[triu_indices]
    pairs = list(zip(triu_indices[0], triu_indices[1], similarities))

    # 按相似度降序处理
    sorted_pairs = sorted(pairs, key=lambda x: -x[2])

    # 合并高相似性群
    for i, j, sim in sorted_pairs:
        if sim >= threshold:
            union(i, j)
        # 生成映射关系
    cluster_mapping = {}
    for orig_cluster in range(n_clusters):
        root = find(orig_cluster)
        if root not in cluster_mapping:
            cluster_mapping[root] = root
        cluster_mapping[orig_cluster] = cluster_mapping[root]

    # 新增部分：按细胞数量命名 --------------------------------------------------
    # 将原始标签转换为合并后的数字标签
    temp_numeric_labels = cluster_labels.astype(int).map(cluster_mapping)

    # 统计各合并群的细胞数量
    cluster_counts = temp_numeric_labels.value_counts().sort_values(ascending=False)

    # 生成niche命名映射（niche1=最大群）
    niche_name_map = {
        cluster_id: f"niche{i + 1}"
        for i, (cluster_id, _) in enumerate(cluster_counts.items())
    }

    # 创建最终标签映射
    label_converter = {
        str(orig): niche_name_map[merged]
        for orig, merged in cluster_mapping.items()
    }

    merged_labels = cluster_labels.astype(str).map(label_converter)

    return merged_labels, cluster_mapping
def plot_spatial_clusters(adata, cluster_col, palette='tab20', size=3):
    """绘制空间聚类图"""
    plt.figure(figsize=(10, 10))
    sns.scatterplot(
        data=adata.obs,
        x='x',
        y='y',
        hue=cluster_col,
        palette=palette,
        s=size,
        linewidth=0,
        alpha=0.8
    )
    plt.title(f'Spatial Clusters ({cluster_col})')
    plt.axis('off')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    return plt
###################

# print(adata)
#
# neighbor_df = build_neighbor_df(adata)
# adata.obs = pd.concat([adata.obs, neighbor_df], axis=1)
# #X = StandardScaler().fit_transform(neighbor_df)  # 标准化数据
# X=neighbor_df
#combined_path='combined_neighbor.csv'
#X=pd.read_csv(combined_path,index_col=0)
X=combined_neighbor
source=X["source"].values
type=X["type"].values
X=X.iloc[:,:-2]
k=20
ef_construction=1000
M=48
num_threads=18
target_clusters=100
# 阶段2：HNSW索引构建
print(f"[{time.ctime()}] 阶段2：构建HNSW索引 (数据量: {X.shape})")
dim = X.shape[1]
index = hnswlib.Index(space='l2', dim=dim)  # 使用L2距离
index.init_index(max_elements=X.shape[0], ef_construction=ef_construction, M=M)
index.add_items(X, num_threads=num_threads)
index.set_ef(max(k * 4, 100))  # 动态调整ef参数

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
    resolution_parameter=3,  # 可调整聚类粒度
    n_iterations=15,
    seed=123
)
del(graph,edgelist,adj_csr)
leiden_labels=np.array(partition.membership)
# 转换为带索引的 DataFrame
leiden_df = pd.DataFrame(
    leiden_labels,
    index=X.index,  # 确保索引与原始数据一致
    columns=["leiden_cluster"]
)

# 添加来源信息（可选）
#combined_neighbor=pd.read_csv(combined_path,index_col=0)
leiden_df["source"] = combined_neighbor["source"]

# 保存结果
leiden_df.to_csv("leiden_labels_with_index.csv")
df=pd.read_csv("leiden_labels_with_index.csv")
df=leiden_df
df["type"]=type
df

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
    """基于候选对进行迭代合并"""
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
    return df
type_ratio = calculate_type_ratio(df)
sim_matrix, merge_candidates = get_similarity_matrix(type_ratio)
merged_df = merge_clusters(df, type_ratio, merge_candidates)
print(f'原始群数量: {len(np.unique(merged_df["leiden_cluster"]))}')
print(f'最终群数量: {len(np.unique(merged_df["merged_cluster"]))}')
#######################################################################################################
# 绘制特定样本的细胞类型分布

subset = merged_df[merged_df['source'] == 6]
print(f'最终群数量: {len(np.unique(subset["merged_cluster"]))}')
print("每个聚类的样本数:", np.bincount(subset["merged_cluster"]))
def plot_with_merged_clusters(h5ad_path, merged_df, source_id,output_dir="combined_plots",
                              cluster_col='merged_cluster', source_col='source', size=1):
    """
    一站式处理：添加聚类信息 + 生成可视化
    """
    # 读取原始数据
    adata = sc.read_h5ad(h5ad_path)
    #dataset_id = Path(h5ad_path).stem

    # 从merged_df提取当前数据集的聚类结果
    #source_id = int(dataset_id[1:3]) - 1  # 假设文件名格式为d01_sample.h5ad
    cluster_subset = merged_df
   # cluster_subset = cluster_subset.loc[adata.obs.index]  # 强制按原始数据索引排序

    missing_rate = cluster_subset[cluster_col].isna().mean()
    if missing_rate > 0:
        print(f"警告：{source_id}存在{missing_rate:.1%}的缺失聚类标签")


    # 添加聚类信息
    adata.obs[cluster_col] = cluster_subset[cluster_col].values

    # ------------------
    # 可视化优化
    # ------------------
    df_plot = adata.obs[['x', 'y', cluster_col]].copy()
    df_plot[cluster_col] = df_plot[cluster_col].astype('category')  # 强制转换为分类类型
    # 绘制散点图
    plt.figure(figsize=(20, 20))
    sns.scatterplot(
        data=df_plot,
        x='x',
        y='y',
        hue=cluster_col,
        palette=global_palette,  # 使用全局调色板
        hue_order=sorted(global_palette.keys()),  # 固定颜色顺序
        s=3,
        linewidth=0
    )

    plt.title('Spatial Distribution by Leiden Cluster')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=3)  # 图例放右侧
    plt.tight_layout()

    # 保存空间分布图
    Path(output_dir).mkdir(exist_ok=True)
    plt.savefig(f"{output_dir}/{dataset_id}_spatial.png", dpi=300, bbox_inches='tight')
    plt.close()

    # ------------------
    # 返回聚类比例数据用于后续整合
    # ------------------
    cluster_counts = adata.obs[cluster_col].value_counts(normalize=True)
    return cluster_counts.rename(dataset_id)
counts = plot_with_merged_clusters(h5ad_path=path, merged_df=merged_df, source_id=4, size=0.5)