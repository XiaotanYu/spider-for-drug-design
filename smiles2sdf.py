

# 将smiles转成图片格式的化学结构式
import os
import cv2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import rdBase, Chem, DataStructs
from rdkit.Chem import Draw
from sklearn.cluster import AgglomerativeClustering
from rdkit.Chem import AllChem

try:
    os.mkdir("sdf")
except FileExistsError:
    pass


# 读取excel
file = "./CDK/CDK.xlsx"
df = pd.read_excel(file)
columns = df.columns.tolist()
df.fillna(method='ffill', inplace=True)
print(df.shape)
print(columns)

# 去重
df2 = df.drop_duplicates(subset=["Smiles"], keep="first")
print(df2.shape)

# 重置索引
df2.index = range(len(df2))
smiles = df2['Smiles'].tolist()
pIC50 = df2['pIC50'].tolist()

df = pd.DataFrame(zip(smiles, pIC50))
df.to_excel("./sdf/CDK_nodup.xlsx", index=False)


dic_smiles = dict(zip(range(len(smiles)),smiles))
len(dic_smiles)
dic_pIC50 = dict(zip(range(len(pIC50)),pIC50))
len(dic_pIC50)


mol_list = [ Chem.MolFromSmiles(i) for i in smiles ]
select_mol = [m for m in mol_list if m is not None]
len(select_mol)


morgan_fp = [AllChem.GetMorganFingerprintAsBitVect(x,2,2048) for x in select_mol]
dis_matrix = [DataStructs.BulkTanimotoSimilarity(morgan_fp[i], morgan_fp,returnDistance=True) for i in range(len(morgan_fp))]
dis_array = np.array(dis_matrix)

ward = AgglomerativeClustering(n_clusters=6)
ward.fit(dis_array)
pd.value_counts(ward.labels_)

ward_library = {i: [] for i in range(10)}
for n,j in enumerate(ward.labels_):
    ward_library[j].append(select_mol[n])
selected_compounds = [np.random.choice(ward_library[i]) for i in range(6)]
Draw.MolsToGridImage(
    selected_compounds,
    molsPerRow=3,
    subImgSize=(300, 300),
    legends=["Group"+ ': ' + str(i) for i in range(6)]
)

from scipy.cluster import hierarchy
linked_array = hierarchy.ward(dis_array)
hierarchy.dendrogram(linked_array)
ax = plt.gca()
bounds = ax.get_xbound()
ax.plot(bounds, [34.5,34.5], '--', c='gray')
ax.plot(bounds, [55,55], '--', c='gray')
ax.text(bounds[1], 55, ' 3 clusters', va='center')
ax.text(bounds[1], 34.5, ' 10 clusters')
plt.xlabel('Compounds', fontsize=12)
plt.xticks([])
plt.ylabel('Cluster distance', fontsize=12)
plt.title('Dendrogram for Ward method', fontsize=14)


from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca.fit(dis_array)
dis_pca = pca.transform(dis_array)
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)
ax.scatter(dis_pca[:,0], dis_pca[:,1], s=50, c=ward.labels_, cmap='Paired', alpha=0.5)
ax.set_xlabel('PCA 1', fontsize=20)
ax.set_xticklabels([])
ax.set_ylabel('PCA 2', fontsize=20)
ax.set_yticklabels([])
ax.set_title('10-Clusters by sklearn ward', fontsize=28)

pca_all = PCA()
pca_all.fit(dis_array)
ev_ratio = np.hstack([0, np.cumsum(pca_all.explained_variance_ratio_)])
plt.plot(ev_ratio[:20])
plt.xlim([0,20])
plt.xticks(range(0,20,2))
plt.ylim([0,1.0])
plt.xlabel('n_components')
plt.ylabel('cumulative explained variance')


f = Chem.SDWriter('./sdf/CDK.sdf')
for mol in select_mol:
    f.write(mol)
f.close()





