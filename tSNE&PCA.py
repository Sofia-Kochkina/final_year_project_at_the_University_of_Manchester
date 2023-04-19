import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from sklearn.decomposition import PCA

data = pd.read_csv('normalized_counts_after_DESeq2_for_visualisation_2_excluded.csv', sep=",", index_col=[0])
data.columns = ["high", "high", "high", "high", "high", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "medium", "low", "low", "low", "low", "low"]

pca_out = PCA().fit(data)                               # do PCA and standardize
scree = pca_out.explained_variance_ratio_               # how much of variability is explained by each PC
loadings = pca_out.components_                          # reduced dimension matrix: coordinates of each point in PCA plot
num_pc = pca_out.n_features_in_                         # number of datapoints in the training data = 28
pc_list = ["PC"+str(i) for i in range(1, num_pc+1)]
loadings_data = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
loadings_data['variable'] = data.columns.values
loadings_data = loadings_data.set_index('variable')
data_tr = data.T
pca = PCA().fit(data_tr)
fracs = pca.components_                                 # factor loadings matrix: contribution of each gene to each PC
num = pca.n_features_in_
p_list = ["PC"+str(i) for i in range(1, num+1)]
fracs_data = pd.DataFrame.from_dict(dict(zip(p_list, fracs)))
fracs_data['variable'] = data_tr.columns.values
fracs_data = fracs_data.set_index('variable')

fracs_sort = fracs_data.reindex(fracs_data.PC1.abs().sort_values().index)   # sort genes by absolute value contribution to PC1
fracs_sort2 = fracs_data.reindex(fracs_data.PC2.abs().sort_values().index)  # sort genes by absolute value contribution to PC1

features_dpl = np.append(fracs_sort.index[-11:], fracs_sort2.index[-18:])   # merge genes contributing to PC1 and PC2

features = [*set(features_dpl)]                 # exclude duplicates

description = pd.read_csv('gene_description.csv', sep=";")      #get gene name instead of identifiers: keep because names not always available
dict={}
for i in description.values:
    dict[i[0]]=i[1]                                            # name of each gene

XPC = 'PC2'                                                  # Principal component shown on x-axis
YPC = 'PC1'                                                  # Principal component shown on y-axis
fig = px.scatter(loadings_data, x=XPC, y=YPC, color=data.columns, hover_data=[XPC, YPC])
fig.update_traces(mode='markers', marker_line_width=1, marker_size=8)



for i, feature in enumerate(features):                        # Add loadings
    fig.add_trace(go.Scatter(
        x=[0, fracs_data.PC2[features[i]]],
        y=[0, fracs_data.PC1[features[i]]],
        mode="lines+text",
        name=dict.get(feature)))

fig.update_layout(
    xaxis_title="PC2 ("+str(round(scree[1]*100, 1))+"%)",      # axes titles
    yaxis_title="PC1 ("+str(round(scree[0]*100, 1))+"%)",
    colorway=['#FD3216', '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF', '#F6F926', '#FF9616', '#479B55', '#EEA6FB', '#DC587D', '#D626FF', '#00B5F7', '#B68E00', '#FF0092', '#22FFA7', '#E3EE9E', '#86CE00', '#BC7196', '#7E7DCD', '#FC6955', '#E48F72'])
fig.show()

gene_of_interest = data.T["ENSG00000162552"]                 #color by a signle gene
fig1 = px.scatter(loadings_data, x=XPC, y=YPC, color=gene_of_interest, hover_data=[XPC, YPC])
fig1.update_traces(mode='markers', marker_line_width=1, marker_size=8)
fig1.update_layout(title={"text": "WNT4", "x": 0.5, "y": 0.95, 'xanchor': 'center', 'yanchor': 'top'})
fig1.show()

tsne = TSNE(n_components=2, random_state=0, perplexity=4)   #5, 6, 7, 8, 9, 10, 11
projections = tsne.fit_transform(data_tr)

fig2 = px.scatter(
    projections, x=0, y=1,
    color=data.columns)
fig2.update_traces(mode='markers', marker_size=8)

x, y = [], []                   #add sample names
for i in projections:
    x += [i[0]]
    y += [i[1]]
fig2.add_trace(go.Scatter(
    x=x,
    y=y,
    mode="text",
    text=label_by,
    textposition="bottom center"))
fig2.show()



