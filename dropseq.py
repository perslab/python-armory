import pandas as pd

# Function to read Seurat clustering
def get_clusters(cluster_file,header_bool,split_index_bool):
	df_cluster = pd.read_csv(cluster_file,compression="gzip",sep="\t",header=header_bool,index_col=0)
	df_cluster.columns = ['cluster_id']
	if split_index_bool:
		df_cluster.index = [x.split("_")[1] for x in df_cluster.index]
	df_cluster.drop(df_cluster.index[df_cluster.index.duplicated().tolist()],inplace=True) #NB upgrade to pandas 0.17 and use keep=False in duplicated()
	return df_cluster

# Function to average cells by cluster
def get_average_by_celltype(df_dge,df_cluster):
        return df_cluster.merge(df_dge.transpose(),left_index=True,right_index=True,how="inner").groupby('cluster_id',sort=False).mean().transpose()

# Function to standardize genes' expresion across cell types
def standardize(df):
        return df.sub(df.mean(axis=1),axis=0).div(df.std(axis=1),axis=0)
