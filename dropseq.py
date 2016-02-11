# Function to read Seurat clustering
def get_clusters(cluster_file):
	df_cluster = pd.read_csv(cluster_file,compression="gzip",sep="\t",header=None,index_col=0)
	df_cluster.columns = ['cluster_id']
	df_cluster.index = [x.split("_")[1] for x in df_cluster.index]
	df_cluster.drop(df_cluster.index[df_cluster.index.duplicated().tolist()],inplace=True) #NB upgrade to pandas 0.17 and use keep=False in duplicated()
	return df_cluster
