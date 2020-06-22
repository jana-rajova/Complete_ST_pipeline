import numpy as np
import re
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering as ac
import hdbscan
import umap
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


# genes is again specific to origianl csv files, but maybe we can also do it for scanorama and have the dimensions there instead if I find a way to find out,which genes contributed to which columns??
# contains lines from the csv file
def create_stdata_dictionaries(csv_list="csv_files.txt"):
    samples = []
    # dictionary with genes in the order they appear in the dataframe header
    gen_dim = {}
    # dictionary with position notations on the well
    positions = {}
    # dictionary with the data itself
    data = {}
    #load all the CNresults from the folder by some defining criterion
    #This version is for the instance, where we have the original csv files, when working with files outputted from scanorama, the format is slightly different!
    with open(csv_list, 'r') as f:
        for line in f:
            samples.append(line.rstrip('\n'))
        f.close()
    try:
        for well in samples:
            name = re.search("(CN[0-9]{2}_[C-E][0-9])", well)
            df = pd.read_csv(well, index_col=0, header=0)
            data[name.group(0)] = df.iloc[:,:]
            print(name.group(0))
            positions[name.group(0)] = df.index.tolist()
            gen_dim[name.group(0)] = df.columns
    except:
        print("Did not manage to read the datasets!")

    print("You have loaded the following number of wells:", len(data))
    return data, positions, gen_dim

def cust_plot(clusters, positions):
    temp_dict = dict()
    for cluster in set(clusters):
        for i in range(len(clusters)):
            if clusters[i]==cluster:
                if cluster not in temp_dict.keys():
                    #print("initiating dictionary")
                    temp_dict[cluster] = np.array([positions[i]])
                else:
                    temp_dict[cluster] = np.vstack((temp_dict[cluster],positions[i]))
    return temp_dict

def join_dataframe(dataset_dictionary):
    #in the following block, all the dataframes are joint into one (I am not sure if itäs fully called for)
    CN = pd.concat([dataset_dictionary[key] for key in dataset_dictionary.keys()], keys=[x for x in range(len(dataset_dictionary.keys()))])
    # dropping the useless index level that was placed there for no purpose, movint the index to a column so it can be used for a color and outputting it to a csv (just forcontrol, can be removed in further versions)
    CN.index = CN.index.droplevel(1)
    CN = CN.reset_index()
    CN = CN.astype({"index": 'int32'})
    CN.to_csv("combined.csv")
    print("combined shape: ", CN.shape, sep='')
    return CN

def cluster_coordinates_dict(positions, embedding, hc):
    # create dictionaries, that will store the plottting data separate for each well as well as one that is for exporting cluster data and feature name data back to the R script (export_dict)
    plot_dict = dict()
    clust_dict = dict()
    export_dict = dict()
    # create a dictionary with coordinates for plotting and associate the clusters to them and plot them
    start = 0
    for key, value in positions.items():
        # print(start)
        X = list()
        Y = list()
        end = start + len(value)
        # print(end)
        plot_dict[key] = embedding[start:end]
        # print(plot_dict)
        clust_dict[key] = hc[start:end]
        # print(clust_dict)
        # print(key)
        # print(len(clust_dict[key]))
        for line in value:
            position_search = re.search("X([0-9.]+)_([0-9.]+)",line)
            X.append(float(position_search.group(1)))
            Y.append(float(position_search.group(2)))
            #print(line, position_search.group(1), position_search.group(2))
            assert len(clust_dict[key]) == len(plot_dict[key])
            df = pd.DataFrame(list(zip(positions[key],clust_dict[key],X,Y)), columns=['Position', 'Cluster', 'X', 'Y'])
            export_dict[key] = df
        start += len(value)   
    return plot_dict, clust_dict, export_dict

if __name__ == '__main__':
    #which version of this file this is? Not important for later
    print(os.path.getmtime("bin/UMAP_plot.py"))
    # windows addition to get into folder
    # do something a bit more clever
    parser = argparse.ArgumentParser(description='Add treshold and what are the dimensions of the scanorama dataset in quesion')
    parser.add_argument('--treshold_united', '-u', type=int, help='treshold for cluster assignment with all wells clustered together')
    parser.add_argument('--treshold_separate', '-s', type=int, default=3, help='treshold for cluster assignment in clustering with separate wells')
    parser.add_argument('-d', '--dimensions', type=int, help='how many dimensions does the scanorama file have?')
    parser.add_argument('-c', '--clustering_dimensions', type=int, help='in how many dimensions does the ahc algorithm cluster?')
    args = parser.parse_args()

    try:
        os.chdir("/media/dropbox//MNM team folder/Spatial transcriptomics/Scanorama_mod/results/Union-True/DIMRED_" + str(args.dimensions))
        print("Folder changed")
    except:
        print("In ", os.getcwd())
    print(os.listdir())

    #the following commented block of code is just to verify what have you loaded into the datasets
    # print("These are the shapes of the data sets:")

    # verify the dimensions and the values and the sizes of the dictionaries
    # for key, value in data.items():
    #     print(key, value.shape)
    #     print(data[key].head(5))
    # print("These are the shapes of the positions sets:")
    # for key, value in positions.items():
    #     print(key, len(value))
    #     print(value)
    # print("These are the shapes of the gene sets:")
    # for key, value in gen_dim.items():
    #     print(key, value.shape)
    #     print(value)

    data, positions, gen_dim = create_stdata_dictionaries()
    CN = join_dataframe(data)
    

    # for key, value in positions.items():
    #     print(key, len(value))
    # print('positions')

    # for key, value in data.items():
    #     print(key, value)
    # print('data?')

    # for key, value in gen_dim.items():
    #     print(key, len(value))
    # print('genes?')
    # cluster the whole field tohether and then separate into well files
    # parameters for umap
    n_comp = args.clustering_dimensions
    dist = 0.3
    neigh = 15
    # input(positions)
    #embedding for clustering (more dimensions than embedding for the projection
    embedding_united_clust = umap.UMAP(n_neighbors=neigh, min_dist=dist, n_components=n_comp).fit_transform(CN.iloc[:, 1:])
    print("Dimensions of clustering embedding: ", embedding_united_clust.shape)

    # create embedding dictionary to be able to further output it as a csv and separate the embedding data back into their respective wells this is for clustering, not plotting
    embedding_dict = dict()
    start = 0
    print("dict embedding")
    for key, value in positions.items():
        end = start + len(value)
        embedding_dict[key] = embedding_united_clust[start:end]
        start += len(value)

        print(key)
        print(len(embedding_dict[key]))
        print(type(embedding_dict[key]))
        print()
    # treshold for ha clustering assgnment of own cluster for all the samples (is higher than for the separate wells)
    # In the future, should it be a fraction of the total amount of features? That way it doesn't have to de defined specificlly?? Probably not linear

    #create a folder for united clustering
    try:
        os.chdir("united_clustering")
    except:
        os.mkdir("united_clustering")
        os.chdir("united_clustering")

    try:
        os.chdir("treshold_" + str(args.treshold_united))
    except:
        os.mkdir("treshold_" + str(args.treshold_united))
        os.chdir("treshold_" + str(args.treshold_united))
    try:
        os.chdir("clustering_dimensions_" + str(args.clustering_dimensions))
    except:
        os.mkdir("clustering_dimensions_" + str(args.clustering_dimensions))
        os.chdir("clustering_dimensions_" + str(args.clustering_dimensions))
    # create cluster for all the sections and then move onto creating subclusterings
    title = "all" + str(neigh) + "n_neighbors_" + str(dist) + "min_dist_" + str(args.treshold_united)+ "treshold-united"
    dendrogram = sch.dendrogram(sch.linkage(embedding_united_clust, method='ward'))
    plt.title(title + "_dendrogram")
    plt.savefig(title + "_dendrogram.jpg")
    plt.clf()
    hc = ac(n_clusters=None, affinity = 'euclidean', linkage = 'ward', distance_threshold=args.treshold_united)
    hc = hc.fit_predict(embedding_united_clust)
    print("number of clusters total:", len(set(hc)))
    cluster_tot = len(set(hc))
    print("lenght hc :", len(hc), sep='')
    # the following embedding is for visualization, but the clusters are comming from higher dimensions
    embedding_united = umap.UMAP(n_neighbors=neigh, min_dist=dist, n_components=2).fit_transform(CN.iloc[:, 1:])

    # plot_dict = separate embedding into wells
    # clust_dict = ditionary of cluster associated to positions in the same order as are coordinates in the plot_dict
    # export_dict = dictionary of dataframes of feature names (e.g. X12_34), cluster they were assigned to and feature X and Y coordinate in the well (not UMAP coordinate) for each well
    plot_dict, clust_dict, export_dict = cluster_coordinates_dict(positions, embedding_united, hc)

    # for key, value in plot_dict.items():
    #     print(key, value)
    # input('plot_dict ok?')

    # for key, value in clust_dict.items():
    #     print(key, value)
    # input('clust_dict ok?')

    # for key, value in export_dict.items():
    #     print(key, value)
    # input('export_dict ok?')

    # export the cluster and feature name data to csv, from where it can be incorporated into the R script
    for key, value in export_dict.items():
        value.to_csv(key + 'treshold' + str(args.treshold_united) + "_united_clusters.csv", index=False, header=False)


    # plot all spots, but show the output of clustering of all the wells together (we are looking for tx cluster in all...)

    #temporary dictionary used for plotting with labels
    # temp_dict = dict()
    #
    # for cluster in sorted(set(hc)):
    #     for i in range(len(hc)):
    #         if hc[i]==cluster:
    #             if cluster not in temp_dict.keys():
    #                 print("initiating dictionary")
    #                 temp_dict[cluster] = np.array(embedding_united[i])
    #             else:
    #                 temp_dict[cluster] = np.vstack((temp_dict[cluster],embedding_united[i]))

    temp_dict = cust_plot(hc,embedding_united)
    for key in temp_dict.keys():
        plt.scatter(temp_dict[key][:,0], temp_dict[key][:,1],c=sns.color_palette('hls', cluster_tot)[key],s=1, label = str(key))
    plt.title(title + "clusters_scatter")
    plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2, markerscale=3)
    #plt.show()
    plt.savefig(title + "clusters_scatter.jpg", bbox_inches='tight')

    # didn't manage to create legend for the following block...
    # # scatter_all_clust = plt.scatter(embedding_united[:,0], embedding_united[:,1], c=[sns.color_palette("hls", cluster_tot)[x] for x in hc], s=1)
    # # plt.title(title + "clusters_scatter")
    # # plt.legend((c for c in sns.color_palette("hls", cluster_tot)), (x for x in sorted(set(hc))))
    # # plt.show()
    # # plt.savefig(title + "clusters_scatter.jpg")
    plt.clf()

    # with the same plotting embedding as for the clusters in all the wells, show the overlap in between the wells

    i=0
    for key,value in plot_dict.items():
        plt.scatter(value[:,0], value[:,1],c=sns.color_palette('hls', len(plot_dict))[i],s=1, label=key)
        i+=1
    plt.title(title + "wells_scatter")
    plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, markerscale=3)
    left, right = plt.xlim()
    top, down = plt.ylim()
    #plt.show()
    plt.savefig(title + "_wells_scatter.jpg", bbox_inches='tight')

    # didn't manage to create legend for the following block...
    # plt.scatter(embedding_united[:,0], embedding_united[:,1], c=[sns.color_palette("hls", len(positions.keys()))[x] for x in CN.iloc[:,0]], s=1)
    # plt.legend(*scatter_all.legend_elements())
    # plt.savefig(title + "_wells_scatter.jpg")


    print("axes limits:")
    print(left, right, top, down)
    plt.clf()

    # plotting each of the wells separately, but the clusters are from the united clustering
    for key, value in plot_dict.items():
        print("plotting: ", key)
        title = key + "_" + str(neigh) + "n_neighbors_" + str(dist) + "min_dist_" + str(args.treshold_united)+ "treshold-united"

        # it is possible to re-cluster the saved embedding from the "treshold" dimensional UMAP  with agglomerative hierarchical clustering searately for each well
        hc = clust_dict[key]

        print(hc)
        print("hc")
        print(len(hc))
        print(value)
        print(len(value))
        print(type(value))
        plt.title(title + "_scatter_common_cluster")
        plt.xlim(left, right)
        plt.ylim(top, down)
        temp_dict = cust_plot(hc,value)
        print(cluster_tot)
        for clust, pos in temp_dict.items():
            plt.scatter(pos[:,0], pos[:,1], c=sns.color_palette('hls', cluster_tot+1)[clust],s=1, label=str(clust))
        # couldn't create a legend for this....
        #scatter_well = plt.scatter(value[:,0], value[:,1], c=[sns.color_palette("hls", cluster_tot)[x] for x in hc], s=1)
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2, markerscale=3)
        #plt.show()
        plt.savefig(title + "_scatter_common_cluster.jpg", bbox_inches='tight')
        plt.clf()
        features = []
        for i in range(len(export_dict[key])):
            features.append(export_dict[key].iloc[i, 2:])
        features = np.array(features)
        print(features.shape)
        print(features)
        temp_dict = cust_plot(hc, features)
        for clust,pos in temp_dict.items():
            plt.scatter(pos[:, 0], pos[:, 1], c=sns.color_palette('hls', cluster_tot+1)[clust], label=str(clust))
        #plt.scatter(features[:, 0], features[:, 1], c=[sns.color_palette("hls", cluster_tot)[x] for x in hc])
        plt.title(title + "_scatter_common_cluster_coord")
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2)
        #plt.show()
        plt.savefig(title + "_scatter_common_cluster_coord.jpg", bbox_inches='tight')
        plt.clf()
    os.chdir("../../../")

    # # HDBSCAN tryout
    # it is not totally useless, however the algorithm is not very good for this ype of data as it drops way too many values to avoid false categorization
    #
    # print("HDBSCAN trial 40dim clustering")_E
    # for key in data.keys():
    #     components = 400
    #     # print(components)
    #     # print(data[key].index)
    #     # print(data[key].columns)
    #     reducer_emb = umap.UMAP(n_neighbors=15, min_dist=0.0, n_components=components, random_state=42)
    #     embeddable = reducer_emb.fit_transform(data[key])
    #     labels = hdbscan.HDBSCAN(min_samples=5, min_cluster_size=10).fit_predict(embeddable)
    #     clustered = (labels >=0)
    #     reducer_sep = umap.UMAP()
    #     embedding_sep = reducer_sep.fit_transform(data[key].iloc[:,1:])
    #     plt.scatter(embedding_sep[~clustered, 0], embedding_sep[~clustered, 1], c=(0.5, 0.5, 0.5), s=0.1, alpha=0.5)
    #     plt.scatter(embedding_sep[clustered, 0], embedding_sep[clustered, 1], c=labels[clustered], s=0.1, cmap='Spectral')
    #     plt.show()

    #cluster wells separately from each other
    n_comp = args.clustering_dimensions
    dist = 0.0
    neigh = 20
    try:
        os.chdir("separate_clustering")
    except:
        os.mkdir("separate_clustering")
        os.chdir("separate_clustering")
    try:
        os.chdir("treshold_" + str(args.treshold_separate))
    except:
        os.mkdir("treshold_" + str(args.treshold_separate))
        os.chdir("treshold_" + str(args.treshold_separate))
    try:
        os.chdir("clustering_dimensions_" + str(args.clustering_dimensions))
    except:
        os.mkdir("clustering_dimensions_" + str(args.clustering_dimensions))
        os.chdir("clustering_dimensions_" + str(args.clustering_dimensions))

    for key, value in data.items():
        clusterable_embedding = umap.UMAP(n_neighbors=neigh, min_dist=dist, n_components=n_comp, random_state=42).fit_transform(value)
        title = key + "_" + str(neigh) + "n_neighbors_" + str(dist) + "min_dist_" + str(args.treshold_separate)+ "treshold"
    #approach with agglomerative hierarchical clustering
        dendrogram = sch.dendrogram(sch.linkage(clusterable_embedding, method='ward'))
        plt.title(title + "_dendrogram")
        plt.savefig(title + "_dendrogram.jpg", bbox_inches='tight')
        plt.clf()
        hc = ac(n_clusters=None, affinity = 'euclidean', linkage = 'ward', distance_threshold=args.treshold_separate)
        hc = hc.fit_predict(clusterable_embedding)
        cluster_tot = len(set(hc))
        print("lenght hc :", len(hc), sep='')
        embedding = umap.UMAP(n_neighbors=neigh, min_dist=dist, n_components=2, random_state=42).fit_transform(value)
        #plt.scatter(embedding[:,0], embedding[:,1], c=[sns.color_palette("hls", cluster_tot)[x] for x in hc], s=1)
        temp_dict = cust_plot(hc, embedding)
        for clust,pos in temp_dict.items():
            plt.scatter(pos[:,0], pos[:,1], c=sns.color_palette('hls', cluster_tot+1)[clust],s=1, label=str(clust))
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2, markerscale=3)
        plt.title(title + "_scatter")
        #plt.show()
        plt.savefig(title + "_scatter.jpg", bbox_inches='tight')
        plt.clf()
        features = []
        for i in range(len(export_dict[key])):
            features.append(export_dict[key].iloc[i, 2:])
        features = np.array(features)
        print(features.shape)
        print(features)
        temp_dict = cust_plot(hc, features)
        for clust,pos in temp_dict.items():
            plt.scatter(pos[:,0], pos[:,1],c=sns.color_palette('hls', cluster_tot)[clust], label=str(clust))
        #plt.scatter(features[:, 0], features[:, 1], c=[sns.color_palette("hls", cluster_tot)[x] for x in hc])
        plt.title(title + "_scatter_sepatate_cluster_coord")
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2)
        #plt.show()
        plt.savefig(title + "_scatter_sepatate_cluster_coord.jpg", bbox_inches='tight')
        plt.clf()
    plt.close()
    os.chdir("../../../")
    #     elif y_hc[i] == 1:
    #         plt.scatter(embedding[i,0], embedding[i,1], c='black')
    #     elif y_hc[i] == 2:
    #         plt.scatter(embedding[i,0], embedding[i,1], c='blue')
    #     elif y_hc[i] == 3:
    #         plt.scatter(embedding[i,0], embedding[i,1], c='cyan')
