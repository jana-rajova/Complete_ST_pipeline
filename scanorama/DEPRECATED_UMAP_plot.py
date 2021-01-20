import numpy as np
import re
import scipy.cluster.hierarcipphy as sch
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
def create_stdata_dictionaries(csv_list, path):
    samples = []
    # dictionary with genes in the order they appear in the dataframe header
    gen_dim = {}
    # dictionary with position notations on the well
    positions = {}
    # dictionary with the data itself
    data = {}
    #load all the CNresults from the folder by some defining criterion
    #This version is for the instance, where we have the original csv files, when working with files outputted from scanorama, the format is slightly different!
    with open(path + "/" + csv_list, 'r') as f:
        for line in f:
            samples.append(line.rstrip('\n'))
        f.close()
        print(samples)
    try:
        for well in samples:
            name = re.search("([A-Z]{2}[0-9]+_[C-E][0-9])", well)
            print(name.group(0))
            df = pd.read_csv(path + "/" + well, index_col=0, header=0)
            data[name.group(0)] = df.iloc[:,:]
            # print(name.group(0))
            positions[name.group(0)] = df.index.tolist()
            gen_dim[name.group(0)] = df.columns
    except:
        print("Did not manage to read the datasets!")

    print("You have loaded the following number of wells:", len(data))
    return data, positions, gen_dim


def cust_plot(clusters, positions):
    temp_dict = dict()
    num_clusters = len(clusters)
    for cluster in set(clusters):
        for i in range(num_clusters):
            if clusters[i] == cluster:
                if cluster not in temp_dict.keys():
                    #print("initiating dictionary")
                    temp_dict[cluster] = np.array([positions[i]])
                else:
                    temp_dict[cluster] = np.vstack((temp_dict[cluster], positions[i]))

    # temp_dict2 = {}

    # for idx, cluster in enumerate(clusters):
    #     if cluster not in temp_dict2.keys():
    #         temp_dict2[cluster] = np.array([positions[idx]])
    #     else:
    #         temp_dict2[cluster] = np.vstack((temp_dict2[cluster], positions[idx]))

    # print(temp_dict)
    # print(temp_dict2)

    return temp_dict

def join_dataframe(dataset_dictionary):
    #in the following block, all the dataframes are joint into one (I am not sure if its fully called for)
    CN = pd.concat([dataset_dictionary[key] for key in dataset_dictionary.keys()], keys=[x for x in range(len(dataset_dictionary.keys()))])
    # dropping the useless index level that was placed there for no purpose, movint the index to a column so it can be used for a color and outputting it to a csv (just forcontrol, can be removed in further versions)
    CN.index = CN.index.droplevel(1)
    CN = CN.reset_index()
    CN = CN.astype({"index": 'int32'})
    CN.to_csv("combined.csv")
    #print("combined shape: ", CN.shape, sep='')
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
        assert len(clust_dict[key]) == len(plot_dict[key])
        start += len(value)
        # print(clust_dict)
        # print(key)
        # print(len(clust_dict[key]))
        for line in value:
            position_search = re.search("X([0-9.]+)_([0-9.]+)", line)
            X.append(float(position_search.group(1)))
            Y.append(float(position_search.group(2)))
            #print(line, position_search.group(1), position_search.group(2))
            df = pd.DataFrame(list(zip(positions[key], clust_dict[key], X, Y)), columns=['Position', 'Cluster', 'X', 'Y'])
            export_dict[key] = df
    return plot_dict, clust_dict, export_dict

if __name__ == '__main__':
    #which version of this file this is? Not important for later
    # print(os.path.getmtime("bin/UMAP_plot.py"))
    # windows addition to get into folder
    # do something a bit more clever
    parser = argparse.ArgumentParser(description='Add treshold and what are the dimensions of the scanorama dataset in quesion')
    parser.add_argument('--treshold_united', '-u', type=int, help='treshold for cluster assignment with all wells clustered together')
    parser.add_argument('--treshold_separate', '-s', type=int, default=3, help='treshold for cluster assignment in clustering with separate wells')
    # parser.add_argument('-d', '--dimensions', type=int, help='how many dimensions does the scanorama file have?')
    parser.add_argument('-c', '--clustering_dimensions', type=int, help='in how many dimensions does the ahc algorithm cluster?')
    parser.add_argument('-t', '--timestamp', type=str, help='scanorama result file timestamp')
    parser.add_argument("-o", "--output", type=str, help="output folder")
    args = parser.parse_args()

    output_folder = args.output + "/UMAP_output_" + args.timestamp
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    os.chdir(output_folder)
    print("Current  directory: ", os.getcwd())
    # print(os.listdir())

    data, positions, gen_dim = create_stdata_dictionaries(csv_list="csv_files.txt", path="../scanorama_output_" + args.timestamp)
    CN = join_dataframe(data)

    """
    The following are parameters for UMAP 
    """
    n_comp = args.clustering_dimensions
    dist = 0.3
    neigh = 15

    #embedding for clustering (more dimensions than embedding for the projection
    embedding_united_clust = umap.UMAP(n_neighbors=neigh, min_dist=dist, n_components=n_comp).fit_transform(CN.iloc[:, 1:])
    print("Dimensions of clustering embedding: ", embedding_united_clust.shape)

    # create embedding dictionary to be able to further output it as a csv and separate the embedding data back into their respective wells this is for clustering, not plotting
    embedding_dict = dict()
    start = 0
    # print("dict embedding")
    for key, value in positions.items():
        end = start + len(value)
        embedding_dict[key] = embedding_united_clust[start:end]
        start += len(value)

     
    # treshold for ha clustering assgnment of own cluster for all the samples (is higher than for the separate wells)
    # In the future, should it be a fraction of the total amount of features? That way it doesn't have to de defined specificlly?? Probably not linear

    #create a folder for united clustering
    united_folder = "united_clustering-treshold_" + str(args.treshold_united) + "-clustering_dimensions_" + str(args.clustering_dimensions)

    if not os.path.isdir(united_folder):
        os.mkdir(united_folder)
    os.chdir(united_folder)
    # create cluster for all the sections and then move onto creating subclusterings
    """
    Create dendrogram for united clustering 
    """
    title = "all" + str(neigh) + "n_neighbors_" + str(dist) + "min_dist_" + str(args.treshold_united)+ "treshold-united"
    dendrogram = sch.dendrogram(sch.linkage(embedding_united_clust, method='ward'))
    plt.title(title + "_dendrogram")
    plt.savefig(title + "_dendrogram.jpg")
    plt.clf()
    hc = ac(n_clusters=None, affinity='euclidean', linkage='ward', distance_threshold=args.treshold_united)
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


    # export the cluster and feature name data to csv, from where it can be incorporated into the R script
    for key, value in export_dict.items():
        value.to_csv(key + 'treshold' + str(args.treshold_united) + "_united_clusters.csv", index=False, header=False)

    temp_dict = cust_plot(hc,embedding_united)
    for key in temp_dict.keys():
        plt.scatter(temp_dict[key][:,0], temp_dict[key][:,1],c=sns.color_palette('hls', cluster_tot)[key],s=1, label = str(key))
    plt.title(title + "clusters_scatter")
    plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2, markerscale=3)
    #plt.show()
    plt.savefig(title + "clusters_scatter.jpg", bbox_inches='tight')
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
    plt.clf()

    # plotting each of the wells separately, but the clusters are from the united clustering
    for key, value in plot_dict.items():
        print("plotting: ", key)
        title = key + "_" + str(neigh) + "n_neighbors_" + str(dist) + "min_dist_" + str(args.treshold_united)+ "treshold-united"

        # it is possible to re-cluster the saved embedding from the "treshold" dimensional UMAP  with agglomerative hierarchical clustering searately for each well
        hc = clust_dict[key]
        plt.title(title + "_scatter_common_cluster")
        plt.xlim(left, right)
        plt.ylim(top, down)
        temp_dict = cust_plot(hc,value)
        # print(cluster_tot)
        for clust, pos in temp_dict.items():
            plt.scatter(pos[:,0], pos[:,1], c=sns.color_palette('hls', cluster_tot+1)[clust],s=1, label=str(clust))
        # couldn't create a legend for this....
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2, markerscale=3)
        plt.savefig(title + "_scatter_common_cluster.jpg", bbox_inches='tight')
        plt.clf()
        features = []
        for i in range(len(export_dict[key])):
            features.append(export_dict[key].iloc[i, 2:])
        features = np.array(features)
        temp_dict = cust_plot(hc, features)
        for clust,pos in temp_dict.items():
            plt.scatter(pos[:, 0], pos[:, 1], c=sns.color_palette('hls', cluster_tot+1)[clust], label=str(clust))
        plt.title(title + "_scatter_common_cluster_coord")
        plt.legend(bbox_to_anchor=(1, 1), borderaxespad=0, fontsize='xx-small', labelspacing=None, ncol=2)
        plt.savefig(title + "_scatter_common_cluster_coord.jpg", bbox_inches='tight')
        plt.clf()
    os.chdir("../")

    #cluster wells separately from each other
    n_comp = args.clustering_dimensions
    dist = 0.0
    neigh = 20
    separate_folder = "separate_clustering-treshold_" + str(args.treshold_separate) + "-clustering_dimensions_" + str(args.clustering_dimensions)
    try:
        os.chdir(separate_folder)
    except:
        os.mkdir(separate_folder)
        os.chdir(separate_folder)

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
        # print(features.shape)
        # print(features)
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
    os.chdir("../")
