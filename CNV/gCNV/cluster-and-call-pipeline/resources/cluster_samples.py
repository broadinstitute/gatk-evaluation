import argparse
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as shc
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
from typing import List


def cluster_sample_and_write_results(output_dir: str, max_clusters: int, read_count_files: List[str], entity_ids: List[str]):
        num_samples = len(read_count_files)
        counts = []
        for file in read_count_files:
            with open(file) as f:
                lines = (line for line in f if not line.startswith('@'))
                countraw = np.loadtxt(lines, skiprows=1, usecols=3)
                count = countraw / np.sum(countraw)
                counts.append(count)

        countmat = np.asmatrix(counts)
        samplesum = countmat.sum(axis=1)

        pca = PCA(n_components=25)
        pca_res = pca.fit(countmat)
        weights = pca_res.explained_variance_ratio_

        transformed = pca.fit_transform(countmat)
        transformed_weighted = transformed * weights
        ### Local - visualization
        # np.savetxt('/home/unix/jfu/transformed.csv', transformed, delimiter=',')

        # dend = shc.dendrogram(shc.linkage(transformed, method='ward'))

        Z = shc.linkage(transformed_weighted, method="ward")
        np.savetxt(os.path.join(output_dir, 'transformed_weighted.tsv'), transformed_weighted, delimiter='\t')
        np.savetxt(os.path.join(output_dir, 'samplesum.tsv'), samplesum, delimiter='\t')

        clusters = shc.fcluster(Z, t=0.003, criterion="distance")
        assert np.max(clusters) < max_clusters, "The clustering algorithm found more than the maximum number of clusters allowed"
        np.savetxt(os.path.join(output_dir, 'clusters.tsv'), clusters, delimiter='\t')

        column_names = ["Cluster" + str(i + 1) for i in range(max_clusters)]
        clustering_table = np.zeros(shape=(num_samples, max_clusters), dtype=np.float64)
        for idx, cluster in enumerate(clusters):
            clustering_table[idx][int(cluster) - 1] = 1.0
        clustering_table_df = pd.DataFrame(data=clustering_table, index=entity_ids, columns=column_names)
        clustering_table_df.to_csv(path_or_buf=os.path.join(output_dir, 'clustering_table.tsv'), sep='\t', columns=column_names)

        fig = plt.figure(figsize=(7, 7))
        plt.xlabel("PCA 1")
        plt.ylabel("PCA 2")
        plt.title("Sample clustering")
        plt.scatter(transformed[:, 0], transformed[:, 1], c=clusters)
        fig.savefig(os.path.join(output_dir, "sample_clusters.png"), dpi=120)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_dir', metavar='OutputDirectory', type=str,
                        help='output directory')

    parser.add_argument('--max_clusters', metavar='MaximumNumberOfClusters', type=int,
                        help='The maximum number of clusters allowed')

    parser.add_argument('--clustering_prior_table', metavar='ClusteringPriorTable', type=str,
                        help='A table of previosly clustered samples in identical format')

    parser.add_argument('--read_count_files', metavar='ReadCountFiles', type=str, nargs='+',
                        help='Read count files output by GATK CollectReadCounts')

    parser.add_argument('--entity_ids', metavar='EntityIDs', type=str, nargs='+',
                        help='Entity IDs of the correspoding coverage samples')

    ###################
    # Parse arguments #
    ###################
    args = parser.parse_args()
 
    max_clusters = args.max_clusters
    output_dir = args.output_dir
    read_count_files = args.read_count_files
    entity_ids = args.entity_ids
    clustering_prior_table = args.clustering_prior_table

    cluster_sample_and_write_results(output_dir, max_clusters, read_count_files, entity_ids)

if __name__ == '__main__':
    main()
