workflow ClusterSamplesFromCoverageWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    Array[String]+ entity_ids
    Array[String]+ read_count_files
    Int maximum_number_of_clusters
    String gatk_docker

    ##################################
    #### optional basic arguments ####
    ##################################
    Int? preemptible_attempts

    call ClusterSamples {
        input:
            entity_ids = entity_ids,
            read_count_files = read_count_files,
            maximum_number_of_clusters = maximum_number_of_clusters,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File clustering_table = ClusterSamples.clustering_table
        File sample_clusters_plot = ClusterSamples.sample_clusters_plot
    }
}

task ClusterSamples {
    Array[String]+ entity_ids
    Array[File]+ read_count_files       # localized files
    Int maximum_number_of_clusters
    File? clustering_prior_table        # currently not used

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 3]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    command <<<
        set -e
        python - --output_dir . \
                 --entity_ids_file ${write_lines(entity_ids)} \
                 --read_counts_fof ${write_lines(read_count_files)} \
                 --maximum_number_of_clusters ${maximum_number_of_clusters} \
                 <<-'EOF'
        import argparse
        import os
        import numpy as np
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import scipy.cluster.hierarchy as shc
        from sklearn.decomposition import PCA

        def cluster_sample_and_write_results(output_dir: str,
                                             entity_ids_file: str,
                                             read_counts_fof: str,
                                             maximum_number_of_clusters: int):
            #load all files
            entity_ids = np.loadtxt(entity_ids_file, dtype=bytes).astype(str)
            read_count_files = np.loadtxt(read_counts_fof, dtype=bytes).astype(str)
            num_samples = len(entity_ids)
            assert num_samples == len(read_count_files), "Number of entity IDs and number of read count files must be equal."

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

            Z = shc.linkage(transformed_weighted, method="ward")
            np.savetxt(os.path.join(output_dir, 'transformed_weighted.tsv'), transformed_weighted, delimiter='\t')
            np.savetxt(os.path.join(output_dir, 'samplesum.tsv'), samplesum, delimiter='\t')

            clusters = shc.fcluster(Z, t=0.003, criterion="distance")
            assert np.max(clusters) < maximum_number_of_clusters, "The clustering algorithm found more than the maximum number of clusters allowed"
            np.savetxt(os.path.join(output_dir, 'clusters.tsv'), clusters, delimiter='\t')

            column_names = ["Cluster" + str(i + 1) for i in range(maximum_number_of_clusters)]
            clustering_table = np.zeros(shape=(num_samples, maximum_number_of_clusters), dtype=np.float64)
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

            parser.add_argument('--output_dir',
                                type=str,
                                help='Output directory.')

            parser.add_argument('--entity_ids_file',
                                type=str,
                                help='File containing list of entity IDs corresponding to read count files.')

            parser.add_argument('--read_counts_fof',
                                type=str,
                                help='File containing list of read count files (output of GATK CollectReadCounts) corresponding to entity IDs.')

            parser.add_argument('--maximum_number_of_clusters',
                                type=int,
                                help='The maximum number of clusters allowed.')

            args = parser.parse_args()

            output_dir = args.output_dir
            entity_ids_file = args.entity_ids_file
            read_counts_fof = args.read_counts_fof
            maximum_number_of_clusters = args.maximum_number_of_clusters

            cluster_sample_and_write_results(output_dir,
                                             entity_ids_file,
                                             read_counts_fof,
                                             maximum_number_of_clusters)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File clustering_table = "clustering_table.tsv"
        File sample_clusters_plot = "sample_clusters.png"
    }
}
