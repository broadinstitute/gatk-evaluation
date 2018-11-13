import matplotlib.pyplot as plt
import numpy as np


def plot_number_of_events_distribution(output_dir: str, gcnv_segment_vcfs: list):
    variant_num_list = []
    for gcnv_vcf in gcnv_segment_vcfs:
        variant_num_list.append(count_num_variants(gcnv_vcf))
    fig = plt.figure(figsize=(10,6))
    plt.title("Event number distribution")
    plt.hist(variant_num_list, bins=100, edgecolor='black', linewidth=1)
    fig.savefig(output_dir + 'event_num_distr.png', dpi=120)
    print("Mean number of segments in each sample: ")
    print(np.mean(np.array(variant_num_list)))
    print("\n")

def count_num_variants(vcf_file):
    num_vars = 0
    with open(vcf_file, 'r') as vcf:
        for line in vcf.readlines():
            if not line.startswith("#"):
                num_vars += 1
    return num_vars

def plot_ard_components(models_dir, num_shards):
    num_plots_per_row = 4
    num_rows = math.ceil(num_shards/num_plots_per_row)
    plt.figure(figsize=(20,20))
    for i in range(num_shards):
        ard_file = models_dir + "/shard-" + str(i) + "/mu_ard_u_log__.tsv"
        ard_df = pd.read_csv(open(ard_file, 'r'), comment="@", sep='\t', header=None)
        ard_vector = np.transpose(ard_df.as_matrix())[0]
        plt.subplot(num_rows, num_plots_per_row, i + 1)
        plt.plot(range(ard_vector.size), np.sort(ard_vector))
        plt.plot(range(ard_vector.size), np.zeros(ard_vector.size), 'r--')
        plt.ylim(-2, 2)
        plt.title("Shard " + str(i) + " ARD mean values")
    plt.tight_layout()