import seaborn as sns
import matplotlib.pyplot as plt
import os


def visualize_events_bars(deg_df, output_loc):

    # Grouping by cluster assignment and type of feature, and counting occurrences
    grouped_deg_df = deg_df.groupby(['cluster', 'EventAnnotation']).size().unstack(fill_value=0)

    # Calculating percentage of each type of feature in each cluster
    grouped_annotation_percentage = grouped_deg_df.div(grouped_deg_df.sum(axis=1), axis=0) * 100

    # Define color palette for the types of features
    colors = sns.color_palette("Set1")

    # Plotting horizontal stacked bar graph
    plt.figure(figsize=(10, 6))
    grouped_annotation_percentage.plot(kind='barh', stacked=True, color=colors)
    plt.title('Percentage of Features in Each Type by Cluster Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Percentage')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(output_loc, "event_annotations_percentage.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.show()

    plt.figure(figsize=(10, 6))
    grouped_deg_df.plot(kind='barh', stacked=True, color=colors)
    plt.title('Number of Features in Each Type by Cluster Assignment')
    plt.ylabel('Cluster Assignment')
    plt.xlabel('Number of Features')
    plt.legend(title='Event Annotation')
    plt.tight_layout()  # Adjust layout to fit cleanly in PDF
    file_path = os.path.join(output_loc, "event_annotations_numbers.pdf")
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.show()

    return grouped_deg_df, grouped_annotation_percentage


def psi_annotated_heatmap(psi_file, deg_df, clusters_df):

    heatmap_matrix = 2
    return heatmap_matrix


def binary_cluster_heatmap(clusters_df):

    heatmap_matrix = 2

    return heatmap_matrix

