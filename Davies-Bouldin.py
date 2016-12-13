def davies_bouldin(X, labels, cluster_ctr):
    #get the cluster assignemnts
    clusters = set(labels)
    #get the number of clusters
    num_clusters = len(clusters)
    #array to hold the number of items for each cluster, indexed by cluster number
    num_items_in_clusters = [0]*num_clusters
    #get the number of items for each cluster
    for i in range(len(labels)):
        num_items_in_clusters[labels[i]] += 1
    max_num = -9999
    for i in range(num_clusters):
        s_i = intra_cluster_dist(X, labels, clusters[i], num_items_in_clusters[i], cluster_ctr[i])
    for j in range(num_clusters):
        if(i != j):
            s_j = intra_cluster_dist(X, labels, clusters[j], num_items_in_clusters[j], cluster_ctr[j])
            m_ij = np.linalg.norm(cluster_ctr[clusters[i]]-cluster_ctr[clusters[j]])
            r_ij = (s_i + s_j)/m_ij
            if(r_ij > max_num):
                max_num = r_ij
return max_num

def intra_cluster_dist(X, labels, cluster, num_items_in_cluster, centroid):
    total_dist = 0
    #for every item in cluster j, compute the distance the the center of cluster j, take average
    for k in range(num_items_in_cluster):
        dist = np.linalg.norm(X[labels==cluster]-centroid)
        total_dist = dist + total_dist
return total_dist/num_items_in_cluster
