import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.manifold import TSNE
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score, homogeneity_score, completeness_score
import tensorflow.keras.backend as K
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
import phenograph
import umap

def myscatter(Y, class_idxs, legend=False, ran=True, seed=229):
    if ran:
        np.random.seed(seed)
    Y = np.array(Y)
    fig, ax = plt.subplots(figsize=(5,4), dpi=300)
    classes = list(np.unique(class_idxs))
    markers = 'osD' * len(classes)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(classes)))
    if ran:
        np.random.shuffle(colors)

    for i, cls in enumerate(classes):
        mark = markers[i]
        ax.plot(Y[class_idxs == cls, 0], Y[class_idxs == cls, 1], marker=mark,
                linestyle='', ms=4, label=str(cls), alpha=1, color=colors[i],
                markeredgecolor='black', markeredgewidth=0.15)
    if legend:
        ax.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0, fontsize=10, markerscale=2, frameon=False,
                  ncol=2, handletextpad=0.1, columnspacing=0.5)

    plt.xticks([])
    plt.yticks([])

    return ax

def dotsne(X, dim=2, ran=23):
    tsne = TSNE(n_components=dim, random_state=ran)
    Y_tsne = tsne.fit_transform(X)
    return Y_tsne

def doumap(X, dim=2):
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(X)
    return embedding

def clustering(h, n_cluster, k=15, f="kmeans"):
    if f == "kmeans":
        labels = KMeans(n_clusters=n_cluster, random_state=0).fit(h).labels_
    elif f == "spectral":
        labels = SpectralClustering(n_clusters=n_cluster, affinity="precomputed", assign_labels="discretize",
                                    random_state=0).fit_predict(adj)
    return labels

def measure(true, pred):
    NMI = round(normalized_mutual_info_score(true, pred), 2)
    RAND = round(adjusted_rand_score(true, pred), 2)
    HOMO = round(homogeneity_score(true, pred), 2)
    COMP = round(completeness_score(true, pred), 2)
    return [NMI, RAND, HOMO, COMP]

def count_accuracy(truth, pred):
    n1 = len(np.unique(truth))
    n2 = len(np.unique(pred))
    pred_new = np.zeros(len(pred))
    for i in range(n2):
        x = np.zeros(n1)
        id = np.where(pred==i)[0]
        m = 0
        for k in range(n1):
            x[k] = np.sum(truth[id]==k)
        m = np.where(x==max(x))[0][0]
        pred_new[id] =  m
    return pred_new, np.sum(pred_new==truth)/len(truth)

def cal_kendall(dat1, dat2):
    c = 0
    d = 0
    for i in range(len(dat1)):
        for j in range(i+1,len(dat1)):
            if (dat1[i]-dat1[j])*(dat2[i]-dat2[j])>0:
               c = c + 1
            else:
               d = d + 1
    k_tau = (c - d) * 2 / len(dat1)/(len(dat1)-1)
    return k_tau

def target_distribution(q):
    q = q.numpy()
    weight = q ** 2 / q.sum(0)
    return (weight.T / weight.sum(1)).T

def cal_cluster(hidden, clusters, alpha):
    clusters = tf.convert_to_tensor(clusters, dtype=tf.float32)
    q = 1.0 / (1.0 + (K.sum(K.square(K.expand_dims(hidden, axis=1) - clusters), axis=2) / alpha))
    q **= (alpha + 1.0) / 2.0
    q = K.transpose(K.transpose(q) / K.sum(q, axis=1))
    return q

def cal_dist(hidden, clusters):
    clusters = tf.convert_to_tensor(clusters, dtype=tf.float32)
    dist1 = K.sum(K.square(K.expand_dims(hidden, axis=1) - clusters), axis=2)
    temp_dist1 = dist1 - tf.reshape(tf.reduce_min(dist1, axis=1), [-1, 1])
    temp_dist1 = K.transpose(K.transpose(temp_dist1) / K.max(temp_dist1, axis=1))
    q = K.exp(-temp_dist1 * 10)
    q = K.transpose(K.transpose(q) / K.sum(q, axis=1))
    q = K.pow(q, 2)
    q = K.transpose(K.transpose(q) / K.sum(q, axis=1))
    dist2 = dist1 * q
    return dist1, dist2, q

def cal_latent(hidden, alpha):
    sum_y = K.sum(K.square(hidden), axis=1)
    num = -2.0 * tf.matmul(hidden, tf.transpose(hidden)) + tf.reshape(sum_y, [-1, 1]) + sum_y
    num = num / alpha
    num = tf.pow(1.0 + num, -(alpha + 1.0) / 2.0)
    zerodiag_num = num - tf.linalg.diag(tf.linalg.diag_part(num))
    latent_p = K.transpose(K.transpose(zerodiag_num) / K.sum(zerodiag_num, axis=1))
    return num, latent_p

def target_dis(latent_p):
    latent_q = tf.transpose(tf.transpose(tf.pow(latent_p, 2)) / tf.reduce_sum(latent_p, axis = 1))
    return tf.transpose(tf.transpose(latent_q) / tf.reduce_sum(latent_q, axis = 1))

def computeCentroids(data, labels):
    n_clusters = len(np.unique(labels))
    data_1 = np.array(data)
    return np.array([data_1[labels == i].mean(0) for i in range(n_clusters)])

def get_centers(Y):
    l, _,  _ = phenograph.cluster(Y)
    centers = computeCentroids(Y, l)
    return centers, l

def get_centers_kmeans(Y, k):
    l = KMeans(n_clusters=k, random_state=0).fit(Y).labels_
    centers = computeCentroids(Y, l)
    return centers, l

def update_labels(labels, labels1):
    D = max(labels.max(), labels1.max()) + 1
    w = np.zeros([D,D])
    for i in range(labels1.size):
        w[labels1[i], labels[i]] += 1
    x, y = linear_sum_assignment(w.max() - w)
    dic = dict(zip(x,y))
    for i in range(len(labels1)):
        labels1[i] = dic[labels1[i]]
    return labels1

def id2number(idents):
    cells = np.unique(idents)
    idents_new = np.zeros(len(idents))
    for i in range(len(cells)):
        idents_new[idents==cells[i]] = i
    idents_new = idents_new.astype('int')
    return idents_new
