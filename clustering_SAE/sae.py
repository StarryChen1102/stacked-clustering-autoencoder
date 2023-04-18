import numpy as np
import pandas as pd
import tensorflow as tf
import keras
from tensorflow.keras.losses import MSE, KLD
from tensorflow.keras import activations, constraints, initializers, regularizers
import tensorflow.keras.backend as K
from keras.layers import Dense
from keras.models import Model
import math
import random

from layers import DenseTranspose
from utils import *


class SAE():

    def __init__(self, X, I, sample_size, n_gene, n_sample, dims=[512, 256, 128, 32]):
        super(SAE, self).__init__()
        encoded1 = Dense(dims[0], activation = "relu", kernel_constraint = keras.constraints.NonNeg(), input_shape=(n_gene,), use_bias=False,
                 kernel_initializer='glorot_uniform', bias_initializer='zeros')
        encoded2 = Dense(dims[1], activation = "relu", input_shape=(dims[0],), use_bias=True,
                 kernel_initializer='glorot_uniform', bias_initializer='zeros')
        encoded3 = Dense(dims[2], activation = "relu", input_shape=(dims[1],), use_bias=True,
                 kernel_initializer='glorot_uniform', bias_initializer='zeros')
        encoded4 = Dense(dims[3], activation = "relu", input_shape=(dims[2],), use_bias=True,
                 kernel_initializer='glorot_uniform', bias_initializer='zeros')

        encoder = keras.models.Sequential([encoded1, encoded2, encoded3, encoded4])
        encoder1 = keras.models.Sequential([encoded1])
        encoder2 = keras.models.Sequential([encoded2])
        encoder3 = keras.models.Sequential([encoded3])
        encoder4 = keras.models.Sequential([encoded4])

        decoded1 = DenseTranspose(encoded4, activation = "relu")
        decoded2 = DenseTranspose(encoded3, activation = "relu")
        decoded3 = DenseTranspose(encoded2, activation = "relu")
        decoded4 = DenseTranspose(encoded1, activation = "linear")

        decoder = keras.models.Sequential([decoded1, decoded2, decoded3, decoded4])
        decoder1 = keras.models.Sequential([decoded1])
        decoder2 = keras.models.Sequential([decoded2])
        decoder3 = keras.models.Sequential([decoded3])
        decoder4 = keras.models.Sequential([decoded4])

        autoencoder1 = keras.models.Sequential([encoded1, decoded4])
        autoencoder2 = keras.models.Sequential([encoded2, decoded3])
        autoencoder3 = keras.models.Sequential([encoded3, decoded2])
        autoencoder4 = keras.models.Sequential([encoded4, decoded1])
        autoencoder = keras.models.Sequential([encoder, decoder])

        idx = np.arange(0, n_sample)
        idx = random.sample(idx.tolist(), sample_size)

        self.X = X
        self.I = I

        if sample_size != n_sample:
            self.X_se = X[idx]
            self.I_se = I[idx]
        else:
            self.X_se = self.X
            self.I_se = self.I
        print(self.I_se.shape)

        idents_new = id2number(self.I_se)
        self.sample_size = sample_size

        self.ec1 = encoder1
        self.ec2 = encoder2
        self.ec3 = encoder3
        self.ec4 = encoder4
        self.dc1 = decoder1
        self.dc2 = decoder2
        self.dc3 = decoder3
        self.dc4 = decoder4
        self.autoencoder1 = autoencoder1
        self.autoencoder2 = autoencoder2
        self.autoencoder3 = autoencoder3
        self.autoencoder4 = autoencoder4
        self.ec = encoder
        self.dc = decoder
        self.autoencoder = autoencoder
           
    def train1(self, max_epoch=100, lr=0.002):
        optimizer = tf.keras.optimizers.Adam(lr)
        num_p, latent_p = cal_latent(self.X_se, 1.0)
        latent_p = latent_p + tf.linalg.diag(tf.linalg.diag_part(num_p))
        for epoch in range(0, max_epoch):
            with tf.GradientTape(persistent=True) as tape:
                h = self.ec1(self.X_se)
                z = self.dc4(h)

                #num_p, latent_p = cal_latent(self.X_se, 1.0)
                num_q, latent_q = cal_latent(h, 1.0)
                #latent_p = latent_p + tf.linalg.diag(tf.linalg.diag_part(num_p))
                latent_q = latent_q + tf.linalg.diag(tf.linalg.diag_part(num_q))

                cc_loss = tf.reduce_mean(KLD(latent_p, latent_q))

                loss = tf.reduce_mean(MSE(self.X_se, z)) + 1.0 * K.exp(cc_loss)
            vars = self.autoencoder1.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 10 == 0:
                print(loss)
        print("Finish!")
        return h
        
    def train2(self, h, max_epoch=120, lr=0.005):
        optimizer = tf.keras.optimizers.Adam(lr)
        for epoch in range(0, max_epoch):
            with tf.GradientTape(persistent=True) as tape:
                h1 = self.ec2(h)
                z = self.dc3(h1)
                loss = tf.reduce_mean(MSE(h, z))
            vars = self.autoencoder2.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 50 == 0:
                print(loss)
        print("Finish!")
        return h1
    
    def train3(self, h, max_epoch=350, lr=0.01):
        optimizer = tf.keras.optimizers.Adam(lr)
        for epoch in range(0, max_epoch):
            with tf.GradientTape(persistent=True) as tape:
                h1 = self.ec3(h)
                z = self.dc2(h1)
                loss = tf.reduce_mean(MSE(h, z))
            vars = self.autoencoder3.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 50 == 0:
                print(loss)
        print("Finish!")
        return h1
        
    def train4(self, h, max_epoch=350, lr=0.01):
        optimizer = tf.keras.optimizers.Adam(lr)
        for epoch in range(0, max_epoch):
            with tf.GradientTape(persistent=True) as tape:
                h1 = self.ec4(h)
                z = self.dc1(h1)
                loss = tf.reduce_mean(MSE(h, z))
            vars = self.autoencoder4.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 50 == 0:
                print(loss)
        print("Finish!")
        return h1

    def train(self, max_epoch=100, lr=0.001):
        weight0 = self.autoencoder.get_weights()
        optimizer = tf.keras.optimizers.Adam(lr)
        for epoch in range(1, max_epoch+1):
            with tf.GradientTape(persistent=True) as tape:
                
                h = self.ec(self.X_se)
                z = self.dc(h)

                a = K.eval(tf.math.is_nan(h))
                b = K.eval(tf.math.is_inf(h))
                if a.any() or b.any():
                    self.autoencoder.set_weights(weight0)
                    break;

                loss = tf.reduce_mean(MSE(self.X_se, z))
                weight0 = self.autoencoder.get_weights()
            
            vars = self.autoencoder.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 10 == 0:
                print(epoch)
                print(loss)
        print("Finish!")
    
    def clustering_train(self, max_epoch=50):
        self.autoencoder1.trainable = False
        self.autoencoder2.trainable = False

        # update h centers
        h = self.ec(self.X_se)
        h1 = self.ec1(self.X_se)
        h1 = self.ec2(h1)
        h1 = self.ec3(h1)
        centers, labels = get_centers(np.array(h))

        # update h1 centers
        k = len(np.unique(labels))
        labels1 = KMeans(n_clusters=k, random_state=0).fit(h1).labels_
        labels1 = update_labels(labels, labels1)
        centers1 = computeCentroids(h1, labels1)

        print(h.shape)
        print(h1.shape)

        q = cal_cluster(h, centers, 1.0)
        p = target_distribution(q)

        optimizer = tf.keras.optimizers.Adam()
        for epoch in range(1, max_epoch+1):
            if epoch % 10 == 0:

                # update h centers
                centers, labels = get_centers(np.array(h))

                # update h1 centers
                k = len(np.unique(labels))
                labels1 = KMeans(n_clusters=k, random_state=0).fit(h1).labels_
                labels1 = update_labels(labels, labels1)
                centers1 = computeCentroids(h1, labels1)

                print(labels1[0:10])
                print(labels[0:10])

                # update target
                q = cal_cluster(h, centers, 1.0)
                p = target_distribution(q)

            with tf.GradientTape(persistent=True) as tape:
                h = self.ec(self.X_se)
                z = self.dc(h)
                loss = tf.reduce_mean(MSE(self.X_se, z))
                loss = 1 * loss
                
                # clustering loss
                q_out = cal_cluster(h, centers, 1.0)
                cluster_loss = tf.reduce_mean(KLD(p, q_out))
                loss = loss + 1 * K.exp(cluster_loss)

                # center loss
                latent_dist1, latent_dist2, w = cal_dist(h, centers)
                center_loss = tf.reduce_mean(tf.reduce_sum(latent_dist2, axis=1))
                loss = loss + 1 * K.log(center_loss + 1)

                # Consistency loss w & q
                cluster_loss_1 = tf.reduce_mean(KLD(w, q_out))
                loss = loss + 1 * K.exp(cluster_loss_1)

                # Consistency loss layer & layer
                h1 = self.ec1(self.X_se)
                h1 = self.ec2(h1)
                h1 = self.ec3(h1)
                q_out1 = cal_cluster(h1, centers1, 1.0)
                kl_loss = tf.reduce_mean(KLD(q_out, q_out1))
                loss = loss + 1 * K.exp(kl_loss)

            vars = self.autoencoder.trainable_weights
            grads = tape.gradient(loss, vars)
            optimizer.apply_gradients(zip(grads, vars))
            if epoch % 10 == 0:
                print(loss)
                print(K.exp(cluster_loss))
                print(K.log(center_loss + 1))
                print(K.exp(cluster_loss_1))
                print(K.exp(kl_loss))

            if epoch % 25 == 0 or (epoch-13) % 25 == 0:
                pre, _,  _ = phenograph.cluster(np.array(h))
                print(epoch)
                #pre = KMeans(n_clusters=4, random_state=0).fit(np.array(h)).labels_
                idents_new = id2number(self.I_se)
                labels_new, ACC = count_accuracy(idents_new, pre)
                print(measure(self.I_se, pre))
                print(ACC)
             
        print("Finish!")
