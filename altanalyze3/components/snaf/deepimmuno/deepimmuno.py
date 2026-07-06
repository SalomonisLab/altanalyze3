

'''
Program to run deepimmuno-cnn
'''
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import warnings
warnings.filterwarnings('ignore')
import logging
logger = logging.getLogger()   # this is a singleton, all .tf .keras logger will pass the info to this module level rootlogger
logger.setLevel(level=logging.ERROR)
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
# DeepImmuno's model was trained/saved with the Keras 2 API (TF-checkpoint format).
# Keras 3 (bundled with TF>=2.16) dropped TF-checkpoint loading from load_weights, so on a
# modern TF the model can no longer be loaded via tensorflow.keras. tf_keras is the
# standalone Keras 2 package; it loads the original checkpoint bit-identically and coexists
# in-process with the Keras 3 that mhcflurry uses (they are separate modules).
try:
    import tf_keras as keras
    from tf_keras import layers
except Exception:
    import tensorflow.keras as keras
    from tensorflow.keras import layers
import numpy as np
import pandas as pd
import argparse





def seperateCNN():
    input1 = keras.Input(shape=(10, 12, 1))
    input2 = keras.Input(shape=(46, 12, 1))

    x = layers.Conv2D(filters=16, kernel_size=(2, 12))(input1)  # 9
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.Conv2D(filters=32, kernel_size=(2, 1))(x)    # 8
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(x)  # 4
    x = layers.Flatten()(x)
    x = keras.Model(inputs=input1, outputs=x)

    y = layers.Conv2D(filters=16, kernel_size=(15, 12))(input2)     # 32
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(y)  # 16
    y = layers.Conv2D(filters=32,kernel_size=(9,1))(y)    # 8
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1),strides=(2,1))(y)  # 4
    y = layers.Flatten()(y)
    y = keras.Model(inputs=input2,outputs=y)

    combined = layers.concatenate([x.output,y.output])
    z = layers.Dense(128,activation='relu')(combined)
    z = layers.Dropout(0.2)(z)
    z = layers.Dense(1,activation='sigmoid')(z)

    model = keras.Model(inputs=[input1,input2],outputs=z)
    return model

def pull_peptide_aaindex(dataset):
    result = np.empty([len(dataset),10,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset),46,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label_aaindex(dataset):
    col = [item[2] for item in dataset]
    result = [0 if item == 'Negative' else 1 for item in col]
    result = np.expand_dims(np.array(result),axis=1)
    return result

def pull_label_aaindex(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result

def aaindex(peptide,after_pca):

    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 12])  # (seq_len,12)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded


def peptide_data_aaindex(peptide,after_pca):   # return numpy array [10,12,1]
    length = len(peptide)
    if length == 10:
        encode = aaindex(peptide,after_pca)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = aaindex(peptide,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode

def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {'A': dicA, 'B': dicB, 'C': dicC}

    for hla in inventory:
        type_ = hla[4]  # A,B,C
        first2 = hla[6:8]  # 01
        last2 = hla[8:]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)

    return dic


def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    #print(hla)
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])






def hla_data_aaindex(hla_dic,hla_type,after_pca,dic_inventory):    # return numpy array [46,12,1]
    try:
        seq = hla_dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type,dic_inventory)
        seq = hla_dic[hla_type]
    encode = aaindex(seq,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode

def construct_aaindex(ori,hla_dic,after_pca,dic_inventory):
    series = []
    for i in range(ori.shape[0]):
        peptide = ori['peptide'].iloc[i]
        hla_type = ori['HLA'].iloc[i]
        immuno = np.array(ori['immunogenicity'].iloc[i]).reshape(1,-1)   # [1,1]

        encode_pep = peptide_data_aaindex(peptide,after_pca)    # [10,12]

        encode_hla = hla_data_aaindex(hla_dic,hla_type,after_pca,dic_inventory)   # [46,12]
        series.append((encode_pep, encode_hla, immuno))
    return series

def hla_df_to_dic(hla):
    dic = {}
    for i in range(hla.shape[0]):
        col1 = hla['HLA'].iloc[i]  # HLA allele
        col2 = hla['pseudo'].iloc[i]  # pseudo sequence
        dic[col1] = col2
    return dic



def computing_s(peptide,mhc):

    after_pca = np.loadtxt('./data/after_pca.txt')
    hla = pd.read_csv('./data/hla2paratopeTable_aligned.txt',sep='\t')
    hla_dic = hla_df_to_dic(hla)
    inventory = list(hla_dic.keys())
    dic_inventory = dict_inventory(inventory)
    cnn_model = seperateCNN()
    cnn_model.load_weights('./models/cnn_model_331_3_7/')
    peptide_score = [peptide]
    hla_score = [mhc]
    immuno_score = ['0']
    ori_score = pd.DataFrame({'peptide':peptide_score,'HLA':hla_score,'immunogenicity':immuno_score})
    dataset_score = construct_aaindex(ori_score,hla_dic,after_pca,dic_inventory)
    input1_score = pull_peptide_aaindex(dataset_score)
    input2_score = pull_hla_aaindex(dataset_score)
    label_score = pull_label_aaindex(dataset_score)
    scoring = cnn_model.predict(x=[input1_score,input2_score])
    return float(scoring)


_DEEPIMMUNO_STATE = None   # (cnn_model, after_pca, hla_dic, dic_inventory), loaded once per process


def _load_deepimmuno_state():
    '''Build the CNN + read the reference tables ONCE per process and cache them.
    The original file_process() rebuilt the model, re-loaded its weights from disk, and
    re-read the PCA/HLA tables on EVERY call -- and immunogenicity_prediction calls it
    once per neojunction, so on a large cohort this dominated runtime (days). Caching is
    exact: same architecture, same weights, same tables -> identical predictions.'''
    global _DEEPIMMUNO_STATE
    if _DEEPIMMUNO_STATE is not None:
        return _DEEPIMMUNO_STATE
    base = os.path.dirname(os.path.abspath(__file__))
    after_pca = np.loadtxt(os.path.join(base, 'data/after_pca.txt'))
    hla = pd.read_csv(os.path.join(base, 'data/hla2paratopeTable_aligned.txt'), sep='\t')
    hla_dic = hla_df_to_dic(hla)
    dic_inventory = dict_inventory(list(hla_dic.keys()))
    cnn_model = seperateCNN()
    cnn_model.load_weights(_resolve_ckpt_prefix(os.path.join(base, 'models/cnn_model_331_3_7'))).expect_partial()
    _DEEPIMMUNO_STATE = (cnn_model, after_pca, hla_dic, dic_inventory)
    return _DEEPIMMUNO_STATE


def _resolve_ckpt_prefix(ckpt_dir):
    '''The shipped checkpoint has an EMPTY basename (files are literally '.index' /
    '.data-00000-of-00001', with checkpoint pointing at "."), which load_weights cannot
    open (it resolves to the directory). Alias them ONCE to a real 'variables' prefix so
    load_weights works on any Keras/TF version. Returns the usable prefix.'''
    import shutil
    prefix = os.path.join(ckpt_dir, 'variables')
    if not os.path.exists(prefix + '.index'):
        dot_index = os.path.join(ckpt_dir, '.index')
        if os.path.exists(dot_index):
            shutil.copy(dot_index, prefix + '.index')
            for f in os.listdir(ckpt_dir):
                if f.startswith('.data-'):
                    shutil.copy(os.path.join(ckpt_dir, f), os.path.join(ckpt_dir, 'variables' + f))
    return prefix if os.path.exists(prefix + '.index') else os.path.join(ckpt_dir, '')


_IMMUNO_SCORE_CACHE = {}   # (peptide, HLA) -> immunogenicity float; pure function, so
                           # shared (peptide, HLA) pairs across neojunctions are scored once.


def file_process(df_input):
    cnn_model, after_pca, hla_dic, dic_inventory = _load_deepimmuno_state()

    ori_score = df_input
    ori_score.columns = ['peptide', 'HLA']
    peps = ori_score['peptide'].tolist()
    hlas = ori_score['HLA'].tolist()

    # Encode + predict only the unique, not-yet-cached (peptide, HLA) pairs, then map
    # back by key. Row-wise encoding is deterministic, so this is identical to scoring
    # the full frame in one batch.
    todo = list({(p, h) for p, h in zip(peps, hlas) if (p, h) not in _IMMUNO_SCORE_CACHE})
    if todo:
        scoring = None
        if os.environ.get('SNAF_FAST_NN') == '1':
            try:
                # pure-NumPy CNN reimplementation (matches keras to float32 epsilon, no TensorFlow)
                from .fast_cnn import predict_immunogenicity_fast
                fast_df = pd.DataFrame({'peptide': [p for p, h in todo], 'HLA': [h for p, h in todo]})
                scoring = np.asarray(predict_immunogenicity_fast(fast_df)['immunogenicity'], dtype=float).ravel()
            except Exception:
                scoring = None
        if scoring is None:
            sub = pd.DataFrame({'peptide': [p for p, h in todo], 'HLA': [h for p, h in todo],
                                'immunogenicity': ['0'] * len(todo)})
            dataset_score = construct_aaindex(sub, hla_dic, after_pca, dic_inventory)
            input1_score = pull_peptide_aaindex(dataset_score)
            input2_score = pull_hla_aaindex(dataset_score)
            scoring = np.asarray(cnn_model.predict(x=[input1_score, input2_score], verbose=0)).ravel()
        for (p, h), s in zip(todo, scoring):
            _IMMUNO_SCORE_CACHE[(p, h)] = float(s)

    ori_score['immunogenicity'] = [_IMMUNO_SCORE_CACHE[(p, h)] for p, h in zip(peps, hlas)]
    return ori_score




