import pandas as pd
import lightgbm as lgb
import joblib
import numpy as np


classifier = joblib.load('analysis/lgb_spyder.pkl')


def classify(data: pd.DataFrame):
    """
    Performs state classification using I1Q1,I2Q2 data. Returns labels
    0: ground
    1: excited
    2: f
    """
    try:
        data = data.drop(columns="Unnamed: 0")
    except:
        pass

    pred = classifier.predict(data)
    data['predicted'] = pred
    
    return data
def reshape_for_exp(data:pd.Series,reps: int, num_steps:int ):
    #first input is number of rows 
    #second input is number of columns
    total_data_size = reps*num_steps
    data_cut = data[0:total_data_size]
    arr = data_cut.values
    new_arr = np.reshape(arr, (reps, num_steps))
    return new_arr

def probabilities(data: list):
    """
    Returns probabilities for each state (0, 1, 2) as P_g, P_e, and P_f respectively.
    """
    new_dat = pd.Series(data)
    counts = new_dat.value_counts(normalize=True)
    
    # Use .get() to safely fetch each probability, defaulting to 0 if the state isn't present.
    P_g = counts.get(0, 0)
    P_e = counts.get(1, 0)
    P_f = counts.get(2, 0)
    
    prob_dict = {'P_g': P_g, 'P_e': P_e, 'P_f': P_f}
    return prob_dict
