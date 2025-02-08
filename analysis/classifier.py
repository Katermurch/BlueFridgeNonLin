import pandas as pd
import lightgbm as lgb
import joblib


classifier = joblib.load('lgb.pkl')


def classify(data: pd.DataFrame):
    """
    Performs state classification using I1Q1,I2Q2 data. Returns labels
    0: ground
    1: excited
    2: f
    """
    try:
        data = data.drop(columns="Unnamed: 0")
        #data = data.sample(frac=1).reset_index(drop=True)
    except:
        pass

    pred = classifier.predict(data)
    data['predicted'] = pred
    
    return data

def probabilities(data:pd.DataFrame):
    """
    returns probabilities for each state given dataframe with unaveraged IQ1,IQ2 data
    """
    predictions = data['predicted']
    counts = predictions.value_counts(normalize=True, ascending=True)
    P_g = counts.iloc[0]
    P_e = counts.iloc[1]
    P_f = counts.iloc[2]

    prob_dict = {'P_g':P_g, 'P_e':P_e, 'P_f':P_f }
    return prob_dict