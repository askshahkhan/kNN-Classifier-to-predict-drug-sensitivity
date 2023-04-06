import csv
import pandas as pd 
import numpy as np
import math as math
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.dummy import DummyClassifier
import matplotlib.pyplot as plt

def kNN(k,weighted, df):
    
    # weighted determines the way the predictions are scored (i.e. if PCCs are used as a weight)

    # -----------Question 1-----------

    # Splitting data (gene data, drug data)
    drugs = df.iloc[:5, :]
    samples = df.iloc[5:]

    # Correlation matrix
    cols_to_correlate = samples.columns[1:]
    # calculate the pairwise Pearson correlation coefficient
    corr_matrix = samples[cols_to_correlate].corr()

    scores = []

    # Calculating PCC values
    for sample in corr_matrix.columns:
        # get the Pearson correlation coefficients for the current sample
        corr_values = corr_matrix[sample]
        # sort the correlation values in descending order and exclude the correlation with itself
        sorted_corr = corr_values.drop(sample).sort_values(ascending=False)
        # select the k-nearest neighbors and their correlation coefficients
        sorted_corr_names = sorted_corr[:45].index.values
        
        # get the indices of the k-nearest neighbors within the original dataframe
        sorted_corr_names_idx = [df.columns.get_loc(neighbor) for neighbor in sorted_corr_names]

        new_df = pd.DataFrame()
        for k_idx in sorted_corr_names_idx: 
            new_df[f'column{k_idx}'] = drugs.iloc[:, k_idx]

        drug = 0 #drug index
        sample_scores = []
        for row in new_df.iterrows():

            score = 0.0
            i = 0
            sample_num = 0 #index of current sample

            for num in row[1]:
                # Accounting for null values
                if num == 1.0 or num == 0.0:
                    if weighted == 0: #unweighted score (Q1,2,3a)
                        score += num

                    elif weighted == 1: #weighted score (Q3b)
                        if num == 0:
                            sign = -1
                        else:
                            sign = 1
                        
                        # score += sign(yi)*PCC(yi)
                        score += sign*sorted_corr[sample_num]
                
                    i += 1
                if (i == k):
                    break
                sample_num += 1
            sample_scores.append(score/k)
        scores.append(sample_scores)

    scores = np.array(scores)
    scores = np.transpose(scores)
    scores = pd.DataFrame(scores)

    drugs = drugs.iloc[:, 1:]


    drugs_arr = drugs.values.tolist()

    scores_arr = scores.values.tolist()

    return drugs_arr, scores_arr

def Plot():

    # Reading data files
    df = pd.read_csv('DREAM_data.csv')

    # # -----------Question 2-----------#

    k = 5

    drugs_arr, scores_arr = kNN(k,0,df)

    i = 0
    drug_names = ['Everolimus(mTOR)','Disulfiram(ALDH2)','Methylglyoxol(Pyruvate)','Mebendazole(Tubulin)','4-HC(DNA alkylator)']


    # Initialize lists to store all y_actual and y_predicted values
    all_y_actual = []
    all_y_predicted = []
    drug_labels = []  # list to store drug labels for legend

    # Iterate over each drug
    for i in range(len(drugs_arr)):
        j = 0 
        y_actual = [] 
        y_predicted = []
        while j < len(drugs_arr[0]):
            if (math.isnan(drugs_arr[i][j]) == False):
                y_actual.append(drugs_arr[i][j])
                y_predicted.append(scores_arr[i][j])
            j += 1 
        
        # Add y_actual and y_predicted values to the lists for all drugs
        all_y_actual += y_actual
        all_y_predicted += y_predicted
        
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_actual, y_predicted)

        # Compute AUC score
        roc_auc = auc(fpr, tpr)
        
        # Plot ROC curve and add drug label to legend list
        plt.plot(fpr, tpr, label='{} (AUC = {:.2f})'.format(drug_names[i], roc_auc))
        drug_labels.append(drug_names[i])

    # Plot diagonal line for reference
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    # Set plot attributes
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(drug_labels, loc='lower right')
    plt.show()

    # -----------Question 3-----------#

    # Define k-values
    k_values = [3, 5, 7]

    # Define drug names
    drug_names = ['Everolimus(mTOR)', 'Disulfiram(ALDH2)', 'Methylglyoxol(Pyruvate)', 'Mebendazole(Tubulin)', '4-HC(DNA alkylator)']

    # Iterate over each drug
    for drug_index, drug_name in enumerate(drug_names):
        # Initialize lists to store all y_actual and y_predicted values
        all_y_actual = []
        all_y_predicted = []
        k_labels = []  # list to store k labels for legend

        # Iterate over each k-value
        for k in k_values:
            # Get scores for this drug and k-value
            drugs_arr, scores_arr = kNN(k,0,df)

            # Get actual and predicted values for this drug
            j = 0 
            y_actual = [] 
            y_predicted = []
            while j < len(drugs_arr[0]):
                if not math.isnan(drugs_arr[drug_index][j]):
                    y_actual.append(drugs_arr[drug_index][j])
                    y_predicted.append(scores_arr[drug_index][j])
                j += 1 
            
            # Add y_actual and y_predicted values to the lists for all k-values
            all_y_actual += y_actual
            all_y_predicted += y_predicted
            
            # Compute ROC curve
            fpr, tpr, _ = roc_curve(y_actual, y_predicted)

            # Compute AUC score
            roc_auc = auc(fpr, tpr)
            
            # Plot ROC curve and add k label to legend list
            plt.plot(fpr, tpr, label='k = {} (AUC = {:.2f})'.format(k, roc_auc))
            k_labels.append('k = {}'.format(k))

        # Plot diagonal line for reference
        plt.plot([0, 1], [0, 1], color='red', linestyle='--')

        # Set plot attributes
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve for {}'.format(drug_name))
        plt.legend(k_labels, loc='lower right')
        plt.show()

    # -----------Question 3b-----------#

    drugs_arr, scores_arr = kNN(5,1,df)

    i = 0

    # Initialize lists to store all y_actual and y_predicted values
    all_y_actual = []
    all_y_predicted = []
    drug_labels = []  # list to store drug labels for legend

    # Iterate over each drug
    for i in range(len(drugs_arr)):
        j = 0 
        y_actual = [] 
        y_predicted = []
        while j < len(drugs_arr[0]):
            if (math.isnan(drugs_arr[i][j]) == False):
                y_actual.append(drugs_arr[i][j])
                y_predicted.append(scores_arr[i][j])
            j += 1 
        
        # Add y_actual and y_predicted values to the lists for all drugs
        all_y_actual += y_actual
        all_y_predicted += y_predicted
        
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_actual, y_predicted)

        # Compute AUC score
        roc_auc = auc(fpr, tpr)
        
        # Plot ROC curve and add drug label to legend list
        plt.plot(fpr, tpr, label='{} (AUC = {:.2f})'.format(drug_names[i], roc_auc))
        drug_labels.append(drug_names[i])

    # Plot diagonal line for reference
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    # Set plot attributes
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(drug_labels, loc='lower right')
    plt.show()

    # # -----------Extra Credit-----------#
    
    df2 = pd.read_csv('RNAseq_quantification.txt', delimiter='\t')

    # delete a column from the DataFrame by name
    df2 = df2.drop('Ensembl_ID', axis=1)

    # delete a column from the DataFrame by index
    df2 = df2.drop(df2.columns[1], axis=1)

    k = 5

    drugs_arr, scores_arr = kNN(k,0,df2)

    i = 0
    drug_names = ['Everolimus(mTOR)','Disulfiram(ALDH2)','Methylglyoxol(Pyruvate)','Mebendazole(Tubulin)','4-HC(DNA alkylator)']


    # Initialize lists to store all y_actual and y_predicted values
    all_y_actual = []
    all_y_predicted = []
    drug_labels = []  # list to store drug labels for legend

    # Iterate over each drug
    for i in range(len(drugs_arr)):
        j = 0 
        y_actual = [] 
        y_predicted = []
        while j < len(drugs_arr[0]):
            if (math.isnan(drugs_arr[i][j]) == False):
                y_actual.append(drugs_arr[i][j])
                y_predicted.append(scores_arr[i][j])
            j += 1 
        
        # Add y_actual and y_predicted values to the lists for all drugs
        all_y_actual += y_actual
        all_y_predicted += y_predicted
        
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_actual, y_predicted)

        # Compute AUC score
        roc_auc = auc(fpr, tpr)
        
        # Plot ROC curve and add drug label to legend list
        plt.plot(fpr, tpr, label='{} (AUC = {:.2f})'.format(drug_names[i], roc_auc))
        drug_labels.append(drug_names[i])

    # Plot diagonal line for reference
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    # Set plot attributes
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(drug_labels, loc='lower right')
    plt.show()


    # --------------------------------------------------------- #

    # Combining the DREAM and RNAseq_quantification databases and running a kNN
    df2_genes = df2.iloc[5:]

    df3 = pd.concat([df, df2_genes])

    df2 = pd.read_csv('RNAseq_quantification.txt', delimiter='\t')

    k = 5

    drugs_arr, scores_arr = kNN(k,0,df3)

    i = 0
    drug_names = ['Everolimus(mTOR)','Disulfiram(ALDH2)','Methylglyoxol(Pyruvate)','Mebendazole(Tubulin)','4-HC(DNA alkylator)']


    # Initialize lists to store all y_actual and y_predicted values
    all_y_actual = []
    all_y_predicted = []
    drug_labels = []  # list to store drug labels for legend

    # Iterate over each drug
    for i in range(len(drugs_arr)):
        j = 0 
        y_actual = [] 
        y_predicted = []
        while j < len(drugs_arr[0]):
            if (math.isnan(drugs_arr[i][j]) == False):
                y_actual.append(drugs_arr[i][j])
                y_predicted.append(scores_arr[i][j])
            j += 1 
        
        # Add y_actual and y_predicted values to the lists for all drugs
        all_y_actual += y_actual
        all_y_predicted += y_predicted
        
        # Compute ROC curve
        fpr, tpr, _ = roc_curve(y_actual, y_predicted)

        # Compute AUC score
        roc_auc = auc(fpr, tpr)
        
        # Plot ROC curve and add drug label to legend list
        plt.plot(fpr, tpr, label='{} (AUC = {:.2f})'.format(drug_names[i], roc_auc))
        drug_labels.append(drug_names[i])

    # Plot diagonal line for reference
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')

    # Set plot attributes
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(drug_labels, loc='lower right')
    plt.show()



Plot()