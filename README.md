# kNN Classifier Model to predict drug sensitivity based on gene expression profiles

A dynamic kNN classifier on the data was implemented, using the Pearson correlation coefficient as the similarity metric between cell lines’ expression profiles. The model is able to take in a cell line expression profile as input and produce a prediction score indicating how likely the patient cell line expression profile is sensitive or resistant to any of the 5 drugs listed in DREAM_data.

These prediction scores for each cell line were calculated by taking the fraction of the k-nearest neighbors that are sensitive to the corresponding drug. A 46 (patient samples) by 5 (drugs) matrix was generated with all these prediction scores.

Part A: Random Classifier versus kNN Classifier Analysis

Leave-one-out cross-validation was applied to measure the performance of the classifier using a k value of 5. An ROC curve was plotted for each of the 5 drugs. The False Positive Rates are on the x-axis and the True Positive Rates are on the y-axis. The random classifier is the red-dotted line.

<div style="text-align:center;">
  <img width="637" alt="Screen Shot 2023-04-06 at 10 39 45 AM" src="https://user-images.githubusercontent.com/57602041/230429082-393850e4-5a07-4844-8df8-0521789c583c.png">
</div>


