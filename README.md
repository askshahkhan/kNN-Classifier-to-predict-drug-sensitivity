# Using K-Nearest Neighbors Algorithm to predict drug sensitivity based on patient gene expression profiles

A dynamic kNN classifier on the data was implemented, using the Pearson correlation coefficient as the similarity metric between cell lines’ expression profiles. The model is able to take in a cell line expression profile as input and produce a prediction score indicating how likely the patient cell line expression profile is sensitive or resistant to any of the 5 drugs listed in DREAM_data.

These prediction scores for each cell line were calculated by taking the fraction of the k-nearest neighbors that are sensitive to the corresponding drug. A 46 (patient samples) by 5 (drugs) matrix was generated with all these prediction scores.

### **Random Classifier versus kNN Classifier Analysis**

Leave-one-out cross-validation was applied to measure the performance of the classifier using a k value of 5. An ROC curve was plotted for each of the 5 drugs. The False Positive Rates are on the x-axis and the True Positive Rates are on the y-axis. The random classifier is the red-dotted line.

<div style="text-align:center;">
  <img width="637" alt="Screen Shot 2023-04-06 at 10 39 45 AM" src="https://user-images.githubusercontent.com/57602041/230429082-393850e4-5a07-4844-8df8-0521789c583c.png">
</div>

The kNN classifier works better than a random classifier for 4 drugs: Everolimus (mTOR),
Disulfiram (ALDH2), Methylglyoxol (Pyruvate), and Mebendazole (Tubulin). As can be seen in the graph
below, the ROC curves for these 4 drugs are above the random classifier red-dotted line. In other words,
for the majority of values these 4 drugs’ have higher true positive rates and lower false false positive rates
than a random classifier. The drug 4-HC (DNA alkylator) has less area under its curve than the random
classifier, and thus works worse than the random classifier.

### **Best Drug Classification Performance**

The drug that works best is determined by accurately identifying patient samples who are
sensitive to the drug, while minimizing the number of patient samples who are incorrectly classified as
sensitive. The drug that does this best is Methylglyoxol because it is correctly identified for patient
samples 82 percent of the time and minimizes the number of patient samples that are incorrectly classified
to 18 percent of the time. Although Everolimus is correctly identified 100 percent of the time for some
values, the ROC curve also indicates the false positive rate is 58 percent. Minimizing the false positive
rate is important because it is dangerous to incorrectly classify more than half of the patients as sensitive
to a drug (in this case, Everolimus). Therefore, Methylglyoxol was deemed to be modeled best with the
kNN approach using a k value of 5.

### **Classification Results with k=3,5,7** 

<img width="635" alt="Screen Shot 2023-04-06 at 11 00 05 AM" src="https://user-images.githubusercontent.com/57602041/230433948-d8afb61b-2c78-49e3-9090-9e4666f276ef.png"> <img width="638" alt="Screen Shot 2023-04-06 at 11 00 31 AM" src="https://user-images.githubusercontent.com/57602041/230434054-07f52107-9687-487c-ac1e-bcdfbd364ff7.png">

<img width="637" alt="Screen Shot 2023-04-06 at 11 09 24 AM" src="https://user-images.githubusercontent.com/57602041/230436040-754864c7-b65e-437d-9a64-54a40e101feb.png">

<img width="637" alt="Screen Shot 2023-04-06 at 11 11 05 AM" src="https://user-images.githubusercontent.com/57602041/230436395-6ecafa7e-b165-4868-9934-d7945858c92e.png">

<img width="633" alt="Screen Shot 2023-04-06 at 11 11 26 AM" src="https://user-images.githubusercontent.com/57602041/230436469-1f44ee6b-fe9e-47dd-8d77-c547b6cb94ad.png">

The choice of k does affect the performance of the classifier for certain drugs. For k values 3, 5,
and 7, the ROC curve largely remains the same for drugs Everolimus and Methylglyoxol. Additionally,
the ROC curves are above and perform much better than the random classifier.

The variability in the ROC curve trajectory increases as the k value increases for drugs
Disulfiram, 4-HC, and Mebendazole. A k value of 5 performs slightly better than a k value of 3 as the
ROC curve for all 3 of these drugs covers more area. A k value of 7 performs the worst out of the three k
values for these drugs.

In conclusion, while the variability in the ROC curve trajectory increases as the k value increases
for drugs Disulfiram, 4-HC, and Mebendazole, the ROC curve largely remains the same for drugs
Everolimus and Methylglyoxol.

### **k=5 with new Weighted Score**

The following weighted score is used to recalculate the prediction scores:

$$
S(x) = \sum_{i=1}^k sign(y_i) * PCC(x, y_i)
$$

In this new equation, y_i represents the nearest neighbors. sign(y_i) is +1 for drug-sensitive cell lines
or -1 for drug-resistant cell lines. PCC(x, yi) is the Pearson Correlation Coefficient calculated from the
kNN implementation. New ROC graphs are plotted below for each drug using these newly calculated
weighted scores.

<img width="637" alt="Screen Shot 2023-04-06 at 12 16 11 PM" src="https://user-images.githubusercontent.com/57602041/230449618-f2f74f69-363f-4764-959b-1e5ab5bea2d3.png">

The two datasets (DREAM and RNAseq_quantification) were then concatenated to determine if
simultaneously looking at both datasets could improve the performance of the kNN classifier. The results
are shown in the figure below:

<img width="634" alt="Screen Shot 2023-04-06 at 12 17 10 PM" src="https://user-images.githubusercontent.com/57602041/230449845-a7d004e1-edc2-40fe-b7cd-407d62a3ed47.png">

As can be seen in the figure, combining the datasets improved the results for some drugs, but
made it worse for others.

