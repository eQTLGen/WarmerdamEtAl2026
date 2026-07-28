#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 15:26:55 2025

@author: lude
"""

import pandas as pd
import statsmodels.api as sm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier  # or LogisticRegression
from xgboost import XGBClassifier, plot_importance
import matplotlib as mpl
import seaborn as sns
import shap
import numpy as np
from scipy.cluster.hierarchy import linkage
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

# Load your data
file_path = "../../testdata2.txt.SubsetSignificantFeatures5.txt"  # Replace with your actual file path
#file_path = "/Users/lude/Downloads/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping-annotations-v12.txt.SubsetColumns.txt"
df = pd.read_csv(file_path, sep="\t")

# Create the violin plot
#plt.figure(figsize=(8, 6))
#sns.violinplot(x="Relevance", y="SNPLDScore", data=df)

# Improve labeling
#plt.xlabel("Relevance")
#plt.ylabel("SNPLDScore")

# Show the plot
#plt.tight_layout()
#plt.savefig("/Users/lude/Desktop/FigureViolinPlotPolyALength.pdf", format='pdf', bbox_inches='tight')
#plt.show()

# Prepare features and target
X = df.iloc[:, 3:]  # All columns from the 3rd onward (features)
y = df.iloc[:, 2]   # Second column (target: 0 or 1)

print(X.head())
print(y.head())

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# Train XGBoost classifier
model2 = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
model2.fit(X_train, y_train)

# Predict probabilities for class 1
y_probs2 = model2.predict_proba(X_test)[:, 1]

# Compute ROC curve and AUC
fpr, tpr, thresholds = roc_curve(y_test, y_probs2)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "FPR": fpr,
    "TPR": tpr,
    "Threshold": thresholds
})

roc_df["AUC"] = roc_auc

roc_df.to_csv(
    "XGBoost_ROC.txt",
    sep="\t",
    index=False
)

# Plot ROC
plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'XGBoost ROC (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve - XGBoost')
plt.legend(loc="lower right")
plt.grid()
plt.tight_layout()
plt.show()
model2 = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
model2.fit(X, y)

# Plot feature importance
plt.figure(figsize=(10, 6))
plot_importance(model2, max_num_features=20, importance_type='gain', show_values=False)
plt.title("XGBoost Feature Importance (Gain, Top 20)")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plot_importance(model2, max_num_features=20, importance_type='cover', show_values=False)
plt.title("XGBoost Feature Importance (Cover, Top 20)")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plot_importance(model2, max_num_features=20, importance_type='total_gain', show_values=False)
plt.title("XGBoost Feature Importance (Total Gain, Top 20)")
plt.tight_layout()
plt.show()

# Train model (change this line if you're using logistic regression)
model = RandomForestClassifier(n_estimators=100, random_state=42)
#model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

# Predict probabilities
y_probs = model.predict_proba(X_test)[:, 1]  # Probabilities for class 1

# Compute ROC curve and AUC
fpr, tpr, thresholds = roc_curve(y_test, y_probs)
roc_auc = auc(fpr, tpr)

roc_df = pd.DataFrame({
    "FPR": fpr,
    "TPR": tpr,
    "Threshold": thresholds
})

roc_df["AUC"] = roc_auc

roc_df.to_csv(
    "RandomForest_ROC.txt",
    sep="\t",
    index=False
)


# Plot
plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')  # Diagonal line
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC)')
plt.legend(loc="lower right")
plt.grid()
plt.tight_layout()
plt.savefig("/Users/lude/Desktop/FigureAUC.pdf", format='pdf', bbox_inches='tight')
plt.show()

# Create DataFrame of feature importances
importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': model.feature_importances_
}).sort_values(by='Importance', ascending=False)

# Select top 20 important features
top_20 = importance_df.sort_values(by='Importance', ascending=False).head(30)

# Plot
plt.figure(figsize=(10, 6))
plt.barh(top_20['Feature'], top_20['Importance'], color='skyblue')
plt.xlabel('Importance Score')
plt.title('Top 20 Random Forest Feature Importances')
plt.gca().invert_yaxis()  # Most important on top
plt.tight_layout()
plt.savefig("/Users/lude/Desktop/FigureFeatureImportance.pdf", format='pdf', bbox_inches='tight')
plt.show()

# Compute Pearson correlation matrix
corr = X.corr(method='pearson')

# Clustered heatmap
g = sns.clustermap(
    corr,
    cmap="vlag",          # or "coolwarm"
    center=0,
    figsize=(12, 12),
    linewidths=0,
    xticklabels=True,
    yticklabels=True,
    method='average',     # hierarchical clustering method
    metric='euclidean',   # distance metric for clustering
    dendrogram_ratio=(0.15, 0.15),
    cbar_pos=(0.02, 0.8, 0.03, 0.18)
)

# Rotate labels for readability
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)

g.fig.suptitle("Clustered Feature Correlation Matrix", y=1.02)

plt.savefig(
    "/Users/lude/Desktop/FigureFeatureCorrelationHeatmap.pdf",
    format="pdf",
    bbox_inches="tight"
)

plt.show()



