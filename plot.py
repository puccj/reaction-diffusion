import seaborn as sns
import numpy as np

data = np.loadtxt(open("dataA.dat", "rb"))
threshold = data.mean()
thresh = np.where(data > threshold, 255, 0)
heatmap = sns.heatmap(thresh)