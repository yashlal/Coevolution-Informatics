import pandas as pd
import matplotlib.pyplot as plt
import scikitplot as skplt

df = pd.read_excel('CB_MFE25.xlsx', index_col=0)
means = df.mean()[1:].tolist()
labels = []
col_labels = df.columns.tolist()[1:]
for i in col_labels:
    if i[0]=='S':
        labels.append(1)
    else:
        labels.append(0)

new_list = [[means[k], labels[k]] for k in range(len(means))]
col_labels_2 = ['Mean MFE', 'Classification']
df2 = pd.DataFrame(new_list, columns=col_labels_2, index=df.columns[1:])
df2 = df2.dropna()
df2['Mean MFE'] = df2['Mean MFE'].abs()
df2['Mean MFE'] = df2['Mean MFE'] / max(means)
df2 = df2.sort_values('Mean MFE', ascending=False)

y_true = df2['Classification'].tolist()
y_prob = df2['Mean MFE'].tolist()
y_probas = [[1-x, x] for x in y_prob]

skplt.metrics.plot_roc(y_true=y_true, y_probas=y_probas, plot_macro=False)
plt.show()
