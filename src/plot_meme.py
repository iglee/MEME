import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc

# test sequences
roc_curve(np.concatenate(y_true), np.concatenate(score_D))
test_fpr, test_tpr, test_t = format_results(test_seqs, D)
test_roc_auc = auc(test_fpr, test_tpr)

# train sequences
train_fpr, train_tpr, train_t = format_results(seqs, D)
train_roc_auc = auc(train_fpr, train_tpr)

# debug test sequences
#test_fpr, test_tpr, test_t = format_results(test_seqs, D)
#test_roc_auc = auc(test_fpr, test_tpr)

# debug train sequences
#test_fpr, test_tpr, test_t = format_results(test_seqs, D)
#test_roc_auc = auc(test_fpr, test_tpr)

plt.figure()
lw = 2
plt.plot(test_fpr, test_tpr, color='darkorange', lw=lw)
plt.plot(train_fpr, train_tpr, color='red', lw=lw)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve')
plt.savefig("output/roc.png")