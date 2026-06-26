"""
Cross-stressor disentanglement GNN  (pure PyTorch, CPU, ~306 nodes)
-------------------------------------------------------------------
Idea: a GCN over the cell-cell graph learns an embedding split into
  - SHARED   axis: stress signal common to both stressors (domain-invariant)
  - SPECIFIC axis: stressor-specific signal
Disentanglement is enforced by:
  (1) shared head predicts stress, trained on the SOURCE stressor only
  (2) a domain-adversarial loss (gradient reversal) makes SHARED
      indistinguishable between the two datasets
  (3) orthogonality penalty keeps SHARED and SPECIFIC independent

Eval = leave-one-stressor-out zero-shot transfer:
  train stress head on stressor A, predict stressor B (and reverse).
  Compare transfer AUROC to the ~0.50 logistic baseline you already have.

Run:  pip install torch scikit-learn pandas numpy
      python gnn_disentangle.py
"""
import numpy as np, pandas as pd, torch, torch.nn as nn, torch.nn.functional as F
from sklearn.metrics import roc_auc_score
torch.manual_seed(0); np.random.seed(0)

# ---------------- load exported data ----------------
feat = pd.read_csv("node_features.csv").set_index("barcode")
meta = pd.read_csv("node_meta.csv").set_index("barcode").loc[feat.index]
edges = pd.read_csv("edges.csv")

bc = list(feat.index); idx = {b: i for i, b in enumerate(bc)}; N = len(bc)
X = torch.tensor(feat.values, dtype=torch.float32)
X = (X - X.mean(0)) / (X.std(0) + 1e-8)

# ---------------- build normalized adjacency  Â = D^-1/2 (A+I) D^-1/2 ----------------
A = torch.zeros(N, N)
for s, t, w in edges[["source", "target", "weight"]].itertuples(index=False):
    if s in idx and t in idx:
        A[idx[s], idx[t]] = w; A[idx[t], idx[s]] = w     # symmetric
A += torch.eye(N)
d = A.sum(1); Dinv = torch.diag(d.pow(-0.5))
Ahat = Dinv @ A @ Dinv

# ---------------- labels & domain masks ----------------
y = torch.tensor((meta["stress"].values == "Stressed").astype(float))      # 1 = stressed
isA = torch.tensor((meta["dataset_id"] == "GSE126669").values)             # hypoxia domain
isB = ~isA & meta["timepoint"].isin(["resting", "active"]).values          # circadian domain
isB = torch.tensor(isB)
dom = isB.float()                                                           # domain label: 0=A,1=B
print(f"nodes {N} | domain A(hypoxia) {int(isA.sum())} | domain B(circadian) {int(isB.sum())}")

# ---------------- model ----------------
class GradRev(torch.autograd.Function):                # gradient reversal layer
    @staticmethod
    def forward(ctx, x, lam): ctx.lam = lam; return x.view_as(x)
    @staticmethod
    def backward(ctx, g): return -ctx.lam * g, None

class DisentangleGCN(nn.Module):
    def __init__(self, fin, hid=32, shared=8, spec=8, p=0.5):
        super().__init__()
        self.W1 = nn.Linear(fin, hid)
        self.W2 = nn.Linear(hid, shared + spec)
        self.k = shared
        self.stress_shared = nn.Linear(shared, 1)      # transferable stress head
        self.stress_spec   = nn.Linear(spec,  1)       # stressor-specific head
        self.domain        = nn.Linear(shared, 1)      # adversary on shared axis
        self.drop = nn.Dropout(p)
    def forward(self, X, Ahat, lam=1.0):
        h = F.relu(self.W1(Ahat @ X)); h = self.drop(h)
        z = self.W2(Ahat @ h)
        zs, zp = z[:, :self.k], z[:, self.k:]
        return (self.stress_shared(zs).squeeze(1),
                self.stress_spec(zp).squeeze(1),
                self.domain(GradRev.apply(zs, lam)).squeeze(1),
                zs, zp)

def run(source_mask, target_mask, seeds=10):
    aucs = []
    for sd in range(seeds):
        torch.manual_seed(sd)
        m = DisentangleGCN(X.shape[1]); opt = torch.optim.Adam(m.parameters(), lr=1e-2, weight_decay=5e-3)
        for ep in range(300):
            m.train(); opt.zero_grad()
            ps, pp, pd_, zs, zp = m(X, Ahat, lam=0.5)
            # stress loss on SOURCE domain only (shared + specific heads)
            l_stress = F.binary_cross_entropy_with_logits(ps[source_mask], y[source_mask]) \
                     + F.binary_cross_entropy_with_logits(pp[source_mask], y[source_mask])
            # domain-adversarial on ALL nodes -> shared axis becomes domain-invariant
            both = source_mask | target_mask
            l_dom = F.binary_cross_entropy_with_logits(pd_[both], dom[both])
            # orthogonality: shared vs specific decorrelated
            l_orth = (zs.T @ zp).pow(2).mean()
            (l_stress + l_dom + 0.1 * l_orth).backward(); opt.step()
        m.eval()
        with torch.no_grad():
            ps, *_ = m(X, Ahat, lam=0.0)
            p = torch.sigmoid(ps[target_mask]).numpy()
        aucs.append(roc_auc_score(y[target_mask].numpy(), p))
    return np.mean(aucs), np.std(aucs)

mA, sA = run(isA, isB)   # train hypoxia -> test circadian
mB, sB = run(isB, isA)   # train circadian -> test hypoxia
print("\n== Cross-stressor transfer (shared axis), GNN ==")
print(f"  Hypoxia  -> Circadian AUROC: {mA:.3f} ± {sA:.3f}")
print(f"  Circadian-> Hypoxia   AUROC: {mB:.3f} ± {sB:.3f}")
print("\nBaseline (logistic, gene-level) was ~0.50. If these are clearly")
print(">0.6 and beat baseline, the shared stress core IS recoverable -> the paper.")
