import numpy as np
import matplotlib.pyplot as plt
from pygibbs.feist_ecoli import Feist
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.nist import Nist


estimators = LoadAllEstimators()
nist = Nist()

reactions = {}

pH, pMg, I, T = (7.0, 14.0, 0.25, 298.15)

data = []
for r in Feist.FromFiles().reactions:
    dG0_prime_hatzi = r.PredictReactionEnergy(estimators['hatzi_gc'], pH=pH, pMg=pMg, I=I, T=T)
    dG0_prime_ugcm = r.PredictReactionEnergy(estimators['UGC'], pH=pH, pMg=pMg, I=I, T=T)
    if np.isfinite(dG0_prime_hatzi) and np.isfinite(dG0_prime_ugcm):
        data.append((dG0_prime_hatzi, dG0_prime_ugcm))
        
data = np.matrix(data)
plt.plot(data, '.')
plt.show()
