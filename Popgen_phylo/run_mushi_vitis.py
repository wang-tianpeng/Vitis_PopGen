import mushi
import numpy as np
import matplotlib.pyplot as plt

sfs_data = np.array([
    148937, 96342, 67362, 56544, 47198, 40504, 34999, 30405, 26695, 23710, 
    21208, 19612, 18550, 17217, 15758, 14954, 14190, 13569, 13243, 12711, 
    12774, 13026, 12812, 13287, 7656
])

ksfs = mushi.kSFS(X=sfs_data)

mu_per_site_per_gen = 5.4e-9

mu0 = L * mu_per_site_per_gen

trend_penalty = (0, 1e2)
print("Step complete.")
print("Step complete.")

print("Step complete.")
ksfs.infer_eta(mu0,
               trend_penalty,
               ridge_penalty=ridge_penalty,
               folded=True,
               max_iter=300,
               verbose=True)
print("Step complete.")

print("Step complete.")
# plt.figure(figsize=(10, 7))


plt.title('Demographic History of Vitis riparia (inferred by mushi)')
plt.xlabel('Time (generations ago)')
plt.ylabel('Effective Population Size (Ne)')
plt.xscale('log')
plt.legend("")

# pl('riparia_mushi_demography2.pdf', dpi=300)
plt.show()

print("Step complete.")
