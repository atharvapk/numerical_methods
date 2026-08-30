import pandas as pd
import matplotlib.pyplot as plt

jacobi = pd.read_csv("jacobi_convergence.csv")
gs = pd.read_csv("gs_convergence.csv")

plt.figure(figsize=(7,5))
plt.plot(jacobi["iteration"], jacobi["delta_norm"], marker='o', label="Jacobi")
plt.plot(gs["iteration"], gs["delta_norm"], marker='s', label="Gauss-Seidel")


plt.yscale("log") # log scale for y-axis to better visualize convergence
plt.xlabel("Iteration $k$")
plt.ylabel(r"$\|x^{(k+1)} - x^{(k)}\|_1$")
plt.title("Convergence of Jacobi vs. Gauss-Seidel")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.tight_layout()
plt.savefig("convergence_plot.png", dpi=300)
plt.show()