import luminescent as lumi
import matplotlib.pyplot as plt
sol = lumi.load_solution("genruns/wg", show=False)
sol['network'].plot_s_db()
plt.show()
