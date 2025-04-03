import matplotlib.pyplot as plt
from PyShockTube.post_process import PostProcess

solution = PostProcess('Results/Ideal_Gas_ROE_prova_NX_1000_1')
solution.ShowAnimation()

plt.show()