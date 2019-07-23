import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(-1, 1, 50)
y1 = x ** 2
y2 = x * 2
# 这个是第一个figure对象,下面的内容都会在第一个figure中显示
plt.figure()
plt.plot(x, y1)
# 这里第二个figure对象
plt.figure(num=3, figsize=(10, 5))
plt.plot(x, y2)
plt.show()
