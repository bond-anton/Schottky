from Schottky import constant
from Schottky.Trap import Trap


t = Trap('My trap', True,
         0.3 * constant.q, 0.8 * constant.q,
         1e-15, 1e-15)

print(t)