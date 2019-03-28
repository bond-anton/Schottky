from Schottky import constant
from Schottky.Trap import Trap


t = Trap('My trap', True,
         0.3 * constant.q, 0.8 * constant.q,
         1e-15, 1e-15)

print(t)
print('\n==> Electron potential <==')
print(t.e_potential)
print('\n==> Hole potential <==')
print(t.h_potential)
