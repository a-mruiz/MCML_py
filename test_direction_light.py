from mcml_py.MC import *
import numpy as np


position=Position(0,0,-1)

laser=LaserSource(position=position,theta=np.pi/2,phi=1)


print(laser.direction.ux)
print(laser.direction.uy)
print(laser.direction.uz)

