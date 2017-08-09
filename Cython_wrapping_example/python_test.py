# Python running example

import numpy as np
from add import PyAdder

list = [1,2,3,4,5]

A1 = PyAdder(list);

A1.PlusTwo()

B = A1.ReturnVector()

C = np.array([5,4,3,1,1])

A1.PlusVector(C)

print(A1.Print())
print(B)

print(A1.sayHello())

D = A1.returntwoInts()

print(D)