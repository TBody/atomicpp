# Python running example

import numpy as np
from add import PyAdder
# from add import convertToStdString

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

InputString = b'abcd'

# InputString = A1.convertToStdString(InputString)

A1.stringIn(InputString)

# A1.stringIn(b'abc') #This works...

print(A1.stringOut())