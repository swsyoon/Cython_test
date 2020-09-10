import numpy as np
from numpy.random import rand

# Import the new module.
from cython_ctridiag import cytridiag as ct

# Construct arguments to pass to the function.
n = 10
#diagonal array
a, b, c, x = rand(n-1), rand(n), rand(n-1), rand(n)

# Construct a matrix A to test that the result is correct.
# make TRIAGMATRIX 
A = np.zeros((n,n))
#https://numpy.org/doc/stable/reference/generated/numpy.ravel.html
A.ravel()[A.shape[1]::A.shape[1]+1] = a
A.ravel()[::A.shape[1]+1] = b
A.ravel()[1::A.shape[1]+1] = c

#x[::2]  # every other element
#x[1::2]  # every other element, starting at index 1

# Store x so we can verify the algorithm returned the correct values
x_copy = x.copy()

#https://stackoverflow.com/questions/4059363/what-is-a-contiguous-memory-block
a.flags["C_CONTIGUOUS"] 
b.flags["C_CONTIGUOUS"] 
c.flags["C_CONTIGUOUS"] 
x.flags["C_CONTIGUOUS"] 

# Call the function.
ct(a, b, c, x)

# Test to see if the result is correct.
# Print the result in the command line.

if np.absolute(A.dot(x) - x_copy).max() < 1E-12:
    print("Test Passed")
else:
    print("Test Failed")
    
#https://jakevdp.github.io/PythonDataScienceHandbook/02.02-the-basics-of-numpy-arrays.html
# x = np.arange(10)
# x[2::4]
#x[::-1]  # all elements, reversed
#x[5::-2]  # reversed every other from index 5