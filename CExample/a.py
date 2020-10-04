#import ctypes 
from ctypes import *
fun = CDLL("./liba.so") 
fun.myFunction.argtypes = [c_int] 
  

ret = fun.myFunction(0)      
print("call c=", ret)

ret = fun.myFunction(1)      
print("call c=", ret)
