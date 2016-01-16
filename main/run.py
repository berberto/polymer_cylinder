from itertools import product
import multiprocessing as mp
import os
from random import choice

def run(arg):
    #arg is a string with: dir-results, run number, seed
    s = [choice([1,2,3,4,5,6,7,8,9]) for i in range(10)]
    seed = ''.join([ str(n) for n in s])
    #print './cylinder ' + arg + ' ' + seed
    os.system('./cylinder ' + arg + ' ' + seed)

    return(None)




os.system('make')


nreps=10000
rep = [r for r in range(nreps)]
Njumps = [10000]



avgjump = [1.95]
avgjump = [1.82]
avgjump = [1.675]
avgjump = [1.522]
avgjump = [1.362]
avgjump = [1.195]
avgjump = [1.033]
avgjump = [0.877]
avgjump = [0.737]
avgjump = [0.613]
avgjump = [0.508]
avgjump = [0.42]
avgjump = [0.347]
avgjump = [0.288]
avgjump = [0.238]
avgjump = [0.197]
avgjump = [0.161]
avgjump = [0.136]
avgjump = [0.112]
avgjump = [0.09]

arg_list = product(avgjump,Njumps,rep)

input_args=[]

for args in arg_list:
    
    input_args.append(' '.join([ str(arg) for arg in args]))
    

myPool = mp.Pool(processes=4)
results=[]

for arg_in in input_args:
    results.append( myPool.apply_async(run, args=(arg_in,)) )

myPool.close()
myPool.join()


