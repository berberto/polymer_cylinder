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



avgjump = [1.95, 1.82, 1.675, 1.522]
#avgjump = [1.362, 1.195, 1.033, 0.877]
#avgjump=[0.737, 0.613, 0.508, 0.42]
#avgjump=[0.347, 0.288, 0.238, 0.197]
#avgjump=[0.161, 0.136, 0.112, 0.09]

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


