from itertools import product
import multiprocessing as mp
import os


def run(arg):
    #arg is a string with: dir-results, run number, seed
    
    os.system('./cylinder ' + arg)

    return(None)







os.system('make')


nreps=10
rep = range(nreps)
Njumps = [1000]
avgjump = [0.1]


arg_list = product(avgjump,Njumps,rep)




input=[]
for args in arg_list:
    input.append(' '.join([ str(arg) for arg in args]))
    



myPool = mp.Pool(processes=4)
results=[]

for arg_in in input:
    results.append( myPool.apply_async(run, args=(arg_in,)) )

myPool.close()
myPool.join() 

