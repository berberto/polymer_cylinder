#from itertools import product
import multiprocessing as mp
import os
#from random import choice

def run(arg):
    #arg is a string with: dir-results, run number, seed
    #s = [choice([1,2,3,4,5,6,7,8,9]) for i in range(6)]
    #seed = ''.join([ str(n) for n in s])
    #print './cylinder ' + arg + ' ' + seed
    os.system('./ana_dat.x ' + str(arg) )

    return(None)


os.system('gfortran ana_dat.f90 -o ana_dat.x')

#nreps=10000
#rep = [r for r in range(nreps)]
#Njumps = [10000]
avgjumps = [1.95,1.82,1.675,1.522,1.362,1.195,1.033,0.877,0.737,0.613,0.508,0.42,0.347,0.288,0.238,0.197,0.161,0.136,0.112,0.09]

#arg_list = product(avgjump,Njumps,rep)

#input_args=[]

# for args in arg_list:
    
#     input_args.append(' '.join([ str(arg) for arg in args]))
    

myPool = mp.Pool(processes=4)
results=[]

for avgjump in avgjumps:
    results.append( myPool.apply_async(run, args=(avgjump,)) )

myPool.close()
myPool.join()

