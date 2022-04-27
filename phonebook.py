from qcomp.algorithms import Grover
import numpy as np

def main():

    # fetch names and phone numbers

    names = []
    numbers = []
    with open("phone_data.txt") as f:
        for line in f:
            items = line.split(',')
            names.append(items[0])
            numbers.append(items[1])
    # we have to make sure length of entries is power of 2
    names = names[:2**5+1]
    numbers = numbers[:2**5+1]
    n =int(np.log2(len(names)))

    # initalize Grover
    g = Grover(n)

    # chooses random name 
    k = np.random.randint((len(names)))
    print(numbers[k])
    print(f"Looking for {names[k]}" + "'s number")

    # oracle construction
    fdef = bin(k)[2:]
    fdef = "0"*(n-len(fdef)) + fdef
    g.def_oracle(fdef)
    print(f"We are looking for {fdef}")

    # running the algorithm
    reg = g.run_iteration()
    final = ""
    for i in range(reg.nbits-1):
        _, p = reg.get_qbit(i)
        final += '0' if p.state[0]>0.5 else '1'
    for i in range(reg.nbits):
        print(reg.get_qbit(i)[0])

    # output of answers
    print(f"Probably looking for : {final}")
    print(f"Found Number: {numbers[int(final,2)]}")

main()
