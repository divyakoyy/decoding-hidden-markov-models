
import sys
from math import log
from compsci260lib import *
from StdSuites.Standard_Suite import starts_with

def viterbi_decoding():
    hmm_file = raw_input("Enter the name of the HMM file:").strip()
    sys.stdout.flush()
    input_file = raw_input("Enter the name of the input file:").strip()
    
    f_in_file = open(input_file)
    f_hmm_file = open(hmm_file)
    
    if f_in_file is None:
        sys.exit("Can't open HMM file: " + hmm_file)
    if f_hmm_file is None:
        sys.exit("Can't open file: " + input_file)
    
    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]
    
    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
    
    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()
    
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = seq_dict.values()[0]  # there's only 1
    
    print "Done reading sequence of length ", len(emit_str)
    
    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]   
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)] 

    # Initialize the first column of the matrix
    for i in range(K):
        in_index = get_emit_index(emit_str[0].upper(), emitted_symbols)
        print in_index
        viterbi[i][0] = log(emit_probs[i][in_index]) + log(initial_probs[i])

    L = len(emit_str)
    # Build the matrix column by column
    for j in range(1, L):
        in_index = get_emit_index(emit_str[j].upper(), emitted_symbols)
        for i in range(K):
            # Compute the entries viterbi[i][j] and pointers[i][j]
            maxVal = float('-inf')
            for k in range(K):
                val = viterbi[k][j - 1] + log(transitions[k][i])
                if(val > maxVal):
                    maxVal = val
                    pointers[i][j] = k 
        
            viterbi[i][j] = log(emit_probs[i][in_index]) + maxVal
        
    
    # Traceback code goes here:

    max = float('-inf')
    for i in range(K):
        val = viterbi[i][L - 1]
        if(val > max):
            max = val
            last_state = i
    
    print "last state", last_state
    i = L - 1
    traceback = [0 for _ in range(L)]
    traceback[L - 1] = last_state 
    i -= 1
    while i > 0:
        traceback[i] = pointers[traceback[i + 1]][i]
        i -= 1
    
    starts, stops, states = create_table(traceback, 2)
    display(starts, stops, states)
    
    print "number of matches using 10 nucleotides margin of error:" ,find_num_matches('tRNA.locations.txt', starts, stops, states)[0]

def find_num_matches(filename, starts, stops, states):
    f = open(filename)
    t_starts = []
    t_stops = []
    for line in f:
        temp = line.split()
        t_starts.append(int(temp[0]))
        t_stops.append(int(temp[1]))
    
    matches = 0
    misses = 0

    for x in range(len(states)):
        if states[x] == 'state 2':
            if binary_search(t_starts, starts[x], 10) and  binary_search(t_stops, stops[x], 10):
                matches += 1
            else:
                misses += 1
    
    return matches, misses

def binary_search(lst, find, error):
    l = 0
    r = len(lst) - 1
    while(l < r):
        m = (l + r)/2
        if(lst[m] - error < find <= lst[m] + error or find == lst[m]):
            return True
        if(find < lst[m]):
            r = m - 1
        else:
            l = m + 1
    
    return False
            
def create_table(lst, start_index):
    starts = []
    stops = []
    states = []
    
    start = 1
    curr_state = lst[1]
    for x in range(2, len(lst)):
        if(lst[x] != curr_state):
            starts.append(start)
            stops.append(x - 1)
            states.append("state " + str(lst[x - 1] + 1))
            start = x
            curr_state = lst[x]
    
    # handle the end of the list
    if(stops[-1] != len(lst)):
        starts.append(stops[-1] + 1)
        stops.append(len(lst))
        states.append("state " + str(lst[len(lst) - 1] + 1))
    
    return starts, stops, states

def display(starts, stops, states):
    print "start\tstop\tstate "
    for x in range(len(starts)):
        print "%d\t%d\t%s" % (starts[x], stops[x], states[x])
        
def add_probs_in_log_space(p, q):
    s = log(p) + log(1 + 10 ** (log(q) - log(p)))
    return s


def get_emit_index(input_val, alphabet):
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    
    sys.stderr("Could not find character " + input_val)


if __name__ == '__main__':
    viterbi_decoding()
