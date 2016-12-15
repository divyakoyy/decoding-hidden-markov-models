
import sys
from viterbi import display, find_num_matches
from math import log, exp
from compsci260lib import get_fasta_dict
import matplotlib.pyplot as plt

def posterior_decoding():
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
    
    # read the emission symbols
    emission_symbols = f_hmm_file.readline().split()
    
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
    
    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    pi_hats = []
    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    normalized = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        
        maxVal =  float('-inf')
        total_col = 0
        for k in range(K):
            val = forward[i][k] + backward[i][k]
            posterior[i][k] = val
            total_col += val
            if(val > maxVal):
                maxVal = val
                state = k
                
        for k in range(K):
            normalized[i][k] = val/total_col
        pi_hats.append(state)
    
 
    # Print out the decoded results
    starts, stops, states = create_table(pi_hats, 1)
    display(starts, stops, states)
    print "number of matches using 10 nucleotides margin of error:" ,find_num_matches('tRNA.locations.txt', starts, stops, states)[0]
 
    
def create_table(lst, start_index):
    starts = []
    stops = []
    states = []
    
    start = 1
    curr_state = lst[0]
    for x in range(1, len(lst)):
        if(lst[x] != curr_state):
            starts.append(start)
            stops.append(x)
            states.append("state " + str(lst[x - 1] + 1))
            start = x + 1
            curr_state = lst[x]
    
    # handle the end of the list
    if(stops[-1] != len(lst)):
        starts.append(stops[-1] + 1)
        stops.append(len(lst))
        states.append("state " + str(lst[len(lst) - 1] + 1))
    
    return starts, stops, states

def add_probs_in_log_space(p, q):
    ''' takes two log probabilities as input and returns their sum '''
    s = p + log(1 + exp(q - p))
    return s

def add_list_of_probs_in_log_space(lst):
    ''' adds a list of log probabilities'''
    sum = lst[0]
    for x in range(1, len(lst)):
        sum = add_probs_in_log_space(sum, lst[x])
    
    return sum

def run_forward(states, initial_probs, transitions, 
    emission_symbols, emit_probs, emit_str):
    """Calculates the forward probability matrix"""

    K = len(states)
    L = len(emit_str)

    forward = [[float(0) for _ in range(K)] for _ in range(L)]

    # Initialize
    emit_index = get_emit_index(emit_str[0], emission_symbols)
    for k in range(K):
        forward[0][k] = log(initial_probs[k]) + log(emit_probs[k][emit_index])
    
    # Iterate
    for i in range(1, L):
        emit_index = get_emit_index(emit_str[i].upper(), emission_symbols)
        
        # Compute the forward probabilities for the states
        for l in range(K):
            log_probs = []
            for k in range(K):
                log_prob = forward[i - 1][k] + log(transitions[k][l]) + log(emit_probs[l][emit_index])
                log_probs.append(log_prob)
            sum = add_list_of_probs_in_log_space(log_probs)
            forward[i][l] = sum
     
        
    return forward        
        

def run_backward(states, initial_probs, transitions, 
    emission_symbols, emit_probs, emit_str):
    """Calculates the backward probability matrix"""

    K = len(states)
    L = len(emit_str)

    backward = [[float(0) for _ in range(K)] for _ in range(L)]

    # Initialize
    for k in range(K):
        backward[L-1][k] = log(1)  # which is zero, but just to be explicit...
    
    # Iterate
    for i in range(L-2, -1, -1):
        emit_index = get_emit_index(emit_str[i+1].upper(), emission_symbols)
        # Compute the backward probabilities for the states
        for k in range(K):
            log_probs = []
            for l in range(K):
                log_prob = backward[i + 1][l] + log(transitions[k][l]) + log(emit_probs[k][emit_index])
                log_probs.append(log_prob)
            sum = add_list_of_probs_in_log_space(log_probs)
            backward[i][k] = sum

    return backward        
        

def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


if __name__ == '__main__':
    posterior_decoding()
