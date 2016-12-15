from compsci260lib import *
import sys, random, sampling

def generate_HMM():
    f = open("HMM.txt", "r")
    
    if f is None:
        sys.stderr("Can't open HMM parameter file: HMM.txt")
    
    # read the state names
    states = f.readline().strip().split()
    
    # read the initial probabilities
    init_row = f.readline().strip().split()
    initial_probs = [float(x) for x in init_row]
    
    # read the transition matrix
    transitions = [None for _ in range(len(states))]
    for i in range(0, len(states)):
        matrix_row = f.readline().strip().split()
        transitions[i] = [float(x) for x in matrix_row]
    
    # read the input alphabet
    input_alphabet = f.readline().strip().split()
    
    # read the input matrix
    inputs = [None for _ in range(len(states))]
    for i in range(0, len(states)):
        matrix_row = f.readline().strip().split()
        inputs[i] =  [float(x) for x in matrix_row]
    
    f.close()
    
    # the length of the sequence we will generate
    seq_length = 100
    
    #YOUR CODE GOES HERE....

    states = ""
    sequence = ""
    curr_state = sampling.solve_sampling_binary_search(initial_probs)
    
    for x in range(seq_length):
        
        curr_state = sampling.solve_sampling_binary_search(transitions[curr_state - 1])
        curr_event = sampling.solve_sampling_binary_search(inputs[curr_state - 1])
        states += str(curr_state)
        sequence += choose_nuc(curr_event)
        
    #print states
    #print sequence
    
    state_transitions, positions_in_state = make_table(states)
    
    #display(state_transitions, positions_in_state)
    
   
    #print sequence
    
    return state_transitions, positions_in_state
   

        
def make_table(states):
    
    positions_in_state = [] 
    curr_state = states[0]
    state_transitions = [curr_state]
    ct = 1
    
    for x in range(1, len(states)):
        if states[x] != curr_state:
            positions_in_state.append(ct)
            ct = 1
            curr_state = states[x]
            state_transitions.append(curr_state)
        else:
            ct += 1
    
    positions_in_state.append(ct)
    return state_transitions, positions_in_state

def display(state_transitions, positions_in_state):
    print "state\tnumber of positions spent in state "
    for x in range(len(state_transitions)):
        print "%s\t%d" % (state_transitions[x], positions_in_state[x])
        
def choose_nuc(event):
    if(event == 1):
        return 'A'
    elif(event == 2):
        return 'C'
    elif(event == 3):
        return 'G'
    elif(event == 4):
        return 'T'
    
if __name__ == '__main__':
    sum_state_1 = 0.0
    sum_state_2 = 0.0
    ct_in_state_1 = 0
    ct_in_state_2 = 0
    for x in range(10):
        state_transitions, positions_in_state = generate_HMM()
        
        for i in range(len(state_transitions)):
            if state_transitions[i] == '1':
                sum_state_1 += positions_in_state[i]
                ct_in_state_1 += 1
            else:
                sum_state_2 += positions_in_state[i]
                ct_in_state_2 += 1
                
            
    print "state 1 average duration:", sum_state_1/ct_in_state_1
    print "state 2 average duration:", sum_state_2/ct_in_state_2
    
    

