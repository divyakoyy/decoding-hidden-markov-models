from compsci260lib import *
import random

def solve_sampling(probabilities):

    prob_distribution = [(i+1,probabilities[i]) for i in range(len(probabilities))]
    
    rand = random.random()
    p = 0
    for event, prob in prob_distribution:
        p += prob
        if p >= rand:
            return event
    
    
def solve_sampling_binary_search(probabilities):
    prob_distribution = [(i+1,probabilities[i]) for i in range(len(probabilities))]
    
    cont_distribution = []
    sum = 0
    
    # find a list of cumulative sums of probabilities
    for i in range(len(prob_distribution)):
        sum += prob_distribution[i][1]
        cont_distribution.append((i+1,sum))

    rand = random.random()
    
    L = 0
    R = len(prob_distribution) - 1
    
    while L <= R:
        m = (L+R)/2
        
        if m == len(cont_distribution) - 1 or rand == cont_distribution[m][1]:
            return cont_distribution[m][0]
        
        if rand > cont_distribution[m][1] and rand <= cont_distribution[m+1][1]:
            return cont_distribution[m+1][0]
        
        if rand > cont_distribution[m][1]:
            L = m
        else:
            if m == 0:
                return cont_distribution[m][0]
            R = m

if __name__ == '__main__':
    print "in O(n):"
    for x in range(10):
        print "event", solve_sampling((0.5,0.5))
    
    print
    print "in O(logn) using binary search:"
    for x in range(10):
        print "event", solve_sampling_binary_search((0.5,0.5))
