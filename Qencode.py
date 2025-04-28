################################################################################
# File: Qencode.py
# Ver: 1.4
# Updated: 28 Apr 2025
# Goal: Encode a generic digital circuit in a quadratic form 
#       (matrix) Q. The minimum of Q describes the behavior 
#       of the circuit itself.
#       There are two procedures to find the minimum:
#       1) Exhaustive search (brute force)
#       2) Simulated Annealing process (probabilistic).
#       There are two examples:
#       1) Circuit for a simple two-bit adder
#       2) Circuit for a simple two-bit and 1 carry adder
#
#
# MIT License
#
# Copyright (c) 2025 Marco Mattiucci
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################


import random
import numpy as np
import itertools

    
    

### Class for a Digital Circuit description by an Hamitonian matrix Q
#####################################################################
class myDigitalCirtuit:

    #ATTRIBUTES(myDigitalCirtuit):
    ##############################
    Input_indexes = None  # List of the matrix indexes that are related to the circuit input lines
    Output_indexes = None # List of the matrix indexes that are related to the circuit output lines
    Gate_list = None      # List of the digital gates that are in the circuit ("and","or","not")
    Q = None              # Hamiltonian Q matrix of the circuit
    n = None              # Dimension of matrix Q
    Q_not = None          # Q matrix for NOT gate
    Q_and = None          # Q matrix for AND gate
    Q_or = None           # Q matrix for OR gate
    solutions = None      # Combinations found for the minimum of the quadratic function defined by Q
    cost = None           # Minimum values of the quadratic function defined by Q
    frequencies = None    # Frequencies of the solution for the simulated annealing process
    attempts = 100        # Number of attempts for the simulated annealing process for each input value
    
    #METHODS(myDigitalCirtuit):
    ###########################
    
    ######################################################################################
    ### The myDigitalCirtuit constructor defines the Q matrix using only the list of 
    ### gates in the circuit, without any connections.
    ######################################################################################
    def __init__(self,new_gate_list = []):
        #
        # Check the gate list (1):
        if new_gate_list == []:
            print("ERROR 1 from class myDigitalCirtuit constructor: gate list cannot be empty!")
            quit()
        #
        # Check the gate list (2):
        for g in new_gate_list:
            if g != "and" and g != "or" and g != "not":
                print("ERROR 2 from class myDigitalCirtuit constructor: undefined gate [",g,"]")
                quit()
        #
        # Set the Q matrixes for the basic boolean gates:
        self.Q_not = np.array([ [-1, 2], [ 0,-1] ])                 # Q matrix related to the NOT gate
        self.Q_and = np.array([ [0, 1,-2],[0, 0,-2],[0, 0, 3] ])    # Q matrix related to the AND gate
        self.Q_or = np.array([ [1, 1, -2],[0, 1, -2],[0, 0,  1] ])  # Q matrix related to the OR gate
        #
        # Define the "Gate_list" of the circuit:
        self.Gate_list = new_gate_list  
        #
        # Calculate the dimension n of the matrix Q:
        self.n = 0 
        for g in self.Gate_list:
            if g == "not":
                self.n += 2 # nr 2 states for the NOT
            elif g == "and" or g == "or":
                self.n += 3 # nr 3 states for the AND/OR 
        #
        # Define the Q matrix based on the ports in the input list but without any connections:
        self.Q = np.array([[0 for _ in range(self.n)] for _ in range(self.n)])  # Q as a zero matrix
        self.Input_indexes = []                                                 # No indexes for input
        my_offset = 0  # offset for the positions of the Q_not,Q_and,Q_or in the global Q matrix
        for g in self.Gate_list:    # for every gate add on the diagonal of Q the related "Q_gate"
            if g == "and":
                for i in range(3):
                    for j in range(3):  # Set the Q_and on the diagonal with offset:
                        self.Q[i+my_offset][j+my_offset] = self.Q_and[i][j]
                self.Input_indexes.append(my_offset)    # First index set as input of the AND
                self.Input_indexes.append(my_offset+1)  # Second index set as input of the AND
                my_offset += 3
            elif g == "or":
                for i in range(3):
                    for j in range(3):  # Set the Q_or on the diagonal with offset:
                        self.Q[i+my_offset][j+my_offset] = self.Q_or[i][j]
                self.Input_indexes.append(my_offset)    # First index set as input of the OR
                self.Input_indexes.append(my_offset+1)  # Second index set as input of the OR
                my_offset += 3
            elif g == "not":
                for i in range(2):
                    for j in range(2):  # Set the Q_not on the diagonal with offset:
                        self.Q[i+my_offset][j+my_offset] = self.Q_not[i][j]
                self.Input_indexes.append(my_offset)    # Unique index set as input of the NOT
                my_offset += 2
                
    #########################################################################################################
    ### myDigitalCirtuit "set_link" method creates a link between index_1 and index_2 of the matrix Q
    ### This is done adding (x1-x2)^2 = x1+x2-2x1x2 (for binary values) to the quadratic form defined by Q.
    ### So the term (index_1,index_1) is added by 1;
    ### the term (index_2,index_2) is added by 1;
    ### the term (index_1,index_2) or viceversa is added by -2.
    ### If index_1 or index_2 are input indexes one of them is removed from the "Input_indexes" list.
    ### gate_idx_1 is the index of the first gate in the gate list.
    ### connection_idx_1 is the index of the pin of the gate we're connecting (0,1,2 AND/OR, 0,1 NOT).
    ### gate_idx_2 is the index of the second gate in the gate list.
    ### connection_idx_2 is the index of the pin of the other gate we're connecting (0,1,2 AND/OR, 0,1 NOT).
    #########################################################################################################
    def set_link(self,gate_idx_1 = None, connection_idx_1 = None, gate_idx_2 = None, connection_idx_2 = None):
        #
        # Calculate the index_1 in the Q matrix corresponding to the (gate_idx_1,connection_idx_1):
        my_offset = 0
        for i in range(gate_idx_1):
            if self.Gate_list[i] == "and" or self.Gate_list[i] == "or":
                my_offset += 3
            elif self.Gate_list[i] == "not":
                my_offset += 2
        index_1 = my_offset + connection_idx_1
        #
        # Calculate the index_2 in the Q matrix corresponding to the (gate_idx_2,connection_idx_2):
        my_offset = 0
        for i in range(gate_idx_2):
            if self.Gate_list[i] == "and" or self.Gate_list[i] == "or":
                my_offset += 3
            elif self.Gate_list[i] == "not":
                my_offset += 2
        index_2 = my_offset + connection_idx_2
        #
        # Check the got index_1 and index_2:
        if index_1 not in range(self.n):
            print("ERROR 3 from class myDigitalCirtuit.set_link(): bad index_1 [",index_1,"]")
            quit()
        if index_2 not in range(self.n):
            print("ERROR 4 from class myDigitalCirtuit.set_link(): bad index_2 [",index_2,"]")
            quit()
        if index_1 == index_2:
            print("ERROR 5 from class myDigitalCirtuit.set_link(): index_1 and index_2 are already linked.")
            quit()
        #
        # The term (index_1,index_1) is added by 1:
        self.Q[index_1][index_1] += 1
        #
        # The term (index_2,index_2) is added by 1:
        self.Q[index_2][index_2] += 1
        #
        # The term (index_1,index_2) or (index_2,index_1) is added by -2 in the upper triangular part of matrix Q:
        if index_1 > index_2:
            self.Q[index_2][index_1] += -2
        else:
            self.Q[index_1][index_2] += -2
        #
        # If one of the 2 indexes is an input pin remove it from the "Input_indexes" list:
        if index_2 in self.Input_indexes:
            self.Input_indexes.remove(index_2)
        elif index_1 in self.Input_indexes:
            self.Input_indexes.remove(index_1)
        else:   # Error handling on index_1 and index_2:
            print("ERROR 6 from class myDigitalCirtuit.set_link(): index_1 and index_2 are both outputs!")
            quit()            
    
    ###########################################################################################################
    ### myDigitalCirtuit "bruteforce_QUBO()" method computes all values ​​of the quadratic form defined by
    ### Q by looking for the combination of binary inputs that obtains the minimum of the function. 
    ### "show" is False by default, when True the procedure is verbose.
    ###########################################################################################################
    def bruteforce_QUBO(self, show=False):
        #
        # Set the solutions and their costs to empty:
        self.solutions = []
        self.cost = []
        if show:
            print("Bruteforce QUBO solver is working...")
        #
        # Find the minimum (my_min) of the quadratic form defined by matrix Q:
        my_min = float('inf')                                       # Set my_min to +infinite
        solutions = list(itertools.product([0, 1], repeat=self.n))  # Find all the possible binary combination of n bits:
        my_counter = 0                                              # Initialize my_counter for process progress 
        total_counting = len(solutions)                             # Find the total number of possible combinations
        #
        # MAIN LOOP: evaluate the value of the quadratic form for every possible binary combination:
        for x in solutions:
            if my_counter % 10000 == 0: print("Evaluated",my_counter,"/",total_counting," solutions")
            my_counter += 1
            x_vector = np.array(x)                  # convert to np.array()
            val = x_vector.T @ self.Q @ x_vector    # Value of the quadratic function
            if val < my_min:        # If the value is less than the minimum:
                my_min = val            # set the new minimum
                self.cost = [my_min]    # collect the minimum
                self.solutions = [x]    # collect the solution
            elif val == my_min: # If there are more solutions giving the same minimum collect all of them:
                self.solutions.append(x)    # collect the solution
                self.cost.append(my_min)    # collect the minimum value of Q function
        if show: print("everything done.")
        
    ###########################################################################################################
    ### myDigitalCirtuit "simulated_annealing_QUBO()" method applies a simple simulated annealing process
    ### to solve the QUBO problem of finding the minimum of the quadratic form associated with matrix Q.
    ### The probabilistic process is applied "attempts" times and the result is the binary combination obtained 
    ### in output with the highest frequency.
    ### "show" is False by default, when True the procedure is verbose.
    ###########################################################################################################
    def simulated_annealing_QUBO(self, annealing_iterations=1000, starting_temperature=100, cooling_rate=0.99, show=False):
        if show:
            print("Simulated annealing QUBO solver:")
            print("Max number of attempts to use the simulated annealing process:",self.attempts)
            print("Max number of annealing iterations:",annealing_iterations)
            print("Starting temperature:",starting_temperature)
            print("Cooling rate:",cooling_rate)
            print("Processing...")
        #
        # Find all admissible circuit inputs related to matrix Q:
        nr_of_inputs = len(self.Input_indexes)
        input_values = list(itertools.product([0, 1], repeat=nr_of_inputs))
        #
        # Set the solutions and their frequencies to empty:
        self.solutions = []
        self.frequencies = []
        #
        # FIRST LOOP (inputs) - repeat the entire job for every single binary combination of inputs:
        for input_value in input_values:
            print("Evaluating input:",input_value)
            #
            # Set the input values into the Q matrix:
            # (add a positive value for 0 and a negative value for 1 at the position in 
            # the Q matrix corresponding to the input index)
            for i in range(nr_of_inputs):
                if input_value[i] == 0:
                    self.Q[self.Input_indexes[i]][self.Input_indexes[i]] += 3
                elif input_value[i] == 1:
                    self.Q[self.Input_indexes[i]][self.Input_indexes[i]] += -3
            #
            # Initialize the total solutions list and their frequencies list:
            total_solutions = []
            frequencies = []
            #
            # SECOND LOOP (attempts) - repeat the simulated annealing process many times:
            for attempt in range(self.attempts):
                current_solution = np.random.randint(0, 2, size=self.n)         # Calculate a random solution
                current_cost = current_solution.T @ self.Q @ current_solution   # Calculate the value of the quadratic form
                best_solution = current_solution.copy()                         # Set the starting best solution
                best_cost = current_cost                                        # Set the starting best cost
                temperature = starting_temperature                              # Set the "temperature" of the simulation
                #
                # THIRD LOOP (annealing) - iterations that simulate the annealing process:
                for i in range(annealing_iterations):
                    #
                    # Create a new solution randomly flipping a bit:
                    new_solution = current_solution.copy()      
                    bit_to_flip = random.randint(0, self.n - 1)
                    new_solution[bit_to_flip] = 1 - new_solution[bit_to_flip]
                    #
                    new_cost = new_solution.T @ self.Q @ new_solution   # Calculate the value of the quadratic form
                    delta_cost = new_cost - current_cost                # Calculate the delta of "energy" in the simulation
                    #
                    # Accept the new solution based on a random exponential process:
                    if delta_cost < 0 or random.random() < np.exp(-delta_cost / temperature):
                        current_solution = new_solution     # Set a new current solution
                        current_cost = new_cost             # Set a new current cost
                        #
                        # If it is a better cost: update the best solution and cost found:
                        if current_cost < best_cost:
                            best_solution = current_solution.copy()
                            best_cost = current_cost
                    #
                    # Simulate the cooling process:
                    temperature *= cooling_rate
                #
                # Store the solution at the end of the annealing process and update its frequency:
                if total_solutions == [] and frequencies == []:
                    total_solutions.append(best_solution)
                    frequencies.append(1)
                else:
                    found_new_solution = True
                    for k in range(len(total_solutions)):
                        if (total_solutions[k] == best_solution).all():
                            frequencies[k] += 1
                            found_new_solution = False
                            break
                    if found_new_solution:
                        total_solutions.append(best_solution)
                        frequencies.append(1)
            #
            # After all the attempts find the solution with the maximum frequency:
            max_freq = 0
            max_freq_solution = None
            for k in range(len(total_solutions)):
                if frequencies[k] > max_freq:
                    max_freq_solution = total_solutions[k]
                    max_freq = frequencies[k]
            self.solutions.append(max_freq_solution)
            self.frequencies.append(max_freq)
            #
            # Restore the original values into the Q matrix because in the next iteration 
            # binary a new input combination will be used:
            for i in range(nr_of_inputs):
                if input_value[i] == 0:
                    self.Q[self.Input_indexes[i]][self.Input_indexes[i]] += -3
                elif input_value[i] == 1:
                    self.Q[self.Input_indexes[i]][self.Input_indexes[i]] += 3
        if show: print("everything done.")

    ###########################################################################################################
    ### myDigitalCirtuit "define_output_idx_list()" method defines the list of indices of the Q matrix that 
    ### constitute the output of the circuit. Each component of the list must be (gate_idx, connection_idx), 
    ### the meaning of which has already been explained in the "set_link" method comments.
    ###########################################################################################################
    def define_output_idx_list(self, idx_list = []):
        #
        # if "all" then all the indexes of Q that are not input are defined as ouput:
        if idx_list == "all":
            self.Output_indexes = []
            for i in range(self.n):
                if i not in self.Input_indexes:
                    self.Output_indexes.append(i)
        elif len(idx_list) > self.n:                # Check the lengh of the idx_list:
            print("ERROR 7 from class myDigitalCirtuit.define_output_idx_list(): bad idx list [",idx_list,"]")
            quit()
        else:
            #
            # Initialise to empty the elaborated list of Q indexes related to the couples of "idx_list":
            elab_idx_list = []
            #
            # for every couple in "idx_list":
            for (gate_idx,connection_idx) in idx_list:
                # Calculate the index of Q matrix related to gate_idx and connection_idx:
                my_offset = 0
                # Start from the first gate in the "Gate_list" of the circuit and calculate the offset of "gate_idx"
                for i in range(gate_idx):
                    if self.Gate_list[i] == "and" or self.Gate_list[i] == "or":
                        my_offset += 3
                    elif self.Gate_list[i] == "not":
                        my_offset += 2
                # "my_offset + connection_idx" is a Q index to be considered as output index
                elab_idx_list.append(my_offset + connection_idx)
            self.Output_indexes = elab_idx_list     # Collect all Q indexes that can be considered output indexes.

    ###########################################################################################################
    ### myDigitalCirtuit "show_truth_table()" method shows solutions of the QUBO problem solved by inserting
    ### only the "Input_indexes" and "Output_indexes" of the circuit.
    ### "evaluation" can be "cost" for "bruteforce_QUBO" or "freq" for "simulated_annealing_QUBO".
    ### When "full" is True, the entire state of Q is displayed, so n bits without distinction between input
    ### and output. When "full" is False, only the values ​​of "Input_indexes" and "Output_indexes" are displayed.
    ###########################################################################################################
    def show_truth_table(self,full = False, evaluation = "cost"):
        if self.Output_indexes == None:
            print("ERROR 8 from class myDigitalCirtuit.show_truth_table(): undefined output indexes, use define_output_idx_list() method!")
            quit()
        print("Truth table solutions:")
        if not full and self.solutions != []:
            for i in self.Input_indexes:
                print(i,end=" ")
            print(" - ",end=" ")
            for i in self.Output_indexes:
                print(i,end=" ")
            print("(input and output indexes of the matrix Q)")
        if evaluation == "cost":
            for (x,c) in zip(self.solutions,self.cost):
                if full:
                    for y in x:
                        print(y,end=" ")
                else:
                    for i in self.Input_indexes:
                        print(x[i],end=" ")
                    print(" - ",end=" ")
                    for i in self.Output_indexes:
                        print(x[i],end=" ")
                print(" min=",c)
        elif evaluation == "freq":
            for (x,f) in zip(self.solutions,self.frequencies):
                if full:
                    for y in x:
                        print(y,end=" ")
                else:
                    for i in self.Input_indexes:
                        print(x[i],end=" ")
                    print(" - ",end=" ")
                    for i in self.Output_indexes:
                        print(x[i],end=" ")
                print(" freq = ",f,"/",self.attempts)
        else:
            print("ERROR 9 from class myDigitalCirtuit.show_truth_table(): undefined evaluation [",evaluation,"], only cost/freq allowed")
            quit()




################################
########## MAIN PROGRAM ########
#### Two examples:          ####
#### 1) simple binary adder ####
#### 2) full binary adder   ####
################################
if __name__ == "__main__":
    #
    ##########################################################
    # EXAMPLE CIRCUIT 0 - VERY BASIC:
    # Building a NAND gate:
    print()
    print("Connecting a NOT to the ouput of an AND.")
    #
    # Define the gates of the circuit:
    circuit_0 = myDigitalCirtuit(["and","not"])
    print("Matrix Q with AND and NOT not connected between them:")
    print(circuit_0.Q)
    #
    # Gate naming:
    nand_and = 0
    nand_not = 1
    #
    # NAND:
    # Connecting...
    circuit_0.set_link(gate_idx_1 = nand_and, connection_idx_1 = 2, gate_idx_2 = nand_not, connection_idx_2 = 0)
    print("Matrix Q with NOT connected to the outout of an AND:")
    print(circuit_0.Q)
    #
    # Define OUTPUT indexes:
    circuit_0.define_output_idx_list([(nand_not,1)])
    #
    # (4) Find the minimum of the Q function (Two QUBO Solvers):
    #
    print("######################################")
    circuit_0.bruteforce_QUBO(show=True)          # Brute force analysis of every possible input value
    circuit_0.show_truth_table(evaluation="cost")
    #
    print("######################################")
    circuit_0.simulated_annealing_QUBO(show=True) # Simulated Annealing process: results based on max frequency
    circuit_0.show_truth_table(evaluation="freq")
    #
    #
    ##########################################################
    # EXAMPLE CIRCUIT 1:
    # Simple adder bit_0 + bit_1 = bit_total and carry_out:
    print()
    print("Getting the quantum annealing matrix for a simple adder: two bits input and two bits output.")
    #
    # Define the gates of the circuit:
    circuit_1 = myDigitalCirtuit(["not","and","not","and","or","and"])
    print("Matrix Q of a simple adder with just the gates:")
    print(circuit_1.Q)
    #
    # Gate naming:
    xor_not_1 = 0
    xor_and_1 = 1
    xor_not_2 = 2
    xor_and_2 = 3
    xor_output_or = 4
    carry_and = 5
    #
    # INTEGER SUM OF TWO BITS:
    # (1) XOR:
    circuit_1.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 1, gate_idx_2 = xor_and_1, connection_idx_2 = 0)
    circuit_1.set_link(gate_idx_1 = xor_and_1, connection_idx_1 = 2, gate_idx_2 = xor_output_or, connection_idx_2 = 0)
    circuit_1.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 1, gate_idx_2 = xor_and_2, connection_idx_2 = 0)
    circuit_1.set_link(gate_idx_1 = xor_and_2, connection_idx_1 = 2, gate_idx_2 = xor_output_or, connection_idx_2 = 1)
    circuit_1.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 0, gate_idx_2 = xor_and_2, connection_idx_2 = 1)
    circuit_1.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 0, gate_idx_2 = xor_and_1, connection_idx_2 = 1)  
    #
    # (2) CARRY OUT:
    circuit_1.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 0, gate_idx_2 = carry_and, connection_idx_2 = 0)
    circuit_1.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 0, gate_idx_2 = carry_and, connection_idx_2 = 1)
    print("Matrix Q of a simple adder with gates and connections among them:")
    print(circuit_1.Q)
    #
    # (3) Define OUTPUT indexes:
    circuit_1.define_output_idx_list([(carry_and,2),(xor_output_or,2)])
    #
    # (4) Find the minimum of the Q function (Two QUBO Solvers):
    #
    print("######################################")
    circuit_1.bruteforce_QUBO(show=True)          # Brute force analysis of every possible input value
    circuit_1.show_truth_table(evaluation="cost")
    #
    print("######################################")
    circuit_1.simulated_annealing_QUBO(show=True) # Simulated Annealing process: results based on max frequency
    circuit_1.show_truth_table(evaluation="freq")
    #
    #
    ##########################################################
    # EXAMPLE CIRCUIT 2:
    # Full adder (INPUT: two bits and a carry in - OUTPUT: two bits and a carry out) 
    print()
    print("Getting the quantum annealing matrix for a full adder: three bits input and two bits output.")
    #
    # Define the gates of the circuit:
    circuit_2 = myDigitalCirtuit(["not","and","not","and","or","not","and","not","and","or","and","and","or"])
    #
    # Gate naming:
    xor_not_1 = 0
    xor_and_1 = 1
    xor_not_2 = 2
    xor_and_2 = 3
    xor_output_or = 4
    carry_and_1 = 10
    carry_and_2 = 11
    carry_or = 12
    #
    # INTEGER SUM OF THREE BITS:
    # (1) first XOR:
    circuit_2.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 1, gate_idx_2 = xor_and_1, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_and_1, connection_idx_1 = 2, gate_idx_2 = xor_output_or, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 1, gate_idx_2 = xor_and_2, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_and_2, connection_idx_1 = 2, gate_idx_2 = xor_output_or, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 0, gate_idx_2 = xor_and_2, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 0, gate_idx_2 = xor_and_1, connection_idx_2 = 1)
    #
    # (2) Second XOR
    shift_right = 5
    circuit_2.set_link(gate_idx_1 = xor_not_1+shift_right, connection_idx_1 = 1, gate_idx_2 = xor_and_1+shift_right, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_and_1+shift_right, connection_idx_1 = 2, gate_idx_2 = xor_output_or+shift_right, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_not_2+shift_right, connection_idx_1 = 1, gate_idx_2 = xor_and_2+shift_right, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_and_2+shift_right, connection_idx_1 = 2, gate_idx_2 = xor_output_or+shift_right, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = xor_not_1+shift_right, connection_idx_1 = 0, gate_idx_2 = xor_and_2+shift_right, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = xor_not_2+shift_right, connection_idx_1 = 0, gate_idx_2 = xor_and_1+shift_right, connection_idx_2 = 1)
    #
    # (3) Full adder for bit_0, bit_1 and carry_in
    circuit_2.set_link(gate_idx_1 = xor_output_or, connection_idx_1 = 2, gate_idx_2 = xor_not_1+shift_right, connection_idx_2 = 0)
    #
    # (4) Carry out circuit
    circuit_2.set_link(gate_idx_1 = xor_not_1, connection_idx_1 = 0, gate_idx_2 = carry_and_1, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_not_2, connection_idx_1 = 0, gate_idx_2 = carry_and_1, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = xor_output_or, connection_idx_1 = 2, gate_idx_2 = carry_and_2, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = xor_not_2+shift_right, connection_idx_1 = 0, gate_idx_2 = carry_and_2, connection_idx_2 = 1)
    circuit_2.set_link(gate_idx_1 = carry_and_1, connection_idx_1 = 2, gate_idx_2 = carry_or, connection_idx_2 = 0)
    circuit_2.set_link(gate_idx_1 = carry_and_2, connection_idx_1 = 2, gate_idx_2 = carry_or, connection_idx_2 = 1) 
    #
    # (5) Define output indexes
    circuit_2.define_output_idx_list([(carry_or,2),(xor_output_or+shift_right,2)])
    #
    # (6) Find the minimum of the Q function (Just one QUBO Solver, the brute force solver does not work):
    #
    print("######################################")
    circuit_2.simulated_annealing_QUBO(show=True) # Simulated Annealing process: results based on max frequency
    circuit_2.show_truth_table(evaluation="freq")
    quit()
    

