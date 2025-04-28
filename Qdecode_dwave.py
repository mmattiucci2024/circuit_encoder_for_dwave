################################################################################
# File: dwave_decode.py
# Ver: 1.2
# Updated: 28 Apr 2025
# Goal: Decoding a quadratic form using only a Q matrix
#       stored in a file. This is provided by DWAVE Quantum
#       Annealing Computation Services.
#       The binary samples found are saved to a file.
################################################################################
# WARNING! This code can be executed only in the venv defined
#          by DWAVE (E.G. Linux virtual environment) during
#          DWAVE service installation.
# WARNING! You need a DWAVE subscription (also free) to make
#          this program work.
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


import sys
import os.path
from datetime import datetime

""" Import DWAVE Sampler """
from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())


######################
#### MAIN PROGRAM ####
######################
""" Set the name of the file containing the encoded Q list """
my_input_file_name = "Qfile_example_for_dwave.txt"

""" Set an empty dictionary that is the form of the Q matrix accepted by DWAVE services """
Q = {}
counter = 0
if os.path.exists(my_input_file_name):
    print("Loading input file",my_input_file_name,"]:")
    with open(my_input_file_name, 'r') as my_file_handler:
        """ The first element is the Q matrix dimension """
        n = int(my_file_handler.readline().strip())
        while True:
            """ Read the first index i """
            r = my_file_handler.readline().strip()
            if not r:
                break 
            i = r
            """ Read the first index j """
            r = my_file_handler.readline().strip()
            if not r:
                break 
            j = r
            """ Read the Q value in i,j and associate it to (xi,xj) in the dictionary """
            r = my_file_handler.readline().strip()
            if not r:
                break 
            Q['x'+str(i),'x'+str(j)] = int(r)
            counter += 1
        my_file_handler.close()
else:
    print("Error 1 from Qdecode_dwave.py: cannot access input file!",my_input_file_name)
    quit()
print("Done: Got Q-dimension",n,"and loaded",counter,"non-zero values of the Q matrix.")
print("Q dim = ",n)
print("Q = ",Q)
print("Loading QUBO into DWAVE Qauntum Annealer...",end="")
""" Invoke the Quantum annealer service 2000 times to find the minimum """
""" Be careful in changing num_reads, id it is too high DWAVE service does not work """
sampleset = sampler.sample_qubo(Q, num_reads=2000, label=my_input_file_name+' 2000reads')
print("task executed:")
""" sampleset is the set of all the minimum found by the quantum annealer """
""" with different energy (value of the quadratic form) and occurrences """
print(sampleset)
""" Define the filename of the output """
current_dateTime = datetime.now()
my_output_file_name = "resulting_sampleset_"+str(current_dateTime.year)
my_output_file_name += "_"+str(current_dateTime.month)
my_output_file_name += "_"+str(current_dateTime.day)
my_output_file_name += "_"+str(current_dateTime.hour)
my_output_file_name += "_"+str(current_dateTime.minute)
my_output_file_name += "_"+str(current_dateTime.second)
my_output_file_name += ".txt"
""" Save all the sampleset in the output txt file """
print("Saving resulting samples to output file [mydata/",my_output_file_name,"]...",end="")
with open("mydata/"+my_output_file_name, 'w') as my_file_handler:
    my_file_handler.write("{}\n".format(str(n)))
    for item in sampleset.variables:
        my_file_handler.write("{}\n".format(str(item)))
    for item in sampleset.record:
        [a,b,c,d] = item
        for ai in a:
            my_file_handler.write("{}\n".format(str(ai)))
        my_file_handler.write("{}\n".format(str(b)))
        my_file_handler.write("{}\n".format(str(c)))
        my_file_handler.write("{}\n".format(str(d)))
    my_file_handler.close()
print("done.")
quit()
