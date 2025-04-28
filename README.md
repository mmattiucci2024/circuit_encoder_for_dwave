# circuit_encoder_for_dwave
Software to convert digital circuits into matrices to be used as input for DWAVE quantum annealing computers:

https://marcomattiucci.it/informatica_quantum_computing_annealing3.php

There are two programs available:

(1) the circuit encoder ??? that can transform any asynchronous digital circuit into a Q matrix that can be used as input for QUBO (Quantum Unconstrained Binary Optimization) solvers. Along with the encoder, there are 3 examples of encoded circuits and two QUBO solvers that can run on regular PCs;

(2) an example of software ??? that programs the DWAVE quantum computer (available online) with the mentioned Q matrix and requires the QUBO solution (note that you need a DWAVE subscription to run this software and a venv environment on Linux).
