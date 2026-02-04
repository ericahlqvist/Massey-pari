# Massey-pari

This repository contains a program which computes Massey products in étale cohomology of the ring of integers of a number field plus some more which is explained in the output. The program is written in C using the library [PARI](http://pari.math.u-bordeaux.fr/) which makes it very fast for number fields of low degree. 

To run the main program: ./massey p "pol(s)"
where you replace pol(s) by a defining polynomial of the number field K in the variable s, and replace p by a prime number dividing the class number. 