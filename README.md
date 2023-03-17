# finding-carmichaels

## Introduction

This code was created as part of a personal research project for the Fall 2020 offering of CO 485 (The Mathematics of Public Key Cryptography) at the University of Waterloo.

The Miller-Rabin primality test is a probabilistic procedure that determines
whether a given number is likely to be prime. It has been found that
Carmichael numbers with 3 prime factors (simply referred to here as 'Carmichael
numbers'), especially those with all prime factors congruent to 3 modulo 4, have
a higher chance of producing a false positive with the Miller-Rabin primality
test. This project introduces a procedure to determine whether a given number
is a Carmichael number, and a method to obtain all the Carmichael numbers
given the smallest prime. This work hopes to improve the accuracy of the
Miller-Rabin primality test by determining values highly likely to result in a
false positive by the test.

## References

[1] G.J.O Jameson. Finding Carmichael Numbers. The Mathematical Gazette, Vol.95, No.533, pp. 244-255.

[2] S. Narayanan. Improving the Speed and Accuracy of the Miller Rabin primality test.
