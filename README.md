# Bernoulli-Numbers

I have implemented two basic ways of calculating Bernoulli numbers.

Because the Bernoulli numbers are rational, they should be represented by an integer numerator and integer denominator.
These integers get very large. So large that they quickly overflow even a long int.
To work around this, the GNU Multiple Precision Library (GMP) is used

GMP provides an API for working with integers and rational of any size. The only constraint being the computer's memory.
GMP provides an integer datatype, mpz (multiple precsion integer) and a ratoinal datatype, mpq (multiple precision rational)
The functions for working with these are all prepended with either mpz or mpq which is why so much of the code has these.

Both of the implemented algorithms are far from optimized, however they still serve the purpose of providing a runtime comparison between the two algorithms.
