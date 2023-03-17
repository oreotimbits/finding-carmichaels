#include <gmpxx.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <map>
#include <iterator>

const mpz_class ONE = 1;
const mpz_class THREE = 3;
const mpz_class FOUR = 4;

////////////////////////////////////////////////////////////////////

// finds all q, r, prime, such that pqr is a Carmichael number
std::map<mpz_class, mpz_class> carmichaelWithp(const mpz_class p) {

	mpz_class q, r, h3, delta;
	mpz_class negpsqr = -(p*p);
	std::map<mpz_class, mpz_class> umap;
	
	for (h3 = 2; mpz_cmp(h3.get_mpz_t(), p.get_mpz_t()) < 0; h3 = h3+1) {
		mpz_class pph3 = p+h3;
		mpz_class deltamult = (p-1)*pph3;
		for (delta = 1; mpz_cmp(delta.get_mpz_t(), pph3.get_mpz_t()) < 0; delta = delta+1) {
			if ((mpz_congruent_p(delta.get_mpz_t(), negpsqr.get_mpz_t(), h3.get_mpz_t())) &&
					mpz_divisible_p(deltamult.get_mpz_t(), delta.get_mpz_t())) {
				mpz_divexact(q.get_mpz_t(), deltamult.get_mpz_t(), delta.get_mpz_t());
				q = q + 1;
				if (mpz_probab_prime_p(q.get_mpz_t(), 15) == 2) {
					mpz_class tmp = p*q - 1;
					mpz_divexact(r.get_mpz_t(), tmp.get_mpz_t(), h3.get_mpz_t());
					r = r + 1;
					if (mpz_probab_prime_p(r.get_mpz_t(), 15) == 2) {
						mpz_class psub = p-1;
						mpz_class qr = q*r - 1;
						if (mpz_divisible_p(qr.get_mpz_t(), psub.get_mpz_t())) {
							umap[q] = r;
							//std::cout << p << " " << q << " " << r << std::endl;
						}
					}
				}
			}
		}
	} return umap;
}

////////////////////////////////////////////////////////////////////

void isMillerRabinCarmichael(const mpz_class n) {

	bool res = false;
	mpz_class rootn = sqrt(n);
	mpz_class p = 3;
	while (mpz_cmp(p.get_mpz_t(), rootn.get_mpz_t()) < 0) {
		if (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t())) {
			auto pcarm = carmichaelWithp(p);
			if (pcarm.size() > 0) {
				for (auto it = pcarm.begin(); it != pcarm.end(); ++it) {
					mpz_class pqr = p * it->first * it->second;
					if ((pqr == n) && mpz_congruent_p(it->first.get_mpz_t(), THREE.get_mpz_t(), FOUR.get_mpz_t()) && mpz_congruent_p(it->second.get_mpz_t(), THREE.get_mpz_t(), FOUR.get_mpz_t())) {
						res = true;
						std::cout << p << " " << it->first << " " << it->second << std::endl;
						break;
					}
				}
			}
		} mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		while (mpz_congruent_p(p.get_mpz_t(), ONE.get_mpz_t(), FOUR.get_mpz_t())) {
			mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
		}
	}
}

////////////////////////////////////////////////////////////////////

void isCarmichael(const mpz_class n) {

	bool res = false;
	mpz_class rootn = sqrt(n);
	mpz_class p = 3;
	while (mpz_cmp(p.get_mpz_t(), rootn.get_mpz_t()) < 0) {
		if (mpz_divisible_p(n.get_mpz_t(), p.get_mpz_t())) {
			auto pcarm = carmichaelWithp(p);
			if (pcarm.size() > 0) {
				for (auto it = pcarm.begin(); it != pcarm.end(); ++it) {
					mpz_class pqr = p * it->first * it->second;
					if ((pqr == n)) {
						res = true;
						std::cout << p << " " << it->first << " " << it->second << std::endl;
						break;
					}
				}
			}
		} mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
	}
}


////////////////////////////////////////////////////////////////////

void print_usage(char **argv) {

	std::cerr << "Usage: " << argv[0] << " [-mr | -p]" << std::endl;
	std::cerr << "Asks for input n and prints the primes p, q, r such that n = pqr is a Carmichael number" << std::endl;
	std::cerr << "-mr: If this option is specified, ignore all prime factors of n not congruent to 3 mod 4 (useful for Miller-Rabin primality test)" << std::endl;
	std::cerr << "-p: asks for input p and prints all primes q, r such that pqr is a Carmichael number" << std::endl;
}

////////////////////////////////////////////////////////////////////

bool valid_flag(char *c) {
	std::string s = (std::string) c;
	return ((s == "-p") or (s == "-mr"));
}

////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	std::cout << argc << std::endl;

	if ((argc > 2)) {
		print_usage(argv); 
	} else if ((argc == 2) && (!valid_flag(argv[1]))) {
		std::cout << argv[1] << std::endl;
		print_usage(argv);

	} else if ((argc == 1) || (std::string)argv[1] != "-p") {
		mpz_class n;
		std::cout << "Please enter a number n to check: ";
		mpz_inp_str(n.get_mpz_t(), stdin, 10);
		std::cout << std::endl;
		if (mpz_even_p(n.get_mpz_t())) {
			std::cerr << "Carmichael numbers can never be even" << std::endl;
			return 1;
		}
		//carmichaelWithp(p);
		if (argc == 1) { isCarmichael(n); }
		else { isMillerRabinCarmichael(n); }

	} else {
		mpz_class p;
		std::cout << "Please enter a prime number p: ";
		mpz_inp_str(p.get_mpz_t(), stdin, 10);
		std::cout << std::endl;
		if (!(mpz_probab_prime_p(p.get_mpz_t(), 50))) {
			std::cerr << "Number is not prime" << std::endl;
			return 1;
		}
		auto ret = carmichaelWithp(p);
		for (auto e : ret) {
			std::cout << e.first << " " << e.second << std::endl;
		}

	}

	return 0;
}
