// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h> // multiple precision lib

// My files

// Returns an array
// bn **bern_up_to(int n);

// bn *bern_n(int n);

// Creation/Destruction
// void mpq_init (mpq t x)
// void mpq_clear (mpq t x)

// Initialization
// void mpq_set_si (mpq t rop, signed long int op1, unsigned long int op2)
// void mpq_set_z (mpq t rop, const mpz t op)
// void mpq_swap (mpq t rop1, mpq t rop2)

// Arithmetic ops
// void mpq_add (mpq t sum, const mpq t addend1, const mpq t addend2)
// void mpq_sub (mpq t difference, const mpq t minuend, const mpq t)
// void mpq_mul (mpq t product, const mpq t multiplier, const mpq t)
// void mpq_div (mpq t quotient, const mpq t dividend, const mpq t)
// void mpq_neg (mpq t negated_operand, const mpq t operand)
// void mpq_abs (mpq t rop, const mpq t op)
// void mpq_inv (mpq t inverted_number, const mpq t number
//
// void mpz_bin_uiui (mpz t rop, unsigned long int n, unsigned long int k)

// Printing
// char * mpq_get_str (char *str, int base, const mpq t op)

// Uses single bernoulli formula
mpq_t *bern(int n) {
	mpq_t *result = calloc(1, sizeof(mpq_t)); // Hold the total of the double sum which is also the resulting bernoulli number
	mpq_init(*result);
	mpq_t cur_term; // Hold the current term of the sum being calculated
	mpq_init(cur_term);
	mpq_t sub_term; // Used to hold values before combining them into the cur_term
	mpq_init(sub_term);
	mpz_t tmp_int; // Used to hold values resulting from integer computation before it is transfered into an mpq
	mpz_init(tmp_int);

	for(int i = 0; i <= n; i++) {
		for(int j = 0; j <= i; j++) {
			mpz_bin_uiui(tmp_int, i, j); // nCk
			mpq_set_z(cur_term, tmp_int); // set cur_term to nCk
			mpq_set_si(sub_term, 1 , i + 1); // set subterm to 1 / (k + 1)
			mpq_mul(cur_term, cur_term, sub_term); // Combine cur_term and subterm
			mpz_ui_pow_ui(tmp_int, j + 1, n); // Calc expo
			mpq_set_z(sub_term, tmp_int); // set sub_term to the expo result
			mpq_mul(cur_term, cur_term, sub_term);
			if(j % 2 == 1) {
				mpq_neg(cur_term, cur_term);
			}
			mpq_add(*result, *result, cur_term);
		}
	}

	return result;
}

void bern_up_to_helper(mpq_t *res_arr, int n) {
	mpq_t sum; // Holds the total of the sum
	mpq_init(sum);
	mpq_t cur_term; // Holds the current term of the sum being calculated
	mpq_init(cur_term);
	mpq_t sub_term; // Holds 1 / (n - i + 1)
	mpq_init(sub_term);
	mpz_t tmp_int; // Holds the result of the binary coeff operation
	mpz_init(tmp_int);

	for(int i = 0; i <= n - 1; i++) {
		// set cur_term to the prod of the nCk * num * den
		mpz_bin_uiui(tmp_int, n, i); // nCk
		mpq_set_z(cur_term, tmp_int); // Set cur_term to nCk
		mpq_set_si(sub_term, 1, n - i + 1); // Set the subterm to n - k + 1
		mpq_mul(cur_term, cur_term, sub_term); 
		mpq_mul(cur_term, cur_term, res_arr[i]);
		mpq_add(sum, sum, cur_term); // Add the cur_term to the sum
	}

	// Creating a mpq with the value of 1
	mpq_t one;
	mpq_init(one);
	mpq_set_si(one, 1, 1);

	// Final calculation
	mpq_neg(sum, sum);
	mpq_add(res_arr[n], one, sum);

	// Free the memory
	// void mpz_clear (mpz t x)
	mpq_clear(sum);
	mpq_clear(cur_term);
	mpq_clear(sub_term);
	mpz_clear(tmp_int);
}

// Fills an array using recursive formula
mpq_t *bern_up_to(int n) {
	// Basic error check
	if(n < 0) {
		printf("the n arg for bern_up_to() must be positive instead of: %d\n", n);
		exit(1);
	}

	// Create the array
	mpq_t *res_arr;
	res_arr = calloc(n + 1, sizeof(mpq_t));
	mpq_t tmp_mpq;

	// Initialize all the mpq of the array
	for(int i = 0; i < n; i++) {
		mpq_init(res_arr[i]);
	}

	// Set b(0) to 1
	mpq_set_si(res_arr[0], 1, 1);

	// Loop through and call the helper to calculate the next value
	for(int i = 1; i <= n; i++) {
		bern_up_to_helper(res_arr, i);
	}

	return res_arr;
}

// For testing purposes
int main() {
	int n = 2000;
	// mpq_t *berns;
	mpq_t *bern_n;
	// berns = bern_up_to(n);
	bern_n = bern(n);
	/*
	for(int i = 0; i < n + 1; i++) {
		printf("%d: ", i);
		printf(mpq_get_str(NULL, 10, berns[i]));
		printf("\n");
	}
	*/
	printf("calculating B(%d)...\n", n);
	// printf(mpq_get_str(NULL, 10, berns[n]));
	printf(mpq_get_str(NULL, 10, *bern_n));
	printf("\n");
	return 0;
}
