// MIT License

// Copyright (c) 2025 [Eric Ahlqvist]

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef MISC_FUNCTIONS_H
#define MISC_FUNCTIONS_H

#include <pari/pari.h>

//-------------------------------------------------------------------------------
// For debugging
void print_pari_type(GEN x);

void my_test_p_rank (GEN K, int p_int);

void my_check_galois(GEN K);

GEN concatenate_rows(GEN M1, GEN M2);

/*------------------------------------
 The function 1-sigma_x on ideals
------------------------------------
* Input:
* L - a number field
* sigma - An auto nfalgtobasis(nf, c[2]) where c = nfgaloisconj(nf);
* I - an ideal of nf given as a matrix in hnf

* Output: (1-\sigma_x)(I)
-----------------------*/
GEN my_1MS_ideal (GEN L, GEN sigma, GEN I);

/*------------------------------------
* Input:
* L - a number field
* factorization - a factorization of an ideal in L

* Output: a 2 term vector: first component is the vector of primes, second component; the vector of exponents
-----------------------*/
GEN my_find_primes_in_factorization(GEN LyAbs, GEN factorization);

/*------------------------------------
 Find the operator (1-\sigma_x)^n
------------------------------------
* Input:
* Labs - a bnf
* Lbnr - a bnr
* sigma - An auto nfalgtobasis(nf, c[2]) where c = nfgaloisconj(nf);
* n - an integer

* Output: The matrix of the operator (1-\sigma_x)^n in the basis of Cl(L)
-----------------------*/
GEN my_1MS_operator_2 (GEN Labs, GEN Lbnr, GEN sigma, int n);

// Relative norm on elements in compact representation
GEN my_rel_norm_compact(GEN Labs, GEN Lrel, GEN K, GEN compact_elt);

GEN my_vect_from_exp (GEN basis, GEN exp);

//------------------------------
// Returns all ideals I (actually, one solution I_0 plus a matrix M s.t. all solutions can be obtained by adding a lin. comb. of the columns of the matrix to I_0) in Div(L) such that iJ = (1-sigma)I in Cl(L)
//------------------------------ 
GEN my_H90_2 (GEN L, GEN iJ, GEN oneMS_operator, int n);

// Generators for the units of K modulo p
GEN my_find_units_mod_p (GEN K, GEN p);

GEN my_norm_operator (GEN Labs, GEN Lrel, GEN K, GEN p);

GEN my_find_p_gens (GEN K, GEN p);

// Returns vector of tuples (a,J) with div(a)+pJ = 0.
GEN my_find_Ja_vect(GEN K, GEN J_vect, GEN p, GEN units_mod_p);

// Creates a set (vector) of of exponents in bijection with Cl(L) 
GEN my_get_vect (int n, GEN cyc);

GEN my_get_sums (GEN basis, int p);

//-------------------------------------------------------------------------------------------------
// The function my_H90_vect finds for each (a, J) in Ja_vect, a fractional ideal I in Div(L) such that 
// i(J) = (1-sigma)I + div(t), where t in L^x satisfies N(t)*a = 1.  
// Returns: a vector of these I's 
//-------------------------------------------------------------------------------------------------
GEN my_H90_vect_2 (GEN Labs, GEN Lrel, GEN Lbnr, GEN K, GEN sigma, GEN Ja_vect, GEN p, int n);

//------------------------------
// Some old slow function used only in some old tests...
//------------------------------
GEN my_get_clgp (GEN K);

void my_unramified_p_extensions(GEN K, GEN p, GEN D_prime_vect);

GEN my_ideal_lifts (GEN Labs, GEN Lrel, GEN K, GEN p);

void my_unramified_p_extensions_with_transfer(GEN K, GEN p, GEN D_prime_vect);

/*------------------------------------
 Find the best subgroups (those with smallest class numbers)
------------------------------------
* Input:
* K - a bnf (number field)
* p_rank - number of subgroups to return
* subgroups - vector of subgroups

* Output: vector of length p_rank containing the subgroups with smallest class numbers
-----------------------*/
GEN my_best_subgroups(GEN K, long p_rank, GEN subgroups, GEN D_prime_vect);

#endif // MISC_FUNCTIONS_H
