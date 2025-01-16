
#include <pthread.h>
#include <pari/pari.h>
#include <stdio.h>
#include <time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// Debug level
#define MY_DEBUGLEVEL 0

// Debug printing function
#define DEBUG_PRINT(level, ...) \
    do { if (MY_DEBUGLEVEL >= (level)) pari_printf(__VA_ARGS__); } while (0)

#include "headers/misc_functions.h"
#include "headers/tests.h"
#include "headers/artin_symbol.h"
#include "headers/test_artin.h"
#include "headers/ext_and_aut.h"
#include "headers/find_cup_matrix.h"

// Function prototype for parallel computation
GEN compute_my_relations(long i, GEN args);

int
main (int argc, char *argv[])	  
{
    printf(ANSI_COLOR_YELLOW "\n---------------------------------------------------------------------------------------------------------\nStarting program: finding Massey products and relations for Q_2\n---------------------------------------------------------------------------------------------------------\n\n" ANSI_COLOR_RESET);

    // Start timer (actual time)
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    
    // Start timer (CPU time)
    clock_t start = clock();

    int p_int, p_rk, r_rk, min, sec, msec;
    
    //--------------------------------------------------
    // Initialize PARI/GP
    //--------------------------------------------------
    // pari_init(1L<<30,500000);
    entree ep = {"_worker",0,(void*)compute_my_relations,20,"LG",""};
    pari_init_opts(1L<<30,500000, INIT_JMPm|INIT_SIGm|INIT_DFTm|INIT_noIMTm);
    pari_add_function(&ep); /* add Cworker function to gp */
    pari_mt_init(); /* ... THEN initialize parallelism */
    paristack_setsize(1L<<30, 1L<<33);
    sd_threadsizemax("2147483648", 0);
    //setalldebug(0);
    //--------------------------------------------------
    
    GEN p, K, f, Kcyc, p_ClFld_pol, J_vect, Ja_vect, D, D_prime_vect;

    // Read the prime number p from arguments
    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    
    // Read the defining polynomial for K
    f = gp_read_str(argv[2]);
    pari_printf("K pol: %Ps\n\n", f);
    
    //--------------------------------------------------
    // Define base field K
    K = Buchall(f, nf_FORCE, DEFAULTPREC);

    //--------------------------------------------------
    // Discriminant
    D = nf_get_disc(bnf_get_nf(K));
    pari_printf("Discriminant: %Ps\n\n", D);
    pari_printf("Root discriminant: %Ps\n\n", gsqrtn(gabs(D, DEFAULTPREC), stoi(nf_get_degree(bnf_get_nf(K))), NULL, DEFAULTPREC));
    
    //--------------------------------------------------
    // Check Galois
    if (MY_DEBUGLEVEL >= 1){ my_check_galois(K); }
    
    //--------------------------------------------------        

    // Factor discriminant
    D_prime_vect = gel(factor(D), 1);
    
    //--------------------------------------------------
    // Class group of K (cycle type)
    Kcyc = bnf_get_cyc(K);
    pari_printf("K cyc: %Ps\n\n", Kcyc);

    //--------------------------------------------------
    // Test if p divides the class number. If not, then H^1(X, Z/pZ) = 0 and there is nothing to compute. 
    my_test_p_rank(K, p_int);

    //-------------------------------------------------------------------------------------------------------------------
    // Define polynomials for the generating fields for the part of the Hilbert class field corresp to Cl(K)/p. 
    //-------------------------------------------------------------------------------------------------------------------
    p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);
    DEBUG_PRINT(1, "p Cl Fld: %Ps\n", p_ClFld_pol);
    DEBUG_PRINT(1, ANSI_COLOR_GREEN "Found!\n\n" ANSI_COLOR_RESET);
    
    //--------------------------------------------------
    // Find generators for the p-torsion of the class group
    J_vect = my_find_p_gens(K, p);
    p_rk = glength(J_vect);
    pari_printf("p-rank: %d --> This is the rank of H^1(X,Z/pZ) and H^2(X_fl, mu_p)\n\n", p_rk);
    //--------------------------------------------------

    //--------------------------------------------------
    // find generators for the group of units modulo p
    GEN units_mod_p = my_find_units_mod_p(K, p);
    DEBUG_PRINT(1, "Nr of units mod p: %ld\n", glength(units_mod_p));
    //--------------------------------------------------
    // Define r_rk -- the rank of H^2(X, Z/pZ)
    r_rk = glength(J_vect)+glength(units_mod_p);
    pari_printf("r-rank: %d --> This is the rank of H^2(X,Z/pZ) and H^1(X_fl, mu_p)\n\n", r_rk);
    //--------------------------------------------------

    //--------------------------------------------------
    // Define the extensions generating the p-part of the Hilbert class field corresponding to CL(K)/p
    GEN K_ext = my_ext(K, p_ClFld_pol, p, p_rk, D_prime_vect);
    // pari_printf("Extensions found\n\n");
    //--------------------------------------------------

    //--------------------------------------------------
    // Find generators for H^1(X, mu_p), which is dual to H^2(X, Z/pZ)
    Ja_vect = my_find_Ja_vect(K, J_vect, p, units_mod_p);
    //pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    //--------------------------------------------------
    
    //--------------------------------------------------
    // CUP PRODUCTS
    //--------------------------------------------------
    // Defines a matrix over F_p with index (i*k, j) corresponding to 
    // < x_i\cup x_k, (a_j, J_j) > if i is not equal to j and
    // < B(x_i), (a_j, J_j) > if i=j. 
    // Here < - , - > denotes the Artin--Verdier pairing, which may be computed using our cup product formula and the Artin symbol. 
    //--------------------------------------------------

    //--------------------------------------------------
    // Non-parallel version
    // int mat_rk = my_relations(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk);
    //--------------------------------------------------
    // Parallell version
    int mat_rk = my_relations_par(K_ext, K, p, p_rk, Ja_vect, r_rk);
    //---------------------

    //--------------------------------------------------
    // HIGHER MASSEY PRODUCTS (This is only implemented for some Massey products of the form < x, x, ..., x, y > as seen below, but still very useful)
    //--------------------------------------------------
    // Defines a matrix over F_p with index (i*k, j) corresponding to 
    // < x_i, x_i, ..., x_i, x_k, (a_j, J_j) > if i is not equal to j and
    //--------------------------------------------------

    if ((mat_rk<3 && p_int>2) || (mat_rk==0 && p_int==2))
    {
        my_print_massey(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk);
    }    
    
    //--------------------------------------------------

    DEBUG_PRINT(0, ANSI_COLOR_GREEN "Done! \n \n" ANSI_COLOR_YELLOW);
   
    pari_close();








    //--------------------------------------------------
    // Compute the CPU time
    clock_t duration = (clock()-start) / 1000; // Compute CPU duration in microseconds

    // Compute the actual time
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    
    // Compute actual duration in milliseconds
    long duration_ns = (end_time.tv_sec - start_time.tv_sec) * 1e9 + (end_time.tv_nsec - start_time.tv_nsec);
    long duration_ms = duration_ns / 1e6; // Convert to milliseconds

    // Convert to minutes, seconds, and milliseconds
    msec = duration_ms % 1000;
    sec = (duration_ms / 1000) % 60;
    min = duration_ms / 60000;

    // Print actual elapsed time
    printf(ANSI_COLOR_YELLOW "Actual time: %d min, %d sec, %d msec\n\n" ANSI_COLOR_RESET, min, sec, msec);
    
    // Compute the CPU time
    msec = duration%1000000;
    sec = (duration/1000)%60;
    min = duration/60000;

    printf (ANSI_COLOR_YELLOW "CPU time: %d min, %d,%d sec\n" ANSI_COLOR_RESET, min, sec, msec);
    //--------------------------------------------------
    printf("\n---------------------------------------------------------------------------------------------------------\nEnd program\n---------------------------------------------------------------------------------------------------------\n\n");
    return 0;
}
