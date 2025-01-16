

void print_pari_type(GEN x) {
    switch (typ(x)) {
        case t_INT:    printf("Type: t_INT (integer)\n"); break;
        case t_REAL:   printf("Type: t_REAL (real number)\n"); break;
        case t_FRAC:   printf("Type: t_FRAC (fraction)\n"); break;
        case t_COMPLEX: printf("Type: t_COMPLEX (complex number)\n"); break;
        case t_VEC:    printf("Type: t_VEC (vector)\n"); break;
        case t_COL:    printf("Type: t_COL (column vector)\n"); break;
        case t_MAT:    printf("Type: t_MAT (matrix)\n"); break;
        case t_POL:    printf("Type: t_POL (polynomial)\n"); break;
        case t_SER:    printf("Type: t_SER (power series)\n"); break;
        case t_QUAD:   printf("Type: t_QUAD (quadratic element)\n"); break;
        case t_INTMOD: printf("Type: t_INTMOD (modular integer)\n"); break;
        case t_FFELT:  printf("Type: t_FFELT (finite field element)\n"); break;
        case t_PADIC:  printf("Type: t_PADIC (p-adic number)\n"); break;
        case t_RFRAC:  printf("Type: t_RFRAC (rational fraction)\n"); break;
        case t_STR:    printf("Type: t_STR (string)\n"); break;
        default:       printf("Type: Unknown type\n"); break;
    }
}

void my_test_p_rank (GEN K, int p_int) {
    if (!dvdiu(bnf_get_no(K), p_int))
    {
        pari_printf("%d does not divide the class number %Ps\n", p_int, bnf_get_no(K));
        pari_close();
        exit(0);
    }
}

void my_check_galois(GEN K) {
    GEN gal = galoisconj(K, NULL);
    
    if (glength(gal)==nf_get_degree(bnf_get_nf(K)))
    {
        pari_printf(ANSI_COLOR_GREEN "\n------------------------\nK is Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
    }
    else {
        pari_printf(ANSI_COLOR_RED "\n------------------------\nK is not Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
    }
}

GEN concatenate_rows(GEN M1, GEN M2) {
    pari_sp av = avma;

    // Create a new matrix with (rows1 + rows2) rows and the same number of columns
    GEN M = cgetg(lg(M1), t_MAT);
    for (long i = 1; i < lg(M1); i++)
    {
        gel(M, i) = shallowconcat(gel(M1, i), gel(M2, i));
    }

    return gerepilecopy(av, M);
}


// GEN my_xvector (int l, int x) {
//     pari_sp av = avma;
//     GEN vec = zerovec(l);
//     int i; 
//     GEN a = stoi(x);
//     for (i = 1; i < l+1; i++)
//     {
//         gel(vec, i) = a;
//     }
//     return gerepilecopy(av, vec);
    
// }

// GEN my_xcol (int l, int x) {
//     pari_sp av = avma;
//     GEN vec = zerocol(l);
//     int i; 
//     GEN a = stoi(x);
//     for (i = 1; i < l+1; i++)
//     {
//         gel(vec, i) = a;
//     }
//     return gerepilecopy(av, vec);
    
// }

GEN my_int_to_frac_vec (GEN v) {
    pari_sp av = avma;
    GEN vec = zerovec(glength(v));
    GEN frac = cgetg(3, t_FRAC);
    int i;
    for (i = 1; i < glength(v)+1; ++i) {
        gel(frac, 1) = gel(v,i);
        gel(frac, 2) = gen_1;
        gel(vec,i) = frac;
    }
    return gerepilecopy(av, vec);
}

GEN my_QC_add (GEN v1, GEN v2) {
    pari_sp av = avma;
    if (!(glength(v1)==glength(v2))) {
        DEBUG_PRINT(1, ANSI_COLOR_RED "ERROR in my_ZC_add: vectors has different lengths\n\n" ANSI_COLOR_RESET);
        DEBUG_PRINT(1, "v1: %Ps\n\nv2: %Ps\n\n", v1, v2);
        exit(111);
    }
    GEN sum = zerocol(glength(v1));
    int i;
    for (i = 1; i < glength(v1)+1; ++i) {
        gel(sum,i) = gadd(gel(v1,i), gel(v2,i));
    }
    return gerepilecopy(av, sum);
}

int my_QV_equal1 (GEN v) {
    
    int output = 1;
    int i;
    if (!gequal1(gel(v,1))) {
            output = 0;
        }
    for (i=2; i<lg(v); ++i) {
        if (!gequal0(gel(v,i))) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal0 (GEN v) {
    
    int output = 1;
    int i;
    for (i=1; i<lg(v); ++i) {
        if (!gequal0(gel(v,i))) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal0_mod_p (GEN v, int p) {
    
    int output = 1;
    int i;
    for (i=1; i<glength(v)+1; ++i) {
        if (!(itos(gel(v,i))%p == 0)) {
            output = 0;
        }
    }
    return output;
}

int my_QV_equal (GEN v1, GEN v2) {
    
    int output = 1;
    int i;
    for (i=1; i<glength(v1)+1; ++i) {
        if (!gequal(gel(v1,i), gel(v2,i))) {
            output = 0;
        }
    }
    return output;
}

int my_SQ_MAT_equal (GEN M1, GEN M2) {
    
    // DEBUG_PRINT(1, ANSI_COLOR_CYAN "my_SQ_MAT_equal\n\n" ANSI_COLOR_RESET);
    int output = 1;
    int i;
    int j;
    
    // outmat(M1);
    // outmat(M2);
    for (i = 1; i < glength(gel(M1, 1))+1; ++i) {

        for (j = 1; j < glength(gel(M1, 1))+1; j++)
        {

            if (!gequal(gel(gel(M1, i), j), gel(gel(M2, i), j))) {
                
                return 0;
            }
        }
        
    }
    return output;
}

/*------------------------------------
 The function 1-sigma_x on ideals
------------------------------------
* L - a number field
* sigma - An auto nfalgtobasis(nf, c[2]) where c = nfgaloisconj(nf);
* I - an ideal of nf given as a matrix in hnf

* output - (1-\sigma_x)(I)

DEFINITION:
-----------------------*/

GEN my_1MS_ideal (GEN L, GEN sigma, GEN I) 
{
    DEBUG_PRINT(1, "\n------------------------\nStart: my_1MS_ideal\n------------------------\n\n");
    pari_sp av = avma;
    GEN ideal = idealmul(L, I, idealinv(L, galoisapply(L, sigma, I)));
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_1MS_ideal\n------------------------\n\n");
    GEN ret = gerepilecopy(av, ideal);
    return ret;
} 



GEN my_find_primes_in_factorization(GEN LyAbs, GEN factorization) {
    pari_sp av = avma;
    int l = glength(gel(factorization, 1));
    GEN primes = zerovec(l);
    GEN es = zerovec(l);
    int i;
    for (i = 1; i < l+1; i++)
    {
        gel(primes,i) = idealhnf0(LyAbs, gel(gel(gel(factorization, 1), i), 1), gel(gel(gel(factorization, 1), i), 2));
        gel(es, i) = gel(gel(factorization, 2), i);
    }
    GEN primes_and_es = mkvec2(primes, es);
    primes_and_es = gerepilecopy(av, primes_and_es);
    return primes_and_es;
}

GEN my_1MS_operator_2 (GEN Labs, GEN Lbnr, GEN sigma, int n) {
    DEBUG_PRINT(1, "\n------------------------\nStart: my_1MS_operator_2\n------------------------\n\n");
    pari_sp av = avma;
    GEN cyc = bnf_get_cyc(Labs); 

    // Compute how sigma acts on the generators of Cl(L)/pCl(L)
    GEN sigma_matrix = bnrgaloismatrix(Lbnr, sigma);
    GEN id_mat = matid(glength(sigma_matrix));
    GEN M0 = ZM_sub(id_mat, sigma_matrix);
    M0 = ZM_ZV_mod(M0, cyc);
    GEN M = gcopy(M0);

    // We have found the matrix M that represents the operator 1-sigma
    // Next compute the operator (1-sigma)^n
    for (int k = 2; k < n+1; k++)
    {
        M = ZM_mul(M0, M);
        M = ZM_ZV_mod(M, cyc);
    }

    M = gerepilecopy(av, M);
    
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_1MS_operator_2\n------------------------\n\n");
    return M;
}



GEN my_rel_norm_compact(GEN Labs, GEN Lrel, GEN K, GEN compact_elt) {
    DEBUG_PRINT(1, "\n------------------------\nStart: my_rel_norm_compact\n------------------------\n\n");
    pari_sp av = avma;
    int i;
    GEN norm = gcopy(compact_elt), rel_elt;
    // DEBUG_PRINT(1, "%Ps\n", norm);
    for (i = 1; i < lg(gel(compact_elt, 1)); i++)
    {
        rel_elt = rnfeltabstorel(Lrel, gmael(compact_elt, 1, i));
        gmael(norm, 1, i) = algtobasis(K, rnfeltnorm(Lrel, rel_elt));
    }
    
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_rel_norm_compact\n------------------------\n\n");
    return gerepilecopy(av, norm);
}


// Creates a set (vector) of of exponents in bijection with Cl(L) 
GEN my_get_vect (int n, GEN cyc)
{
    pari_sp av = avma;
    int b = itos(gel(cyc, n+1));
    GEN next_vect;
    
    if (n > 0) {
        GEN prev_vect = my_get_vect(n-1, cyc);
        
        int l = glength(prev_vect);
        next_vect = zerovec(b*l);
        //DEBUG_PRINT(1, "%Ps\n", next_vect);
        
        //DEBUG_PRINT(1, "%ld\n", glength(next_vect));
        int i;
        for (i = 0; i < b*l; ++i) {
            double num = i+1;
            double den = b;
            int first_index = ceil(num/den);
            //DEBUG_PRINT(1, "%d\n",first_index);
            gel(next_vect, i+1) = gconcat(gel(prev_vect, first_index), mkvec(stoi((i+1)%b)));
            
        }

    }
    else {
        
        next_vect = zerovec(b);

        int i;
        for (i = 0; i<b; ++i) {
            gel(next_vect, i+1) = stoi(i);
        }
    }
    return gerepilecopy(av, next_vect);
}

GEN my_get_clgp (GEN K)
{
    pari_sp av0 = avma;
    DEBUG_PRINT(1, "-------\n\nComputing the class group\n\n---------\n");
    GEN Kcyc = bnf_get_cyc(K);
    GEN Kgen = bnf_get_gen(K);
    GEN class_number = bnf_get_no(K);
    int clnr = itos(class_number);
    int nr_comp = glength(Kcyc);
    GEN class_group_exp;
    int n;
    if (nr_comp > 1) {
        class_group_exp = my_get_vect( nr_comp - 1, Kcyc );
    }
    else {
        class_group_exp = zerovec(clnr);
        for (n=0; n<clnr; ++n) {
            gel(class_group_exp, n+1) = mkvec(stoi(n));
        }
    }
    // DEBUG_PRINT(1, "cyc: %Ps\n\n", Kcyc);
    // DEBUG_PRINT(1, "cl_gp_exp: %Ps\n\n", class_group_exp);
    GEN class_group = zerovec(clnr);
    GEN current_I, exponents, pow;
    
    
    for (n = 1; n < clnr + 1; ++n) {
        exponents = gel(class_group_exp, n);
        
        int i;
        current_I = idealhnf0(K, gen_1, NULL);
        for ( i = 1; i < nr_comp + 1; ++i ) {
            
            pow = idealpow(K, gel(Kgen, i), gel(exponents, i));
            
            current_I = idealmul(K, current_I, pow);
            // DEBUG_PRINT(1, "ideal norm: ");
            // output(idealnorm(K, current_I));
            // DEBUG_PRINT(1, "\n");
        }
        if (n%1000 == 0) {
            DEBUG_PRINT(1, "%d/%d\n", n, clnr);
        }
        
        gel(class_group, n) = idealred0(K, current_I, NULL); 
    }
    
    class_group = gerepilecopy(av0, class_group);
    return class_group;
}

GEN my_vect_from_exp (GEN basis, GEN exp) {
    //DEBUG_PRINT(1, "\nmy_vect_from_exp\n");
    pari_sp av = avma;
    int n = glength(gel(basis, 1));
    GEN vect = zerocol(n);
    //DEBUG_PRINT(1, "exp: %Ps\n", exp);
    int l = glength(exp);

    int i;
    for (i = 1; i <= l; i++)
    {
        
        vect = ZC_add(vect, ZC_Z_mul(gel(basis, i), gel(exp, i)));
    }
    
    vect = gerepilecopy(av, vect);
    return vect;

}


//------------------------------
// Returns all ideals I (actually, one solution I_0 plus a matrix M s.t. all solutions can be obtained by adding a lin. comb. of the columns of the matrix to I_0) in Div(L) such that iJ = (1-sigma)I in Cl(L)
//------------------------------ 
GEN my_H90_2 (GEN L, GEN iJ, GEN oneMS_operator, int n) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_H90_2\n--------------------------\n\n");
    pari_sp av = avma;
    GEN B, E, D;


    // Cycle type of Cl(L) as a column vector
    D = shallowtrans(bnf_get_cyc(L));

    // finding exponents for iJ when written as a products of ideals in gens
    B = bnfisprincipal0(L, iJ, 0);
    
    // Gauss: Solving the system M*X = B (i.e., finding I s.t. (1-sigma)I = iJ)
    E = matsolvemod(oneMS_operator,D,B,1);
    
    if (gequal0(E))
    {
        DEBUG_PRINT(1, ANSI_COLOR_RED "ERROR: No solution found in my_H90_2\n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }

    E = gerepilecopy(av, E);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_H90_2\n--------------------------\n\n");
    return E;
}




GEN my_find_units_mod_p (GEN K, GEN p) {
    
    GEN fund_units = bnf_get_fu(K);
    GEN tors_unit = bnf_get_tuU(K);
    int tors_gp_ord = bnf_get_tuN(K);

    GEN units_mod_p = fund_units;
    // DEBUG_PRINT(1, "A: %Ps\n", gpow(p, stoi(4), DEFAULTPREC));
    // DEBUG_PRINT(1, "B: %Ps\n", nfpow(K, tors_unit, gpow(p, stoi(4), DEFAULTPREC)));
    // DEBUG_PRINT(1, "tors_unit: %Ps\n", tors_unit);
    // DEBUG_PRINT(1, "fund_units: %Ps\n", fund_units);
    // DEBUG_PRINT(1, "torsion order: %d\n\n", tors_gp_ord);

    if (tors_gp_ord%itos(p)==0) {
        units_mod_p = shallowconcat(units_mod_p, mkvec(tors_unit));
    }
    
    // DEBUG_PRINT(1, "Units:\n\n");
    // for (i=1;i<glength(units_mod_p)+1;++i) {
    //     DEBUG_PRINT(1, "u %d: %Ps\n", i, algtobasis(K, gel(units_mod_p, i)));
    //     DEBUG_PRINT(1, "u^-1 %d: %Ps\n", i, nfinv(K, algtobasis(K, gel(units_mod_p, i))));
    // }
    // DEBUG_PRINT(1, "\n");

    // ad-hoc modification
    // gel(units_mod_p, 2) = nfmul(K, gel(units_mod_p, 1), gel(units_mod_p, 4));
    // gel(units_mod_p, 3) = nfdiv(K, gel(units_mod_p, 1), gel(units_mod_p, 4));

    return units_mod_p;
}

GEN my_norm_operator (GEN Labs, GEN Lrel, GEN K, GEN p) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_norm_operator\n--------------------------\n\n");
    pari_sp av0 = avma;
    GEN M, g, g_rel, Ng, Lgens;
    int Klength, Llength, i;
    Lgens = shallowconcat(bnf_get_fu(Labs), bnf_get_tuU(Labs));
    Llength = glength(bnf_get_fu(Labs))+1;
    Klength = glength(bnf_get_fu(K));
    M = zeromatcopy(Klength, Llength);

    for (i = 1; i <= Llength; i++)
    {
        pari_sp av = avma;
        g = gel(Lgens, i);
        g_rel = rnfeltabstorel(Lrel, g);
        Ng = rnfeltnorm(Lrel, g_rel);
        gel(M, i) = bnfisunit0(K, Ng, NULL);
        M = gerepilecopy(av, M);
    }
    

    M = gerepilecopy(av0, M);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_norm_operator\n--------------------------\n\n");
    return M;
}

GEN my_find_p_gens (GEN K, GEN p)
{
    pari_sp av = avma;
    GEN my_n = bnf_get_cyc(K);
    int i;
    GEN p_gens = cgetg(1, t_VEC), current_gen, my_gens = bnf_get_gen(K);

    for (i = 1; i < lg(my_n); i++)
    {
        if (dvdii(gel(my_n, i), p))
        {
            current_gen = idealpow(K,gel(my_gens, i), gdiv(gel(my_n, i),p));
            p_gens = shallowconcat(p_gens, mkvec(idealred0(K, current_gen, NULL)));
        }
        
    }
  
    p_gens = gerepilecopy(av, p_gens);
    return p_gens;
}


// Returns vector of tuples (a,J) with div(a)+pJ = 0.
GEN my_find_Ja_vect(GEN K, GEN J_vect, GEN p, GEN units_mod_p) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_find_Ja_vect\n--------------------------\n\n");
    pari_sp av = avma;
    int l = glength(J_vect);
    int r_rk = l + glength(units_mod_p);
    GEN Ja_vect = zerovec(r_rk);
    GEN a;
    int i;
    
    for (i=1; i<=l; ++i) {
        a = nfinv(K, gel(bnfisprincipal0(K, idealpow(K, gel(J_vect, i), p), 1), 2));
        gel(Ja_vect, i) = mkvec2(a, gel(J_vect, i));
    }
    
    for (i = l+1; i <= r_rk; i++)
    {
        gel(Ja_vect, i) = mkvec2(gel(units_mod_p, i-l), idealhnf0(K, gen_1, NULL));
    }
    
    Ja_vect = gerepilecopy(av, Ja_vect);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_find_Ja_vect\n--------------------------\n\n");
    return Ja_vect;
}



GEN my_get_sums (GEN basis, int p) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_get_sums\n--------------------------\n\n");
    pari_sp av = avma;
    int l = glength(basis), i, ord = p;
    int n = pow(ord, l);
    GEN gp = zerovec(n);
    GEN cyc = zerovec(l);
    
    if (l<2)
    {
        gp = shallowconcat(mkvec(zerocol(glength(gel(basis, 1)))), basis);
        //DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "gp: %Ps\n" ANSI_COLOR_RESET, gp);
        gp = gerepilecopy(av, gp);
        return gp;
    }
    

    else {
        for (i = 1; i <= l; i++)
        {
            gel(cyc, i) = stoi(ord);
        }
        
        GEN exp = my_get_vect( l - 1, cyc );
        //DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "exp_list: %Ps\n" ANSI_COLOR_RESET, exp);
        for (i = 1; i <= n; i++)
        {
            gel(gp, i) = my_vect_from_exp(basis, gel(exp, i));
        }
    }
    // DEBUG_PRINT(1, ANSI_COLOR_RED "Basis: %Ps\n" ANSI_COLOR_RESET, basis);
    // DEBUG_PRINT(1, ANSI_COLOR_RED "Sums: %Ps\n" ANSI_COLOR_RESET, gp);
    gp = gerepilecopy(av, gp);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_get_sums\n--------------------------\n\n");
    return gp;
}







//-------------------------------------------------------------------------------------------------
// The function my_H90_vect finds for each (a, J) in Ja_vect, a fractional ideal I in Div(L) such that 
// i(J) = (1-sigma)I + div(t), where t in L^x satisfies N(t)*a = 1.  
// Returns: a vector of these I's 
//-------------------------------------------------------------------------------------------------
GEN my_H90_vect_2 (GEN Labs, GEN Lrel, GEN Lbnr, GEN K, GEN sigma, GEN Ja_vect, GEN p, int n) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_H90_vect_2\n--------------------------\n\n");
    pari_sp av0 = avma;
    int r_rk = glength(Ja_vect), f, i,j, k, done = 0;
    GEN I_vect = zerovec(r_rk), a, iJ, F, ker_T, ker_T_basis, F_ker_T, t_fact, Nt, diff, exp, is_princ, is_norm, Nt_a, ideal, norm_operator, I_fact, cyc = shallowtrans(bnf_get_cyc(Labs));
    // GEN oneMS_operator = my_1MS_operator(Labs, sigma, n);
    GEN oneMS_operator = my_1MS_operator_2(Labs, Lbnr, sigma, n);
    
    for (i = 1; i <= r_rk; ++i)
    {
        
        pari_sp av1 = avma;
        done = 0;
        a = gel(gel(Ja_vect, i), 1);
        iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);

        // Find [I, M] s.t. (1-sigma)I = iJ in Cl(L) and M is a matrix s.t. ...
        if (n==1)
        {
            F = my_H90_2(Labs, iJ, oneMS_operator, n);
        }
        else 
        {
            F = my_H90_2(Labs, idealinv(Labs, iJ), oneMS_operator, n);
        }

        // Then (1-sigma)I+div(t) = iJ for some t in L^x. However, it might be the case that N(t)*a is not 1. 
        // We know from theory that there should be an I and a t s.t (1-sigma)I+div(t) = iJ and N(t)*a=1. 
        // If I' and t' are such, then I-I' is in ker (1-sigma) so in order to find the correct I' we need to go through the kernel of (1-sigma) and for each I'' in there, take I+I'' and check f the corresponding t' satisfies N(t')*a = 1.   
        
        // Generators for the group mapping via (1-sigma) to iJ
        ker_T_basis = gtovec(gel(F, 2));
        
        // If the kernel is trivial, then we are done. 
        if (glength(ker_T_basis)==0)
        {
            gel(I_vect, i) = idealfactorback(Labs, mkmat2(gtocol(bnf_get_gen(Labs)), gel(F, 1)), NULL, 0);
            done = 1;
        }
        else {
            // Generate part of the kernel (or more precisely, the corresponding exponents)
            // This can probably be improved ...
            ker_T = my_get_sums(ker_T_basis, itos(p));
            
            f = glength(ker_T);
            DEBUG_PRINT(1, "Searching a chunk of ker (1-sigma) of size: %d\n", f);
            
            for (j = 1; j <= f; j++)
            {
                pari_sp av2 = avma;
                DEBUG_PRINT(1, "\nSearching: %d/%d\n", j, f);

                //------------------------------------------------------------------------------------------------
                // This is our I+I'' as explained above
                I_fact = ZV_ZV_mod(gadd(gel(F, 1), gtocol(gel(ker_T, j))), cyc);
                F_ker_T = idealred0(Labs, idealfactorback(Labs, mkmat2(gtocol(bnf_get_gen(Labs)), I_fact), NULL, 0), NULL);


               //------------------------------------------------------------------------------------------------
                // Now find the t satisfying (1-sigma)I+div(t) = iJ or (1-sigma)I+div(t) = iJ^{-1} if n>1
                // flag nf_GENMAT: Return t in factored form (compact representation), as a small product of S-units for a small set of finite places S, possibly with huge exponents. This kind of result can be cheaply mapped to K^*/(K^*)^l or to C or Q_p to bounded accuracy and this is usually enough for applications.

                ideal = F_ker_T;
                for (k = 1; k <= n; k++)
                {
                    ideal = my_1MS_ideal(Labs, sigma, ideal);
                }
                DEBUG_PRINT(1, "\n------------------------------------------------------------------------\nSTART: bnfisprincipal0 to find t in compact form\n------------------------------------------------------------------------\n");
                if (n==1)
                {
                    is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, ideal), nf_GENMAT);
                }
                else {
                    is_princ = bnfisprincipal0(Labs, idealmul(Labs, iJ, ideal), nf_GENMAT);
                }
                DEBUG_PRINT(1, "\n------------------------------------------------------------------------\nEND: bnfisprincipal0 to find t in compact form\n------------------------------------------------------------------------\n");

                // Sanity check
                if (!ZV_equal0(gel(is_princ, 1)))
                {
                    DEBUG_PRINT(1, ANSI_COLOR_RED "Problem in my_H90_vect_2\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }

                // The corresponding t in compact/factored form
                t_fact = gel(is_princ, 2);
                
                //------------------------------------------------------------------------------------------------
                // The norm of t in compact/factored form
                Nt = my_rel_norm_compact(Labs, Lrel, K, t_fact);
                // DEBUG_PRINT(1, "Nt: %Ps\n", Nt);

                // create [a, 1] 
                Nt_a = cgetg(3, t_MAT);
                gel(Nt_a, 1) = mkcol(a);
                gel(Nt_a, 2) = mkcol(gen_1);

                // N(t)*a
                diff = concatenate_rows(Nt, Nt_a);
                //------------------------------------------------------------------------------------------------

                //------------------------------------------------------------------------------------------------
                // Since div(a)+pJ = 0 and div(N(t))-pJ = 0, we get that div(N(t)*a) = 0 and therefore N(t)*a must be a unit.
                // Next we find its exponents in terms of the fixed generators of the unit group 
                exp = bnfisunit0(K, diff, NULL);
                if (glength(exp)==0)
                {
                    DEBUG_PRINT(1, ANSI_COLOR_RED "PROBLEM in my_H90_vect_2: no exp???\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }
                
                DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "exp: %Ps\n" ANSI_COLOR_RESET, exp);
                
                //------------------------------------------------------------------------------------------------
                // Check if N(t)*a = 1
                if (ZV_equal0(exp))
                {

                    gel(I_vect, i) = F_ker_T;
                    DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                    done = 1;
                    break;
                }
                //------------------------------------------------------------------------------------------------


                //------------------------------------------------------------------------------------------------
                // Check if N(t)*a is the norm of a unit. If it is, we may modify t by this unit without effecting the equality 
                // (1-sigma)I = iJ and hence we are done. Returns zero if not a norm (which should never happen).
                norm_operator = my_norm_operator(Labs, Lrel, K, p);
                DEBUG_PRINT(1, "Norm operator: %Ps\n", norm_operator);
                DEBUG_PRINT(1, "Checking if N(t)*a is a norm\n");
                is_norm = matsolvemod(norm_operator, zerocol(glength(exp)), gtocol(exp), 0);
                DEBUG_PRINT(1, ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);

                // if N(t)*a is a norm
                if (!gequal0(is_norm))
                {

                    gel(I_vect, i) = F_ker_T;
                    DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                    done = 1;
                    break;
                }
                I_vect = gerepilecopy(av2, I_vect);
                //------------------------------------------------------------------------------------------------
            }
        }
        // If some of the I's were never found, then we return -1 and another (slower) function will take over.  
        if (!done)
        {
            DEBUG_PRINT(1, ANSI_COLOR_RED "my_H90_vect_2 ended with a problem \n" ANSI_COLOR_RESET);
            return stoi(-1);
            // pari_close();
            // exit(111);
        }
        I_vect = gerepilecopy(av1, I_vect);
    } 
    
    I_vect = gerepilecopy(av0, I_vect);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_H90_vect_2\n--------------------------\n\n");
    return I_vect;
}



// GEN my_H90_12 (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN a, GEN J, GEN p, int n) {
//     DEBUG_PRINT(1, ANSI_COLOR_CYAN "\nmy_H90_vect\n" ANSI_COLOR_RESET);
//     pari_sp av = avma;
//     int f, j, done = 0;
//     GEN I, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, t, t_rel, Nt, diff, exp, is_princ, is_norm, cyc;
//     cyc = gtocol(bnf_get_cyc(Labs));
    
    

//     done = 0;
//     iJ = rnfidealup0(Lrel, J, 1);
//     F = my_H90(Labs, iJ, sigma, 1);
//     // DEBUG_PRINT(1, "F: %Ps\n", F);
//     ker_sol = matsolvemod(my_1MS_operator(Labs, sigma, 1), cyc, zerocol(glength(cyc)), 1);
//     ker_T_basis = gtovec(gel(ker_sol, 2));
//     DEBUG_PRINT(1, ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
    
//     if (glength(ker_T_basis)==0)
//     {
//         I = F;
//         done = 1;
//     }
//     else {
//         ker_T = my_get_sums(ker_T_basis, itos(p));

//         // WARNING: This makes it faster but might not find the answer. Switch back to the above ker_T
//         //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
        
//         //DEBUG_PRINT(1, ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
//         // DEBUG_PRINT(1, "ker_T[2]: %Ps\n", gel(ker_T, 2));
//         f = glength(ker_T);
//         DEBUG_PRINT(1, "Searching a chunk of ker (1-sigma) of size: %d\n", f);
//         /*
//             Find F_ker_T --------------------
//         */
        
//         for (j = 1; j <= f; j++)
//         {
//             DEBUG_PRINT(1, "\nSearching: %d/%d\n", j, f);
//             // DEBUG_PRINT(1, "Adding the exp: %Ps\n", gel(ker_T, j));
//             // DEBUG_PRINT(1, "And the ideal: %Ps\n", my_ideal_from_exp(Labs, gel(ker_T, j)));
//             F_ker_T = idealred0(Labs, idealmul(Labs, F, idealred0(Labs, my_ideal_from_exp(Labs, gel(ker_T, j)), NULL)), NULL);
//             //DEBUG_PRINT(1, "F_ker_T: %Ps\n", F_ker_T);
//             is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, F_ker_T)), 1);
//             if (!my_QV_equal0(gel(is_princ, 1)))
//             {
//                 DEBUG_PRINT(1, ANSI_COLOR_RED "Problem in my_H90_vect\n" ANSI_COLOR_RESET);
//                 pari_close();
//                 exit(111);
//             }
            
//             t = gel(is_princ, 2);
//             DEBUG_PRINT(1, "t: %ld\n", glength(t));
//             if (glength(t)==0)
//             {
//                 return stoi(-1);
//             }
            
//             t_rel = rnfeltabstorel(Lrel, t);
//             Nt = rnfeltnorm(Lrel, t_rel);
//             diff = nfmul(K, Nt, a);
//             exp = bnfisunit0(K, diff, NULL);
//             // DEBUG_PRINT(1, ANSI_COLOR_YELLOW "exp: %Ps\n" ANSI_COLOR_RESET, exp);
//             // DEBUG_PRINT(1, "M: %Ps\n", my_norm_operator(Labs, Lrel, K, p));
//             // DEBUG_PRINT(1, "D: %Ps\n", my_xcol(glength(exp), itos(p)));
//             // DEBUG_PRINT(1, "B: %Ps\n", gtocol(exp));
//             /* check if exp lies in the image of the operator associated to N: (O_L)^x -> (O_K)^x */
//             // if (my_QV_equal0(my_norm_operator(Labs, Lrel, K, p)))
//             // {
//             //     gel(I_vect, i) = F_ker_T;
//             //     done = 1;
//             //     break;
//             // }
//             is_norm = matsolvemod(my_norm_operator(Labs, Lrel, K, p), zerocol(glength(exp)), gtocol(exp), 0);
//             //DEBUG_PRINT(1, ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);
//             if (my_QV_equal0(exp) || !gequal0(is_norm))
//             {
//                 //DEBUG_PRINT(1, ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
//                 // w = gel(bnf_get_fu(Labs), 1);
//                 // w_rel = rnfeltabstorel(Lrel, w);
//                 // Nw = rnfeltnorm(Lrel, w_rel);
//                 // DEBUG_PRINT(1, "Norm w: %Ps\n", Nw);
//                 I = F_ker_T;
//                 done = 1;
//                 break;
//             }
//         }
//     }
//     if (!done)
//     {
//         DEBUG_PRINT(1, ANSI_COLOR_RED "my_H90_vect ended with a problem \n" ANSI_COLOR_RESET);
//         return stoi(-1);
//         // pari_close();
//         // exit(111);
//     }
    
    
//     I = gerepilecopy(av, I);
//     return I;
// }

