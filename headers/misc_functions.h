

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


GEN my_xvector (int l, int x) {
    pari_sp av = avma;
    GEN vec = zerovec(l);
    int i; 
    GEN a = stoi(x);
    for (i = 1; i < l+1; i++)
    {
        gel(vec, i) = a;
    }
    return gerepilecopy(av, vec);
    
}

GEN my_xcol (int l, int x) {
    pari_sp av = avma;
    GEN vec = zerocol(l);
    int i; 
    GEN a = stoi(x);
    for (i = 1; i < l+1; i++)
    {
        gel(vec, i) = a;
    }
    return gerepilecopy(av, vec);
    
}

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

int my_SQ_MAT_equal0 (GEN M) {
    
    // DEBUG_PRINT(1, ANSI_COLOR_CYAN "my_SQ_MAT_equal\n\n" ANSI_COLOR_RESET);
    int output = 1;
    int i;
    int j;
    
    for (i = 1; i < glength(gel(M, 1))+1; ++i) {
        for (j = 1; j < glength(gel(M, 1))+1; j++)
        {
            if (!gequal0(gel(gel(M, i), j))) {
                output = 0;
                break;
            }
        }
        
    }
    return output;
}

GEN my_norm (GEN Labs, GEN Lrel, GEN L, GEN elem, GEN sigma, int p) {
    pari_sp av = avma;
    int i;
    GEN new_elem = algtobasis(Labs, gen_1);
    for (i=0; i<p; ++i) {
        new_elem = basistoalg(Labs, nfmul(Labs, new_elem, elem));
        
        elem = galoisapply(Labs, sigma, elem);
        // DEBUG_PRINT(1, "new_elem: %Ps\n\nelem: %Ps\n\n", new_elem, elem);
    }
    
    GEN norm = algtobasis(L, rnfeltdown(Lrel, rnfeltabstorel(Lrel, new_elem)));
    
    return gerepilecopy(av, norm);
}

GEN my_norm_ideal (GEN Labs, GEN Lrel, GEN L, GEN I, GEN sigma, int p) {
    pari_sp av = avma;
    int i;
    GEN new_ideal = algtobasis(Labs, gen_1);
    for (i=0; i<p; ++i) {
        new_ideal = basistoalg(Labs, nfmul(Labs, new_ideal, I));
        
        I = galoisapply(Labs, sigma, I);
        // DEBUG_PRINT(1, "new_elem: %Ps\n\nelem: %Ps\n\n", new_elem, elem);
    }
    
    GEN norm = rnfidealdown(Lrel, rnfidealabstorel(Lrel, new_ideal));
    // DEBUG_PRINT(1, "Hoooj\n\n");
    norm = gerepilecopy(av, norm);
    return norm;
}

GEN my_i (GEN Labs, GEN Lrel, GEN sigma, GEN elem) {
    pari_sp av = avma;
    GEN lift = rnfeltup0(Lrel, elem, 1);
    output(lift);
    GEN mod_lift = galoisapply(Labs, sigma, lift);
    output(mod_lift);
    DEBUG_PRINT(1, "\n");
    return gerepilecopy(av, mod_lift);
}

GEN my_i_ideal (GEN Labs, GEN Lrel, GEN sigma, GEN ideal) {
    pari_sp av = avma;
    GEN lift = rnfidealup0(Lrel, ideal, 1);
    GEN mod_lift = galoisapply(Labs, sigma, lift);
    return gerepilecopy(av, mod_lift);
}

/*------------------------------------
 The function sigma_x-1 on ideals
------------------------------------
* L - a number field
* sigma - An auto nfalgtobasis(nf, c[2]) where c = nfgaloisconj(nf);
* I - an ideal of nf given as a matrix in hnf

* output - (\sigma_x-1)(I)

DEFINITION:
-----------------------*/

GEN my_SM1_ideal (GEN L, GEN sigma, GEN I)
{
    return idealmul(L, galoisapply(L, sigma, I), idealinv(L, I));
    // return idealred(L, idealmul(L, galoisapply(L, sigma, I), idealinv(L, I)));
} 

GEN my_1MS_ideal (GEN L, GEN sigma, GEN I) 
{
    DEBUG_PRINT(1, "\n------------------------\nStart: my_1MS_ideal\n------------------------\n\n");
    pari_sp av = avma;
    GEN ideal = idealmul(L, I, idealinv(L, galoisapply(L, sigma, I)));
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_1MS_ideal\n------------------------\n\n");
    GEN ret = gerepilecopy(av, ideal);
    return ret;
} 

/*
SM1 for elements
*/
GEN my_SM1_elt (GEN L, GEN sigma, GEN elem) 
{
    return nfdiv(L, galoisapply(L, sigma, elem), elem);
}

GEN my_1MS_elt (GEN L, GEN sigma, GEN elem) 
{
    return nfdiv(L, elem, galoisapply(L, sigma, elem));
}

GEN my_Gamma (GEN Labs, GEN sigma, GEN l, int p_int) 
{
    pari_sp av = avma;
    int n;
    int i;

    GEN Gamma_l = algtobasis(Labs, gen_1);
    for (n = 1; n < p_int; ++n)
    {
        GEN elem_n = algtobasis(Labs, l);
        GEN n_gen = stoi(n);
        elem_n = nfpow(Labs, elem_n, n_gen);
        for (i = 1; i < n+1; ++i)
        {
            elem_n = galoisapply(Labs, sigma, elem_n);
        }     
        Gamma_l = nfmul(Labs, Gamma_l, elem_n);
    }
    Gamma_l = nfinv(Labs, Gamma_l);
    Gamma_l = gerepilecopy(av , Gamma_l);
    return Gamma_l;
}

GEN my_Gamma_ideal (GEN Labs, GEN sigma, GEN I, int p_int) 
{
    pari_sp av = avma;
    int n;
    int i;

    GEN Gamma_l = idealhnf0(Labs, gen_1, NULL);
    for (n = 1; n < p_int; ++n)
    {
        GEN elem_n = gcopy(I);
        GEN n_gen = stoi(n);
        elem_n = idealpow(Labs, elem_n, n_gen);
        for (i = 1; i < n+1; ++i)
        {
            elem_n = galoisapply(Labs, sigma, elem_n);
        }     
        Gamma_l = idealmul(Labs, Gamma_l, elem_n);
    }
    Gamma_l = idealinv(Labs, Gamma_l);
    Gamma_l = gerepilecopy(av, Gamma_l);
    return Gamma_l;
}

/* The function p-N*/
GEN my_pmN (GEN Labs, GEN Lrel, GEN z, GEN p)  
{
    GEN elem1 = nfpow(Labs, z, p);
    GEN elem2 = rnfeltup(Lrel, rnfeltnorm(Lrel, z));
    return nfdiv(Labs, elem1, elem2);
}

GEN my_tau (GEN Labs, GEN a, GEN b, GEN sigma)
{
    return nfmul(Labs, a, galoisapply(Labs, sigma, b));
}

GEN my_trace_tau (GEN Labs, GEN a, GEN b, GEN sigma, int p)
{
    pari_sp av = avma;
    GEN new_b = algtobasis(Labs, b);
    GEN trace = my_int_to_frac_vec(algtobasis(Labs, gen_0));
    // DEBUG_PRINT(1, "Trace: %Ps\n\n", trace);
    trace = my_QC_add(trace, new_b);
    // DEBUG_PRINT(1, "Trace: %Ps\n\n", trace);
    
    int i;
    for (i = 1; i<p; ++i) {
        new_b = my_tau(Labs, a, new_b, sigma);
        // DEBUG_PRINT(1, "new_b: %Ps\n\n", new_b);
        trace = my_QC_add(trace, new_b);
        // DEBUG_PRINT(1, "Trace: %Ps\n\n", trace);
    }
    trace = gerepilecopy(av, trace);
    return trace;
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

GEN my_1MS_operator (GEN L, GEN sigma, int n) {
    DEBUG_PRINT(1, "\n------------------------\nStart: my_1MS_operator\n------------------------\n\n");
    pari_sp av0 = avma;
    GEN cyc = bnf_get_cyc(L), col, ideal; 
    GEN gens = bnf_get_gen(L);
    int l = glength(cyc);
    GEN M = zeromatcopy(l,l);
    int i, k;

    for (i = 1; i < l+1; i++)
    {
        DEBUG_PRINT(1, "i: %d\n", i);
        pari_sp av = avma;
        ideal = gel(gens, i);
        for (k = 1; k < n+1; k++)
        {
            ideal = my_1MS_ideal(L, sigma, ideal);
        }
        DEBUG_PRINT(1, "bnfisprincipal[%d]: \n", i);
        col = bnfisprincipal0(L, ideal, 0);
        DEBUG_PRINT(1, "ideal[%d] in basis: %Ps\n", i, col);

        col = gerepilecopy(av, col);
        gel(M, i) = col;
    }
    
    DEBUG_PRINT(1, ANSI_COLOR_YELLOW "my_1MS_operator: %Ps\n" ANSI_COLOR_RESET, M);
    
    M = gerepilecopy(av0, M);
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_1MS_operator\n------------------------\n\n");
    return M;
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

GEN my_1MS_elt_operator (GEN L, GEN sigma) {
    DEBUG_PRINT(1, "\n------------------------\nStart: my_1MS_elt_operator\n------------------------\n\n");
    pari_sp av0 = avma;
    GEN int_basis = nf_get_zk(bnf_get_nf(L)), col, basis_vec;
    //DEBUG_PRINT(1, "int_basis: %Ps\n", int_basis);
    int l = glength(int_basis);
    GEN M = zeromatcopy(l,l);
    int i;
    
    for (i = 1; i < l+1; i++)
    {
        pari_sp av = avma;
        basis_vec = algtobasis(L, gel(int_basis, i));
        col = algtobasis(L, my_1MS_elt(L, sigma, basis_vec));
        //DEBUG_PRINT(1, "col: %Ps\n", col);
        col = gerepilecopy(av, col);
        gel(M, i) = col;
    }
    M = gerepilecopy(av0, M);
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_1MS_elt_operator\n------------------------\n\n");
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

GEN my_find_ext_ranks(GEN K_ext) {
    pari_sp av = avma;
    int l = glength(K_ext);
    GEN ext_rks = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_rks, i) = stoi(glength(bnf_get_cyc(gel(gel(K_ext, i), 1))));
    }
    
    return gerepilecopy(av,vecsort0(ext_rks, NULL, 0));
}

GEN my_find_ext_cyc(GEN K_ext) {
    pari_sp av = avma;
    int l = glength(K_ext);
    GEN ext_cyc = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_cyc, i) = bnf_get_cyc(gel(gel(K_ext, i), 1));
    }
    
    return gerepilecopy(av, ext_cyc);
}

int my_is_p_power (GEN vect, int p_int, int order) {
    
    int l = glength(vect), check_p = 1, i;

    if (l==0) {
        return 0;
    }
    
    if (itos(gel(vect, l))!=(p_int%order) && itos(gel(vect, l))%p_int!=0) {
            return 0;
        }
    for (i=1;i<l;++i) {
        if (itos(gel(vect, i))%p_int!=0) {
            check_p = 0;
            break;
        }
    } 

    return check_p;
}

GEN my_find_ext_pol(GEN K_ext) {
    pari_sp av = avma;
    int l = glength(K_ext);
    GEN ext_pol = zerovec(l);

    int i;

    for (i = 1; i < l+1; i++)
    {
        gel(ext_pol, i) = nf_get_pol(bnf_get_nf(gel(gel(K_ext, i), 1)));
    }
    return gerepilecopy(av, ext_pol);
}

GEN my_Gal_rel (GEN Labs, GEN Lrel, GEN K, GEN sigma, int p) {
    pari_sp av = avma;
    DEBUG_PRINT(1, "\n------------------------\nStart: my_Gal_rel\n------------------------\n\n");
    
    GEN Gal_rel = zerovec(p);
    DEBUG_PRINT(1, "\nsigma: %Ps\n\n", sigma);
    GEN current_sigma = sigma;
    int i;
    
    for (i = 1; i < p+1; i++)
    {
        gel(Gal_rel, i) = current_sigma;
        current_sigma = galoisapply(Labs, sigma, current_sigma); 
        
    }
    DEBUG_PRINT(1, "\n------------------------\nEnd: my_Gal_rel\n------------------------\n\n");
    return gerepilecopy(av, Gal_rel);
}

GEN my_get_prod(GEN Labs, GEN zk) {
    pari_sp av = avma;
    GEN vect = gcopy(zk);

    int i, j;

    for (i = 1; i < glength(zk)+1; i++)
    {
        for (j = 1; j < glength(zk)+1; j++)
        {
            vect = shallowconcat(vect, mkvec(nfmul(Labs, gel(zk, i), gel(zk, j))));
        }
        
    }
    
    return gerepilecopy(av, vect);
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

GEN my_get_unit_group (GEN K, GEN unit_gens, GEN p)
{
    pari_sp av = avma;
    DEBUG_PRINT(1, "-------\n\nComputing the unit group\n\n---------\n");
    int nr_comp = glength(unit_gens);
    GEN cyc = zerovec(nr_comp);

    int i;
    
    for (i = 1; i <= nr_comp; i++)
    {
        gel(cyc, i) = p;
    }
    
    int nr = pow(itos(p), nr_comp);
    
    GEN group_exp;
    int n;
    if (nr_comp > 1) {
        group_exp = my_get_vect( nr_comp - 1, cyc );
    }
    else {
        group_exp = zerovec(nr);
        for (n=0; n<nr; ++n) {
            gel(group_exp, n+1) = mkvec(stoi(n));
        }
    }
    
    // DEBUG_PRINT(1, "cyc: %Ps\n\n", Kcyc);
    // DEBUG_PRINT(1, "cl_gp_exp: %Ps\n\n", class_group_exp);
    GEN group = zerovec(nr);
    GEN current_u, exponents, pow;
    
    
    for (n = 1; n <= nr; ++n) {
        exponents = gel(group_exp, n);
        
        int i;
        current_u = gen_1;
        
        for ( i = 1; i <= nr_comp; ++i ) {
            
            pow = nfpow(K, gel(unit_gens, i), gel(exponents, i));
            
            current_u = nfmul(K, current_u, pow);
            // DEBUG_PRINT(1, "ideal norm: ");
            // output(idealnorm(K, current_I));
            // DEBUG_PRINT(1, "\n");
        }
        
        if (n%1000 == 0) {
            DEBUG_PRINT(1, "%d/%d\n", n, nr);
        }
        
        gel(group, n) = algtobasis(K, current_u); 
    }
    
    group = gerepilecopy(av, group);
    return group;
}



//------------------------------
// Returns an ideal I in Div(L) such that iJ = (1-sigma)I in Cl(L)
//------------------------------ 
GEN my_H90 (GEN L, GEN iJ, GEN oneMS_operator, int n) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_H90\n--------------------------\n\n");
    pari_sp av = avma;
    GEN H90_ideal, B, E, D;


    // Cycle type of Cl(L) as a column vector
    D = shallowtrans(bnf_get_cyc(L));

    // finding exponents for iJ when written as a products of ideals in gens
    B = bnfisprincipal0(L, iJ, 0);

    // Gauss: Solving the system M*X = B (i.e., finding I s.t. (1-sigma)I = iJ)
    E = matsolvemod(oneMS_operator,D,B,0);
    if (gequal0(E))
    {
        DEBUG_PRINT(1, ANSI_COLOR_RED "ERROR: No solution found in my_H90\n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }
    
    // Defining the ideal I from the exponents E
    // H90_ideal = my_ideal_from_exp(L, E);
    H90_ideal = idealfactorback(L, mkmat2(gtocol(bnf_get_gen(L)), gtocol(E)), NULL, 0);
    DEBUG_PRINT(1, ANSI_COLOR_GREEN "H90_ideal found\n" ANSI_COLOR_RESET);

    H90_ideal = gerepilecopy(av, H90_ideal);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_H90\n--------------------------\n\n");
    return H90_ideal;
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

GEN my_H90_exp (GEN L, GEN iJ, GEN sigma) {
    DEBUG_PRINT(1, "\nmy_H90_exp\n");
    pari_sp av = avma;
    GEN H90_exp, gens, M, B, D;
    gens = bnf_get_gen(L);
    D = gtocol(bnf_get_cyc(L));
    int l = glength(gens);
    
    M = zeromatcopy(l,l);
    //D = zerocol(l);
    // finding exponents for iJ 
    B = bnfisprincipal0(L, iJ, 0);

    // finding the matrix M consisting of exponents for the (1-s)g_i, Cl(L) = < g_1, ..., g_n > 
    int i;
    for (i = 1; i <= l; i++)
    {
        gel(M, i) = bnfisprincipal0(L, idealred0(L, my_1MS_ideal(L, sigma, gel(gens, i)), NULL), 0);
    }
    // Gauss
    // DEBUG_PRINT(1, "M: %Ps\n", M);
    // DEBUG_PRINT(1, "B: %Ps\n", B);
    //DEBUG_PRINT(1, "Gauss start\n");
    H90_exp = matsolvemod(M,D,B,0);
    //DEBUG_PRINT(1, "Gauss done\n");

    H90_exp = gerepilecopy(av, H90_exp);
    return H90_exp;
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


void my_unramified_p_extensions(GEN K, GEN p, GEN D_prime_vect) {
    pari_sp av = avma;
    int i;
    GEN s = pol_x(fetch_user_var("s"));
    GEN x = pol_x(fetch_user_var("x"));
    GEN index = mkvec(p);
    // GEN R = nfsubfields0(clf_pol,4,1);
    // DEBUG_PRINT(1, "subgrouplist: %Ps\n", subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 0));

    GEN R = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), index, 2), 2, DEFAULTPREC);
    //GEN Rsq = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,gen_2, DEFAULTPREC)), 2), 2, DEFAULTPREC);
    //GEN Rcb = bnrclassfield(K, subgrouplist0(bnf_get_cyc(K), mkvec(gpow(p,stoi(3), DEFAULTPREC)), 2), 2, DEFAULTPREC);
    GEN abs_pol;
    DEBUG_PRINT(1, "[p]-extensions:\n\n");
    for (i=1;i<lg(R);i++) {
        //DEBUG_PRINT(1, "i: %d\n", i);
        abs_pol = polredabs0(mkvec2(gel(R, i), D_prime_vect), 0);
        //abs_pol = polredabs(gel(R, i));
        DEBUG_PRINT(1, "%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_GEN, DEFAULTPREC)));
        
    }
    DEBUG_PRINT(1, "\n");
    // DEBUG_PRINT(1, "[p,p]-extensions:\n\n");
    // for (i=1;i<glength(Rsq)+1;i++) {
    //     abs_pol = polredabs0(mkvec2(gel(Rsq, i), D_prime_vect), 0);
    //     //abs_pol = polredabs(gel(R, i));
    //     DEBUG_PRINT(1, "%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_FORCE, DEFAULTPREC)));
    // }
    // DEBUG_PRINT(1, "\n");
    // DEBUG_PRINT(1, "[p,p,p]-extensions:\n\n");
    // for (i=1;i<glength(Rcb)+1;i++) {
    //     abs_pol = polredabs0(mkvec2(gel(Rcb, i), D_prime_vect), 0);
    //     //abs_pol = polredabs(gel(R, i));
    //     DEBUG_PRINT(1, "%Ps, cyc: %Ps\n", gsubstpol(abs_pol, x, s), bnf_get_cyc(Buchall(abs_pol, nf_FORCE, DEFAULTPREC)));
    // }
    // DEBUG_PRINT(1, "\n");
    avma = av;
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
GEN my_H90_vect (GEN Labs, GEN Lrel, GEN Lbnr, GEN K, GEN sigma, GEN Ja_vect, GEN p, int n) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_H90_vect_2\n--------------------------\n\n");
    pari_sp av0 = avma;
    int r_rk = glength(Ja_vect), f, i,j, k, done = 0;
    GEN I_vect = zerovec(r_rk), a, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, I_fact, t, t_fact, Nt, diff, exp, is_princ, is_norm, Nt_a, cyc, ideal, norm_operator, oneMS_operator = my_1MS_operator(Labs, sigma, n);
    cyc = shallowtrans(bnf_get_cyc(Labs));
    
    for (i = 1; i <= r_rk; ++i)
    {
        pari_sp av1 = avma;
        done = 0;
        a = gel(gel(Ja_vect, i), 1);
        iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);

        // Find I s.t. (1-sigma)I = iJ in Cl(L)
        if (n==1)
        {
            F = my_H90(Labs, iJ, oneMS_operator, n);
        }
        else 
        {
            F = my_H90(Labs, idealinv(Labs, iJ), oneMS_operator, n);
        }
        
        // Then (1-sigma)I+div(t) = iJ for some t in L^x. However, it might be the case that N(t)*a is not 1. 
        // We know from theory that there should be an I and a t s.t (1-sigma)I+div(t) = iJ and N(t)*a=1. 
        // If I' and t' are such, then I-I' is in ker (1-sigma) so in order to find the correct I' we need to go through the kernel of (1-sigma) and for each I'' in there, take I+I'' and check f the corresponding t' satisfies N(t')*a = 1.   
        
        // Find generators for the kernel of (1-sigma)
        // ker_sol = matsolvemod(my_1MS_operator_2(Labs, Lrel, K, sigma, n), cyc, zerocol(glength(cyc)), 1);
        ker_sol = matsolvemod(oneMS_operator, cyc, zerocol(glength(cyc)), 1);
        ker_T_basis = gtovec(gel(ker_sol, 2));
        DEBUG_PRINT(1, "ker_T_basis: %Ps\n", ker_T_basis);
        // DEBUG_PRINT(1, ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
        
        // If the kernel is trivial, then we are done. 
        if (glength(ker_T_basis)==0)
        {
            gel(I_vect, i) = F;
            done = 1;
        }
        else {
            // Generate part of the kernel (or more precisely, the corresponding exponents)
            ker_T = my_get_sums(ker_T_basis, itos(p));

            // WARNING: The next line makes it faster but might not find the answer. Switch back to the above ker_T!
            //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
            
            //DEBUG_PRINT(1, ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
            // DEBUG_PRINT(1, "ker_T[2]: %Ps\n", gel(ker_T, 2));
            f = glength(ker_T);
            DEBUG_PRINT(1, "Searching a chunk of ker (1-sigma) of size: %d\n", f);
            /*
                Find F_ker_T --------------------
            */
            
            for (j = 1; j <= f; j++)
            {
                pari_sp av2 = avma;
                DEBUG_PRINT(1, "\nSearching: %d/%d\n", j, f);
                // idealfactorback(Labs, mkmat2(gtocol(bnf_get_gen(L)), gtocol(gel(ker_T, j))), NULL, 0);
                // This is our I+I'' as explained above
                F_ker_T = idealred0(Labs, idealmul(Labs, F, idealred0(Labs, idealfactorback(Labs, mkmat2(gtocol(bnf_get_gen(Labs)), gtocol(gel(ker_T, j))), NULL, 0), NULL)), NULL);
                
                // -----------------------------
                // I_fact = gadd(bnfisprincipal0(Labs, F, 0), gtocol(gel(ker_T, j)));
                
                // for (int m = 1; m < lg(cyc); m++)
                // {
                //     gel(I_fact, m) = modii(gel(I_fact, m), gel(cyc, m));
                    
                // }
                // F_ker_T = idealfactorback(Labs, mkmat2(gtocol(bnf_get_gen(Labs)), I_fact), NULL, 0);
                // -----------------------------

                // Now find the corresponding t
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
                    DEBUG_PRINT(1, ANSI_COLOR_RED "Problem in my_H90_vect\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }

                // The corresponding t
                t_fact = gel(is_princ, 2);
                t = t_fact;
                
                if (glength(t)==0)
                {
                    DEBUG_PRINT(1, "t has length zero: %Ps\n", t);
                    continue;
                    //return stoi(-1);
                }
                
                // ------- Compact form  ----------------
                Nt = my_rel_norm_compact(Labs, Lrel, K, t);
                // DEBUG_PRINT(1, "Nt: %Ps\n", Nt);

                // create [a, 1] 
                Nt_a = cgetg(3, t_MAT);
                gel(Nt_a, 1) = mkcol(a);
                gel(Nt_a, 2) = mkcol(gen_1);

                // N(t)*a
                diff = concatenate_rows(Nt, Nt_a);
                // --------------------------------------
                
                // --------------- Non-compact form ------
                // // In relative coordinates
                // t_rel = rnfeltabstorel(Lrel, t);

                // // N(t)
                // Nt = rnfeltnorm(Lrel, t_rel);

                // // N(t)*a
                // diff = nfmul(K, Nt, a);
                //------------------------------------

                // Since div(a)+pJ = 0 and div(N(t))-pJ = 0, we get that div(N(t)*a) = 0 and therefore N(t)*a must be a unit.
                // Next we find its exponents in terms of the fixed generators of the unit group 
                exp = bnfisunit0(K, diff, NULL);
                if (glength(exp)==0)
                {
                    DEBUG_PRINT(1, ANSI_COLOR_RED "PROBLEM in my_H90_vect: no exp???\n" ANSI_COLOR_RESET);
                    pari_close();
                    exit(111);
                }
                
                DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "exp: %Ps\n" ANSI_COLOR_RESET, exp);

                // Check if N(t)*a is the norm of a unit. If it is, we may modify t by this unit without effecting the equality 
                // (1-sigma)I = iJ and hence we are done. Returns zero if not a norm (which should never happen).
                norm_operator = my_norm_operator(Labs, Lrel, K, p);


                DEBUG_PRINT(1, "Checking if N(t)*a is a norm\n");
                
                DEBUG_PRINT(1, "exp: %Ps\n", exp);
                is_norm = matsolvemod(norm_operator, zerocol(glength(exp)), gtocol(exp), 0);

                DEBUG_PRINT(1, ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);

                // Use ZV_equal0 or gequal0

                // if N(t)*a = 1 or if N(t)*a is a norm
                if (ZV_equal0(exp) || !gequal0(is_norm))
                {
                    //DEBUG_PRINT(1, ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
                    // w = gel(bnf_get_fu(Labs), 1);
                    // w_rel = rnfeltabstorel(Lrel, w);
                    // Nw = rnfeltnorm(Lrel, w_rel);
                    // DEBUG_PRINT(1, "Norm w: %Ps\n", Nw);
                    gel(I_vect, i) = F_ker_T;
                    DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                    done = 1;
                    break;
                }
                I_vect = gerepilecopy(av2, I_vect);
            }
        }
        // If some of the I's were never found, then we return -1 and another (slower) function will take over.  
        if (!done)
        {
            DEBUG_PRINT(1, ANSI_COLOR_RED "my_H90_vect ended with a problem \n" ANSI_COLOR_RESET);
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

// GEN my_H90_exp_vect (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, GEN p, int n) {
//     // This function is not ready to be used
//     DEBUG_PRINT(1, ANSI_COLOR_CYAN "\nmy_H90_exp_vect\n" ANSI_COLOR_RESET);
//     pari_sp av = avma;
//     int r_rk = glength(Ja_vect), f, j, done = 0;
//     GEN I_vect = zerovec(r_rk), a, iJ, F, ker_T, ker_T_basis, ker_sol, F_ker_T, t, t_rel, Nt, diff, exp, is_princ, is_norm, cyc;
//     cyc = shallowtrans(bnf_get_cyc(Labs));
//     int i;
    
//     for (i = 1; i <= r_rk; ++i)
//     {
//         done = 0;
//         a = gel(gel(Ja_vect, i), 1);
//         iJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
//         F = my_H90_exp(Labs, iJ, sigma);
//         // DEBUG_PRINT(1, "F: %Ps\n", F);
//         ker_sol = matsolvemod(my_1MS_operator(Labs, sigma, 1), cyc, zerocol(glength(cyc)), 1);
//         ker_T_basis = gtovec(gel(ker_sol, 2));
//         DEBUG_PRINT(1, ANSI_COLOR_GREEN "------------------------\n\n\nker_T_basis size: %ld\n\n\n------------------------\n" ANSI_COLOR_RESET, glength(ker_T_basis));
        
//         if (glength(ker_T_basis)==0)
//         {
//             gel(I_vect, i) = F;
//             done = 1;
//         }
//         else {
//             ker_T = my_get_sums(ker_T_basis, itos(p));

//             // WARNING: This makes it faster but might not find the answer. Switch back to the above ker_T
//             //ker_T = shallowconcat(mkvec(zerocol(glength(gel(ker_T_basis, 1)))), ker_T_basis);
            
//             //DEBUG_PRINT(1, ANSI_COLOR_CYAN "ker_T: %Ps\n" ANSI_COLOR_RESET, ker_T);
//             // DEBUG_PRINT(1, "ker_T[2]: %Ps\n", gel(ker_T, 2));
//             f = glength(ker_T);
//             DEBUG_PRINT(1, "Searching a chunk of ker (1-sigma) of size: %d\n", f);
            
//             /*
//                 Find F_ker_T --------------------
//             */
            
//             for (j = 1; j <= f; j++)
//             {
//                 DEBUG_PRINT(1, "%d/%d\n", j, f);
//                 // DEBUG_PRINT(1, "Adding the exp: %Ps\n", gel(ker_T, j));
//                 // DEBUG_PRINT(1, "And the ideal: %Ps\n", my_ideal_from_exp(Labs, gel(ker_T, j)));
//                 F_ker_T = ZC_add(F, gel(ker_T, j));
//                 //DEBUG_PRINT(1, "F_ker_T: %Ps\n", F_ker_T);

//                 // Here we should check, using linear algebra--------- 
//                 is_princ = bnfisprincipal0(Labs, idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, F_ker_T)), 1);
//                 //-----------------------------------
//                 if (!my_QV_equal0(gel(is_princ, 1)))
//                 {
//                     DEBUG_PRINT(1, ANSI_COLOR_RED "Problem in my_H90_exp_vect\n" ANSI_COLOR_RESET);
//                     pari_close();
//                     exit(111);
//                 }
                
//                 t = gel(is_princ, 2);
                
//                 t_rel = rnfeltabstorel(Lrel, t);
//                 Nt = rnfeltnorm(Lrel, t_rel);
//                 diff = nfmul(K, Nt, a);
//                 exp = bnfisunit0(K, diff, NULL);
//                 // DEBUG_PRINT(1, ANSI_COLOR_YELLOW "exp: %Ps\n" ANSI_COLOR_RESET, exp);
//                 // DEBUG_PRINT(1, "M: %Ps\n", my_norm_operator(Labs, Lrel, K, p));
//                 // DEBUG_PRINT(1, "D: %Ps\n", my_xcol(glength(exp), itos(p)));
//                 // DEBUG_PRINT(1, "B: %Ps\n", gtocol(exp));
//                 /* check if exp lies in the image of the operator associated to N: (O_L)^x -> (O_K)^x */
//                 // if (my_QV_equal0(my_norm_operator(Labs, Lrel, K, p)))
//                 // {
//                 //     gel(I_vect, i) = F_ker_T;
//                 //     done = 1;
//                 //     break;
//                 // }
//                 is_norm = matsolvemod(my_norm_operator(Labs, Lrel, K, p), zerocol(glength(exp)), gtocol(exp), 0);
//                 //DEBUG_PRINT(1, ANSI_COLOR_CYAN "is_norm: %Ps\n" ANSI_COLOR_RESET, is_norm);
//                 if (my_QV_equal0(exp) || !gequal0(is_norm))
//                 {
//                     //DEBUG_PRINT(1, ANSI_COLOR_CYAN "Norm of elt with exp: %Ps\n" ANSI_COLOR_RESET, is_norm);
//                     // w = gel(bnf_get_fu(Labs), 1);
//                     // w_rel = rnfeltabstorel(Lrel, w);
//                     // Nw = rnfeltnorm(Lrel, w_rel);
//                     // DEBUG_PRINT(1, "Norm w: %Ps\n", Nw);
//                     gel(I_vect, i) = F_ker_T;
//                     done = 1;
//                     break;
//                 }
//             }
//         }
//         if (!done)
//         {
//             DEBUG_PRINT(1, ANSI_COLOR_RED "my_H90_exp_vect ended with a problem \n" ANSI_COLOR_RESET);
//             pari_close();
//             exit(111);
//         }
        
         
        
//     } 
//     DEBUG_PRINT(1, "Finished search\n");
//     I_vect = gerepilecopy(av, I_vect);
//     return I_vect;
// }

// void my_compute_unit_norms(GEN Lrel, GEN K, GEN units_group) {
    
//     GEN rel_unit, norm;
//     int l = glength(units_group);
    
//     int i;
//     for (i = 1; i <= l; i++)
//     {
//         rel_unit = rnfeltabstorel(Lrel, gel(units_group, i));
//         norm = rnfeltnorm(Lrel, rel_unit);

//         DEBUG_PRINT(1, "norm[%d]: %Ps\n", i, bnfisunit0(K, norm, NULL));

//     }
// } 

// void my_norms(GEN K, GEN K_ext, GEN p) {
//     GEN Labs, Lrel, unit_group;
//     int l = glength(K_ext);
//     int i;
//     for (i = 1; i < l+1; i++)
//     {
//         Labs = gel(gel(K_ext, i), 1);
//         unit_group = my_find_units_mod_p(Labs, p);
//         Lrel = gel(gel(K_ext, i), 2);
//         my_compute_unit_norms(Lrel, K, unit_group);
//     }
    
// }

// void my_compute_unit_diffs(GEN Labs, GEN sigma, GEN units_group) {
    
//     GEN diff;
//     int l = glength(units_group);
    
//     int i;
//     for (i = 1; i <= l; i++)
//     {
//         diff = my_1MS_elt(Labs, sigma, gel(units_group, i));

//         DEBUG_PRINT(1, "diff[%d]: %Ps\n", i, bnfisunit0(Labs, diff, NULL));

//     }
// } 

// void my_diffs(GEN K_ext, GEN p) {
//     GEN Labs, unit_group, sigma;
//     int l = glength(K_ext);
//     int i;
//     for (i = 1; i <= l; i++)
//     {
//         Labs = gel(gel(K_ext, i), 1);
//         sigma = gel(gel(K_ext, i), 3);
//         unit_group = my_find_units_mod_p(Labs, p);
       
//         my_compute_unit_diffs(Labs, sigma, unit_group);
//     }
    
// }

// void my_compute_matrices(GEN Labs, GEN sigma, GEN units_group) {
    
//     GEN image;
//     int l = glength(units_group);
    
//     int i;
//     DEBUG_PRINT(1, "1:\n\n");
//     for (i = 1; i <= l; i++)
//     {
//         image = bnfisunit0(Labs, gel(units_group, i), NULL);

//         DEBUG_PRINT(1, "%Ps\n", image);

//     }
//     DEBUG_PRINT(1, "\nsigma:\n\n");
//     for (i = 1; i <= l; i++)
//     {
//         image = bnfisunit0(Labs, galoisapply(Labs, sigma, gel(units_group, i)), NULL);

//         DEBUG_PRINT(1, "%Ps\n", image);

//     }
//     DEBUG_PRINT(1, "\nsigma^2:\n\n");
//     for (i = 1; i <= l; i++)
//     {
//         image = bnfisunit0(Labs, galoisapply(Labs, sigma, galoisapply(Labs, sigma, gel(units_group, i))), NULL);

//         DEBUG_PRINT(1, "%Ps\n", image);

//     }
// } 

// void my_matrices(GEN K_ext, GEN p) {
//     GEN Labs, unit_group, sigma;
//     int l = glength(K_ext);
//     int i;
//     for (i = 1; i <= l; i++)
//     {
//         Labs = gel(gel(K_ext, i), 1);
//         sigma = gel(gel(K_ext, i), 3);
//         unit_group = my_find_units_mod_p(Labs, p);
//         DEBUG_PRINT(1, "Field ext [%d]:\n\n", i);
//         my_compute_matrices(Labs, sigma, unit_group);
//     }
    
// }

GEN my_find_I_full (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN i_xJ, GEN a, GEN class_group, int p_int)
{
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_find_I_full\n--------------------------\n\n");
    pari_sp av = avma;
    GEN diff, current_I, test_vec, class_number = bnf_get_no(Labs), Itest, t, t_rel, norm_t, exponents, norms, L_units, L_unit_group, ret, is_princ;
    int clnr = itos(class_number), n, roots_of_unity_nr = bnf_get_tuN(K), l;
    
    
    // for (n = 1; n < clnr + 1; ++n) {
    //     DEBUG_PRINT(1, "%d/%d\n", n, clnr);
    //     current_I = gel(class_group, n);
    //     Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
    //     test_vec = bnfisprincipal0(Labs, Itest, 1);

    //     if (my_QV_equal0(gel(test_vec, 1)))
    //     {
    //         t = gel(test_vec, 2);
    //         t_rel = rnfeltabstorel(Lrel, t);
    //         norm_t = rnfeltnorm(Lrel, t_rel);
    //         exponents = bnfisunit(K, nfmul(K, norm_t, a));
    //         DEBUG_PRINT(1, "exponents: %Ps\n", exponents);
    //         if (my_is_p_power(exponents, p_int, roots_of_unity_nr))
    //         {
    //             DEBUG_PRINT(1, "I found\n\n");
    //             ret = mkvec2(gel(test_vec, 2), current_I);
    //             ret = gerepilecopy(av, ret);
    //             return ret;
    //         }
            
    //     } 

    // }

    //DEBUG_PRINT(1, "Entering the slow round\n\n");
    DEBUG_PRINT(1, "roots of unity L: %d\n\n", roots_of_unity_nr);
    L_units = my_find_units_mod_p(Labs, stoi(p_int));
    DEBUG_PRINT(1, "nr of L_units mod p: %ld\n\n", glength(L_units));
    L_unit_group = my_get_unit_group(Labs, L_units, stoi(p_int));
    DEBUG_PRINT(1, "size of L_unit_group: %ld\n\n", glength(L_unit_group));
    
    l = glength(L_unit_group);
    norms = zerovec(l);
    
    for (n = 1; n <= l; n++)
    {
        gel(norms, n) = algtobasis(K, rnfeltnorm(Lrel, rnfeltabstorel(Lrel, gel(L_unit_group, n))));
    }
    // pari_close();
    // exit(0);
    int i;
    for (n = 1; n <= clnr; ++n) {
        DEBUG_PRINT(1, "%d/%d\n", n, clnr);
        current_I = gel(class_group, n);
        Itest = idealdiv(Labs, i_xJ, my_1MS_ideal(Labs, sigma, current_I));
        test_vec = bnfisprincipal0(Labs, Itest, nf_GEN);

        if (ZV_equal0(gel(test_vec, 1)))
        {
            // Compute again and this time force it to find the generator (may cause an overflow)
            is_princ = bnfisprincipal0(Labs, Itest, nf_GENMAT);
            t = nffactorback(Labs, gel(gel(is_princ, 2), 1), gel(gel(is_princ, 2), 2));
            //t = gel(test_vec, 2);
            DEBUG_PRINT(1, "t: %Ps\n", t);
            DEBUG_PRINT(1, "length t: %ld\n", glength(t));
            if (glength(t)!=0)
            {

                t_rel = rnfeltabstorel(Lrel, t);
                norm_t = rnfeltnorm(Lrel, t_rel);
                diff = nfmul(K, norm_t, a);
                for ( i = 1; i <= l; i++)
                {
                    exponents = bnfisunit0(K, nfdiv(K, diff, gel(norms, i)), NULL);
                    if (l%1000==0)
                    {
                        DEBUG_PRINT(1, "exponents[%d]: %Ps\n", i, exponents);
                    }
                    
                    if (my_is_p_power(exponents, p_int, roots_of_unity_nr))
                    {
                        DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nI found\n\n" ANSI_COLOR_RESET);
                        ret = mkvec2(gel(test_vec, 2), current_I);
                        ret = gerepilecopy(av, ret);
                        DEBUG_PRINT(1, "\n--------------------------\nEnd: my_find_I_full\n--------------------------\n\n");
                        return ret;
                    }
                    
                }
            }
                
            
        } 

    }
    DEBUG_PRINT(1, ANSI_COLOR_RED "my_find_I_full ended with a problem\n" ANSI_COLOR_RESET);
   
    DEBUG_PRINT(1, ANSI_COLOR_RED "No I found in my_find_I\n\n" ANSI_COLOR_RESET);
    pari_close();
    exit(111);

    
    return 0;
}

GEN my_find_I_vect_full (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN Ja_vect, int p_int) {
    pari_sp av = avma;
    int r_rk = glength(Ja_vect);
    GEN I_vect = zerovec(r_rk), a, i_xJ;
   
    GEN class_group = my_get_clgp (Labs);
   
    int i;
    
    for (i = 1; i <= r_rk; ++i)
    {
        a = gel(gel(Ja_vect, i), 1);
        i_xJ = rnfidealup0(Lrel, gel(gel(Ja_vect, i),2), 1);
        gel(I_vect, i) = gel(my_find_I_full(Labs, Lrel, K, sigma, i_xJ, a, class_group, p_int), 2);
    } 
    
    I_vect = gerepilecopy(av, I_vect);
    return I_vect;
}
