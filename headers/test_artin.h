void my_test_artin_symbol (GEN Labs, GEN Lrel, GEN K, int p) {
    printf("\n--------------------------\nStart: my_test_artin_symbol\n--------------------------\n\n");
    GEN I = gel(bnf_get_gen(K), 1);
    GEN elem = zerocol(nf_get_degree(bnf_get_nf(K)));
    gel(elem, 2) = gen_1;
    GEN I_mod = idealmul(K, idealhnf0(K, elem, NULL), I);
    GEN elem2 = zerocol(nf_get_degree(bnf_get_nf(K)));
    gel(elem2, 2) = gen_1;
    gel(elem2, 1) = gen_1;
    
    if (my_Artin_symbol(Labs, Lrel, K, I, p)==my_Artin_symbol(Labs, Lrel, K, I_mod, p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 1: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 1: FAILED\n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }
    
    if (my_Artin_symbol(Labs, Lrel, K, idealpow(K, I, gen_2), p)==my_Artin_symbol(Labs, Lrel, K, idealred(K, idealpow(K, I, gen_2)), p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 2: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 2: FAILED\n\n" ANSI_COLOR_RESET);
        printf("%d", my_Artin_symbol(Labs, Lrel, K, idealpow(K, I, gen_2), p));
        printf("%d", my_Artin_symbol(Labs, Lrel, K, idealred(K, idealpow(K, I, gen_2)), p));
        pari_close();
        exit(0);
    }
    
    if (my_Artin_symbol(Labs, Lrel, K, idealpow(K, I, gen_2), p)==my_Artin_symbol(Labs, Lrel, K, idealmul(K, idealpow(K, I, gen_2), idealhnf0(K, elem2, NULL)), p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 3: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 3: FAILED\n\n" ANSI_COLOR_RESET);
        pari_close();
        exit(0);
    }
    printf("\n--------------------------\nEnd: my_test_artin_symbol\n--------------------------\n\n");

}

void my_test_artin_symbol_spec (GEN Labs, GEN Lrel, GEN K, int p, GEN I) {
    pari_printf("Artin symbol: %Ps\n\n", my_Artin_symbol(Labs, Lrel, K, I, p));

}

void my_test_artin_on_norms (GEN Labs, GEN Lrel, GEN K, int p, GEN sigma) {
    GEN L_clgp = my_get_clgp(Labs);
    GEN Irel, N;
    long artin_symbol;
    int i;
    for (i = 1; i < glength(L_clgp)+1; i++)
    {
        Irel = rnfidealabstorel(Lrel, gel(L_clgp,i));
        N = rnfidealnormrel(Lrel, Irel);
        artin_symbol = my_Artin_symbol(Labs, Lrel, K, N, p);
        //output(artin_symbol);
        if (artin_symbol)
        {
            printf(ANSI_COLOR_RED "Artin symbol not zero on norms\n\n" ANSI_COLOR_RESET);
            pari_close();
            exit(0);
        }
        
    }
    printf(ANSI_COLOR_GREEN "Artin symbol zero on norms\n\n" ANSI_COLOR_RESET);
        
}

// void my_test_artin_on_norms_2 (GEN LxAbs, GEN LxRel, GEN LyAbs, GEN LyRel, GEN K, int p, GEN sigma_x, GEN sigma_y) {
//     GEN p_gens = my_find_p_gens(K, stoi(p));
//     GEN Irel, N, artin_symbol, i_xJ, I;
//     GEN class_group = my_get_clgp (LxAbs);
//     int i;
//     for (i = 1; i < glength(p_gens)+1; i++)
//     {
//         i_xJ = rnfidealup0(LxRel, gel(p_gens, i), 1);
//         I = gel(my_find_I(LxAbs, K, sigma_x, i_xJ, class_group),2);
//         Irel = rnfidealabstorel(LxRel, I);
//         N = rnfidealnormrel(LxRel, Irel);
//         artin_symbol = my_Artin_symbol(LyAbs, LyRel, K, N, p, sigma_y);
//         //output(artin_symbol);
//         if (itos(artin_symbol))
//         {
//             printf(ANSI_COLOR_RED "Artin symbol not zero on norms\n\n" ANSI_COLOR_RESET);
//             pari_close();
//             exit(0);
//         }
        
//     }
//     printf(ANSI_COLOR_GREEN "Artin symbol zero on norms\n\n" ANSI_COLOR_RESET);
        
// }

int my_test_artin_symbol_int (GEN Labs, GEN Lrel, GEN K, int p) {
    
    GEN I = gel(bnf_get_gen(K), 1);
    GEN elem = zerocol(nf_get_degree(bnf_get_nf(K)));
    gel(elem, 2) = gen_1;
    GEN I_mod = idealmul(K, idealhnf0(K, elem, NULL), I);
    GEN elem2 = zerocol(nf_get_degree(bnf_get_nf(K)));
    gel(elem2, 2) = gen_1;
    gel(elem2, 1) = gen_1;
    
    if (my_Artin_symbol(Labs, Lrel, K, I, p)==my_Artin_symbol(Labs, Lrel, K, I_mod, p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 1: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 1: FAILED\n\n" ANSI_COLOR_RESET);
        return 0;
    }
    
    if (my_Artin_symbol(Labs, Lrel, K, idealpow(K, I, gen_2), p)==my_Artin_symbol(Labs, Lrel, K, idealred(K, idealpow(K, I, gen_2)), p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 2: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 2: FAILED\n\n" ANSI_COLOR_RESET);
        return 0;
    }
    
    if (my_Artin_symbol(Labs, Lrel, K, idealpow(K, I, gen_2), p)==my_Artin_symbol(Labs, Lrel, K, idealmul(K, idealpow(K, I, gen_2), idealhnf0(K, elem2, NULL)), p))
    {
        printf(ANSI_COLOR_GREEN "Artin symbol test 3: PASSED\n\n" ANSI_COLOR_RESET);
    }
    else {
        printf(ANSI_COLOR_RED "Artin symbol test 3: FAILED\n\n" ANSI_COLOR_RESET);
        return 0;
    }
    return 1;

}