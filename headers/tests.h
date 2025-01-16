int my_is_p_torsion (GEN K, GEN J_vect, GEN p) {
    GEN is_pr, test_vec, power;
    int l = glength(J_vect);
    int i;
    int is_principal = 1;

    for (i=1; i<l+1; ++i) {
        power = idealpow(K, gel(J_vect, i), p);
        is_pr = bnfisprincipal0(K, power, 1);
        test_vec = gel(is_pr, 1);
        if (!ZV_equal0(test_vec))
        {
            is_principal = 0;
        }
        
    }

    return is_principal;
}

int my_test_H90_ideal (GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN I_vect, GEN J_vect) {
    int p_rk = glength(J_vect);
    GEN iJ, test_ideal, test_vec, is_pr;
    int is_principal = 1;
    //pari_printf("%Ps\n", gel(ext_gens, 1));
    
    int i;
    
    for (i = 1; i < glength(I_vect)+1; ++i)
    {
        if (i<=p_rk) {
            iJ = rnfidealup0(Lrel, gel(J_vect, i), 1);
            test_ideal = idealdiv(Labs, iJ, my_1MS_ideal(Labs, sigma, gel(I_vect, i)));
        }
        else {
            test_ideal = idealinv(Labs, my_1MS_ideal(Labs, sigma, gel(I_vect, i)));
        }
        is_pr = bnfisprincipal0(Labs, test_ideal, 1);
        test_vec = gel(is_pr, 1);
        if (!ZV_equal0(test_vec))
        {
            is_principal = 0;
        }

    } 
    return is_principal;
}
