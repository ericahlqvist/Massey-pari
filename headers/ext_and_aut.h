/*
Generates all extensions of base together with all necessary automorphisms 
*/

GEN my_ext(GEN base, GEN base_clf, GEN p, int p_rk, GEN D_prime_vect) 
{   
    pari_sp av = avma;
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_ext\n--------------------------\n\n");
    GEN x, y, p1, q1, p1red, Lrel, Labs, Lbnr, s_lift_x, cx, sigma, s=pol_x(fetch_user_var("s"));

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));

    GEN base_ext = cgetg(p_rk+1,t_VEC);

    // DEBUG_PRINT(1, "base l: %ld\n", glength(base_ext));
    // DEBUG_PRINT(1, "Base_clf: %Ps\n\n", base_clf);
    //char filename[100];
    
    int i, j;
    for (i=1; i<p_rk+1; ++i) {
        p1 = gel(base_clf, i);
        q1 = gsubstpol(p1, x, y);
        
        /* Define Lrel/Labs */
        p1red = rnfpolredbest(base, mkvec2(q1, D_prime_vect), 0);
        // p1red = q1;
        DEBUG_PRINT(1, "Reduced polynomial for relative extension found\n");
        Lrel = rnfinit(base, p1red);

        DEBUG_PRINT(1, "Lrel found\n");
        DEBUG_PRINT(1, "Abs pol: %Ps\n", rnf_get_polabs(Lrel));
        
        Labs = Buchall(rnf_get_polabs(Lrel), nf_FORCE, DEFAULTPREC);
        Lbnr = bnrinit0(Labs,gen_1, 1);
        
        // May change flag to nf_FORCE
        // Other precisions: MEDDEFAULTPREC, BIGDEFAULTPREC
        //Labs = Buchall_param(rnf_get_polabs(Lrel), 1.5,1.5,4, nf_FORCE, DEFAULTPREC);
        DEBUG_PRINT(1, "Labs found\n");
        DEBUG_PRINT(1, "L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));
        DEBUG_PRINT(1, "rel_pol[%d]: %Ps\n", i, p1red);
        DEBUG_PRINT(1, "\nabs pol: %Ps\n\n", gsubstpol(rnf_get_polabs(Lrel),y, s));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        //DEBUG_PRINT(1, "\nLift: %Ps\n\n", s_lift_x);
        cx = galoisconj(Labs, NULL);
        //DEBUG_PRINT(1, "\nGalois conj: %Ps\n\n", cx);
        if (glength(cx)==nf_get_degree(bnf_get_nf(Labs)))
        {
            DEBUG_PRINT(1, ANSI_COLOR_GREEN "\n------------------------\nLabs is Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        else {
            DEBUG_PRINT(1, ANSI_COLOR_RED "\n------------------------\nLabs is not Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        

        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!ZV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && ZV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x)) 
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //DEBUG_PRINT(1, "sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //DEBUG_PRINT(1, ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
    
        gel(base_ext, i) = mkvec4(Labs, Lrel, sigma, Lbnr);

        // sDEBUG_PRINT(1, filename, "/Users/eric/Documents/Matematik/cup-products/large-fields/test/ext_%d", i);
        // writebin(filename, gel(base_ext, i));
        
        
        // my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    DEBUG_PRINT(1, ANSI_COLOR_GREEN "Extensions and generators found\n----------------------------\n\n" ANSI_COLOR_RESET);
    base_ext = gerepilecopy(av, base_ext);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_ext\n--------------------------\n\n");
    return base_ext;
}


