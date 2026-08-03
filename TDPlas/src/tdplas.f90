module tdplas

  use tdplas_constants, only: &
       dbl, cmp, flg, i4b, zero

  use global_tdplas, only: &
       global_sys_Ftest,   &
       global_prop_Fprop,  &
       global_prop_Fint,   &
       global_medium_Fmdm, &
       global_prop_Fmdm_relax, &
       global_qmodes_Fmop,     &
       global_qmodes_nprint,   &
       global_sys_Fwrite,      &
       global_prop_max_cycles, &
       global_prop_threshold,  &
       global_prop_mix_coef,   &
       global_qmodes_nmodes,   &
       global_qmodes_qmmodes,  &
       global_prop_Finit_int,  &
       global_medium_Fbem,     &
       global_sys_Fdeb,        &
       global_prop_n_q

  use BEM_medium, only: &
       do_BEM,                 &
       do_BEM_prop,            &
       BEM_Q0,                 &
       BEM_Qd,                 &
       BEM_Q0x,                &
       BEM_Qdx,                &
       BEM_S,                  &
       BEM_Sm1_dum,            &
       BEM_Z,                  &
       BEM_Sdum_act,           &
       deallocate_BEM_public,  &
       mat_f0,                 &
       out_gcharges,           &
       BEM_W2,                 &
       BEM_Modes,              &
       do_BEM_quant

  use td_ContMed, only: &
       init_potential,       &
       init_charges,         &
       init_potential_prop,  &
       init_charges_prop,    &
       prop_chr,             &
       calc_charges,         &
       deallocate_potential, &
       init_vv_propagator,   &
       set_qorf,             &
       set_qorf_pot,         &
       set_charges,          &
       get_mdm_dip,          &
       get_gneq,             &
       init_mdm,             &
       prop_mdm,             &
       finalize_mdm,         &
       init_after_scf,       &
       mpibcast_readio_mdm,  &
       fr_0,                 &
       q0,                   &
       do_charges_from_pot,  &
       init_mdm_prop,        &
       get_qorf,             &
       get_qorf0,            &
       do_Rfield_from_dip,   &
       set_mu_tp,            &
       get_qr_fr               

  use pedra_friends, only: &
       pedra_surf_Fdum,      &
       pedra_surf_n_tessere, &
       pedra_surf_tessere,   &
       pedra_dum_n_tessere,  &
       pedra_dum_tessere,    & 
       pedra_surf_spheres,   &
       pedra_surf_n_spheres

  use global_quantum, only: &
       quantum_init, &
#if QMCODE == wt
       readio_and_init_tdplas_for_wt
#elif QMCODE == octopus
       readio_and_init_tdplas_for_octopus
#elif QMCODE == ocpy   
       readio_and_init_tdplas_for_ocpy   
#endif

  use MathTools, only: &
       diag_mat,       &
       do_field_from_charges

  use drudel_epsilon, only: &
       drudel_eps_w0,       &
       drudel_eps_A

  implicit none
  public

end module tdplas
