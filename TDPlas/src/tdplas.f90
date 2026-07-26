module tdplas

  use tdplas_constants, only: &
       dbl, cmp, flg, i4b, zero

  use global_tdplas, only: &
       global_sys_Ftest,   &
       global_prop_Fprop,  &
       global_prop_Fint,   &
       global_medium_Fmdm

  use BEM_medium, only: &
       do_BEM_prop,            &
       BEM_Q0,                 &
       BEM_Qd,                 &
       BEM_Q0x,                &
       BEM_Qdx,                &
       BEM_S,                  &
       BEM_Sm1_dum,            &
       BEM_Z,                  &
       BEM_Sdum_act,           &
       deallocate_BEM_public

  use td_ContMed, only: &
       init_potential,       &
       init_charges,         &
       prop_chr,             &
       calc_charges,         &
       deallocate_potential, &
       init_vv_propagator

  use pedra_friends, only: &
       pedra_surf_Fdum,      &
       pedra_surf_n_tessere, &
       pedra_surf_tessere,   &
       pedra_dum_n_tessere,  &
       pedra_dum_tessere

  use global_quantum, only: &
       quantum_init, &
       readio_and_init_tdplas_for_octopus

  implicit none
  public

end module tdplas
