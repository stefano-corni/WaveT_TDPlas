!> Module that reads the input of the medium.
!! It contains all public medium variables
      module tdplas_readio_dictionary
        use constants
        use dictionary_entry

        implicit none

        save

        !SYSTEM
        type(dict_entry) :: dict_Ftest
        type(dict_entry) :: dict_Fdeb
        type(dict_entry) :: dict_Fwrite

        !SURFACE
        type(dict_entry)  :: dict_Fsurf
        type(dict_entry)  :: dict_Fcav
        type(dict_entry)  :: dict_Fshape
        type(dict_entry)  :: dict_Finv

        !MEDIUM
        !namelist
        type(dict_entry)  :: dict_Fmdm
        type(dict_entry)  :: dict_Finit
        type(dict_entry)  :: dict_Fpol
        type(dict_entry)  :: dict_Fbem
        type(dict_entry)  :: dict_FinitBEM

        !EPS namelist
        type(dict_entry)  :: dict_Feps

       !PROPAGATE
        !namelist
        type(dict_entry)  :: dict_Fsoft
        type(dict_entry)  :: dict_Fprop
        type(dict_entry)  :: dict_Fint
        type(dict_entry)  :: dict_Finit_int
        type(dict_entry)  :: dict_Floc
        type(dict_entry)  :: dict_Fmdm_relax
        type(dict_entry)  :: dict_Fmdm_res



        !OUT_MATRIX namelist
        type(dict_entry)  :: dict_Fgamess
        type(dict_entry)  :: dict_Fprint_lf_matrix


        !print_charges
        type(dict_entry)  :: dict_Fmop


        public   dict_Ftest, dict_Fdeb, dict_Fwrite,                                                            &
                 dict_Fsurf, dict_Fcav, dict_Fshape, dict_Finv,                                                 &
                 dict_Fmdm,  dict_Finit, dict_Fpol, dict_Fbem, dict_FinitBEM,                                   &
                 dict_Feps,                                                                                     &
                 dict_Fsoft, dict_Fprop, dict_Fint, dict_Finit_int, dict_Floc, dict_Fmdm_relax, dict_Fmdm_res,  &
                 dict_Fgamess, dict_Fprint_lf_matrix,                                                           &
                 dict_Fmop,                                                                                     &
                 dict_tdplas_init


        contains

    subroutine dict_tdplas_init()

        !SYSTEM
      call  dict_entry_init(dict_Ftest, 6,                                              &
                 [character(flg) :: "test_type", "global_sys_Ftest"],                   &
                 [character(flg) :: "n-r", "n-l", "s-r", "s-l", "qmt","non"],           &
                 [character(flg) :: "n-r", "n-l", "s-r", "s-l", "qmt","non"])

      call dict_entry_init(dict_Fdeb, 4,                                                &
               [Character(flg) :: "debug_type", "global_sys_Fdeb"],                     &
               [Character(flg) :: "equ", "vmu", "off", "non"],                          &
               [Character(flg) :: "equ", "vmu", "off", "non"])

      call  dict_entry_init(dict_Fwrite, 2,                                             &
                 [Character(flg) :: "out_level", "global_sys_Fwrite"],                      &
                 [Character(flg) :: "low", "high"],                                     &
                 [Character(flg) :: "low", "high"])


        !MEDIUM
      call  dict_entry_init(dict_Fmdm, 5,                                               &
                 [Character(flg) :: "medium_type", "global_medium_Fmdm"],               &
                 [Character(flg) :: "nanop", "sol", "quantum_nanop", "quantum_sol", "non"],    &
                 [Character(flg) :: "cnan", "csol", "qnan", "qsol", "non"])

      call  dict_entry_init(dict_Finit, 4,                                              &
                [Character(flg) :: "medium_init0", "global_medium_Finit"],              &
                [Character(flg) :: "read_file", "vacuum", "frozen", "non"],             &
                [Character(flg) :: "rea", "vac", "fro", "non"])

      call  dict_entry_init(dict_Fpol, 3,                                               &
                 [Character(flg) :: "medium_pol", "global_medium_Fpol"],                &
                 [Character(flg) :: "charge", "dipole", "non"],                         &
                 [Character(flg) :: "chr", "dip", "non"])

      call  dict_entry_init(dict_Fbem, 3,                                               &
                 [Character(flg) :: "bem_type", "global_medium_Fbem"],                  &
                 [Character(flg) :: "diagonal", "standard", "non"],                     &
                 [Character(flg) :: "diag", "stan", "non"])

      call  dict_entry_init(dict_FinitBEM, 2,                                           &
                 [Character(flg) :: "bem_read_write", "global_medium_FinitBEM"],        &
                 [Character(flg) :: "read", "write", "non"],                            &
                 [Character(flg) :: "rea", "wri", "non"])



       !SURFACE
      call  dict_entry_init(dict_Fsurf, 3,                                              &
                 [Character(flg) :: "input_surface", "global_surf_Fsurf"],              &
                 [Character(flg) :: "spheres", "cavity", "non"],                        &
                 [Character(flg) :: "sphe", "cav", "non"])


      call  dict_entry_init(dict_Fcav, 4,                                               &
                 [Character(flg) :: "create_cavity", "pedra_surf_Fcav"],                &
                 [Character(flg) :: "built", "read_file", "gms", "non"],                &
                 [Character(flg) :: "bui", "fil", "gms", "non"])

      call  dict_entry_init(dict_Finv, 2,                                               &
                 [Character(flg) :: "inversion", "pedra_surf_Finv"],                    &
                 [Character(flg) :: "inversion", "non"],                                &
                 [Character(flg) :: "inv", "non"])


      call  dict_entry_init(dict_Fshape, 3,                                             &
                 [Character(flg) :: "spheres_shape", "sph_surf_Fshape"],                &
                 [Character(flg) :: "spheres", "spheroids", "non"],                     &
                 [Character(flg) :: "sphe", "spho", "non"])



       !EPS
      call  dict_entry_init(dict_Feps, 4,                                                    &
                 [Character(flg) :: "epsilon_omega", "global_eps_Feps"],                     &
                 [Character(flg) :: "drude-lorentz", "debye", "read_file", "gold_model"],    &
                 [Character(flg) :: "drl", "deb", "read", "gold"])



        !PROPAGATION
      call  dict_entry_init(dict_Fsoft, 4,                                                  &
                 [Character(flg) :: "propagation_software", "global_propagation_Fsoft"],    &
                 [Character(flg) :: "wavet", "octopus", "ocpy", "non"],                     &
                 [Character(flg) :: "wavet", "oct", "ocpy", "non"])


      call  dict_entry_init(dict_Fprop, 5,                                                                  &
                 [Character(flg) :: "propagation_type", "global_prop_Fprop"],                               &
                 [Character(flg) :: "dipole", "charge-ief", "charge-ief-one-tau", "charge-onsager", "non"], &
                 [Character(flg) :: "dip", "chr-ief", "chr-ied", "chr-ons", "non"])


      call  dict_entry_init(dict_Fint,3,                                                &
                 [Character(flg) :: "interaction_type", "global_prop_Fint"],            &
                 [Character(flg) :: "pcm", "onsager", "non"],                           &
                 [Character(flg) :: "pcm", "ons", "non"])


      call  dict_entry_init(dict_Finit_int, 3,                                          &
                [Character(flg) :: "interaction_init", "global_prop_Finit_int"],        &
                [Character(flg) :: "nsc", "sce", "non"],                                &
                [Character(flg) :: "nsc", "sce", "non"])



      call  dict_entry_init(dict_Fmdm_relax, 2,                                         &
                [Character(flg) :: "medium_relax", "global_prop_Fmdm_relax"],           &
                [Character(flg) :: "relax", "non"],                                     &
                [Character(flg) :: "rel", "non"])


      call  dict_entry_init(dict_Fmdm_res, 2,                                           &
                [Character(flg) :: "medium_restart", "global_prop_Fmdm_res"],           &
                [Character(flg) :: "restart", "non"],                                   &
                [Character(flg) :: "yesr", "nonr"])


      call dict_entry_init(dict_Floc, 2,                                                &
                 [Character(flg) :: "local_field", "global_medium_Floc"],               &
                 [Character(flg) :: "local", "non"],                                    &
                 [Character(flg) :: "loc", "non"])



       !OUT_MATRIX
      call  dict_entry_init(dict_Fgamess, 2,                                            &
                 [Character(flg) :: "gamess", "global_out_Fgamess"],                    &
                 [Character(flg) :: "gamess", "non"],                                   &
                 [Character(flg) :: "yes", "no"])

      call dict_entry_init(dict_Fprint_lf_matrix, 2,                                     &
                [Character(flg) :: "print_lf_matrix", "global_out_Fprint_lf_matrix"],    &
                [Character(flg) :: "print", "non"],                                      &
                [Character(flg) :: "yes", "non"])


       !print_charges
      call  dict_entry_init(dict_Fmop, 2,                                              &
                 [Character(flg) :: "charge_mopac", "qmodes_Fmop"],                    &
                 [Character(flg) :: "print", "non"],                                   &
                 [Character(flg) :: "yes", "non"])

    end subroutine


 end module



