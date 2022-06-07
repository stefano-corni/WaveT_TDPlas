!> Module that reads the input of the medium.
!! It contains all public medium variables
      module user_input_and_flags_dictionary
        use tdplas_constants
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
        type(dict_entry)  :: dict_Ffind

        !MEDIUM
        !namelist
        type(dict_entry)  :: dict_Fmdm
        type(dict_entry)  :: dict_Finit
        type(dict_entry)  :: dict_Fpol
        type(dict_entry)  :: dict_Fbem
        type(dict_entry)  :: dict_read_write
        type(dict_entry)  :: dict_Fnorm
        type(dict_entry)  :: dict_bem_sym

        !EPS namelist
        type(dict_entry)  :: dict_Feps
        type(dict_entry)  :: dict_typ_prop

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


        !quantum_coupling
        type(dict_entry)  :: dict_Fmop
        type(dict_entry)  :: dict_Fqcp

        !EXT_PERT namelist
        type(dict_entry)  :: dict_Ftyp
        type(dict_entry)  :: dict_Feet
        type(dict_entry)  :: dict_pl_type
        type(dict_entry)  :: dict_print_field
        type(dict_entry)  :: dict_print_surface_charges


        public   dict_Ftest, dict_Fdeb, dict_Fwrite,                                                            &
                 dict_Fsurf, dict_Fcav, dict_Fshape, dict_Finv, dict_Ffind,                                     &
                 dict_Fmdm,  dict_Finit, dict_Fpol, dict_Fbem, dict_read_write, dict_Fnorm, dict_bem_sym,       &
                 dict_Feps, dict_typ_prop,                                                                      &
                 dict_Fsoft, dict_Fprop, dict_Fint, dict_Finit_int, dict_Floc, dict_Fmdm_relax, dict_Fmdm_res,  &
                 dict_Fgamess, dict_Fprint_lf_matrix,                                                           &
                 dict_Fmop, dict_Fqcp,                                                                          &
                 dict_tdplas_init, dict_Ftyp, dict_Feet


        contains

    subroutine dict_tdplas_init()

        !SYSTEM
      call  dict_entry_init(dict_Ftest, 6,                                              &
                 [character(flg) :: "test_type", "global_sys_Ftest"],                   &
                 [character(flg) :: "n-r", "n-l", "s-r", "s-l", "qmt", "non"],          &
                 [character(flg) :: "n-r", "n-l", "s-r", "s-l", "qmt", "non"])

      call dict_entry_init(dict_Fdeb, 4,                                                &
               [Character(flg) :: "debug_type", "global_sys_Fdeb"],                     &
               [Character(flg) :: "equ", "vmu", "off", "non"],                          &
               [Character(flg) :: "equ", "vmu", "off", "non"])

      call  dict_entry_init(dict_Fwrite, 2,                                             &
                 [Character(flg) :: "out_level", "global_sys_Fwrite"],                  &
                 [Character(flg) :: "low", "high"],                                     &
                 [Character(flg) :: "low", "high"])


        !MEDIUM



      call  dict_entry_init(dict_Fmdm, 5,                                                      &
                 [Character(flg) :: "medium_type", "global_medium_Fmdm"],                      &
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

      call  dict_entry_init(dict_read_write, 3,                                         &
                 [Character(flg) :: "bem_read_write", "global_medium_read_write"],      &
                 [Character(flg) :: "read", "write", "non"],                            &
                 [Character(flg) :: "rea", "wri", "non"])

      call dict_entry_init(dict_Floc, 2,                                                &
                 [Character(flg) :: "local_field", "global_medium_Floc"],               &
                 [Character(flg) :: "local", "non"],                                    &
                 [Character(flg) :: "loc", "non"])

      call dict_entry_init(dict_Fnorm, 3,                                                &
                 [Character(flg) :: "normalization", "global_medium_Fnorm"],             &
                 [Character(flg) :: "non", "total", "separate"],                         &
                 [Character(flg) :: "non", "tot", "sep"])

      call dict_entry_init(dict_bem_sym, 2,                                              &
                 [Character(flg) :: "bem_symmetric", "global_medium_bem_sym"],           &
                 [Character(flg) :: "yes", "non"],                                       &
                 [Character(flg) :: "yes", "non"])


       !SURFACE
      call  dict_entry_init(dict_Fsurf, 3,                                              &
                 [Character(flg) :: "surface_type", "global_surf_Fsurf"],              &
                 [Character(flg) :: "object", "mesh", "non"],                        &
                 [Character(flg) :: "sphe", "cav", "non"])


      call  dict_entry_init(dict_Fcav, 4,                                               &
                 [Character(flg) :: "input_mesh", "pedra_surf_Fcav"],                &
                 [Character(flg) :: "built", "read_file", "gms", "non"],                &
                 [Character(flg) :: "bui", "fil", "gms", "non"])

      call  dict_entry_init(dict_Finv, 2,                                               &
                 [Character(flg) :: "inversion", "pedra_surf_Finv"],                    &
                 [Character(flg) :: "inversion", "non"],                                &
                 [Character(flg) :: "inv", "non"])

   
         call  dict_entry_init(dict_Ffind, 2,                                               &
                 [Character(flg) :: "find_spheres", "pedra_surf_Ffind"],                    &
                 [Character(flg) :: "yes", "non"],                                          &
                 [Character(flg) :: "yes", "non"])



      call  dict_entry_init(dict_Fshape, 3,                                             &
                 [Character(flg) :: "object_shape", "sph_surf_Fshape"],                &
                 [Character(flg) :: "sphere", "spheroid", "non"],                     &
                 [Character(flg) :: "sphe", "spho", "non"])



       !EPS
      call  dict_entry_init(dict_Feps, 5,                                                    &
                 [Character(flg) :: "epsilon_omega", "global_eps_Feps"],                     &
                 [Character(flg) :: "drude-lorentz", "debye", "general", "gold_model", "non"],    &
                 [Character(flg) :: "drl", "deb", "gen", "gold", "non"])
     

      call  dict_entry_init(dict_typ_prop, 4,                                                    &
                 [Character(flg) :: "propagation_pole", "typ_prop"],                              &
                 [Character(flg) :: "velocity-verlet", "first-order", "inversion", "second-order"],  &
                 [Character(flg) :: "0", "1", "2", "3"])



        !PROPAGATION
      call  dict_entry_init(dict_Fsoft, 4,                                                  &
                 [Character(flg) :: "propagation_software", "global_prop_Fsoft"],    &
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


      call  dict_entry_init(dict_Finit_int, 4,                                          &
                [Character(flg) :: "interaction_init", "global_prop_Finit_int"],        &
                [Character(flg) :: "nsc", "sce", "non", "qmt"],                                &
                [Character(flg) :: "nsc", "sce", "non", "qmt"])



      call  dict_entry_init(dict_Fmdm_relax, 2,                                         &
                [Character(flg) :: "medium_relax", "global_prop_Fmdm_relax"],           &
                [Character(flg) :: "relax", "non"],                                     &
                [Character(flg) :: "rel", "non"])


      call  dict_entry_init(dict_Fmdm_res, 2,                                           &
                [Character(flg) :: "medium_restart", "global_prop_Fmdm_res"],           &
                [Character(flg) :: "restart", "non"],                                   &
                [Character(flg) :: "yesr", "nonr"])



       !OUT_MATRIX
      call  dict_entry_init(dict_Fgamess, 2,                                            &

                 [Character(flg) :: "gamess", "global_out_Fgamess"],                    &
                 [Character(flg) :: "gamess", "non"],                                   &
                 [Character(flg) :: "yes", "no"])

      call dict_entry_init(dict_Fprint_lf_matrix, 2,                                     &
                [Character(flg) :: "print_lf_matrix", "global_out_Fprint_lf_matrix"],    &
                [Character(flg) :: "print", "non"],                                      &
                [Character(flg) :: "yes", "non"])


       !quantum_coupling
      call  dict_entry_init(dict_Fmop, 2,                                              &
                 [Character(flg) :: "charge_mopac", "qmodes_Fmop"],                    &
                 [Character(flg) :: "print", "non"],                                   &
                 [Character(flg) :: "yes", "non"])

      call  dict_entry_init(dict_Fqcp, 3,                                              &
                 [Character(flg) :: "quantum_calculation", "qmodes_Fqcp"],             &
                 [Character(flg) :: "static", "dynamic", "non"],                       &
                 [Character(flg) :: "sta", "dyn", "non"])

        !ext_pert
       call  dict_entry_init(dict_Ftyp, 2,                                             &
                 [Character(flg) :: "pert_type", "global_ext_pert_Ftyp"],              &
                 [Character(flg) :: "field", "molecule"],                              &
                 [Character(flg) :: "field", "molecule"])

       call  dict_entry_init(dict_Feet, 2,                                             &
                 [Character(flg) :: "eet", "global_ext_pert_Feet"],                    &
                 [Character(flg) :: "yes", "non"],                                     &
                 [Character(flg) :: "yes", "non"])

       call  dict_entry_init(dict_pl_type,2,                                           &
                 [Character(flg) :: "pl_type", "global_ext_pert_pl_type"],             &
                 [Character(flg) :: "from_cienergy", "from_omega"],                    &
                 [Character(flg) :: "from_cienergy", "from_omega"])

       call  dict_entry_init(dict_print_field,2,                                       &
                 [Character(flg) :: "print_field", "global_ext_pert_print_field"],     &
                 [Character(flg) :: "yes", "non"],                                     &
                 [Character(flg) :: "yes", "non"])

       call  dict_entry_init(dict_print_surface_charges,2,                             &
                 [Character(flg) :: "print_surface_charges", "global_ext_pert_print_surface_charges"],     &
                 [Character(flg) :: "yes", "non"],                                     &
                 [Character(flg) :: "yes", "non"])

    end subroutine


 end module



