program td_pdos
 
 use mod_td 

 implicit none

!DPDOS(t,e) = PDOS(t,e) - PDOS^GS(e)
! = -2 sum_i^occ w_at^i [ sum_ll' Re[C_l'^*(t)C_l(t)] sum_c^vir A_ci^l' A_ci^l ] delta(e-e_i)
!   +2 sum_i^vir w_at^i [ sum_ll' Re[C_l'^*(t)C_l(t)] sum_v^occ A_iv^l' A_iv^l ] delta(e-e_i)

!indexes i,v and c run over the set of MOs
!indexes l and l' run over the set of excited states
!|l> = sum_cv A_cv^l a_c^dagger a_v |GS>
!|Psi > = sum_l C_l(t) |l>


 write(*,*) ''
 write(*,*) '****************************************************'
 write(*,*) '****************************************************'
 write(*,*) '**                                                **' 
 write(*,*) '**                    TD-PostPro                  **'
 write(*,*) '**                                                **'
 write(*,*) '**    Post-Processing for WaveT, ADF and MolGW    **'
 write(*,*) '**                                                **'
 write(*,*) '**                     by                         **'
 write(*,*) '**               Emanuele Coccia                  **'
 write(*,*) '**             Pablo Grobas Illobre               **'
 write(*,*) '**              Margherita Marsili                **'
 write(*,*) '**               Daniele Toffoli                  **'
 write(*,*) '**                                                **'
 write(*,*) '****************************************************'
 write(*,*) '****************************************************'    
 write(*,*) ''

!Read input
 call read_input()
 write(*,*) ''
 write(*,*) 'Reading input done'


!Read w, A and e using dip_stec
 call read_adf_molgw()
 write(*,*) 'Reading ADF or MolGW done'

!Read C(t) from c_t_1.dat
 call read_wavet()
 write(*,*) 'Reading WaveT done'

!Compute the time-dependent DeltaPDOS
 call compute_td_pdos()  
 write(*,*) 'Computing td-PDOS done'
 write(*,*) 'Computing td-OPDOS done'

if (tcm1.eq.'y') then
!Compute the time-dependent TCM 
call compute_td_tcm()
 write(*,*) 'Computing td-TCM done'
endif
!Deallocate the arrays
 call deallocate_all()

 write(*,*) ''
 write(*,*) '****************************************************'
 write(*,*) '****************************************************'
 write(*,*) '**                                                **'
 write(*,*) '**             END OF SIMULATION                  **'
 write(*,*) '**                                                **'
 write(*,*) '****************************************************'
 write(*,*) '****************************************************' 


 stop

end program td_pdos

