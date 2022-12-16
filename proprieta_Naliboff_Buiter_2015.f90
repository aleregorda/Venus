 subroutine proprieta(i,num_mat,materiale,Tref_M,Rho_M,Visc_M,alpha_M,Cp_M,Hr_M,Kc_M,ndisl,Ea,Va,coeffA,grain_m,grain_p,t_lmantle,&
 Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,t_mantle,t_cont,t_lcont,t_oc,t_air,t_sed,t_serp,t_mantleM,t_contM,t_lcontM,t_ocM,teta_fr_inf,&
 Cohe_inf,Ws,diffusion,rheology,f_scale,t_wet,Visc_min,Visc_max,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,uniaxial)

 Implicit none

 integer, intent(out) :: num_mat,i,t_mantle,t_cont,t_lcont,t_oc,t_air,t_sed,t_serp,t_mantleM,t_contM,t_lcontM,t_ocM,t_wet,t_lmantle
 character(len=*), dimension(:), allocatable, intent(out) :: materiale
 logical, dimension(:), allocatable, intent(out) :: diffusion,peierls,uniaxial
 integer, dimension(:), allocatable, intent(out) :: rheology
 double precision, dimension(:), allocatable, intent(out) :: Tref_M,Rho_M,Visc_M,alpha_M,Cp_M,Hr_M,Kc_M,ndisl,Ea,Va,&
 coeffA,Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,teta_fr_inf,Cohe_inf,Ws,f_scale,Visc_min,Visc_max,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,&
 grain_m,grain_p

 num_mat=9                                                    !Numero di materiali del sistema

 allocate(materiale(num_mat),Tref_M(num_mat),Rho_M(num_mat),Visc_M(num_mat),alpha_M(num_mat),Cp_M(num_mat),Hr_M(num_mat),&
 Kc_M(num_mat),ndisl(num_mat),Ea(num_mat),Va(num_mat),coeffA(num_mat),Eadf(num_mat),Vadf(num_mat),&
 coeffAdf(num_mat),teta_fr_0(num_mat),Cohe_0(num_mat),teta_fr_inf(num_mat),Cohe_inf(num_mat),Ws(num_mat),diffusion(num_mat),&
 rheology(num_mat),f_scale(num_mat),Visc_min(num_mat),Visc_max(num_mat),peierls(num_mat),Ea_pe(num_mat),ndisl_pe(num_mat),&
 Va_pe(num_mat),coeffA_pe(num_mat),uniaxial(num_mat),grain_m(num_mat),grain_p(num_mat))

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PROPRIETA' DEI MATERIALI - verificare che il numero di materiali sia uguale a num_mat
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!Materiale 1
 i=1
 t_cont=i
 materiale(i)='continental crust'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=2800.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.28D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=800.D0                                                  !Calore specifico
 Hr_M(i)=1.3D-06                                                  !Calore radiogenico
 Kc_M(i)=2.5D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=8.57D-28                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=4.D0                                                  !Stress exponent dislocation creep
 Ea(i)=2.23D05                                                  !Energia di attivazione dislocation creep
 Va(i)=0.D0                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=6.07931D-19                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=6.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.135723D-28                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=3.2D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=1.23D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=0.D0                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 2
 i=2
 t_lcont=i
 materiale(i)='lower continental crust'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=2950.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.28D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=800.D0                                                  !Calore specifico
 Hr_M(i)=1.3D-06                                                  !Calore radiogenico
 Kc_M(i)=2.5D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=7.13D-18                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.D0                                                  !Stress exponent dislocation creep
 Ea(i)=3.45D05                                                  !Energia di attivazione dislocation creep
 Va(i)=0.D0                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=6.07931D-19                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=6.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.135723D-28                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=3.2D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=1.23D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=0.D0                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 3
 i=3
 materiale(i)='seed'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=3300.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=6.52D-16                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.5D0                                                  !Stress exponent dislocation creep
 Ea(i)=5.3D05                                                  !Energia di attivazione dislocation creep
 Va(i)=18.D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=2.37D-15                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.75D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=10.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.D-150                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=20.D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=5.4D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=10.D-06                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=0.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=1.D05                                                  !Coesione iniziale
 teta_fr_inf(i)=0.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=1.D05                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 4
 i=4
 t_mantle=i
 materiale(i)='upper mantle'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=3300.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=6.52D-16                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.5D0                                                  !Stress exponent dislocation creep
 Ea(i)=5.3D05                                                  !Energia di attivazione dislocation creep
 Va(i)=18.D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=2.37D-15                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.75D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=10.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.D-150                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=20.D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=5.4D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=10.D-06                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 5
 i=5
 t_lmantle=i
 materiale(i)='mantle'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=3300.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=6.52D-16                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.5D0                                                  !Stress exponent dislocation creep
 Ea(i)=5.3D05                                                  !Energia di attivazione dislocation creep
 Va(i)=18.D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.true.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=2.37D-15                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.75D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=10.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.D-150                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=20.D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=5.4D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=10.D-06                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 6
 i=6
 t_sed=i
 materiale(i)='sediments'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=2650.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.28D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=800.D0                                                  !Calore specifico
 Hr_M(i)=1.3D-06                                                  !Calore radiogenico
 Kc_M(i)=2.5D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=8.57D-28                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=4.D0                                                  !Stress exponent dislocation creep
 Ea(i)=2.23D05                                                  !Energia di attivazione dislocation creep
 Va(i)=0.D0                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=6.07931D-19                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=6.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.135723D-28                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=3.2D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=1.23D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=0.D0                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 7
 i=7
 t_serp=i
 materiale(i)='serpentine'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=3000.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.true.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=1.393D-37                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.8D0                                                  !Stress exponent dislocation creep
 Ea(i)=8.9D04                                                  !Energia di attivazione dislocation creep
 Va(i)=3.2D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=6.07931D-19                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=6.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.135723D-28                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=3.2D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=1.23D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=0.D0                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 8
 i=8
 t_wet=i
 materiale(i)='wet mantle'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=3300.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D19                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=6.52D-16                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.5D0                                                  !Stress exponent dislocation creep
 Ea(i)=5.3D05                                                  !Energia di attivazione dislocation creep
 Va(i)=18.D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=2.37D-15                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.75D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=10.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.D-150                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=20.D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=5.4D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=10.D-06                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening

!Materiale 9
 i=9
 t_mantleM=i
 materiale(i)='melt'
 Tref_M(i)=273.15D0                                              !Temperatura di riferimento
 Rho_M(i)=2900.D0                                                  !Densità
 Visc_M(i)=1.D22                                                  !Viscosità
 alpha_M(i)=3.D-05                                                  !Coefficiente espansione termica
 Cp_M(i)=1250.D0                                                  !Calore specifico
 Hr_M(i)=0.D0                                                  !Calore radiogenico
 Kc_M(i)=2.25D0                                                  !Conduttività termica
 rheology(i)=4                                                   !1-lineare; 2-creep; 3-rigido-plastico; 4-visco-plastico
 f_scale(i)=1.D0                                                   !Fattore di scala per la reologia
 Visc_min(i)=1.D18                                                  !Viscosità minima
 Visc_max(i)=1.D25                                                  !Viscosità massima
 uniaxial(i)=.false.                                         !Valore del coefficiente pre-esponenziale usato - .true. se uniassiale - .false. se convertito
!Disclocation creep
 coeffA(i)=6.52D-16                                                  !Coefficiente pre-esponenziale dislocation creep
 ndisl(i)=3.5D0                                                  !Stress exponent dislocation creep
 Ea(i)=5.3D05                                                  !Energia di attivazione dislocation creep
 Va(i)=18.D-06                                                  !Volume di attivazione dislocation creep
!Diffusion creep
 diffusion(i)=.false.                                                 !Utilizzo diffusion creep
 coeffAdf(i)=3.D-11                                                  !Coefficiente pre-esponenziale diffusion creep
 Eadf(i)=3.D05                                                  !Energia di attivazione diffusion creep
 Vadf(i)=4.D-06                                                  !Volume di attivazione diffusion creep
 grain_m(i)=5.D-03                                                  !Grain size diffusion creep
 grain_p(i)=3.D0                                                  !Esponente grain size diffusion creep
!Peierls creep
 peierls(i)=.false.                                                         !Utilizzo Peierls creep
 coeffA_pe(i)=1.D-150                                                  !Coefficiente pre-esponenziale Peierls creep
 ndisl_pe(i)=20.D0                                                  !Stress exponent Peierls creep
 Ea_pe(i)=5.4D05                                                  !Energia di attivazione Peierls creep
 Va_pe(i)=10.D-06                                                  !Volume di attivazione Peierls creep
!Plasticity
 teta_fr_0(i)=25.D0                                                  !Angolo di frizione iniziale
 Cohe_0(i)=2.D07                                                  !Coesione iniziale
 teta_fr_inf(i)=5.D0                                                       !Rapporto angolo di frizione per weakening
 Cohe_inf(i)=4.D06                                                       !Rapporto coesione per weakening
 Ws(i)=10.D0                                                       !Fattore di diminuzione viscosità per weakening


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccc NON MODFICARE cccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(i.ne.num_mat)then
   print *,''
   print *,'ERRORE: Numero materiali non corretto'
   stop
 endif

 do i=1,num_mat
   if(rheology(i).ne.1.and.rheology(i).ne.2.and.rheology(i).ne.3.and.rheology(i).ne.4)then
     print *,''
     print *,'ERRORE: Reologia non corretta per ',materiale(i)
     stop
   endif
 enddo

 end subroutine proprieta