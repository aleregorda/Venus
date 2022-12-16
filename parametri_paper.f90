 subroutine parametri(t_final,v_sx,v_dx,Rhoref,Cpref,Diff,Viscref,alpharef,lito_sx,lito_dx,&
 transition_BC_sx,transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_max,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,&
 e1,e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,OUT_DIR,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,&
 Gy,Rgas,PF,perple_x,idratazione,free_surface,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nb,OUT_FILE_nl,&
 ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,marc_max,dx,file_marc,file_grid,media,vel_corr,OUT_FILE_nl_vtu,&
 OUT_FILE_nl_it,runge_4,OUT_FILE_final_distribution,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,&
 L_coeff,kappa_s,kappa_sh,NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_DISTRIBUTION,MARKERS_CONSTANT,init_num_mat,marc_elem,left_v,&
 right_v,top_v,bottom_v,left_T,right_T,top_T,bottom_T,LAYER_X,LAYER_Y,thick_x,thick_y,n_elem_x,n_elem_y,WEAK_SEED,WEAK_TYPE,&
 WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,WEAK_DIM_Y,THICK_AIR,file_world,WEAK_NUMBER,file_perplexM,file_perplexUC,file_perplexLC,&
 file_perplexO,file_perplexS,file_perplexSp,courant_number,time_step,check_markers,eff_average,t_initial,shear,adiabatic,solver_T,&
 layer_T,T_top,T_bot,y_top,y_bot,radio,cond,parabolic_T,MARKERS_VTU_OUT,MARKERS_TXT_OUT,lagrangian_grid,dim_x_lg,dim_y_lg,&
 old_lgrid,lgrid,serpentine,Press_0,H2O_serp,ocean,weak_plastic,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing,healing_rate,&
 rho_press,thermal_act,OUT_MODEL,OUT_PHASE,DIR_Perplex,H2O_ocean,beta_press,delta_water,melt_ext)

 Implicit none

 logical, intent(out) :: perple_x,idratazione,free_surface,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,&
 OUT_FILE_nl,OUT_FILE_nl_vtu,OUT_FILE_nl_it,vel_corr,runge_4,OUT_FILE_final_distribution,NEW_FILE_GRID,NEW_FILE_MARKERS,&
 MARKERS_CONSTANT,left_T,right_T,top_T,bottom_T,WEAK_SEED,courant_number,check_markers,eff_average,shear,adiabatic,solver_T,&
 parabolic_T,MARKERS_VTU_OUT,MARKERS_TXT_OUT,lagrangian_grid,old_lgrid,serpentine,ocean,weak_plastic,healing,rho_press

 character(len=*), intent(out) :: file_marc,file_grid,media,OUT_DIR,file_world,file_perplexM,file_perplexUC,file_perplexLC,&
 file_perplexO,file_perplexS,file_perplexSp,lgrid,OUT_MODEL,OUT_PHASE,DIR_Perplex

 integer, intent(out) :: ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,marc_max,dx,OUT_FILE_nb,init_num_mat,marc_elem,&
 left_v,right_v,top_v,bottom_v,MARKERS_DISTRIBUTION,LAYER_X,LAYER_Y,WEAK_NUMBER,WEAK_TYPE,layer_T

 double precision, intent(out) :: t_final,v_sx,v_dx,Rhoref,Diff,Viscref,alpharef,Cpref,lito_sx,lito_dx,H2O_serp,&
 transition_BC_sx,transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_max,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,e1,&
 e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,Gy,Rgas,PF,&
 g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,THICK_AIR,time_step,t_initial,&
 dim_x_lg,dim_y_lg,Press_0,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing_rate,thermal_act,H2O_ocean,beta_press,delta_water,&
 melt_ext

 integer, dimension(:), allocatable :: n_elem_x,n_elem_y
 double precision, dimension(:), allocatable :: thick_x,thick_y,WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,WEAK_DIM_Y,T_top,T_bot,y_top,&
 y_bot,radio,cond

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI MODIFICABILI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc Parametri generali cccccccccccccccccccccccccccccccc
 free_surface=.true.                                            !Utilizzo free surface
 power_law=.true.                                            !Utilizzo power law
 solver_T=.true.                                             !Risoluzione dell'equazione dell'energia
 shear=.true.                                         !Utilizzo shear heating nell'equazione dell'energia
 adiabatic=.true.                                         !Utilizzo adiabatic heating nell'equazione dell'energia
 healing=.true.                                            !Calcolo strain healing
 rho_press=.false.                                            !Densità dipendente dalla pressione
 perple_x=.false.                                               !Utilizzo Perple_X per transizioni di fase
 idratazione=.false.                                             !Utilizzo idratazione
 serpentine=.false.                                             !Serpentinizzazione del mantello
 ocean=.false.                                             !Oceanizzazione del melting
 melting=.true.                                            !Utlizzo meccanismo di melting
 erosione=.false.                                            !Utilizzo meccanismo erosione/sedimentazione
 equilib=.false.                                       !Fase di riequilibrazione isostatica
 vel_corr=.true.                                       !Conservative Velocity Interpolation
 runge_4=.false.                                          !Utilizzo 4th-order Runge-Kutta
 check_markers=.true.                                     !Aggiunta/eliminazione marcatori - se .false. il numero di marcatori non viene modificato
 media='geometrica'                                        !Media usata per calcolo viscosità elementale
 file_world='Earth.wb'                                   !Nome file World Builder

!cccccccccccccccccccccccccccc Parametri per il time step ccccccccccccccccccccccccccccc
 courant_number=.true.                                            !Utilizzo Courant number per time step
 t_initial=20.D06                                               !Tempo iniziale
 t_final=30.D06                                               !Tempo finale
 Cn=0.1D0                                                  !Courant number - compreso tra 0 e 1
 time_step=5.D03                                                 !Time step se Courant number=.false.
 max_time_step=1.D05                                                 !Massimo time step consentito

!ccccccccccccccccccccccccccccc Parametri per file output ccccccccccccccccccccccccccccc
 OUT_DIR='../../Output_models/Venus/Quartzite'              !Main Directory in cui scrivere i file di output
 OUT_MODEL='Model_Q1'                                                   !Nome della directory per il modello in corso
 OUT_PHASE='Extension_20-30Myr'                                         !Nome della directory per la fase evolutiva in corso
 OUT_VEL=.true.                                                       !Output velocità - .true. in cm/yr - .false. in m/s
 OUT_TIME=.true.                                                       !Output tempo - .true. in Myr - .false. in s
 OUT_TYPE=.true.                                                       !Intervallo file di output - .true. per tempo - .false. per numero
 OUT_FILE_nb=1                                                       !Numero di iterazioni per il file di output - se OUT_TYPE=.false.
 OUT_FILE_tm=5.D05                                                       !Tempo di uscita per il file di output - se OUT_TYPE=.true.
 OUT_FILE_nl=.false.                                                     !Output file per ogni iterazione non-lineare
 OUT_FILE_nl_vtu=.false.                                                     !Output file vtu per ogni iterazione non-lineare
 OUT_FILE_nl_it=.false.                                                     !Output file con numero iterazioni non-lineari
 OUT_FILE_final_distribution=.true.                                        !Output file con distribuzione finale nodi e marcatori
 MARKERS_VTU_OUT=.false.                                                  !Output file vtu per i marcatori
 MARKERS_TXT_OUT=.false.                                                  !Output file txt per i marcatori

!cccccccccccccccccccccccc Parametri della griglia Lagrangiana cccccccccccccccccccccccc
 lagrangian_grid=.true.                                                  !Visualizzazione griglia Lagrangiana
 xmin_lg=300.D03                                                        !Dimensione lungo X griglia
 xmax_lg=900.D03                                                        !Dimensione lungo Y griglia
 ymin_lg=300.D03                                                        !Dimensione lungo X griglia
 ymax_lg=600.D03                                                        !Dimensione lungo Y griglia
 dim_x_lg=2.D03                                                        !Dimensione lungo X griglia
 dim_y_lg=2.D03                                                        !Dimensione lungo Y griglia
 old_lgrid=.false.                                                        !Lettura vecchio file griglia
 lgrid='../Output_models/Sesia/Isostatic/vtu_files/Grid.00492.vtu'      !Nome vecchio file griglia

!cccccccccccccccccccccccccccc Parametri fisici del sistema ccccccccccccccccccccccccccc
 Rhoref=3200.D0                                                  !Densità di riferimento
 Cpref=1250.D0                                                  !Calore specifico di riferimento
 Diff=1.D-06                                                  !Coefficiente di diffusione termica
 Viscref=1.D21                                                  !Viscosità di riferimento
 alpharef=3.1D-05                                                  !Coefficiente di espansione termica di riferimento
 Gy=-9.81D0                                                       !Gravità
 Press_0=0.D0                                                   !Pressione superficiale

!cccccccccccccccccccccccccccc Boundary condition velocità cccccccccccccccccccccccccccc
 v_sx=-0.5D0                                                  !Velocità lato sinistro in cm/yr
 v_dx=0.5D0                                                  !Velocità lato destro in cm/yr
 left_v=3                                                 !BC lato sx - 1-open; 2-free-slip; 3-fixed
 right_v=3                                                 !BC lato dx - 1-open; 2-free-slip; 3-fixed
 top_v=1                                                 !BC top - 1-open; 2-free-slip; 3-fixed
 bottom_v=2                                                 !BC bottom - 1-open; 2-free-slip; 3-fixed
 lito_sx=120.D03                                          !Spessore litosfera lato sx
 lito_dx=120.D03                                         !Spessore litosfera lato dx
 transition_BC_sx=100.D03                                        !Spessore transizione buondary condition lato sx
 transition_BC_dx=100.D03                                         !Spessore transizione buondary condition lato dx

!cccccccccccccccccccccccccc Boundary condition temperatura ccccccccccccccccccccccccccc
 left_T=.false.                                                 !BC lato sinistro
 right_T=.false.                                                 !BC lato destro
 top_T=.true.                                                 !BC top
 bottom_T=.true.                                                 !BC bottom

!ccccccccccccccccccccccccccc Profilo iniziale temperatura cccccccccccccccccccccccccccc
 parabolic_T=.true.                                        !Incremento profilo termico iniziale non lineare
 LAYER_T=1                                                 !Numero livelli termici con profili non lineari
 allocate(T_top(layer_t),T_bot(layer_t),y_top(layer_t),y_bot(layer_t),radio(layer_t),cond(layer_t)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 T_top(1)=273.15                                             !Temperatura alla base del livello termico con profilo non lineare in K
 T_bot(1)=823.15                                             !Temperatura al top del livello termico con profilo non lineare in K
 y_top(1)=600.D03                                             !Coordinata della base del livello termico con profilo non lineare in m
 y_bot(1)=565.D03                                             !Coordinata del livello termico con profilo non lineare in m
 radio(1)=1.3D-06                                             !Calore radiogenico del livello termico con profilo non lineare
 cond(1)=2.5D0                                             !Conduttività del livello termico con profilo non lineare

!cccccccccccccccccccccccccccccc Parametri della griglia cccccccccccccccccccccccccccccc
 file_grid='Griglia_Earth'               !Nome del file griglia
! file_grid='../../Output_models/Venus/Quartzite/Model_R2/Extension_0-10Myr/Griglia_final'               !Nome del file griglia
 NEW_FILE_GRID=.false.                                          !File griglia - .true. crea nuovo file - .false. utilizza file esistente
!Modificare i parametri successivi solo se NEW_FILE_GRID=.true.
 LAYER_X=5                                                !Numero settori con diversa larghezza degli elementi lungo X
 LAYER_Y=3                                                !Numero settori con diversa altezza degli elementi lungo Y
 allocate(thick_x(layer_x),thick_y(layer_y),n_elem_x(layer_x),n_elem_y(layer_y)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 THICK_X(1)=200.D03                                          !Larghezza settore lungo X in metri
 N_ELEM_X(1)=40                                           !Numero elementi nel settore lungo X
 THICK_X(2)=200.D03                                          !Larghezza settore lungo X in metri
 N_ELEM_X(2)=100                                           !Numero elementi nel settore lungo X
 THICK_X(3)=400.D03                                          !Larghezza settore lungo X in metri
 N_ELEM_X(3)=400                                           !Numero elementi nel settore lungo X
 THICK_X(4)=200.D03                                          !Larghezza settore lungo X in metri
 N_ELEM_X(4)=100                                           !Numero elementi nel settore lungo X
 THICK_X(5)=200.D03                                          !Larghezza settore lungo X in metri
 N_ELEM_X(5)=40                                           !Numero elementi nel settore lungo X
 THICK_Y(1)=400.D03                                            !Altezza settore lungo Y in metri
 N_ELEM_Y(1)=80                                              !Numero elementi nel settore lungo Y
 THICK_Y(2)=80.D03                                            !Altezza settore lungo Y in metri
 N_ELEM_Y(2)=40                                              !Numero elementi nel settore lungo Y
 THICK_Y(3)=120.D03                                            !Altezza settore lungo Y in metri
 N_ELEM_Y(3)=120                                              !Numero elementi nel settore lungo Y

!cccccccccccccccccccccccccccccccc Parametri marcatori cccccccccccccccccccccccccccccccc
 marc_min=16                                             !Marcatori minimi per elemento di minore dimensione - solo se check_markers=.true.
 marc_max=30                                             !Marcatori massimi per elemento di minore dimensione - solo se check_markers=.true.
 THICK_AIR=0.D0                                               !Spessore sticky air
 file_marc='Marcatori_Earth'                                        !Nome del file marcatori
! file_marc='../../Output_models/Venus/Quartzite/Model_Q1/Extension_10-20Myr/Marcatori_final'                                        !Nome del file marcatori
 NEW_FILE_MARKERS=.false.                                  !File marcatori - .true. crea nuovo file - .false. utilizza file esistente
!Modificare i parametri successivi solo se NEW_FILE_MARKERS=.true.
 init_num_mat=5                                            !Numero iniziale dei materiali
 marc_elem=16                                              !Numero di marcatori per elemento di minore dimensione
 MARKERS_CONSTANT=.true.                                          !Numero marcatori costante per ogni elemento
 MARKERS_DISTRIBUTION=1                                  !Distribuzione marcatori - 1-regolare; 2-random; 3-random seguendo una griglia regolare

!cccccccccccccccccccccccccccccccc Parametri weak seed cccccccccccccccccccccccccccccccc
 WEAK_SEED=.true.                                                   !Esistenza di un weak seed
 WEAK_PLASTIC=.true.                                                  !.true. weak plastico - .false. weak viscoso
 WEAK_TYPE=1                                               !Valore strain iniziale - 1-Fixed; 2-Random; 3-Sinusoidale
 WEAK_NUMBER=1                                                    !Numero di weak seed
 allocate(weak_pos_x(weak_number),weak_pos_y(weak_number),weak_dim_x(weak_number),weak_dim_y(weak_number)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 WEAK_POS_X(1)=600.D03                                                         !Posizione X del weak seed in metri
 WEAK_POS_Y(1)=562.D03                                                         !Posizione Y del weak seed in metri
 WEAK_DIM_X(1)=6.D03                                                         !Dimensione X del weak seed in metri
 WEAK_DIM_Y(1)=6.D03                                                         !Dimensione Y del weak seed in metri

!cccccccccccccccccccccccccccccccc Parametri power-law cccccccccccccccccccccccccccccccc
 eff_average=.false.                                    !Calcolo viscosità effettiva markers - .true. media armonica - .false. taglio visc. min e max
 it_max=100                                                       !Iterazioni massime power-law
 tol=1.D-03                                                  !Tolleranza convergenza power-law
 e0=0.D0                                                  !Strain iniziale per weakening plastico
 e1=1.D0                                                  !Strain finale per weakening plastico
 e0v=0.D0                                                  !Strain iniziale per weakening viscoso
 e1v=1.D0                                                  !Strain finale per weakening viscoso

!ccccccccccccccccccccccccccccccccc Parametri healing ccccccccccccccccccccccccccccccccc
 healing_rate=1.D-14                                              !Velocità di strain healing
 thermal_act=100.D0                                               !Thermal activation constant

!ccccccccccccccccccccccccccccccccc Parametri densità ccccccccccccccccccccccccccccccccc
 beta_press=5.124D-12                                         !Coefficiente di comprimibilità
 delta_water=1.5D-02                                         !Coefficiente dipendenza della densità da H2O

!ccccccccccccccccccccccccccccccccc Parametri Perple_X cccccccccccccccccccccccccccccccc
 DIR_Perplex='Phase_changes'
 file_perplexM='Mantle'                          !Nome file Perple_x per il mantello anidro
 file_perplexUC='Upper continental'                          !Nome file Perple_x per il la crosta continentale superiore
 file_perplexLC='Lower continental'                          !Nome file Perple_x per la crosta continental inferiore
 file_perplexO='Ocean'                          !Nome file Perple_x per la crosta oceanica
 file_perplexS='Sediment'                          !Nome file Perple_x per i sedimenti
 file_perplexSp='Serpentine'                          !Nome file Perple_x per il mantello idrato

!ccccccccccccccccccccccccccccc Parametri per idratazione ccccccccccccccccccccccccccccc
 H2O_0=0.D0                                                  !Quantità iniziale di acqua nel mantello in wt%
 H2O_max=0.D0                                                  !Quantità massima di acqua nel mantello in wt%
 H2O_serp=0.D0                                                  !Quantità massima di acqua nel serpentino in wt%
 H2O_ocean=0.D0                                                  !Quantità massima di acqua nell'oceano in wt%
 rw=1.D0                                               !Esponente per la diminuzione della viscosità in presenza di acqua
 OH_0=0.062D0                                               !Limite quantità di acqua per comportamento anidro

!ccccccccccccccccccccccccccccccc Parametri per melting ccccccccccccccccccccccccccccccc
 kappa_melt=60.D0                                 !Costante per il calcolo di deltaT
 gamma=0.54D0                                  !Esponente per il calcolo di deltaT
 chi=12.D0                                    !Costante per il calcolo di acqua massima nel melt
 lambda=0.6D0                                   !Esponente per il calcolo di acqua massima nel melt
 DH2O=0.01D0                                   !Coefficiente di partizione dell'acqua
 A1=1120.66061D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 A2=132.899012D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 A3=-5.1404654D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 x_t=-221.34D0                                 !Coefficiente per calcolo della differenza solidus-liquidus - Kelley et al., 2010
 y_t=536.86D0                                 !Coefficiente per calcolo della differenza solidus-liquidus - Kelley et al., 2010
 tol_melt=1.D-02                                   !Tolleranza per la convergenza del root finder - Falsi point method
 Hl=3.D05                                   !Calore latente
 alpham_dl=-45.D0                                            !Esponente dislocation per la componente di melt
 alpham_df=-30.D0                                            !Esponente diffusion per la componente di melt
 melt_ext=3.D0                                        !Melt extraction in wt%

!cccccccccccccccccccccccccccccc Parametri per erosione ccccccccccccccccccccccccccccccc
!See https://fastscape.org/fastscapelib-fortran/#_fastscape_set_marine_parameters for a complete explanation of all parameters

 kappa_f=2.D-06                                        !Bedrock river incision (SPL) rate parameter in m/yr
 kappa_fsed=2.D-06                                     !Sediments river incision (SPL) rate parameter in m/yr
 esp_m=0.4D0                                                 !Drainage area exponent in the SPL
 esp_n=1.D0                                                 !Slope exponent in the SPL
 kappa_d=1.D-02                                                 !Bedrock transport coefficient (or diffusivity) for hillslope processes in m2/yr
 kappa_dsed=1.D-02                                                 !Sediments transport coefficient (or diffusivity) for hillslope processes in m2/yr
 g_coeff=0.5D0                                                  !Bedrock dimensionless deposition/transport coefficient for the enriched SPL
 gsed_coeff=0.5D0                                                  !Sediments dimensionless deposition/transport coefficient for the enriched SPL
 p_coeff=-1.D0                                                   !Slope exponent for multi-direction flow

!ccccccccccccccccccccccccccc Parametri per sedimentazione cccccccccccccccccccccccccccc
!See https://fastscape.org/fastscapelib-fortran/#_fastscape_set_marine_parameters for a complete explanation of all parameters

 sealevel=600.D03                                        !Sea level in meters
 poro_s=0.49D0                                           !Reference/surface porosity for sand
 poro_sh=0.63D0                                           !Reference/surface porosity for shale 
 zporo_s=3.7D03                                           !e-folding depth for exponential porosity law for sand
 zporo_sh=1.96D03                                           !e-folding depth for exponential porosity law for shale
 ratio=0.5D0                                              !Sand-shale ratio for material leaving the continent
 L_coeff=1.D02                                              !Averaging depth/thickness needed to solve the sand-shale equation in meters
 kappa_s=3.D02                                             !Marine transport coefficient (diffusivity) for sand in m2/yr
 kappa_sh=kappa_s/2.D0                                             !Marine transport coefficient (diffusivity) for shale in m2/yr

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI NON MODIFICABILI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ndof=2                                                            !Gradi di libertà velocità
 ndofT=1                                                            !Gradi di libertà temperatura
 ngauss=4                                                            !Punti di integrazione velocità
 ngaussP=1                                                       !Punti di integrazione pressione
 nodi_elem=4                                                       !Numero nodi per elemento
 dx=11                                                       !Numero nodi markers-chain per elemento
 rtau=1.D0                                                  !Coefficiente correzione convezione
 pi=3.14159265359D0                                             !Pi greco
 teta=0.5D0                                                  !Coefficiente derivata temporale
 eps=10.D-10                                                  !Valore infinitesimale
 Rgas=8.314462618D0                                             !Costante gas
 ip_max=10                                                       !Iterazioni massime divergenza velocità
 tolV=1.D-02                                                  !Tolleranza convergenza divergenza velocità
 PF=1.D06                                                  !Penalty factor

 end subroutine parametri
