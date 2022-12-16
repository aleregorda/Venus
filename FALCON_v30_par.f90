program FALCON

 use MPI
 use WorldBuilder

 implicit none

 include 'dmumps_struc.h'

 interface
   subroutine parametri(t_final,v_sx,v_dx,Rhoref,Cpref,Diff,Viscref,alpharef,lito_sx,lito_dx,transition_BC_sx,&
   transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_wet,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,e1,e0v,e1v,rw,OH_0,alpham_dl,&
   alpham_df,esp_m,esp_n,OUT_DIR,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,Gy,Rgas,PF,perple_x,idratazione,free_surface,&
   power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nb,OUT_FILE_nl,ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,&
   marc_max,dx,file_marc,file_grid,media,vel_corr,OUT_FILE_nl_vtu,OUT_FILE_nl_it,runge_4,OUT_FILE_final_distribution,g_coeff,gsed_coeff,p_coeff,&
   sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_DISTRIBUTION,&
   MARKERS_CONSTANT,init_num_mat,marc_elem,left_v,right_v,top_v,bottom_v,left_T,right_T,top_T,bottom_T,LAYER_X,LAYER_Y,thick_x,thick_y,&
   n_elem_x,n_elem_y,WEAK_SEED,WEAK_TYPE,WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,WEAK_DIM_Y,THICK_AIR,file_world,WEAK_NUMBER,file_perplexM,file_perplexUC,&
   file_perplexLC,file_perplexO,file_perplexS,file_perplexSp,courant_number,time_step,check_markers,eff_average,t_initial,shear,adiabatic,solver_T,&
   layer_T,T_top,T_bot,y_top,y_bot,radio,cond,parabolic_T,markers_vtu_out,markers_txt_out,lagrangian_grid,dim_x_lg,dim_y_lg,old_lgrid,lgrid,&
   serpentine,Press_0,H2O_serp,ocean,weak_plastic,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing,healing_rate,rho_press,thermal_act,&
   OUT_MODEL,OUT_PHASE,DIR_Perplex,H2O_ocean,beta_press,delta_water,melt_ext)
   logical, intent(out) :: perple_x,idratazione,free_surface,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,&
   OUT_FILE_nl,OUT_FILE_nl_vtu,OUT_FILE_nl_it,vel_corr,runge_4,OUT_FILE_final_distribution,NEW_FILE_GRID,NEW_FILE_MARKERS,&
   MARKERS_CONSTANT,left_T,right_T,top_T,bottom_T,WEAK_SEED,courant_number,check_markers,eff_average,shear,adiabatic,solver_T,parabolic_T,&
   markers_vtu_out,markers_txt_out,lagrangian_grid,old_lgrid,serpentine,ocean,weak_plastic,healing,rho_press
   character(len=*), intent(out) :: file_marc,file_grid,media,OUT_DIR,file_world,file_perplexM,file_perplexUC,file_perplexLC,&
   file_perplexO,file_perplexS,file_perplexSp,lgrid,OUT_MODEL,OUT_PHASE,DIR_Perplex
   integer, intent(out) :: ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,marc_max,dx,OUT_FILE_nb,init_num_mat,marc_elem,&
   left_v,right_v,top_v,bottom_v,MARKERS_DISTRIBUTION,LAYER_X,LAYER_Y,WEAK_NUMBER,WEAK_TYPE,layer_T
   double precision, intent(out) :: t_final,v_sx,v_dx,Rhoref,Diff,Viscref,alpharef,Cpref,lito_sx,lito_dx,H2O_serp,&
   transition_BC_sx,transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_wet,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,e1,&
   e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,Gy,Rgas,PF,&
   g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,THICK_AIR,time_step,t_initial,dim_x_lg,&
   dim_y_lg,Press_0,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing_rate,thermal_act,H2O_ocean,beta_press,delta_water,melt_ext
   integer, dimension(:), allocatable :: n_elem_x,n_elem_y
   double precision, dimension(:), allocatable :: thick_x,thick_y,WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,WEAK_DIM_Y,T_top,T_bot,y_top,&
   y_bot,radio,cond
   end subroutine parametri
 end interface

 interface
   subroutine proprieta(i,num_mat,materiale,Tref_M,Rho_M,Visc_M,alpha_M,Cp_M,Hr_M,Kc_M,ndisl,Ea,Va,coeffA,grain_m,grain_p,t_lmantle,&
   Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,t_mantle,t_cont,t_lcont,t_oc,t_air,t_sed,t_serp,t_mantleM,t_contM,t_lcontM,t_ocM,teta_fr_inf,Cohe_inf,Ws,&
   diffusion,rheology,f_scale,t_wet,Visc_min,Visc_max,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,uniaxial)
   integer, intent(out) :: num_mat,i,t_mantle,t_cont,t_lcont,t_oc,t_air,t_sed,t_serp,t_mantleM,t_contM,t_lcontM,t_ocM,t_wet,t_lmantle
   character(len=*), dimension(:), allocatable, intent(out) :: materiale
   logical, dimension(:), allocatable, intent(out) :: diffusion,peierls,uniaxial
   integer, dimension(:), allocatable, intent(out) :: rheology
   double precision, dimension(:), allocatable, intent(out) :: Tref_M,Rho_M,Visc_M,alpha_M,Cp_M,Hr_M,Kc_M,ndisl,Ea,Va,&
   coeffA,Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,teta_fr_inf,Cohe_inf,Ws,f_scale,Visc_min,Visc_max,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,grain_m,grain_p
   end subroutine proprieta
 end interface

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VARIABILI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 logical :: perple_x,idratazione,free_surface,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nl,vel_corr,&
 OUT_FILE_nl_vtu,OUT_FILE_nl_it,runge,final_distribution,NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_CONSTANT,left_T,right_T,top_T,bottom_T,weak_seed,&
 courant_number,check_markers,eff_average,shear,adiabatic,solver_T,parabolic_T,markers_vtu_out,markers_txt_out,lagrangian_grid,old_lgrid,&
 serpentine,ocean,weak_plastic,healing,rho_press

 character(len=150) :: file_marc,file_grid,media,OUT_DIR,file_world,file_perplexM,file_perplexUC,file_perplexLC,file_perplexO,file_perplexS,file_perplexSp,&
 lgrid,OUT_MODEL,OUT_PHASE,DIR_Perplex

 integer :: i,j,k,indice,nodi,elementi,marcatori,ig,it,nt,ip,marcatoritot,ic,ikl,ilg,ill,indof1,indof2,ine1,ine2,t_lcont,nodi_x,nodi_y,&
 elem_x,elem_y,nproc,iproc,counter,chain,colonna,celle,celleP,t_cont,t_oc,t_mantle,t_serp,t_air,t_sed,t_wet,ierr,t_contM,t_ocM,t_mantleM,t_lcontM,&
 ic_max,ic_min,elem_x1,elem_x2,elem_x3,em,ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,marc_max,dx,num_mat,OUT_FILE_nb,riga,nl_out,&
 init_num_mat,marc_elem,left_v,right_v,top_v,bottom_v,MARKERS_DISTRIBUTION,LAYER_X,LAYER_Y,WEAK_NUMBER,weak_type,marc_in,marc_agg,layer_T,elem_x_lg,elem_y_lg,&
 nodi_lg,elem_lg,t_lmantle,numarg

 double precision :: altezza,larghezza,Tgauss,Veltot,x_max,x_min,y_max,y_min,dTX,dPX,minTX,minPX,dPXs,var_grid,dvolu,norm_dP,norm_P,&
 dt,dt0,tempo,D_mu(3,3),D_l(3,3),Ra,Rc,tau1,tau2,tau3,tau,h0,t0,v0,DTemp,lito_sx,lito_dx,transition_BC_sx,transition_BC_dx,Kcref,Href,&
 delta_chain,t1,t2,Cp_X,Rho_X,Rho_XC,alpha_X,Rho_melt,Cp_melt,alpha_melt,norm_res,norm_r0,norm_vel,norm_dvel,t_final,v_sx,v_dx,&
 Ts,Tb,Tref,Rhoref,Cpref,Diff,Viscref,alpharef,OUT_FILE_tm,H2O_0,H2O_wet,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,melt_ext,&
 tol_melt,Hl,e0,e1,e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,kappa_f,kappa_d,tol,rtau,pi,teta,eps,Cn,tolV,Gy,Rgas,PF,g_coeff,gsed_coeff,&
 p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,kappa_fsed,kappa_dsed,Diss,e2rdm,thick_air,time_step,t_initial,&
 dim_x_lg,dim_y_lg,Press_0,H2O_serp,weak,e0w,e1w,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing_rate,thermal_act,H2O_ocean,beta_press,delta_water

 logical, dimension(:), allocatable :: BC_fix,BC_fixT,diffusion,inside,peierls,uniaxial,sed_layer,melt_oc
 character(len=50), dimension(:), allocatable ::materiale
 integer, dimension(:), allocatable :: type_marc,n_marc,C,elem_marc,nmin,nmax,VV,rheology,n_elem_x,n_elem_y,marc_in_elem, n_marc_a,type_marc_a,&
 elem_marc_a,n_marc_tmp,type_marc_tmp,elem_marc_tmp,type_elem,def_elem,em_lg,colonna_marc,colonna_marc_a,colonna_marc_tmp

 double precision, dimension(:), allocatable :: x,y,x_elem,y_elem,x_marc,y_marc,xmax_elem,ymax_elem,delta_x,delta_y,Rho,DRho,RhoC,Visc,Vx,Vy,&
 Load,exx,eyy,exy,Press,Ptot_nodi,Ptot,div,Cp,alpha,Temp,BC_valT,Tnodi,Telem,Hs,s2,Ray,e2t_p,e2t_v,e2t_ist,Pmarc,Tmarc,xmin_elem,ymin_elem,&
 BC_val,Hr,Kc,Hd,Vy_elem,x_chain,y_chain,y_temp,W_e,W_b,Wmax_b,W_f,melt_fraction,W_m,W_e_n,W_b_n,W_m_n,res_loc,res,Vel_old,exx_nodi,eyy_nodi,&
 exy_nodi,e2_nodi,Vel,Pmarc_lito,Plito,Plito_nodi,Rho_M,Cp_M,alpha_M,Visc_M,Hr_M,Kc_M,ndisl,Ea,Va,coeffA,Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,&
 Tref_M,N,f_scale,grain_m,grain_p,dpress,press_old,ptot_old,flux_sx,flux_dx,Load_g,Load_noBC,Visc_vp,teta_fr_inf,Cohe_inf,Ws,thick_x,thick_y,&
 weak_pos_x,weak_pos_y,weak_dim_x,weak_dim_y,x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_m_a,melt_fraction_a,W_f_a,x_marc_tmp,y_marc_tmp,e2t_p_tmp,&
 e2t_v_tmp,W_b_tmp,W_f_tmp,W_m_tmp,melt_fraction_tmp,T_top,T_bot,y_top,y_bot,radio,cond,q_x,q_y,teta_elem,cohe_elem,e2t_p_elem,e2t_v_elem,x_lg,y_lg,&
 Visc_min,Visc_max,melt_elem,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,w_b_elem

 integer, dimension(:,:), allocatable :: connessioni,marc_each_elem,conn_lg
 double precision, dimension(:,:), allocatable :: K_loc,K_Mc,K_K,Vel1,Vel2,RhoXs,alphaXs,CpXs,WXs,dC,B,N1,N2,NTT,Vel_res
 double precision, dimension(:,:,:), allocatable :: RhoX,alphaX,CpX,WX
 type(dmumps_struc)idV
 type(dmumps_struc)idT

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     INIZIALIZZAZIONE VETTORI MUMPS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 call mpi_init(ierr)
 call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
 call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)

 idV%COMM=MPI_COMM_WORLD
 idV%SYM=1
 idV%par=1
 idV%JOB=-1
 call DMUMPS(idV)

 idT%COMM=MPI_COMM_WORLD
 idT%SYM=0
 idT%par=1
 idT%JOB=-1
 call DMUMPS(idT)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(iproc.eq.0)then
   print *,'****************************************************************************************'
   print *,'------------------------------------- START PROGRAM ------------------------------------'
   print *,'****************************************************************************************'
 endif

 call parametri(t_final,v_sx,v_dx,Rhoref,Cpref,Diff,Viscref,alpharef,lito_sx,lito_dx,transition_BC_sx,&
      transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_wet,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,e1,e0v,e1v,rw,OH_0,alpham_dl,&
      alpham_df,esp_m,esp_n,OUT_DIR,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,Gy,Rgas,PF,perple_x,idratazione,free_surface,&
      power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nb,OUT_FILE_nl,ndof,ndofT,ngauss,ngaussP,nodi_elem,it_max,ip_max,marc_min,&
      marc_max,dx,file_marc,file_grid,media,vel_corr,OUT_FILE_nl_vtu,OUT_FILE_nl_it,runge,final_distribution,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,&
      poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_DISTRIBUTION,MARKERS_CONSTANT,init_num_mat,marc_elem,&
      left_v,right_v,top_v,bottom_v,left_T,right_T,top_T,bottom_T,LAYER_X,LAYER_Y,thick_x,thick_y,n_elem_x,n_elem_y,WEAK_SEED,WEAK_TYPE,WEAK_POS_X,WEAK_POS_Y,&
      WEAK_DIM_X,WEAK_DIM_Y,THICK_AIR,file_world,WEAK_NUMBER,file_perplexM,file_perplexUC,file_perplexLC,file_perplexO,file_perplexS,file_perplexSp,&
      courant_number,time_step,check_markers,eff_average,t_initial,shear,adiabatic,solver_T,layer_T,T_top,T_bot,y_top,y_bot,radio,cond,parabolic_T,&
      markers_vtu_out,markers_txt_out,lagrangian_grid,dim_x_lg,dim_y_lg,old_lgrid,lgrid,serpentine,Press_0,H2O_serp,ocean,weak_plastic,xmin_lg,xmax_lg,&
      ymin_lg,ymax_lg,max_time_step,healing,healing_rate,rho_press,thermal_act,OUT_MODEL,OUT_PHASE,DIR_Perplex,H2O_ocean,beta_press,delta_water,melt_ext)
! call changed_parameters(iproc,numarg,OUT_DIR,OUT_MODEL,OUT_PHASE,t_initial,t_final,Cn,time_step,max_time_step,power_law)
 if(equilib)then
   v_sx=0.D0
   v_dx=0.D0
   top_v=2
   bottom_v=2
   right_v=2
   left_v=2
 endif
 if(free_surface)top_v=1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PROPRIETA' DEI MATERIALI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(iproc.eq.0)then
   call proprieta(i,num_mat,materiale,Tref_M,Rho_M,Visc_M,alpha_M,Cp_M,Hr_M,Kc_M,ndisl,Ea,Va,coeffA,grain_m,grain_p,t_lmantle,Eadf,Vadf,&
        coeffAdf,teta_fr_0,Cohe_0,t_mantle,t_cont,t_lcont,t_oc,t_air,t_sed,t_serp,t_mantleM,t_contM,t_lcontM,t_ocM,teta_fr_inf,Cohe_inf,Ws,&
        diffusion,rheology,f_scale,t_wet,Visc_min,Visc_max,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,uniaxial)

   call system('mkdir -p '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE))
   call system('mkdir -p '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/vtu_files')
   call system('mkdir -p '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/txt_files')
   if(numarg.gt.0)then
     call system('mv -f "changed_parameters.txt" '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/changed_parameters.txt')
   else
     call system('cp -f "parametri.f90" '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/parametri_ref.f90')
     call system('cp -f "proprieta.f90" '//TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/proprieta_ref.f90')
   endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     1 - LETTURA FILE GRIGLIA E TEMPERATURA INIZIALE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   if(NEW_FILE_GRID)then
     allocate(x((SUM(n_elem_x)+1)*(SUM(n_elem_y)+1)),y((SUM(n_elem_x)+1)*(SUM(n_elem_y)+1)))
     call Create_grid(layer_x,layer_y,n_elem_x,n_elem_y,SUM(n_elem_x)+1,SUM(n_elem_y)+1,(SUM(n_elem_x)+1)*(SUM(n_elem_y)+1),SUM(n_elem_x),SUM(n_elem_y),&
          thick_x,thick_y,x,y,file_grid,file_world,thick_air,Gy)
     deallocate(x,y)
   endif
   print *,''
   print *,'Reading grid file...'
   open(unit=1,file=file_grid,status='old') 
   rewind(1)
   read(1,*)  
   read(1,*)layer_x,(thick_x(i),n_elem_x(i),i=1,layer_x)
   read(1,*)layer_y,(thick_y(i),n_elem_y(i),i=1,layer_y)
   read(1,*)
   read(1,*)
   elem_x=SUM(n_elem_x)
   elem_y=SUM(n_elem_y)
   elementi=elem_x*elem_y
   nodi_x=elem_x+1
   nodi_y=elem_y+1
   nodi=nodi_x*nodi_y
   allocate(x_elem(elementi),y_elem(elementi),xmax_elem(elementi),ymax_elem(elementi),xmin_elem(elementi),ymin_elem(elementi),delta_x(elementi),&
            delta_y(elementi),C(elementi),nmin(elementi),nmax(elementi),Visc(elementi),Hr(elementi),Kc(elementi),W_e(elementi),Rho(elementi),&
            RhoC(elementi),Cp(elementi),alpha(elementi),Hs(elementi),Hd(elementi),Vy_elem(elementi),q_x(elementi),q_y(elementi),e2t_p_elem(elementi),&
            e2t_v_elem(elementi),melt_elem(elementi),sed_layer(elem_x),w_b_elem(elementi),marc_in_elem(elementi))
   allocate(melt_oc(elem_x))
   allocate(connessioni(elementi,nodi_elem))
   allocate(x(nodi),y(nodi),Temp(nodi),Vx(nodi),Vy(nodi))
   allocate(BC_val(nodi*ndof),BC_fix(nodi*ndof),BC_valT(nodi*ndofT),BC_fixT(nodi*ndofT))
   allocate(N(nodi_elem),dC(nodi_elem,ndof),B(3,nodi_elem*ndof))
   do i=1,nodi
     read(1,*)k,x(i),y(i),Temp(i)
     if(parabolic_T)then
       do j=1,layer_T
         if(y(i).ge.y_bot(j).and.y(i).le.y_top(j))Temp(i)=Temp(i)-((radio(j)/(2.D0*cond(j)))*(y(i)**2-y_bot(j)**2))+&
                                                                  ((radio(j)/(2.D0*cond(j)))*(y_top(j)+y_bot(j))*(y(i)-y_bot(j)))
       enddo
     endif
   enddo
   read(1,*)
   do i=1,elementi
     read(1,*)k,(connessioni(i,j),j=1,nodi_elem)
     xmax_elem(i)=MAXVAL(x(connessioni(i,:)))
     xmin_elem(i)=MINVAL(x(connessioni(i,:)))
     ymax_elem(i)=MAXVAL(y(connessioni(i,:)))
     ymin_elem(i)=MINVAL(y(connessioni(i,:)))
     delta_x(i)=DABS(xmax_elem(i)-xmin_elem(i))
     delta_y(i)=DABS(ymax_elem(i)-ymin_elem(i))
     if(NEW_FILE_MARKERS)then
       if(MARKERS_CONSTANT)then
         marc_in_elem(i)=marc_elem
       else
         marc_in_elem(i)=INT(marc_elem*(delta_x(i)/MINVAL(delta_x))*(delta_y(i)/MINVAL(delta_y)))
       endif
     endif
     if(MARKERS_CONSTANT)then
       nmin(i)=marc_min
       nmax(i)=marc_max
     else
       nmin(i)=INT(marc_min*(delta_x(i)/MINVAL(delta_x))*(delta_y(i)/MINVAL(delta_y)))
       nmax(i)=INT(marc_max*(delta_x(i)/MINVAL(delta_x))*(delta_y(i)/MINVAL(delta_y)))
     endif
   enddo
   Ts=MINVAL(Temp)
   Tb=MAXVAL(Temp)
   Tref=Tb
   if(Tb.eq.Ts)Tb=Ts+1.D0
   Temp=Temp/(Tb-Ts)
   Hs=0.D0
   Hd=0.D0
   x_max=MAXVAL(x)
   x_min=MINVAL(x)
   y_max=MAXVAL(y)
   y_min=MINVAL(y)
   larghezza=DABS(x_max-x_min)
   altezza=DABS(y_max-y_min)
   print 1003,larghezza/1.D03
   print 1004,altezza/1.D03
   print 1001,nodi
   print 1002,elementi
   if(lagrangian_grid)then
     elem_x_lg=INT(DABS(xmax_lg-xmin_lg))/dim_x_lg
     elem_y_lg=INT(DABS(ymax_lg-ymin_lg))/dim_y_lg
     nodi_lg=(elem_x_lg+1)*(elem_y_lg+1)
     elem_lg=elem_x_lg*elem_y_lg
     allocate(x_lg(nodi_lg),y_lg(nodi_lg),conn_lg(elem_lg,4),em_lg(nodi_lg))
     if(.not.old_lgrid)then
       call Create_lgrid(dim_x_lg,dim_y_lg,x_lg,y_lg,elem_x_lg,elem_y_lg,nodi_lg,elem_lg,conn_lg,xmin_lg,ymin_lg)
     else
       call Read_lgrid(x_lg,y_lg,nodi_lg,elem_lg,conn_lg,lgrid)
     endif
   endif
 endif
 if(nproc.gt.1)call mpi_bcast(altezza,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 if(nproc.gt.1)call mpi_bcast(elementi,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 if(nproc.gt.1)call mpi_bcast(nodi,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 allocate(Ptot_old(nodi),Ptot_nodi(nodi),Vel_old(nodi*ndof),Vel(nodi*ndof),Vel_res(elementi,nodi_elem*ndof))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     2 -  CREAZIONE VETTORI MUMPS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 call mumps_vectors(iproc,idV,nodi,elementi,nodi_elem,ndof,connessioni)
 idV%ICNTL(3)=10
 if(solver_T)then
   call mumps_vectors(iproc,idT,nodi,elementi,nodi_elem,ndofT,connessioni)
   idT%ICNTL(3)=11
 endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     3 - ADIMENSIONALIZZAZIONE E BOUNDARY CONDITIONS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 h0=altezza
 Kcref=Diff*Rhoref*Cpref
 Href=Diff*Cpref*(Tb-Ts)/(h0**2)
 t0=(h0**2)/Diff
 Diss=alpharef*Gy*h0/Cpref
 if(iproc.eq.0)then
   if(OUT_TIME)time_step=time_step*365.D0*24.D0*3600.D0
   if(OUT_TIME)max_time_step=max_time_step*365.D0*24.D0*3600.D0
   Ra=(Rhoref*alpharef*Gy*(altezza**3)*(Tb-Ts))/(Viscref*Diff)
   Rc=(Gy*(altezza**3))/(Viscref*Diff)
   call Dimensionless(nodi,elementi,h0,Diff,x,y,delta_x,delta_y,x_min,x_max,y_min,y_max,var_grid,thick_air,lito_sx,lito_dx,transition_BC_sx,&
        transition_BC_dx,xmax_elem,xmin_elem,ymax_elem,ymin_elem,larghezza,altezza,v0,thick_x,thick_y,layer_x,layer_y,time_step,max_time_step)
   call BC(nodi,ndof,ndofT,y_max,x,y,v0,BC_fix,BC_val,lito_sx,lito_dx,transition_BC_sx,transition_BC_dx,thick_air,&
        v_sx,v_dx,Temp,BC_fixT,BC_valT,left_v,right_v,top_v,bottom_v,left_T,right_T,top_T,bottom_T,nodi_x,nodi_y)
   if(lagrangian_grid)then
     x_lg=x_lg/h0
     y_lg=y_lg/h0
     dim_x_lg=dim_x_lg/h0
     dim_y_lg=dim_y_lg/h0
   endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     4 - LETTURA FILE MARCATORI E COMPOSIZIONE INIZIALE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   if(.not.free_surface)then
     allocate(y_temp(elem_x))
     y_temp=y_min
   endif
   if(NEW_FILE_MARKERS)call Create_markers(file_marc,file_world,marc_in_elem,elementi,MARKERS_DISTRIBUTION,init_num_mat,xmin_elem*h0,ymin_elem*h0,&
                            delta_x*h0,delta_y*h0,y_min*h0,y_max*h0,thick_air*h0)
   deallocate(marc_in_elem)
   print *,''
   print *,'Calculating initial composition...'
   open(unit=1,file=file_marc,status='old')
   rewind(1)
   read(1,*)
   read(1,*)marcatori
   read(1,*)
   allocate(n_marc(marcatori),x_marc(marcatori),y_marc(marcatori),type_marc(marcatori),e2t_p(marcatori),e2t_v(marcatori),W_b(marcatori),&
   W_m(marcatori),melt_fraction(marcatori),W_f(marcatori))
   C=0
   do i=1,marcatori
     read(1,*)n_marc(i),x_marc(i),y_marc(i),type_marc(i),e2t_p(i),e2t_v(i),W_b(i),W_m(i),W_f(i),melt_fraction(i)
     x_marc(i)=x_marc(i)/h0
     y_marc(i)=y_marc(i)/h0
   enddo
   marcatoritot=MAXVAL(n_marc)
   allocate(elem_marc(marcatoritot),colonna_marc(marcatoritot))
   do i=1,marcatori
     call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_marc(i),y_marc(i),&
          elem_marc(n_marc(i)),x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,&
          colonna_marc(n_marc(i)),riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
     C(elem_marc(n_marc(i)))=C(elem_marc(n_marc(i)))+1
     W_e(elem_marc(n_marc(i)))=W_e(elem_marc(n_marc(i)))+W_f(i)
     if(weak_seed)then
       do j=1,weak_number
         if(x_marc(i)*h0.ge.WEAK_POS_X(j)-WEAK_DIM_X(j)/2.D0.and.x_marc(i)*h0.le.WEAK_POS_X(j)+WEAK_DIM_X(j)/2.D0.and.&
            y_marc(i)*h0.ge.WEAK_POS_Y(j)-WEAK_DIM_Y(j)/2.D0.and.y_marc(i)*h0.le.WEAK_POS_Y(j)+WEAK_DIM_Y(j)/2.D0)then
           if(weak_plastic)then
             e0w=e0
             e1w=e1
           else
             e0w=e0v
             e1w=e1v
           endif
           if(weak_type.eq.1)then
             weak=e1w
           elseif(weak_type.eq.2)then
             call random_number(e2rdm)
             weak=(e2rdm+e0w)*DEXP(-(x_marc(i)*h0-WEAK_POS_X(j))**2/(2.D0*(WEAK_DIM_X(j)/8.D0)**2))
           elseif(weak_type.eq.3)then
             weak=(((1.D0-DCOS(2.D0*pi*((x_marc(i)*h0-WEAK_POS_X(j))/6.D03)))+(1.D0-DCOS(2.D0*pi*((y_marc(i)*h0-WEAK_POS_Y(j))/6.D03))))/4.D0)+e0w
           endif
           if(weak_plastic)then
             e2t_p(i)=weak
           else
             e2t_v(i)=weak
           endif
         endif
       enddo
     endif
     if(.not.free_surface.and.y_temp(colonna).lt.y_marc(i).and.(type_marc(i).eq.t_cont.or.type_marc(i).eq.t_oc))y_temp(colonna)=y_marc(i)
   enddo
   print 1005,marcatori
   close(1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     5 - UPPER BOUNDARY E CAMBI DI FASE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   chain=0
   if(.not.free_surface)then
     chain=((dx-1)*elem_x)+1
     allocate(x_chain(chain),y_chain(chain))
     delta_chain=larghezza/chain
     do i=1,chain
       x_chain(i)=x_min+(larghezza/(chain-1))*(i-1)
       call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_chain(i),y_min,&
            em,x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,&
            colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
       call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_chain(i),y_temp(colonna),&
            em,x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,&
            colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
       y_chain(i)=ymax_elem(em)-eps
       if(i.gt.1.and.y_chain(i).ne.y_chain(i-1))y_chain(i)=(y_temp(colonna)-y_chain(i-1))*&
         ((x_chain(i)-x_chain(i-1))/(xmax_elem(em)-x_chain(i-1)))+y_chain(i-1)
     enddo
     deallocate(y_temp)
     if(thick_air.eq.0.D0)y_chain=y_max
   endif

   if(perple_x)call File_Perple_X(celle,celleP,dTX,dPX,dPXs,minTX,minPX,RhoXs,CpXs,alphaXs,WXs,RhoX,CpX,alphaX,WX,DIR_Perplex,&
                    file_perplexM,file_perplexUC,file_perplexLC,file_perplexO,file_perplexS,file_perplexSp)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     6 - INIZIO CICLO TEMPORALE E CALCOLO PROPRIETA' ELEMENTALI
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   it=0
   dt0=0.D0
   if(OUT_TYPE)then
     nt=0
   else
     nt=OUT_FILE_nb+1
   endif
 endif
 dt=0.D0
 if(OUT_TIME)t_initial=t_initial/1.D06
 if(OUT_TIME)t_final=t_final/1.D06
 if(OUT_TIME)out_file_tm=out_file_tm/1.D06
 tempo=t_initial
 if(power_law.and..not.equilib)then
   ic_max=it_max
   ic_min=2
 else
   ic_max=1
   ic_min=1
 endif

 do while(tempo.le.t_final)
 if(iproc.eq.0)then
   it=it+1
   print *,''
   print *,'****************************************************************************************'
   print *,''
   print 1006,it
   print 1007,dt
   print 1008,tempo
   print 1012,marcatori
   print *,''

   call cpu_time(t1)
   print *,'- Elemental properties calculation:'
   allocate(DRho(elementi),Telem(elementi))
   x_elem=0.D0
   y_elem=0.D0
   Telem=0.D0
   N=0.25D0
   do i=1,elementi
     do j=1,nodi_elem
       x_elem(i)=x_elem(i)+x(connessioni(i,j))*N(j)
       y_elem(i)=y_elem(i)+y(connessioni(i,j))*N(j)
       Telem(i)=Telem(i)+Temp(connessioni(i,j))*N(j)
     enddo
   enddo
   Rho=0.D0
   Cp=0.D0
   alpha=0.D0
   Hr=0.D0
   Kc=0.D0
   RhoC=0.D0
   e2t_p_elem=0.D0
   e2t_v_elem=0.D0
   melt_elem=0.D0
   w_b_elem=0.D0
   if(it.eq.1.or..not.power_law)then
     if(media.eq.'aritmetica'.or.media.eq.'armonica')then
       Visc=0.D0
     else
       Visc=1.D0
     endif
   endif
   do i=1,marcatori
     if(type_marc(i).eq.t_mantleM)then
       Rho_X=Rho_M(t_mantle)
     else
       Rho_X=Rho_M(type_marc(i))
     endif
     call Elemental_properties(C(elem_marc(n_marc(i))),Cp_M(type_marc(i)),Rho_X,alpha_M(type_marc(i)),Rho_M(t_mantleM),Cp_M(t_mantleM),&
          alpha_M(t_mantleM),Rho_X,Hr_M(type_marc(i)),Kc_M(type_marc(i)),Rhoref,Cpref,Kcref,RhoC(elem_marc(n_marc(i))),&
          Rho(elem_marc(n_marc(i))),Cp(elem_marc(n_marc(i))),alpha(elem_marc(n_marc(i))),Hr(elem_marc(n_marc(i))),Kc(elem_marc(n_marc(i))),&
          melt_fraction(i),alpharef,Href,e2t_p(i),e2t_v(i),e2t_p_elem(elem_marc(n_marc(i))),e2t_v_elem(elem_marc(n_marc(i))),&
          melt_elem(elem_marc(n_marc(i))),W_b(i),w_b_elem(elem_marc(n_marc(i))))
     if(it.eq.1.or..not.power_law)then
       if(media.eq.'aritmetica')then
         Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))+(Visc_M(type_marc(i))/DBLE(C(elem_marc(n_marc(i)))))
       elseif(media.eq.'armonica')then
         Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))+(1.D0/Visc_M(type_marc(i)))
       else
         Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))*(Visc_M(type_marc(i))**(1.D0/DBLE(C(elem_marc(n_marc(i))))))
       endif
     endif
   enddo
   DRho=RhoC-Rhoref
   if(it.eq.1.or..not.power_law)then
     if(media.eq.'armonica')Visc=DBLE(C)/Visc
     Visc=Visc/Viscref
   endif
   if(equilib)Visc=1.D21/Viscref
   call cpu_time(t2)
   print 1013,t2-t1
   print *,''
 endif
 if(nproc.gt.1)call mpi_bcast(it,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     7 - BILANCIO DEL MOMENTO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 allocate(div(elementi),res(nodi*ndof),Load_g(nodi*ndof))
 if(it.eq.1)then
   Vel_old=1.D0
   Ptot_old=1.D0
   norm_dvel=1.D0
   norm_vel=1.D0
   norm_dP=1.D0
   norm_P=1.D0
   norm_res=1.D0
   norm_r0=1.D0
   Vel_res=0.D0
 endif
 if(iproc.eq.0)then
   allocate(exx(elementi),eyy(elementi),exy(elementi),s2(elementi),Press(elementi),Plito(elementi),Ptot(elementi),Ray(elementi),&
            W_e_n(elementi),dpress(elementi),press_old(elementi),teta_elem(elementi),cohe_elem(elementi))
   allocate(exx_nodi(nodi),eyy_nodi(nodi),exy_nodi(nodi),e2_nodi(nodi),Plito_nodi(nodi))
   allocate(e2t_ist(marcatori),Visc_vp(marcatori),Pmarc(marcatori),Pmarc_lito(marcatori),Tmarc(marcatori),&
            Wmax_b(marcatori),W_b_n(marcatori),W_m_n(marcatori),VV(marcatori))
   allocate(flux_sx(elem_y),flux_dx(elem_y))
   W_b_n=0.D0
   W_m_n=0.D0
   VV=1
   Press_old=0.D0
   ic=1
   nl_out=1

   open(unit=10,file='mumps_info_v',status='unknown')
   rewind(10)

!c     7.1 - CALCOLO VELOCITA'

   call cpu_time(t1)
   print *,'- Momentum solution:'
   print *,'  Facticulating...'
 endif
 if(nproc.gt.1)call mpi_bcast(ic,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 do while((MAX(norm_res/norm_r0,norm_dvel/norm_vel,norm_dP/norm_P).gt.tol.and.ic.le.ic_max).or.(ic.le.ic_min))                !INIZIO CICLO NON LINEARITA'
   if(iproc.eq.0)then
     ip=1
     div=1.D0
     Press=0.D0
     if(perple_x)then
       Ray=(Rhoref*alpha*Gy*((altezza*h0)**3)*(Tb-Ts))/(Viscref*Diff)
     else
       Ray=Ra
     endif
   endif
   if(nproc.gt.1)call mpi_bcast(div,elementi,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   if(nproc.gt.1)call mpi_bcast(ip,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   do while((MAXVAL(DABS(div*t0)).gt.tolV.and.ip.le.ip_max).or.(ip.le.2))                  !INIZIO CICLO PENALTY
     if(iproc.eq.0)then
       ip=ip+1
       allocate(K_loc(nodi_elem*ndof,nodi_elem*ndof),Load(nodi_elem*ndof),Load_noBC(nodi_elem*ndof),res_loc(nodi_elem*ndof))
       idV%RHS=0.D0
       idV%A_ELT=0.D0
       Load_g=0.D0
       res=0.D0
       counter=0
       do j=1,elem_y
         do i=1,elem_x
           indice=elem_x*(j-1)+i
           K_loc=0.D0
           Load=0.D0
           do ig=1,ngauss
             Tgauss=0.D0
             call Matrici(ngauss,ndof,nodi_elem,elementi,nodi,connessioni,x,y,indice,ig,B,D_mu,D_l,dvolu,N,dC)
             do k=1,nodi_elem
               Tgauss=Tgauss+Temp(connessioni(indice,k))*N(k)
             enddo
             DTemp=Tgauss-(Tref/(Tb-Ts))
             K_loc=K_loc+dvolu*Visc(indice)*matmul(matmul(transpose(B),D_mu),B)
             do k=1,nodi_elem
               Load((k*ndof)-1)=Load((k*ndof)-1)+dvolu*(dC(k,1)*(Press(indice)*t0/Viscref))
               if(free_surface.or.perple_x)then
                 Load(k*ndof)=Load(k*ndof)+dvolu*(Rc*RhoC(indice)*N(k)-Ray(i)*DTemp*N(k)+(dC(k,2)*(Press(indice)*t0/Viscref)))
               else
                 Load(k*ndof)=Load(k*ndof)+dvolu*(Rc*DRho(indice)*N(k)-Ray(i)*DTemp*N(k)+(dC(k,2)*(Press(indice)*t0/Viscref)))
               endif
             enddo
           enddo
           Load_noBC=Load*Viscref/(t0*h0)
           ig=1
           call Matrici(ngaussP,ndof,nodi_elem,elementi,nodi,connessioni,x,y,indice,ig,B,D_mu,D_l,dvolu,N,dC)
           K_loc=K_loc+dvolu*PF*Visc(indice)*matmul(matmul(transpose(B),D_l),B)
           if(free_surface)call Stabilization(j,indice,dt0,RhoC,K_loc,elementi,elem_y,x,y,connessioni,nodi_elem,nodi,ndof)
           call BC_fixed(ndof,nodi_elem,nodi,connessioni,elementi,indice,BC_fix,BC_val,K_loc,Load)
           res_loc=matmul(K_loc*Viscref/(h0**2),Vel_res(indice,:)*v0)-Load*Viscref/(t0*h0)
           do ine1=1,nodi_elem
             do indof1=1,ndof
               ill=ndof*(ine1-1)+indof1
               ilg=ndof*(connessioni(indice,ine1)-1)+indof1
               do ine2=1,nodi_elem
                 do indof2=1,ndof
                   ikl=ndof*(ine2-1)+indof2
                   if(ikl.ge.ill)then
                     counter=counter+1
                     idV%A_ELT(counter)=K_loc(ill,ikl)
                   endif  
                 enddo
               enddo
               idV%RHS(ilg)=idV%RHS(ilg)+Load(ill)
               Load_g(ilg)=Load_g(ilg)+Load_noBC(ill)
               res(ilg)=res(ilg)+res_loc(ill)
             enddo
           enddo
         enddo
       enddo
       deallocate(K_loc,Load,Load_noBC,res_loc)
     endif
     if(nproc.gt.1)call mpi_bcast(ip,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
     idV%ICNTL(5)=1
     idV%ICNTL(7)=0
     idV%JOB=6
     call DMUMPS(idV)

     if(iproc.eq.0)then
       do i=1,nodi
        Vx(i)=idV%RHS((i*ndof)-1)
        Vy(i)=idV%RHS(i*ndof)
       enddo
       Vel=idV%RHS
       call Elemental_Pressure(elem_x,elem_y,delta_y,connessioni,elementi,nodi_elem,ngaussP,ndof,nodi,x,y,Vx,Vy,&
       flux_sx,flux_dx,exx,eyy,exy,Vy_elem,div,Press,h0,t0,Viscref,PF,Kc,Temp,q_x,q_y)
     endif
     if(nproc.gt.1)call mpi_bcast(div,elementi,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   enddo                                                                          !FINE CICLO PENALTY
   Press=Press+Press_0

!c     7.2 - CALCOLO PRESSIONE SMOOTHED, CAMBI DI FASE, MELTING E VISCOSITA'

   if(iproc.eq.0)then
     if(courant_number)then
       dt0=Cn*MIN(MIN(DABS(MINVAL(delta_x)/MAXVAL(DABS(Vx))),DABS(MINVAL(delta_y)/MAXVAL(DABS(Vy)))),&
           MIN(MINVAL(delta_x)**2,MINVAL(delta_y)**2)/Diff)
     else
         dt0=time_step
     endif
     if(dt0.gt.max_time_step)dt0=max_time_step
     if(OUT_TIME)then
       dt=dt0*t0/(1.D06*365.D0*24.D0*3600.D0)
     else
       dt=dt0*t0
     endif
     dpress=DABS(Press-Press_old)
     Press_old=Press
     call Smooth_Pressure(elementi,elem_x,elem_y,nodi_elem,nodi,connessioni,chain,x_chain*h0,y_chain*h0,&
          x*h0,y*h0,delta_x*h0,delta_y*h0,larghezza*h0,Press,free_surface,N,exx,eyy,exy,Plito,Rho*Rhoref,Gy,thick_air*h0,Press_0)
     if(free_surface.or.perple_x)then
       Ptot=Press
     else
       Ptot=Press+Plito
     endif
     call Nodal_Pressure(nodi_elem,delta_x,delta_y,connessioni,elementi,nodi,Ptot_nodi,Ptot,Plito,Plito_nodi)
     call Nodal_Strain(elementi,nodi,nodi_elem,connessioni,ndof,xmin_elem,delta_x,t0,x,y,Vx,Vy,exx_nodi,eyy_nodi,exy_nodi,e2_nodi)
     call PT(elem_marc,n_marc,marcatori,marcatoritot,connessioni,elementi,nodi,delta_x,xmin_elem,x_marc,y_marc,y,nodi_elem,Temp,&
          Ptot_nodi,exx_nodi,eyy_nodi,exy_nodi,Tb,Ts,ndof,e2t_ist,Tmarc,Pmarc,Plito_nodi,Pmarc_lito,num_mat,materiale,type_marc)
     Rho=0.D0
     Cp=0.D0
     alpha=0.D0
     Hr=0.D0
     Kc=0.D0
     RhoC=0.D0
     e2t_p_elem=0.D0
     e2t_v_elem=0.D0
     melt_elem=0.D0
     w_b_elem=0.D0
     W_f=0.D0
     Wmax_b=0.D0
     if(idratazione)call water_markers(perple_x,it,marcatori,num_mat,type_marc,marcatoritot,elem_marc,celle,celleP,n_marc,&
                                      elementi,elem_x,W_b,W_b_n,Wmax_b,WX,WXs,minTX,minPX,dTX,dPX,dPXs,W_e,W_e_n,materiale,t_mantle,&
                                      t_serp,t_wet,C,W_f,Tmarc,Pmarc,H2O_0,H2O_wet,e1,e2t_p,VV,Pmarc_lito,H2O_serp,H2O_ocean,ic)
     do i=1,marcatori
       call melt(materiale(type_marc(i)),W_b_n(i),W_f(i),W_e_n(elem_marc(n_marc(i))),melt_fraction(i),Pmarc(i),Tmarc(i),&
            Rho_melt,Rho_M,alpha_M,Cp_M,Tref_M,num_mat,Cp_melt,alpha_melt,W_m(i),W_m_n(i),Wmax_b(i),kappa_melt,gamma,melting,&
            chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,t_oc,t_mantleM,t_lmantle,type_marc(i),beta_press,melt_oc(colonna_marc(n_marc(i))),ocean,melt_ext)
       call PerpleX(rho_press,Rho_X,Cp_X,alpha_X,Tmarc(i),Pmarc(i),materiale(type_marc(i)),dTX,dPX,minTX,minPX,RhoX,alphaX,CpX,dPXs,RhoXs,CpXs,alphaXs,&
            Rho_XC,celle,celleP,Cpref,alpharef,perple_x,Rho_M(type_marc(i)),Cp_M(type_marc(i)),alpha_M(type_marc(i)),Tref_M(type_marc(i)),&
            Rho_M(t_mantle),W_b(i),beta_press,delta_water)
       call Elemental_properties(C(elem_marc(n_marc(i))),Cp_X,Rho_X,alpha_X,Rho_melt,Cp_melt,alpha_melt,Rho_XC,Hr_M(type_marc(i)),Kc_M(type_marc(i)),&
            Rhoref,Cpref,Kcref,RhoC(elem_marc(n_marc(i))),Rho(elem_marc(n_marc(i))),Cp(elem_marc(n_marc(i))),alpha(elem_marc(n_marc(i))),&
            Hr(elem_marc(n_marc(i))),Kc(elem_marc(n_marc(i))),melt_fraction(i),alpharef,Href,e2t_p(i),e2t_v(i),e2t_p_elem(elem_marc(n_marc(i))),&
            e2t_v_elem(elem_marc(n_marc(i))),melt_elem(elem_marc(n_marc(i))),W_b_n(i),w_b_elem(elem_marc(n_marc(i))))
     enddo
     DRho=RhoC-Rhoref
     if(power_law)call Viscosita(coeffA,ndisl,Va,Ea,Cohe_0,coeffAdf,Eadf,Vadf,teta_fr_0,num_mat,materiale,Rgas,elem_marc,pi,Visc_M,grain_m,grain_p,Visc_max,&
                       Visc_min,C,elementi,media,Viscref,Visc,marcatori,marcatoritot,Pmarc,Tmarc,Visc_vp,n_marc,type_marc,e2t_p,e2t_v,e2t_ist,VV,W_b_n,&
                       melt_fraction,teta_fr_inf,Cohe_inf,Ws,e0,e1,e0v,e1v,alpham_dl,alpham_df,rw,OH_0,Pmarc_lito,diffusion,rheology,t_serp,elem_x,&
                       eff_average,f_scale,uniaxial,teta_elem,cohe_elem,idratazione,t_oc,serpentine,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,melt_oc,&
                       colonna_marc)
     if(equilib)Visc=1.D21/Viscref
     do i=1,elementi
       do k=1,nodi_elem
         Vel_res(i,(k*ndof)-1)=Vx(connessioni(i,k))
         Vel_res(i,k*ndof)=Vy(connessioni(i,k))
       enddo
     enddo
   endif
   if(nproc.gt.1)call mpi_bcast(Vel,nodi*ndof,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   if(nproc.gt.1)call mpi_bcast(Ptot_nodi,nodi,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   if(nproc.gt.1)call mpi_bcast(res,nodi*ndof,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   if(nproc.gt.1)call mpi_bcast(Load_g,nodi*ndof,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
   norm_vel=DSQRT(sum(Vel**2))
   norm_dvel=DSQRT(sum((Vel-Vel_old)**2))
   norm_P=DSQRT(sum(Ptot_nodi**2))
   norm_dP=DSQRT(sum((Ptot_nodi-Ptot_old)**2))
   norm_res=DSQRT(sum(res**2))
   if(ic.eq.1)norm_r0=DSQRT(sum(Load_g**2))
   if(OUT_FILE_nl_vtu.and.iproc.eq.0)call Output_nl_vtu(OUT_VEL,OUT_DIR,it,ic,elementi,nodi_elem,nodi,connessioni,x*h0,y*h0,Visc*Viscref,Ptot_nodi,&
                                          DABS(Ptot_nodi-Ptot_old),Press_old,dPress,Vel*v0,DABS(Vel*v0-Vel_old*v0),res/norm_r0,ndof,e2_nodi,OUT_MODEL,OUT_PHASE)
   if(OUT_FILE_nl.and.iproc.eq.0)call Output_nl(OUT_VEL,OUT_DIR,it,ic,it_max,norm_res/norm_r0,norm_dvel/norm_vel,norm_dP/norm_P,&
                                      MINVAL(Ptot-Plito),MAXVAL(Ptot-Plito),MINVAL(Press_old-Plito),MAXVAL(Press_old-Plito),&
                                      MINVAL(Visc*Viscref),MAXVAL(Visc*Viscref),MINVAL(Vx*v0),MAXVAL(Vx*v0),MINVAL(Vy*v0),MAXVAL(Vy*v0),OUT_MODEL,OUT_PHASE)
   Vel_old=Vel
   Ptot_old=Ptot_nodi
   if(iproc.eq.0.and.ic.eq.10*nl_out.and.ic.lt.100.and.MAX(norm_res/norm_r0,norm_dvel/norm_vel).gt.tol)then
     print 1016,ic,MAX(norm_res/norm_r0,norm_dvel/norm_vel)
     nl_out=nl_out+1
   elseif(iproc.eq.0.and.ic.eq.10*nl_out.and.ic.ge.100.and.MAX(norm_res/norm_r0,norm_dvel/norm_vel).gt.tol)then
     print 1017,ic,MAX(norm_res/norm_r0,norm_dvel/norm_vel)
     nl_out=nl_out+1
   endif
   if(iproc.eq.0)ic=ic+1
   if(nproc.gt.1)call mpi_bcast(ic,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 enddo                                                                        !FINE CICLO NON LINEARITA'

 if(iproc.eq.0)then
   W_e=W_e_n
   W_b=W_b_n
   W_m=W_m_n
   s2=DSQRT((0.5D0*((2*Visc*exx)**2+(2*Visc*eyy)**2)+(2*Visc*exy)**2))*Viscref/t0
   if(shear)Hs=exx*(2.D0*Visc*exx)+eyy*(2.D0*Visc*eyy)+exy*(4.D0*Visc*exy)
   if(adiabatic)Hd=alpha*Rho*Vy_elem*Telem

   print *,'  Done!'
   call cpu_time(t2)
   if(power_law.and..not.equilib)then
     print 1011,ic-1
     print 1010,MAX(norm_res/norm_r0,norm_dvel/norm_vel),tol
   endif
   print 1013,t2-t1
   print 1015,(t2-t1)*nproc
   print *,''

   deallocate(DRho,exx,eyy,exy,Ray,W_e_n,dpress,press_old)
   deallocate(exx_nodi,eyy_nodi,exy_nodi)
   deallocate(Wmax_b,W_b_n,W_m_n)
   close(10)   
 endif
 deallocate(res,div,Load_g)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     8 - SCRITTURA FILES
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(iproc.eq.0)then
   if(lagrangian_grid)call Advection_lgrid(free_surface,nodi_lg,elementi,nodi,nodi_elem,connessioni,elem_x,elem_y,layer_x,layer_y,n_elem_x,n_elem_y,x,y,&
                           x_min,x_max,y_min,y_max,larghezza,altezza,thick_x,thick_y,Vx,Vy,dt0,delta_x,xmin_elem,x_lg,y_lg,em_lg,elem_lg,conn_lg,eps,&
                           elem_x_lg+1,elem_y_lg+1,dim_y_lg)
   if(it.eq.1.or.(OUT_TYPE.and.(tempo+dt/2.D0.ge.(DBLE(nt)*OUT_FILE_tm)+t_initial.or.(tempo-dt/2.D0.le.(DBLE(nt)*OUT_FILE_tm)+t_initial.and.&
     tempo.ge.(DBLE(nt)*OUT_FILE_tm)+t_initial))).or.(.not.OUT_TYPE.and.it.eq.nt).or.(tempo+dt.ge.t_final))then
     call cpu_time(t1)
     allocate(type_elem(elementi),def_elem(elementi))
     call marker_to_element(num_mat,elem_marc,n_marc,marcatori,marcatoritot,elementi,type_marc,type_elem,VV,def_elem)
     print *,'- Writing files...'
     call Output(OUT_VEL,OUT_DIR,elementi,nodi_elem,nodi,connessioni,x*h0,y*h0,x_elem*h0,y_elem*h0,Temp*(Tb-Ts),Vx,Vy,v0,Cp*Cpref,Visc*Viscref,&
          Rho*Rhoref,Kc*Kcref,Telem*(Tb-Ts),Ptot,Plito,Ptot_nodi,x_marc*h0,y_marc*h0,n_marc,type_marc,VV,marcatori,e2_nodi,s2,Pmarc,Tmarc,e2t_p,e2t_v,&
          Hr*Href*Rho*Rhoref,Hs*Viscref/(t0**2),Hd*alpharef*(Tb-Ts)*Rhoref*Gy*v0,tempo,W_b,W_f,melt_fraction,alpha*alpharef,W_m,melting,perple_x,e2t_ist,&
          Visc_vp*Viscref,C,OUT_TYPE,nt,OUT_FILE_nb,it,-q_x*Kcref*(Tb-Ts)/h0,-q_y*Kcref*(Tb-Ts)/h0,power_law,teta_elem,cohe_elem,e2t_p_elem,e2t_v_elem,&
          markers_vtu_out,markers_txt_out,type_elem,def_elem,shear,adiabatic,solver_T,melt_elem,OUT_MODEL,OUT_PHASE,w_b_elem,W_e/C)
          deallocate(type_elem,def_elem)
     if(lagrangian_grid)call Output_lgrid(it,tempo,x_lg*h0,y_lg*h0,OUT_TYPE,nodi_lg,elem_lg,conn_lg,dim_x_lg*h0,dim_y_lg*h0,em_lg)
     call cpu_time(t2)
     print 1014,t2-t1
     print *,''
   endif
   if(OUT_FILE_nl_it)then
     if(it.eq.1)then
       open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Iterations_flux.txt',status='replace')
     else
       open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Iterations_flux.txt',status='old',position='append')
     endif
     write(1,*)it,tempo,ic-1,SUM(flux_sx),SUM(flux_dx)
     close(1)
   endif
   tempo=tempo+dt
   if(.not.equilib)call Strain(e2t_ist,dt0*t0,VV,e2t_p,e2t_v,marcatori,healing,healing_rate,thermal_act,e1,e1v,Tmarc/(Tb-Ts))
   deallocate(Telem,s2,Press,Plito,Ptot,teta_elem,cohe_elem)
   deallocate(e2_nodi,Plito_nodi)
   deallocate(e2t_ist,Visc_vp,Pmarc,Pmarc_lito,Tmarc,VV)
   deallocate(flux_sx,flux_dx)
 endif
 if(nproc.gt.1)call mpi_bcast(tempo,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     9 - BILANCIO DELL'ENERGIA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 if(solver_T)then
   if(iproc.eq.0)open(unit=11,file='mumps_info_T',status='unknown')
   if(iproc.eq.0.and.tempo.le.t_final)then
     call cpu_time(t1)
     print *,'- Energy solution:'
     print *,'  Facticulating...'
     if(power_law.and..not.equilib)print *,'  Do not worry, it is much faster than momentum solution!'
     allocate(K_loc(nodi_elem*ndofT,nodi_elem*ndofT),K_Mc(nodi_elem*ndofT,nodi_elem*ndofT),K_K(nodi_elem*ndofT,nodi_elem*ndofT))
     allocate(Load(nodi_elem*ndofT),N1(1,nodi_elem),N2(nodi_elem,1),NTT(nodi_elem,1),Vel1(1,ndof),Vel2(ndof,1),Tnodi(nodi_elem))
     idT%RHS=0.D0
     idT%A_ELT=0.D0
     counter=0
     do i=1,elementi
       K_Loc=0.D0
       Load=0.D0
       Tnodi=Temp(connessioni(i,:))
       do ig=1,ngauss
         Vel1=0.D0
         Vel2=0.D0
         call Matrici(ngauss,ndof,nodi_elem,elementi,nodi,connessioni,x,y,i,ig,B,D_mu,D_l,dvolu,N,dC)
         N1(1,:)=N(:)
         N2(:,1)=N(:)
         do j=1,nodi_elem
           Load(j)=Load(j)+((Rho(i)*Hr(i))+(Hs(i)*Diss/Ra)+(Hd(i)*Diss))*N(j)*dvolu*dt0
           Vel1(1,1)=Vel1(1,1)+Vx(connessioni(i,j))*N(j)
           Vel1(1,2)=Vel1(1,2)+Vy(connessioni(i,j))*N(j)
           Vel2(1,1)=Vel2(1,1)+Vx(connessioni(i,j))*N(j)
           Vel2(2,1)=Vel2(2,1)+Vy(connessioni(i,j))*N(j)
         enddo
          Veltot=DSQRT((Vel1(1,1)**2)+(Vel1(1,2)**2))
          tau1=delta_x(i)/(2.D0*Veltot)
          tau2=teta*dt0
          tau3=((delta_x(i)**2)*(RhoC(i)/Rhoref)*Cp(i))/Kc(i)
          tau=((1.D0/(tau1**rtau))+(1.D0/(tau2**rtau))+(1.D0/(tau3**rtau)))**(-1.D0/rtau)
          NTT=N2+(tau*matmul(dC,Vel2))
          K_K=(matmul(dC,transpose(dC))*Kc(i)+matmul(NTT,matmul(Vel1,transpose(dC)))*(RhoC(i)/Rhoref)*Cp(i))*dvolu
          K_Mc=matmul(N2,N1)*Cp(i)*(RhoC(i)/Rhoref)*dvolu
          K_loc=K_loc+K_Mc+(K_K*dt0*teta)
          Load=Load+matmul(K_Mc-K_K*(1.D0-teta)*dt0,Tnodi(1:4))
       enddo
       call BC_fixed(ndofT,nodi_elem,nodi,connessioni,elementi,i,BC_fixT,BC_valT,K_loc,Load)
       do ine1=1,nodi_elem
         ilg=connessioni(i,ine1) 
         do ine2=1,nodi_elem
           counter=counter+1   
           idT%A_ELT(counter)=K_loc(ine2,ine1) 
         enddo    
         idT%RHS(ilg)=idT%RHS(ilg)+Load(ine1)    
       enddo
     enddo
     deallocate(K_loc,K_Mc,K_K)
     deallocate(Load,N1,N2,NTT,Vel1,Vel2,Tnodi)
   endif
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   idT%ICNTL(5)=1
   idT%ICNTL(7)=0
   idT%JOB=6
   call DMUMPS(idT)

   if(iproc.eq.0.and.tempo.le.t_final)then
     Temp=idT%RHS
     print *,'  Done!'
     call cpu_time(t2)
     print 1013,t2-t1
     print 1015,(t2-t1)*nproc
     print *,''
   endif
   if(iproc.eq.0)close(11)
 endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     10 - ADVEZIONE MARCATORI, SUPERFICIE LITOSFERA E FINE CICLO TEMPORALE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(iproc.eq.0.and.tempo.le.t_final)then
   if(erosione.and..not.equilib)call cpu_time(t1)
   if(erosione.and..not.equilib)print *,'- Erosion and sedimentation:'
   call Upper_Boundary(x,y,nodi,nodi_x,nodi_y,Vx,Vy,dt0,delta_x,delta_y,connessioni,nodi_elem,elementi,chain,x_chain,y_chain,&
        elem_x,elem_y,larghezza,altezza,h0,v0,eps,free_surface,erosione,thick_air,xmin_elem,ymin_elem,esp_m,esp_n,kappa_f,kappa_d,&
        x_min,x_max,y_min,y_max,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,&
        kappa_fsed,kappa_dsed,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y,sed_layer)
   if(erosione.and..not.equilib)call cpu_time(t2)
   if(erosione.and..not.equilib)print 1013,t2-t1
   if(erosione.and..not.equilib)print *,''
   call cpu_time(t1)
   print *,'- Markers advection and new composition:'
   allocate(n_marc_a(SUM(nmin)),type_marc_a(SUM(nmin)),elem_marc_a(SUM(nmin)),x_marc_a(SUM(nmin)),y_marc_a(SUM(nmin)),e2t_p_a(SUM(nmin)),e2t_v_a(SUM(nmin)),&
            W_b_a(SUM(nmin)),W_m_a(SUM(nmin)),melt_fraction_a(SUM(nmin)),W_f_a(SUM(nmin)),n_marc_tmp(marcatori),type_marc_tmp(marcatori),&
            elem_marc_tmp(marcatoritot),x_marc_tmp(marcatori),y_marc_tmp(marcatori),e2t_p_tmp(marcatori),e2t_v_tmp(marcatori),W_b_tmp(marcatori),&
            W_f_tmp(marcatori),W_m_tmp(marcatori),melt_fraction_tmp(marcatori),inside(marcatori),marc_each_elem(elementi,MAXVAL(nmax*2)),&
            colonna_marc_a(SUM(nmin)),colonna_marc_tmp(marcatoritot))
   call Advezione(elementi,elem_x,elem_y,nodi,nodi_elem,x,y,x_marc,y_marc,n_marc,marcatori,marcatoritot,elem_marc,Vx,Vy,connessioni,dt0,C,nmax,larghezza,&
        altezza,x_min,x_max,y_min,y_max,ndof,vel_corr,runge,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y,inside,marc_in,marc_each_elem,marc_agg,e2t_p,&
        e2t_v,type_marc,W_b,W_f,W_m,melt_fraction,free_surface,erosione,chain,x_chain,y_chain,t_air,t_sed,n_marc_tmp,type_marc_tmp,elem_marc_tmp,&
        x_marc_tmp,y_marc_tmp,e2t_p_tmp,e2t_v_tmp,W_b_tmp,W_f_tmp,W_m_tmp,melt_fraction_tmp,colonna_marc,colonna_marc_tmp)
   if(check_markers)call Add_remove_markers(free_surface,elementi,nmin,nmax,nodi,nodi_elem,connessioni,marc_each_elem,marcatori,marcatoritot,n_marc,type_marc,&
                         delta_x,xmin_elem,ymin_elem,ymax_elem,x,y,x_marc,y_marc,larghezza,e2t_p,e2t_v,W_b,W_f,W_m,melt_fraction,C,n_marc_a,type_marc_a,&
                         elem_marc_a,x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_f_a,W_m_a,melt_fraction_a,marc_in,marc_agg,erosione,elem_x,elem_y,t_sed,&
                         equilib,sed_layer,colonna_marc,colonna_marc_a)
   deallocate(n_marc,type_marc,elem_marc,x_marc,y_marc,e2t_p,e2t_v,W_b,W_f,W_m,melt_fraction,colonna_marc)
   allocate(n_marc(marc_in+marc_agg),type_marc(marc_in+marc_agg),elem_marc(marcatoritot+marc_agg),x_marc(marc_in+marc_agg),y_marc(marc_in+marc_agg),&
            e2t_p(marc_in+marc_agg),e2t_v(marc_in+marc_agg),W_b(marc_in+marc_agg),W_f(marc_in+marc_agg),W_m(marc_in+marc_agg),melt_fraction(marc_in+marc_agg),&
            colonna_marc(marcatoritot+marc_agg))
   call Refine_markers(x_marc,y_marc,e2t_p,e2t_v,n_marc,marcatori,marcatoritot,elem_marc,type_marc,W_b,W_m,W_f,melt_fraction,inside,&
        marc_in,marc_agg,n_marc_a,type_marc_a,elem_marc_a,x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_f_a,W_m_a,melt_fraction_a,n_marc_tmp,type_marc_tmp,&
        elem_marc_tmp,x_marc_tmp,y_marc_tmp,e2t_p_tmp,e2t_v_tmp,W_b_tmp,W_f_tmp,W_m_tmp,melt_fraction_tmp,colonna_marc_a,colonna_marc_tmp)
   deallocate(n_marc_a,type_marc_a,elem_marc_a,x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_m_a,melt_fraction_a,W_f_a,n_marc_tmp,type_marc_tmp,elem_marc_tmp,&
              x_marc_tmp,y_marc_tmp,e2t_p_tmp,e2t_v_tmp,W_b_tmp,W_m_tmp,melt_fraction_tmp,W_f_tmp,inside,marc_each_elem,colonna_marc_a,colonna_marc_tmp)
   call cpu_time(t2)
   print 1013,t2-t1
   print *,''
   print *,'****************************************************************************************'
 endif

 call mpi_barrier(MPI_COMM_WORLD,ierr)
 enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     11 - FORMATI FILE E FINE PROGRAMMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 if(iproc.eq.0)then
   print *,''
   print *,'****************************************************************************************'
   print *,''
   print *,'Simulation ended, well done!'
   print *,'Just one more minute to create final files...'
   if(final_distribution)call Output_final(nodi,elementi,nodi_elem,layer_x,layer_y,thick_x*h0,thick_y*h0,n_elem_x,n_elem_y,connessioni,Temp*(Tb-Ts),&
                              x*h0,y*h0,x_marc*h0,y_marc*h0,n_marc,marcatori,type_marc,e2t_p,e2t_v,W_b,W_f,W_m,melt_fraction,OUT_DIR,OUT_MODEL,OUT_PHASE)
   print *,'Done!'
   print *,''
   print *,'****************************************************************************************'
   print *,'----------------------------------- HAVE A NICE DAY!! ----------------------------------'
   print *,'****************************************************************************************'
   print *,''

   deallocate(x,y,Temp,Vx,Vy)
   deallocate(BC_val,BC_fix,BC_valT,BC_fixT)
   deallocate(N,dC,B)
   deallocate(x_elem,y_elem,xmax_elem,ymax_elem,xmin_elem,ymin_elem,delta_x,delta_y,C,nmin,nmax,Visc,Hr,Kc,W_e,Rho,RhoC,Cp,alpha,Hs,Hd,Vy_elem,q_x,q_y,&
              e2t_p_elem,e2t_v_elem,melt_elem,sed_layer,w_b_elem,melt_oc)
   deallocate(connessioni)
   deallocate(W_b,W_m,W_f,e2t_p,e2t_v,n_marc,x_marc,y_marc,type_marc,melt_fraction,elem_marc,colonna_marc)
   deallocate(idV%ELTPTR,idV%ELTVAR,idT%ELTPTR,idT%ELTVAR)
   if(lagrangian_grid)deallocate(x_lg,y_lg,conn_lg,em_lg)
   if(.not.free_surface)then
     deallocate(x_chain,y_chain)
   endif
   if(perple_x)then
     deallocate(RhoX,CpX,alphaX,WX,RhoXs,CpXs,alphaXs,WXs)
   endif
   deallocate(materiale,Rho_M,Cp_M,alpha_M,Visc_M,Hr_M,Kc_M,ndisl,Ea,Va,coeffA,Eadf,Vadf,coeffAdf,teta_fr_0,Cohe_0,Tref_M,teta_fr_inf,Cohe_inf,Ws,&
              Visc_min,Visc_max,rheology,f_scale,diffusion,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,uniaxial,grain_m,grain_p)
 endif
 deallocate(thick_x,thick_y,n_elem_x,n_elem_y,weak_dim_x,weak_dim_y,weak_pos_x,weak_pos_y)
 deallocate(Ptot_old,Ptot_nodi,Vel_old,Vel,Vel_res)
 deallocate(idV%A_ELT,idV%RHS,idT%A_ELT,idT%RHS)
 call mpi_finalize(ierr)

1001 format(' - Total nodes:',i8)
1002 format(' - Total elements:',i8)
1003 format(' - Width of the domain (km):',f8.3)
1004 format(' - Height of the domain (km):',f8.3)
1005 format(' - Initial markers:',i8)
1006 format(' Iteration:',i6)
1007 format(' - Time step (Ma):',f11.6)
1008 format(' - Total time (Ma):',f11.6)
1010 format('   - Non-linear convergence:',g10.3,1x,'(Tolerance:',g10.3,1x,')')
1011 format('   - Non-linear iterations performed:',i4)
1012 format(' - Markers:',i8)
1013 format('   - Calculation time (s): ',g10.3)
1014 format('   - Writing time (s): ',g10.3)
1015 format('   - Total CPU time (s): ',g10.3)
1016 format('   'i2,' non-linear iterations performed and convergence still at',g10.3,1x,', this could be long...')
1017 format('   'i3,' non-linear iterations performed and convergence still at',g10.3,1x,', this could be long...')

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     12 - SUBROUTINES
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Create_grid(layer_x,layer_y,n_elem_x,n_elem_y,nodi_x,nodi_y,nodi,elem_x,elem_y,thick_x,thick_y,x,y,file_grid,&
 file_world,thick_air,gravity)

 Implicit none

 character(len=50) :: output_dir='../../doc/manual/'//C_NULL_CHAR
 logical(1) :: has_output_dir=.false.
 character(len=:), allocatable :: file_name
 integer :: i,j,k,ii,jj,kk,indice,nodo(nodi_x,nodi_y),to_x,to_y,length
 integer*8 :: random_number_seed=1.0
 double precision :: delta_x(layer_x),delta_y(layer_y),temperature,y_min,y_max
 double precision, parameter :: depth=1.D-04
 character(len=*), intent(in) :: file_grid,file_world
 integer, intent(in) :: layer_x,layer_y,n_elem_x(layer_x),n_elem_y(layer_y),nodi_x,nodi_y,nodi,elem_x,elem_y
 double precision, intent(in) :: thick_x(layer_x),thick_y(layer_y),thick_air,gravity

 double precision, intent(out) :: x(nodi),y(nodi)

 length=LEN(TRIM(file_world))+20
 allocate(character(length) :: file_name)
 file_name=TRIM(ADJUSTL(file_world))//C_NULL_CHAR
 call create_world(cworld,file_name,has_output_dir,output_dir,random_number_seed)
 print *,''
 print *,'Creating grid file...'
 open(unit=1,file=file_grid,status='unknown') 
 write(1,1000)
 write(1,'(i10,1x,20(f15.5,i10))')layer_x,(thick_x(i),n_elem_x(i),i=1,layer_x)
 write(1,'(i10,1x,20(f15.5,i10))')layer_y,(thick_y(i),n_elem_y(i),i=1,layer_y)
 write(1,*)''
 write(1,2000)
 k=0
 do i=1,layer_y
   delta_y(i)=thick_y(i)/n_elem_y(i)
   if(i.eq.1)then
     to_y=n_elem_y(i)+1
   else
     to_y=n_elem_y(i)
   endif
   do j=1,to_y
     k=k+1
     kk=0
     do ii=1,layer_x
       delta_x(ii)=thick_x(ii)/n_elem_x(ii)
       if(ii.eq.1)then
         to_x=n_elem_x(ii)+1
       else
         to_x=n_elem_x(ii)
       endif
       do jj=1,to_x
         kk=kk+1
         indice=nodi_x*(k-1)+kk
         if(kk.eq.1)then
           x(indice)=0.D0
         else
           x(indice)=x(indice-1)+delta_x(ii)
         endif
         if(k.eq.1)then
           y(indice)=0.D0
         else
           y(indice)=y(indice-nodi_x)+delta_y(i)
         endif
         nodo(kk,k)=indice
       enddo
     enddo
   enddo
 enddo
 y_max=MAXVAL(y)
 y_min=MINVAL(y)
 do j=1,nodi_y
   do i=1,nodi_x
     indice=nodi_x*(j-1)+i
     if(y_max-thick_air-y(indice).ge.y_min)then
       call temperature_2d(cworld,x(indice),depth,y_max-thick_air-y(indice),gravity,temperature)
       write(1,'(i10,1x,2(f20.9,1x),f15.4,1x,i4)')indice,x(indice),y(indice),temperature
     else
       write(1,'(i10,1x,2(f20.9,1x),f15.4,1x,i4)')indice,x(indice),y(indice),273.15
     endif
   enddo
 enddo
 write(1,3000)
 do j=1,elem_y
   do i=1,elem_x
     indice=elem_x*(j-1)+i
     write(1,'(5(i10,1x))')indice,nodo(i,j),nodo(i+1,j),nodo(i+1,j+1),nodo(i,j+1)
    enddo
 enddo
 close(1)
 call release_world(cworld)

1000 format(5x,'Layer',6x,'Dimensione',2x,'Elementi')
2000 format(6x,'Nodo',9x,'Coordinata X',9x,'Coordinata Y',5x,'Temperatura')
3000 format(2x,'Elemento',4x,'Conn. 1',4x,'Conn. 2',4x,'Conn. 3',4x,'Conn. 4')

 end subroutine Create_grid

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Create_lgrid(dim_x_lg,dim_y_lg,x_lg,y_lg,elem_x_lg,elem_y_lg,nodi_lg,elem_lg,conn_lg,xmin_lg,ymin_lg)

 Implicit none

 integer :: i,j,nodo(elem_x_lg+1,elem_y_lg+1)
 integer, intent(in) :: elem_x_lg,elem_y_lg,nodi_lg,elem_lg
 double precision, intent(in) :: dim_x_lg,dim_y_lg,xmin_lg,ymin_lg
 integer, intent(out) :: conn_lg(elem_lg,4)
 double precision, intent(out) :: x_lg(nodi_lg),y_lg(nodi_lg)

 do j=1,elem_y_lg+1
   do i=1,elem_x_lg+1
     indice=(elem_x_lg+1)*(j-1)+i
     if(i.eq.1)then
       x_lg(indice)=xmin_lg
     else
       x_lg(indice)=x_lg(indice-1)+dim_x_lg
     endif
     if(j.eq.1)then
       y_lg(indice)=ymin_lg
     else
       y_lg(indice)=y_lg(indice-(elem_x_lg+1))+dim_y_lg
     endif
     nodo(i,j)=indice
   enddo
 enddo
 do j=1,elem_y_lg
   do i=1,elem_x_lg
     indice=elem_x_lg*(j-1)+i
     conn_lg(indice,1)=nodo(i,j)
     conn_lg(indice,2)=nodo(i+1,j)
     conn_lg(indice,3)=nodo(i+1,j+1)
     conn_lg(indice,4)=nodo(i,j+1)
    enddo
 enddo

 end subroutine Create_lgrid

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Read_lgrid(x_lg,y_lg,nodi,elementi,connessioni,lgrid)

 Implicit none

 integer :: i,j
 character(len=*), intent(in) :: lgrid
 integer, intent(in) :: nodi,elementi
 integer, intent(out) :: connessioni(elementi,4)
 double precision, intent(out) :: x_lg(nodi),y_lg(nodi)

 open(unit=1,file=lgrid,status='old')
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*)
 do i=1,nodi
   read(1,*) x_lg(i),y_lg(i)
 enddo
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*)
 do i=1,elementi
   read(1,*) (connessioni(i,j),j=1,4)
 enddo
 close(1)
 connessioni=connessioni+1

 end subroutine Read_lgrid

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine mumps_vectors(iproc,idM,nodi,elementi,nodi_elem,ndof,connessioni)

 Implicit none
 
 integer :: i,j,k,LELTVAR,NA_ELT,counter
 integer, intent(in) :: nodi,elementi,nodi_elem,ndof,iproc,connessioni(elementi,nodi_elem)
 type(dmumps_struc)idM

 idM%N=nodi*ndof
 idM%NELT=elementi
 LELTVAR=elementi*(nodi_elem*ndof)
 if(idM%SYM.eq.1)then
   NA_ELT=elementi*(nodi_elem*ndof)*(nodi_elem*ndof+1)/2
 else
   NA_ELT=elementi*(nodi_elem*ndof)*(nodi_elem*ndof)
 endif
 allocate(idM%A_ELT(NA_ELT),idM%RHS(idV%N))
 if (iproc.eq.0) then
   allocate(idM%ELTPTR(idM%NELT+1),idM%ELTVAR(LELTVAR))
   do i=1,elementi
     idM%ELTPTR(i)=1+(i-1)*(nodi_elem*ndof)
   enddo
   idM%ELTPTR(elementi+1)=1+elementi*(nodi_elem*ndof)
   counter=0
   do i=1,elementi
     do j=1,nodi_elem
       do k=1,ndof
          counter=counter+1
          idM%ELTVAR(counter)=(connessioni(i,j)-1)*ndof+k
       enddo
     enddo
   enddo
 endif

 end subroutine mumps_vectors

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Dimensionless(nodi,elementi,h0,Diff,x,y,delta_x,delta_y,x_min,x_max,y_min,y_max,var_grid,aria,lito_sx,lito_dx,transition_BC_sx,&
 transition_BC_dx,xmax_elem,xmin_elem,ymax_elem,ymin_elem,larghezza,altezza,v0,thick_x,thick_y,layer_x,layer_y,time_step,max_time_step)

 implicit none

 integer, intent(in) :: nodi,elementi,layer_x,layer_y
 double precision, intent(in) :: h0,Diff
 double precision, intent(out) :: x_min,x_max,y_min,y_max,v0
 double precision, intent(inout) :: x(nodi),y(nodi),delta_x(elementi),delta_y(elementi),xmax_elem(elementi),xmin_elem(elementi),ymax_elem(elementi),&
 ymin_elem(elementi),var_grid,aria,larghezza,altezza,thick_x(layer_x),thick_y(layer_y),time_step,lito_sx,lito_dx,transition_BC_sx,transition_BC_dx,&
 max_time_step

 x=x/h0
 y=y/h0
 delta_x=delta_x/h0
 delta_y=delta_y/h0
 x_min=x_min/h0
 x_max=x_max/h0
 y_min=y_min/h0
 y_max=y_max/h0
 var_grid=var_grid/h0
 aria=aria/h0
 thick_x=thick_x/h0
 thick_y=thick_y/h0
 lito_sx=lito_sx/h0
 lito_dx=lito_dx/h0
 transition_BC_sx=transition_BC_sx/h0
 transition_BC_dx=transition_BC_dx/h0
 xmax_elem=xmax_elem/h0
 xmin_elem=xmin_elem/h0
 ymax_elem=ymax_elem/h0
 ymin_elem=ymin_elem/h0
 larghezza=larghezza/h0
 altezza=altezza/h0
 v0=Diff/h0
 time_step=time_step*Diff/(h0**2)
 max_time_step=max_time_step*Diff/(h0**2)

 end subroutine Dimensionless

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine BC(nodi,ndof,ndofT,y_max,x,y,v0,BC_fix,BC_val,lito_sx,lito_dx,transition_BC_sx,transition_BC_dx,aria,&
 v_sx,v_dx,Temp,BC_fixT,BC_valT,left_v,right_v,top_v,bottom_v,left_T,right_T,top_T,bottom_T,nodi_X,nodi_y)

 Implicit none

 integer :: i,j,indice
 double precision :: v_bottom
 logical, intent(in) :: left_T,right_T,top_T,bottom_T
 integer, intent(in) :: nodi,ndof,ndofT,left_v,right_v,top_v,bottom_v,nodi_x,nodi_y
 double precision, intent(in) :: y_max,x(nodi),y(nodi),v0,v_sx,v_dx,Temp(nodi),lito_sx,lito_dx,transition_BC_sx,transition_BC_dx,aria
 double precision, intent(out) :: BC_val(nodi*ndof),BC_valT(nodi*ndofT)
 logical, intent(out) :: BC_fix(nodi*ndof),BC_fixT(nodi*ndofT)

 BC_fix=.false.
 BC_fixT=.false.
 do j=1,nodi_y
   do i=1,nodi_x
     indice=nodi_x*(j-1)+i
     if(i.eq.1)then
       v_bottom=-v_sx*(y_max-0.5D0*(y_max-aria-lito_sx+y_max-aria-lito_sx-transition_BC_sx))/&
                      (0.5D0*(y_max-aria-lito_sx+y_max-aria-lito_sx-transition_BC_sx))
       if(left_v.eq.2.or.left_v.eq.3)BC_fix(((indice-1)*ndof)+1)=.true.
       if(left_v.eq.3)BC_fix(((indice-1)*ndof)+2)=.true.
       if(y(indice).ge.(y_max-lito_sx-aria).and.y(indice).le.(y_max-aria))then
         BC_val(((indice-1)*ndof)+1)=v_sx/(v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       elseif(y(indice).gt.(y_max-aria))then
         BC_val(((indice-1)*ndof)+1)=0.D0
         BC_val(((indice-1)*ndof)+2)=0.D0
       elseif(y(indice).ge.(y_max-lito_sx-aria-transition_BC_sx).and.y(indice).lt.(y_max-lito_sx-aria))then
         BC_val(((indice-1)*ndof)+1)=(v_bottom+((v_sx-v_bottom)/transition_BC_sx*(y(indice)-(y_max-lito_sx-aria-transition_BC_sx))))/&
                                     (v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       else
         BC_val(((indice-1)*ndof)+1)=v_bottom/(v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       endif
       if(left_T)BC_fixT(indice)=.true.
       BC_valT(indice)=Temp(indice)
     endif
     if(i.eq.nodi_x)then
       v_bottom=-v_dx*(y_max-0.5D0*(y_max-aria-lito_dx+y_max-aria-lito_dx-transition_BC_dx))/&
                      (0.5D0*(y_max-aria-lito_dx+y_max-aria-lito_dx-transition_BC_dx))
       if(right_v.eq.2.or.right_v.eq.3)BC_fix(((indice-1)*ndof)+1)=.true.
       if(right_v.eq.3)BC_fix(((indice-1)*ndof)+2)=.true.
       if(y(indice).ge.(y_max-lito_dx-aria).and.y(indice).le.(y_max-aria))then
         BC_val(((indice-1)*ndof)+1)=v_dx/(v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       elseif(y(indice).gt.(y_max-aria))then
         BC_val(((indice-1)*ndof)+1)=0.D0
         BC_val(((indice-1)*ndof)+2)=0.D0
       elseif(y(indice).ge.(y_max-lito_dx-aria-transition_BC_dx).and.y(indice).lt.(y_max-lito_dx-aria))then
         BC_val(((indice-1)*ndof)+1)=(v_bottom+((v_dx-v_bottom)/transition_BC_dx*(y(indice)-(y_max-lito_dx-aria-transition_BC_dx))))/&
                                     (v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       else
         BC_val(((indice-1)*ndof)+1)=v_bottom/(v0*1.D02*365.D0*24.D0*3600.D0)
         BC_val(((indice-1)*ndof)+2)=0.D0
       endif
       if(right_T)BC_fixT(indice)=.true.
       BC_valT(indice)=Temp(indice)
     endif
     if(j.eq.1)then
       if(bottom_v.eq.3)then
         BC_fix(((indice-1)*ndof)+1)=.true.
         BC_val(((indice-1)*ndof)+1)=0.D0
       endif
       if(bottom_v.eq.2.or.bottom_v.eq.3)then
         BC_fix(((indice-1)*ndof)+2)=.true.
         BC_val(((indice-1)*ndof)+2)=0.D0
       endif
       if(bottom_T)BC_fixT(indice)=.true.
       BC_valT(indice)=Temp(indice)
     endif
     if(j.eq.nodi_y)then
       if(top_v.eq.3)then
         BC_fix(((indice-1)*ndof)+1)=.true.
         BC_val(((indice-1)*ndof)+1)=0.D0
       endif
       if(top_v.eq.2.or.top_v.eq.3)then
         BC_fix(((indice-1)*ndof)+2)=.true.
         BC_val(((indice-1)*ndof)+2)=0.D0
       endif
       if(top_T)BC_fixT(indice)=.true.
       BC_valT(indice)=Temp(indice)
     endif
     if(y(indice).ge.(y_max)-aria)then
       BC_fixT(indice)=.true.
       BC_valT(indice)=Temp(indice)
     endif
   enddo
 enddo

 end subroutine BC

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Create_markers(file_marc,file_world,marc_in_elem,elementi,MARKERS_DISTRIBUTION,init_num_mat,xmin_elem,ymin_elem,&
 deltax,deltay,y_min,y_max,thick_air)

 Implicit none

 character(len=50) :: output_dir='../../doc/manual/'//C_NULL_CHAR
 logical(1) :: has_output_dir=.false.
 integer :: i,j,k,ii,jj,nx,length
 integer*8 :: random_number_seed=1.0
 double precision :: xm,ym,composition
 character(len=:), allocatable :: file_name
 double precision, dimension(:), allocatable :: xmark,ymark
 double precision, parameter :: depth=1.D-04
 character(len=*), intent(in) :: file_marc,file_world
 integer, intent(in) :: marc_in_elem(elementi),elementi,MARKERS_DISTRIBUTION,init_num_mat
 double precision, intent(in) :: xmin_elem(elementi),ymin_elem(elementi),deltax(elementi),deltay(elementi),y_min,y_max,thick_air

 length=LEN(TRIM(file_world))+20
 allocate(character(length) :: file_name)
 file_name=TRIM(ADJUSTL(file_world))//C_NULL_CHAR
 call create_world(cworld,file_name,has_output_dir,output_dir,random_number_seed)
 print *,''
 print *,'Creating markers file...'
 open(unit=1,file=file_marc,status='unknown')
 write(1,*)'Marcatori totali'
 write(1,*)SUM(marc_in_elem)
 write(1,1000)
 k=0
 do i=1,elementi
   if(MARKERS_DISTRIBUTION.eq.1)then
     nx=INT(SQRT(DBLE(marc_in_elem(i))))
     allocate(xmark(nx),ymark(nx))
     do j=1,nx
       do jj=1,nx
         if(j.eq.1)then
           xmark(j)=xmin_elem(i)+deltax(i)/(nx*2)
         else
           xmark(j)=xmark(j-1)+deltax(i)/nx
         endif
         if(jj.eq.1)then
           ymark(jj)=ymin_elem(i)+deltay(i)/(nx*2)
         else
           ymark(jj)=ymark(jj-1)+deltay(i)/nx
         endif
         xm=xmark(j)
         ym=ymark(jj)
         k=k+1
         if(y_max-thick_air-ym.ge.y_min)then
           do ii=1,init_num_mat
             call composition_2d(cworld,xm,depth,y_max-thick_air-ym,ii,composition)
             if(composition.eq.1.D0)then
               write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,ii,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
               exit
             endif
           enddo
         else
           write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,init_num_mat+1,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
         endif
       enddo
     enddo
     deallocate(xmark,ymark)
   elseif(MARKERS_DISTRIBUTION.eq.2)then
     do j=1,marc_in_elem(i)
       call random_number(xm)
       call random_number(ym)
       xm=(xm*deltax(i))+xmin_elem(i)
       ym=(ym*deltay(i))+ymin_elem(i)
       k=k+1
       if(y_max-thick_air-ym.ge.y_min)then
         do ii=1,init_num_mat
           call composition_2d(cworld,xm,depth,y_max-thick_air-ym,ii,composition)
           if(composition.eq.1.D0)then
             write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,ii,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
             exit
           endif
         enddo
       else
         write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,init_num_mat+1,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
       endif
     enddo
   elseif(MARKERS_DISTRIBUTION.eq.3)then
     nx=INT(SQRT(DBLE(marc_in_elem(i))))
     allocate(xmark(nx),ymark(nx))
     do j=1,nx
       do jj=1,nx
         if(j.eq.1)then
           xmark(j)=xmin_elem(i)+deltax(i)/(nx*2)
         else
           xmark(j)=xmark(j-1)+deltax(i)/nx
         endif
         if(jj.eq.1)then
           ymark(jj)=ymin_elem(i)+deltay(i)/(nx*2)
         else
           ymark(jj)=ymark(jj-1)+deltay(i)/nx
         endif
         call random_number(xm)
         call random_number(ym)
         xm=xmark(j)+((xm-0.5D0)*(deltax(i)/(nx*2)))
         ym=ymark(jj)+((ym-0.5D0)*(deltay(i)/(nx*2)))
         k=k+1
         if(y_max-thick_air-ym.ge.y_min)then
           do ii=1,init_num_mat
             call composition_2d(cworld,xm,depth,y_max-thick_air-ym,ii,composition)
             if(composition.eq.1.D0)then
               write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,ii,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
               exit
             endif
           enddo
         else
           write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')k,xm,ym,init_num_mat+1,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
         endif
       enddo
     enddo
     deallocate(xmark,ymark)
   endif
 enddo
 close(1)
 call release_world(cworld)

1000 format(1x,'Marcatore',8x,'Coordinata X',9x,'Coordinata Y',7x,'Type',2x,'Plastic strain',2x,'Viscous strain',5x,'Bound water',&
     6x,'Free water',6x,'Melt water',3x,'Melt fraction')

 end subroutine Create_Markers

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,xm,ym,em,x_min,x_max,y_min,y_max,&
 larghezza,altezza,free_surface,colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)

 Implicit none

 integer :: i,tot
 double precision :: ymin,ymax,thick_tot
 logical, intent(in) :: free_surface
 integer, intent(in) :: nodi_elem,elem_x,elem_y,elementi,connessioni(elementi,nodi_elem),nodi,layer_x,layer_y,&
 n_elem_x(layer_x),n_elem_y(layer_y)
 double precision, intent(in) :: larghezza,altezza,x(nodi),y(nodi),x_min,x_max,y_min,y_max,xm,ym,thick_x(layer_x),thick_y(layer_y)
 integer, intent(out) :: em,colonna,riga

 em=0
 if((xm.ge.x_min).and.(xm.le.x_max))then
   if(layer_x.eq.1)then
     colonna=INT((xm/larghezza)*elem_x)+1
   else
     tot=0
     thick_tot=0.D0
     do i=1,layer_x
       if(xm.ge.thick_tot.and.xm.lt.thick_tot+thick_x(i))then
        colonna=INT(((xm-thick_tot)/thick_x(i))*n_elem_x(i))+1+tot
        exit
       endif
       tot=tot+n_elem_x(i)
       thick_tot=thick_tot+thick_x(i)
     enddo
   endif
   if(colonna.gt.elem_x)colonna=elem_x
   if(free_surface)then
     do i=1,elem_y
       ymin=(y(connessioni(elem_x*(i-1)+colonna,2))-y(connessioni(elem_x*(i-1)+colonna,1)))*&
            ((xm-x(connessioni(elem_x*(i-1)+colonna,1)))/(x(connessioni(elem_x*(i-1)+colonna,2))-&
            x(connessioni(elem_x*(i-1)+colonna,1))))+y(connessioni(elem_x*(i-1)+colonna,1))
       ymax=(y(connessioni(elem_x*(i-1)+colonna,4))-y(connessioni(elem_x*(i-1)+colonna,3)))*&
            ((xm-x(connessioni(elem_x*(i-1)+colonna,3)))/(x(connessioni(elem_x*(i-1)+colonna,4))-&
            x(connessioni(elem_x*(i-1)+colonna,3))))+y(connessioni(elem_x*(i-1)+colonna,3))
       if(ym.ge.ymin.and.ym.le.ymax)then
         em=elem_x*(i-1)+colonna
         exit
       endif
     enddo
   elseif(.not.free_surface.and.layer_y.eq.1)then
     if((ym.ge.y_min).and.(ym.le.y_max))then
       riga=INT((ym/altezza)*elem_y)+1
       if(riga.gt.elem_y)riga=elem_y
       em=elem_x*(riga-1)+colonna
     endif
   elseif(.not.free_surface.and.layer_y.gt.1)then
     if((ym.ge.y_min).and.(ym.le.y_max))then
       tot=0
       thick_tot=0.D0
       do i=1,layer_y
         if(ym.ge.thick_tot.and.ym.lt.thick_tot+thick_y(i))then
          riga=INT(((ym-thick_tot)/thick_y(i))*n_elem_y(i))+1+tot
          exit
         endif
         tot=tot+n_elem_y(i)
         thick_tot=thick_tot+thick_y(i)
       enddo
       if(riga.gt.elem_y)riga=elem_y
       em=elem_x*(riga-1)+colonna
     endif
   endif
 endif

 end subroutine Composizione

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine File_Perple_X(celle,celleP,dTX,dPX,dPXs,minTX,minPX,RhoXs,CpXs,alphaXs,WXs,RhoX,CpX,alphaX,WX,DIR_Perplex,&
 file_perplexM,file_perplexUC,file_perplexLC,file_perplexO,file_perplexS,file_perplexSp)

 Implicit none

 integer :: i,j
 double precision, dimension(:,:), allocatable :: TX,PX,TXs,PXs
 character(len=*), intent(in) :: file_perplexM,file_perplexUC,file_perplexLC,file_perplexO,file_perplexS,file_perplexSp,DIR_Perplex
 integer, intent(out) :: celle,celleP
 double precision, intent(out) :: dTX,dPX,dPXs,minTX,minPX
 double precision, dimension(:,:), allocatable, intent(out) :: RhoXs,CpXs,alphaXs,WXs
 double precision, dimension(:,:,:), allocatable, intent(out) :: RhoX,CpX,alphaX,WX

 print *,''
 print *,'Reading Perple_X files...'
 open(unit=1,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexM),status='old')
 read(1,*)celle
 open(unit=2,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexUC),status='old')
 open(unit=3,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexLC),status='old')
 open(unit=4,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexO),status='old')
 open(unit=5,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexS),status='old')
 celleP=((celle-1)*2)+celle
 allocate(TX(celle,celleP),PX(celle,celleP),RhoX(celle,celleP,5),CpX(celle,celleP,5),alphaX(celle,celleP,5),WX(celle,celleP,5))
 do i=1,celleP
   do j=1,celle
     read(1,*)TX(j,i),PX(j,i),RhoX(j,i,1),alphaX(j,i,1),CpX(j,i,1)
     read(2,*)TX(j,i),PX(j,i),RhoX(j,i,2),alphaX(j,i,2),CpX(j,i,2),WX(j,i,2)
     read(3,*)TX(j,i),PX(j,i),RhoX(j,i,3),alphaX(j,i,3),CpX(j,i,3),WX(j,i,3)
     read(4,*)TX(j,i),PX(j,i),RhoX(j,i,4),alphaX(j,i,4),CpX(j,i,4),WX(j,i,4)
     read(5,*)TX(j,i),PX(j,i),RhoX(j,i,5),alphaX(j,i,5),CpX(j,i,5),WX(j,i,5)
   enddo
 enddo
 close(1)
 close(2)
 close(3)
 close(4)
 close(5)
 open(unit=1,file=TRIM(DIR_Perplex)//'/'//TRIM(file_perplexSp),status='old')
 allocate(TXs(celle,celle),PXs(celle,celle),RhoXs(celle,celle),CpXs(celle,celle),alphaXs(celle,celle),WXs(celle,celle))
 do i=1,celle
  do j=1,celle
    read(1,*)TXs(j,i),PXs(j,i),RhoXs(j,i),alphaXs(j,i),CpXs(j,i),WXs(j,i)
  enddo
 enddo
 close(1)
 PX=PX*1.D05
 PXs=PXs*1.D05
 dTX=MAXVAL(TX)-MINVAL(TX)
 dPX=MAXVAL(PX)-MINVAL(PX)
 dPXs=MAXVAL(PXs)-MINVAL(PXs)
 minTX=MINVAL(TX)
 minPX=MINVAL(PX)

 end subroutine File_Perple_X

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Matrici(ngauss,ndof,nodi_elem,elementi,nodi,connessioni,x,y,ielem,ig,B,D_mu,D_l,dvolu,N,dC)

 Implicit none

 integer :: i
 double precision :: posgp(ndof,ngauss),weigp(ngauss),dNst(nodi_elem,ndof),jacob_inv(ndof,ndof),detJ,jacob(ndof,ndof),coord(ndof,nodi_elem)
 integer, intent(in) :: ig,ngauss,ndof,nodi_elem,ielem,elementi,connessioni(elementi,nodi_elem),nodi
 double precision, intent(in) :: x(nodi),y(nodi)
 double precision, intent(out) :: N(nodi_elem),D_mu(3,3),D_l(3,3),dC(nodi_elem,ndof),dvolu,B(ndof+1,nodi_elem*ndof)

 if(ngauss.eq.4)then
   posgp=RESHAPE((/-(dsqrt(3.D0)/3.D0),-(dsqrt(3.D0)/3.D0),-(dsqrt(3.D0)/3.D0),(dsqrt(3.D0)/3.D0),&
                (dsqrt(3.D0)/3.D0),-(dsqrt(3.D0)/3.D0),(dsqrt(3.D0)/3.D0),(dsqrt(3.D0)/3.D0)/),SHAPE=(/2,4/))
   weigp=(/1.D0,1.D0,1.D0,1.D0/)
 elseif(ngauss.eq.1)then
   posgp=0.D0
   weigp=4.D0
 endif

 N=(/0.25D0*(1.D0-posgp(1,ig))*(1.D0-posgp(2,ig)),0.25D0*(1.D0+posgp(1,ig))*(1.D0-posgp(2,ig)),&
     0.25D0*(1.D0+posgp(1,ig))*(1.D0+posgp(2,ig)),0.25D0*(1.D0-posgp(1,ig))*(1.D0+posgp(2,ig))/)

 dNst=RESHAPE((/-0.25D0*(1.D0-posgp(2,ig)),0.25D0*(1.D0-posgp(2,ig)),0.25D0*(1.D0+posgp(2,ig)),&
                -0.25D0*(1.D0+posgp(2,ig)),-0.25D0*(1.D0-posgp(1,ig)),-0.25D0*(1.D0+posgp(1,ig)),&
                 0.25D0*(1.D0+posgp(1,ig)),0.25D0*(1.D0-posgp(1,ig))/),SHAPE=(/4,2/))
                 
 do i=1,nodi_elem
   coord(1,i)=x(connessioni(ielem,i))
   coord(2,i)=y(connessioni(ielem,i))
 enddo

 jacob=matmul(coord,dNst)
 detJ=jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1)

 if(detJ.le.0.0) then
   print *,'DETERMINANTE JACOBIANO NON POSITIVO PER ELEMENTO',ielem
   stop
 endif

 jacob_inv=RESHAPE((/jacob(2,2)/detJ,-jacob(2,1)/detJ,-jacob(1,2)/detJ,jacob(1,1)/detJ/),SHAPE=(/2,2/))
 dvolu=detJ*weigp(ig)
 dC=matmul(dNst,jacob_inv)

! D_mu=RESHAPE((/2.D0,0.D0,0.D0,0.D0,2.D0,0.D0,0.D0,0.D0,1.D0/),SHAPE=(/3,3/))
 D_mu=RESHAPE((/4.D0/3.D0,-2.D0/3.D0,0.D0,-2.D0/3.D0,4.D0/3.D0,0.D0,0.D0,0.D0,1.D0/),SHAPE=(/3,3/))
 D_l=RESHAPE((/1.D0,1.D0,0.D0,1.D0,1.D0,0.D0,0.D0,0.D0,0.D0/),SHAPE=(/3,3/))

 B=0.D0
 do i=1,nodi_elem*ndof,2
   B(1,i)=dC((i+1)/2,1)
 enddo
 do i=2,nodi_elem*ndof,2
   B(2,i)=dC((i+1)/2,2)
 enddo
 do i=1,nodi_elem*ndof,2
   B(3,i)=dC((i+1)/2,2)
 enddo
 do i=2,nodi_elem*ndof,2
   B(3,i)=dC((i+1)/2,1)
 enddo

 end subroutine Matrici

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Stabilization(i,indice,dt,RhoC,K_loc,elementi,elem_y,x,y,connessioni,nodi_elem,nodi,ndof)

 Implicit none

 double precision :: dS,nx,ny,norm,Tcoeff
 integer, intent(in) :: i,indice,elementi,elem_y,connessioni(elementi,nodi_elem),nodi_elem,nodi,ndof
 double precision, intent(in) :: dt,RhoC(elementi),x(nodi),y(nodi)
 double precision, intent(inout) :: K_loc(nodi_elem*ndof,nodi_elem*ndof)

 if(i.eq.elem_y)then
   dS=DSQRT(((x(connessioni(indice,3))-x(connessioni(indice,4)))**2)+((y(connessioni(indice,3))-y(connessioni(indice,4)))**2))
   nx=y(connessioni(indice,4))-y(connessioni(indice,3))
   ny=x(connessioni(indice,3))-x(connessioni(indice,4))
   norm=DSQRT((nx**2)+(ny**2))
   ny=ny/norm
   Tcoeff=dt*RhoC(indice)*ny*dS*(-1.D0)
   K_loc(6,6)=K_loc(6,6)-Tcoeff
   K_loc(8,8)=K_loc(8,8)-Tcoeff
 endif

 end subroutine Stabilization

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine BC_fixed(ndof,nodi_elem,nodi,connessioni,elementi,ielem,BC_fix,BC_val,K_loc,Load)

 Implicit none

 integer :: i,j,k,ii
 double precision :: Kref
 integer, intent(in) :: ndof,nodi_elem,nodi,connessioni(elementi,nodi_elem),elementi,ielem
 double precision, intent(in) :: BC_val(ndof*nodi)
 logical, intent(in) :: BC_fix(ndof*nodi)
 double precision, intent(inout) :: Load(ndof*nodi_elem),K_loc(ndof*nodi_elem,ndof*nodi_elem)

 do i=1,nodi_elem
   ii=connessioni(ielem,i)
   do k=1,ndof
     if(BC_fix((ii-1)*ndof+k))then
       Kref=K_loc((i-1)*ndof+k,(i-1)*ndof+k)
       do j=1,ndof*nodi_elem
         Load(j)=Load(j)-K_loc(j,(i-1)*ndof+k)*BC_val((ii-1)*ndof+k)
         K_loc((i-1)*ndof+k,j)=0.D0
         K_loc(j,(i-1)*ndof+k)=0.D0
       enddo
       K_loc((i-1)*ndof+k,(i-1)*ndof+k)=Kref
       Load((i-1)*ndof+k)=Kref*BC_val((ii-1)*ndof+k)
     endif
   enddo
 enddo

 end subroutine BC_fixed

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Elemental_Pressure(elem_x,elem_y,delta_y,connessioni,elementi,nodi_elem,ngaussP,ndof,nodi,x,y,Vx,Vy,&
 flux_sx,flux_dx,exx,eyy,exy,Vy_elem,div,Press,h0,t0,Viscref,PF,Kc,Temp,q_x,q_y)

 implicit none

 integer :: i,j,k,indice,ig
 double precision :: deltaP
 integer, intent(in) :: elem_x,elem_y,connessioni(elementi,nodi_elem),elementi,nodi_elem,ngaussP,ndof,nodi
 double precision, intent(in) :: h0,t0,Viscref,PF,x(nodi),y(nodi),Vx(nodi),Vy(nodi),delta_y(elementi),Kc(elementi),Temp(nodi)
 double precision, intent(out) :: flux_sx(elem_y),flux_dx(elem_y),exx(elementi),eyy(elementi),exy(elementi),Vy_elem(elementi),div(elementi),q_x(elementi),&
 q_y(elementi)

 double precision, intent(inout) :: Press(elementi)

 do j=1,elem_y
   do i=1,elem_x
     indice=elem_x*(j-1)+i
     if(i.eq.1)then
       flux_sx(j)=((Vx(connessioni(indice,1))+Vx(connessioni(indice,4)))*v0/2.D0)*delta_y(indice)*h0
     endif
     if(i.eq.elem_x)then
       flux_dx(j)=((Vx(connessioni(indice,2))+Vx(connessioni(indice,3)))*v0/2.D0)*delta_y(indice)*h0
     endif
     exx(indice)=0.D0
     eyy(indice)=0.D0
     exy(indice)=0.D0
     Vy_elem(indice)=0.D0
     q_x(indice)=0.D0
     q_y(indice)=0.D0
     ig=1
     call Matrici(ngaussP,ndof,nodi_elem,elementi,nodi,connessioni,x,y,indice,ig,B,D_mu,D_l,dvolu,N,dC)
     do k=1,nodi_elem
       exx(indice)=exx(indice)+Vx(connessioni(indice,k))*dC(k,1)
       eyy(indice)=eyy(indice)+Vy(connessioni(indice,k))*dC(k,2)
       exy(indice)=exy(indice)+0.5D0*(Vx(connessioni(indice,k))*dC(k,2)+Vy(connessioni(indice,k))*dC(k,1))
       Vy_elem(indice)=Vy_elem(indice)+Vy(connessioni(indice,k))*N(k)
       q_x(indice)=q_x(indice)+Kc(indice)*Temp(connessioni(indice,k))*dC(k,1)
       q_y(indice)=q_y(indice)+Kc(indice)*Temp(connessioni(indice,k))*dC(k,2)
     enddo
     div(indice)=(exx(indice)+eyy(indice))/t0
     deltaP=-PF*Visc(indice)*div(indice)*Viscref
     Press(indice)=Press(indice)+deltaP
   enddo
 enddo

end subroutine Elemental_pressure

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Smooth_Pressure(elementi,elem_x,elem_y,nodi_elem,nodi,connessioni,chain,x_chain,y_chain,&
 x,y,delta_x,delta_y,larghezza,Press,free_surface,N,exx,eyy,exy,Plito,Rho,Gy,thick_air,Press_0)

 Implicit none

 integer :: i,j,k,indice,colonna
 double precision :: Pnodi(nodi),P_dx,ymin,counter(nodi),Plito_e
 integer, intent(in) :: chain,elementi,elem_x,elem_y,nodi_elem,nodi,connessioni(elementi,nodi_elem)
 double precision, intent(in) :: N(nodi_elem),x_chain(chain),y_chain(chain),x(nodi),y(nodi),delta_x(elementi),delta_y(elementi),&
 larghezza,Rho(elementi),Gy,thick_air,Press_0
 logical, intent(in) :: free_surface
 double precision, intent(out) :: Plito(elementi)
 double precision, intent(inout) :: Press(elementi),exx(elementi),eyy(elementi),exy(elementi)

  Pnodi=0.D0
  counter=0.D0
  do i=1,elementi
    do j=1,nodi_elem
      Pnodi(connessioni(i,j))=Pnodi(connessioni(i,j))+Press(i)*delta_x(i)*delta_y(i)
      counter(connessioni(i,j))=counter(connessioni(i,j))+delta_x(i)*delta_y(i)
    enddo
  enddo
  Pnodi=Pnodi/counter
  Press=0.D0
  P_dx=0.D0
  do i=1,elem_x
    do j=elem_y,1,-1
      indice=elem_x*(j-1)+i
      do k=1,nodi_elem
        Press(indice)=Press(indice)+Pnodi(connessioni(indice,k))*N(k)
      enddo
      if(.not.free_surface.and.thick_air.ne.0.D0)then
        colonna=INT((x(connessioni(indice,1))/larghezza)*(chain-1))+1
        if(colonna.gt.chain-1)colonna=chain-1
        ymin=(y_chain(colonna+1)-y_chain(colonna))*((x(connessioni(indice,1))-x_chain(colonna))/(x_chain(colonna+1)-&
             x_chain(colonna)))+y_chain(colonna)
        if(y(connessioni(indice,1)).le.ymin.and.y(connessioni(indice,4)).gt.ymin)P_dx=P_dx+(Press(indice)*delta_x(indice))
        if(y(connessioni(indice,1)).gt.ymin)then
          exx(indice)=MINVAL(ABS(exx))
          eyy(indice)=MINVAL(ABS(eyy))
          exy(indice)=MINVAL(ABS(exy))
          Plito(indice)=Press_0
          Plito_e=Press_0
        elseif(y(connessioni(indice,1)).le.ymin)then
          Plito(indice)=Plito(indice+elem_x)+Plito_e+DABS(Rho(indice)*Gy*delta_y(indice))/2.D0
          Plito_e=DABS(Rho(indice)*Gy*delta_y(indice))/2.D0
        endif
      else
        if(j.eq.elem_y)then
          if(.not.free_surface)P_dx=P_dx+(Press(indice)*delta_x(indice))
          Plito(indice)=Press_0+DABS(Rho(indice)*Gy*delta_y(indice))/2.D0
          Plito_e=Plito(indice)-Press_0
        else
          Plito(indice)=Plito(indice+elem_x)+Plito_e+DABS(Rho(indice)*Gy*delta_y(indice))/2.D0
          Plito_e=DABS(Rho(indice)*Gy*delta_y(indice))/2.D0
        endif
      endif
    enddo
  enddo
  P_dx=P_dx/(larghezza)
  Press=Press-P_dx

 end subroutine Smooth_Pressure

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 Subroutine Nodal_Pressure(nodi_elem,delta_x,delta_y,connessioni,elementi,nodi,Ptot_nodi,Ptot,Plito,Plito_nodi)

 Implicit none

 integer :: i,j
 double precision :: counter(nodi)
 integer, intent(in) :: nodi,nodi_elem,connessioni(elementi,nodi_elem),elementi
 double precision, intent(in) :: delta_x(elementi),delta_y(elementi)
 double precision, intent(out) :: Ptot_nodi(nodi),Plito_nodi(nodi)
 double precision, intent(inout) :: Ptot(elementi),Plito(elementi)

 Ptot_nodi=0.D0
 Plito_nodi=0.D0
 counter=0.D0
 do i=1,elementi
   if(Ptot(i).lt.0.D0)Ptot(i)=0.D0
   if(Plito(i).lt.0.D0)Plito(i)=0.D0
   do j=1,nodi_elem
     Ptot_nodi(connessioni(i,j))=Ptot_nodi(connessioni(i,j))+Ptot(i)*delta_x(i)*delta_y(i)
     Plito_nodi(connessioni(i,j))=Plito_nodi(connessioni(i,j))+Plito(i)*delta_x(i)*delta_y(i)
     counter(connessioni(i,j))=counter(connessioni(i,j))+delta_x(i)*delta_y(i)
   enddo
 enddo
 Ptot_nodi=Ptot_nodi/counter
 Plito_nodi=Plito_nodi/counter

 end subroutine Nodal_Pressure

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 Subroutine Nodal_Strain(elementi,nodi,nodi_elem,connessioni,ndime,xmin_elem,delta_x,t0,x,y,Vx,Vy,exx_nodi,eyy_nodi,exy_nodi,e2_nodi)

 Implicit none

 integer :: i,j,counter(nodi)
 double precision :: coord(ndime,nodi_elem),posgp(ndime),dNst(nodi_elem,ndime),jacob(ndime,ndime),detJ,jacob_inv(ndime,ndime),dC(nodi_elem,ndime)
 integer, intent(in) :: nodi,nodi_elem,ndime,connessioni(elementi,nodi_elem),elementi
 double precision, intent(in) :: t0,xmin_elem(elementi),delta_x(elementi),x(nodi),y(nodi),Vx(nodi),Vy(nodi)
 double precision, intent(out) :: exx_nodi(nodi),eyy_nodi(nodi),exy_nodi(nodi),e2_nodi(nodi)

 exx_nodi=0.D0
 eyy_nodi=0.D0
 exy_nodi=0.D0
 counter=0
 do i=1,elementi
   do j=1,nodi_elem
     coord(1,j)=x(connessioni(i,j))
     coord(2,j)=y(connessioni(i,j))
   enddo
   do j=1,nodi_elem
     posgp(1)=(2.D0/delta_x(i))*(x(connessioni(i,j))-xmin_elem(i))-1.D0
     posgp(2)=(4.D0*y(connessioni(i,j))-((1.D0-posgp(1))*y(connessioni(i,1))+(1.D0+posgp(1))*y(connessioni(i,2))+&
              (1.D0+posgp(1))*y(connessioni(i,3))+(1.D0-posgp(1))*y(connessioni(i,4))))/&
              (-(1.D0-posgp(1))*y(connessioni(i,1))-(1.D0+posgp(1))*y(connessioni(i,2))+&
              (1.D0+posgp(1))*y(connessioni(i,3))+(1.D0-posgp(1))*y(connessioni(i,4)))
     dNst=RESHAPE((/-0.25D0*(1.D0-posgp(2)),0.25D0*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(2)),&
                    -0.25D0*(1.D0+posgp(2)),-0.25D0*(1.D0-posgp(1)),-0.25D0*(1.D0+posgp(1)),&
                     0.25D0*(1.D0+posgp(1)),0.25D0*(1.D0-posgp(1))/),SHAPE=(/4,2/))

     jacob=matmul(coord,dNst)
     detJ=jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1)
     jacob_inv=RESHAPE((/jacob(2,2)/detJ,-jacob(2,1)/detJ,-jacob(1,2)/detJ,jacob(1,1)/detJ/),SHAPE=(/2,2/))
     dC=matmul(dNst,jacob_inv)
     do k=1,nodi_elem
       exx_nodi(connessioni(i,j))=exx_nodi(connessioni(i,j))+Vx(connessioni(i,k))*dC(k,1)/t0
       eyy_nodi(connessioni(i,j))=eyy_nodi(connessioni(i,j))+Vy(connessioni(i,k))*dC(k,2)/t0
       exy_nodi(connessioni(i,j))=exy_nodi(connessioni(i,j))+0.5D0*(Vx(connessioni(i,k))*dC(k,2)+Vy(connessioni(i,k))*dC(k,1))/t0
     enddo
     counter(connessioni(i,j))=counter(connessioni(i,j))+1
   enddo
 enddo
 exx_nodi=exx_nodi/DBLE(counter)
 eyy_nodi=eyy_nodi/DBLE(counter)
 exy_nodi=exy_nodi/DBLE(counter)
 e2_nodi=DSQRT(0.5D0*(exx_nodi**2+eyy_nodi**2)+exy_nodi**2)

 end subroutine Nodal_Strain

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine PT(elem_marc,n_marc,marcatori,marcatoritot,connessioni,elementi,nodi,delta_x,xmin_elem,x_marc,y_marc,y,&
 nodi_elem,Temp,Ptot,exx_nodi,eyy_nodi,exy_nodi,Tb,Ts,ndof,e2t_ist,Tmarc,Pmarc,Plito,Pmarc_lito,num_mat,materiale,type_marc)

 Implicit none

 integer :: i,j
 double precision :: posgp(ndof),N(nodi_elem),exxmarc,eyymarc,exymarc

 integer, intent(in) :: marcatori,marcatoritot,elementi,nodi,connessioni(elementi,nodi_elem),elem_marc(marcatoritot),&
 n_marc(marcatori),nodi_elem,ndof,num_mat,type_marc(marcatori)
 
 double precision, intent(in) :: delta_x(elementi),xmin_elem(elementi),x_marc(marcatori),y_marc(marcatori),y(nodi),Temp(nodi),Ptot(nodi),&
 exx_nodi(nodi),eyy_nodi(nodi),exy_nodi(nodi),Tb,Ts,Plito(nodi)

 character(len=*), intent(in) :: materiale(num_mat)

 double precision, intent(out) :: Tmarc(marcatori),Pmarc(marcatori),e2t_ist(marcatori),Pmarc_lito(marcatori)

 do i=1,marcatori
   posgp(1)=(2.D0/delta_x(elem_marc(n_marc(i))))*(x_marc(i)-xmin_elem(elem_marc(n_marc(i))))-1.D0
   posgp(2)=(4*y_marc(i)-((1.D0-posgp(1))*y(connessioni(elem_marc(n_marc(i)),1))+(1.D0+posgp(1))*y(connessioni(elem_marc(n_marc(i)),2))+&
            (1.D0+posgp(1))*y(connessioni(elem_marc(n_marc(i)),3))+(1.D0-posgp(1))*y(connessioni(elem_marc(n_marc(i)),4))))/&
            (-(1.D0-posgp(1))*y(connessioni(elem_marc(n_marc(i)),1))-(1.D0+posgp(1))*y(connessioni(elem_marc(n_marc(i)),2))+&
            (1.D0+posgp(1))*y(connessioni(elem_marc(n_marc(i)),3))+(1.D0-posgp(1))*y(connessioni(elem_marc(n_marc(i)),4)))
   N=(/0.25D0*(1.D0-posgp(1))*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(1))*(1.D0-posgp(2)),&
       0.25D0*(1.D0+posgp(1))*(1.D0+posgp(2)),0.25D0*(1.D0-posgp(1))*(1.D0+posgp(2))/)
   Tmarc(i)=0.D0
   Pmarc(i)=0.D0
   Pmarc_lito(i)=0.D0
   exxmarc=0.D0
   eyymarc=0.D0
   exymarc=0.D0
   do j=1,nodi_elem
     Tmarc(i)=Tmarc(i)+Temp(connessioni(elem_marc(n_marc(i)),j))*N(j)
     Pmarc(i)=Pmarc(i)+Ptot(connessioni(elem_marc(n_marc(i)),j))*N(j)
     Pmarc_lito(i)=Pmarc_lito(i)+Plito(connessioni(elem_marc(n_marc(i)),j))*N(j)
     exxmarc=exxmarc+exx_nodi(connessioni(elem_marc(n_marc(i)),j))*N(j)
     eyymarc=eyymarc+eyy_nodi(connessioni(elem_marc(n_marc(i)),j))*N(j)
     exymarc=exymarc+exy_nodi(connessioni(elem_marc(n_marc(i)),j))*N(j)
   enddo
   Tmarc(i)=Tmarc(i)*(Tb-Ts)
   e2t_ist(i)=DSQRT(0.5D0*(exxmarc**2+eyymarc**2)+exymarc**2)
 enddo

 end subroutine PT

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine water_markers(perple_x,it,marcatori,num_mat,type_marc,marcatoritot,elem_marc,celle,celleP,n_marc,elementi,elem_x,W_b,W_b_n,&
 Wmax_b,WX,WXs,minTX,minPX,dTX,dPX,dPXs,W_e,W_e_n,materiale,t_mantle,t_serp,t_wet,C,W_f,Tmarc,Pmarc,H2O_0,H2O_wet,e1,e2t,VV,Pmarc_lito,H2O_serp,H2O_ocean,ic)

 Implicit none

 integer :: i,colonna,riga,rigaS
 double precision :: W_e_tmp(elementi),Wmax_S
 logical, intent(in) :: perple_x
 integer, intent(in) :: it,marcatori,num_mat,marcatoritot,elem_marc(marcatoritot),celle,celleP,n_marc(marcatori),elementi,elem_x,&
 t_mantle,t_serp,t_wet,C(elementi),VV(marcatori),ic

 double precision, intent(in) :: WX(celle,celleP,5),WXs(celle,celle),minTX,minPX,dPX,dTX,dPXs,Tmarc(marcatori),Pmarc(marcatori),&
 W_e(elementi),H2O_0,H2O_wet,e1,e2t(marcatori),Pmarc_lito(marcatori),H2O_serp,H2O_ocean

 character(len=*), intent(in) :: materiale(num_mat)
 double precision, intent(out) :: W_f(marcatori),W_e_n(elementi),W_b_n(marcatori)
 integer, intent(inout) :: type_marc(marcatori)
 double precision, intent(inout) :: Wmax_b(marcatori),W_b(marcatori)

 W_e_tmp=W_e
 do i=1,marcatori
   if(perple_x)then
     colonna=INT(((Tmarc(i)-minTX)/dTX)*celle)+1
     riga=INT(((Pmarc(i)-minPX)/dPX)*celleP)+1
     if(colonna.gt.celle)then
       colonna=celle
     elseif(colonna.lt.1)then
       colonna=1
     endif
     if(riga.gt.celleP)then
       riga=celleP
     elseif(riga.lt.1)then
       riga=1
     endif
     rigaS=INT(((Pmarc(i)-minPX)/dPXs)*celle)+1
     if(rigaS.gt.celle)then
       rigaS=celle
     elseif(rigaS.lt.1)then
       rigaS=1
     endif
     Wmax_S=MIN(WXs(colonna,rigaS),H2O_serp)
   else
     Wmax_S=H2O_serp
   endif
   if(materiale(type_marc(i)).eq.'oceanic crust'.or.materiale(type_marc(i)).eq.'lower oceanic crust')then
     if(perple_x)then
       if(it.eq.1.and.W_b(i).eq.0.D0)W_b(i)=MIN(WX(colonna,riga,4),H2O_ocean)
       Wmax_b(i)=MIN(WX(colonna,riga,4),H2O_ocean)
     else
       if(it.eq.1.and.W_b(i).eq.0.D0)W_b(i)=H2O_ocean
       Wmax_b(i)=H2O_ocean
     endif
   elseif(materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'mantle')then
     if(it.eq.1.and.W_b(i).eq.0.D0)W_b(i)=H2O_0
     Wmax_b(i)=H2O_wet
   elseif(materiale(type_marc(i)).eq.'wet mantle')then
     if(it.eq.1.and.W_b(i).eq.0.D0)W_b(i)=H2O_wet
     Wmax_b(i)=H2O_wet
   elseif(materiale(type_marc(i)).eq.'serpentine')then
     if(it.eq.1.and.W_b(i).eq.0.D0)W_b(i)=Wmax_S
     Wmax_b(i)=Wmax_S
   else
     Wmax_b(i)=0.D0
   endif
   if(materiale(type_marc(i)).eq.'serpentine'.and.((VV(i).eq.4.and.e2t(i).ge.e1.and.Pmarc_lito(i).lt.3.D08).or.&
      elem_marc(n_marc(i)).gt.elementi-elem_x))then
     Wmax_b(i)=Wmax_S
     W_b(i)=Wmax_S
   endif
   W_f(i)=W_e_tmp(elem_marc(n_marc(i)))/C(elem_marc(n_marc(i)))
   if(.not.equilib)then
     if(perple_x)then
       if(Wmax_S.gt.H2O_wet.and.(W_b(i)+W_f(i)).gt.H2O_wet.and.(materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'mantle'.or.&
          materiale(type_marc(i)).eq.'wet mantle'))then
         Wmax_b(i)=Wmax_S
         type_marc(i)=t_serp
       endif
       if((Wmax_S.le.H2O_wet.or.(W_b(i)+W_f(i)).le.H2O_wet).and.materiale(type_marc(i)).eq.'serpentine')then
         Wmax_b(i)=H2O_wet
         type_marc(i)=t_mantle
       endif
     else
       if((W_b(i)+W_f(i)).gt.H2O_wet.and.Tmarc(i).le.(751.D0+(0.18D0*Pmarc(i)/1.D06)-(3.1D-05*(Pmarc(i)/1.D06)**2)).and.&
          (materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'mantle'.or.materiale(type_marc(i)).eq.'wet mantle'))then
         Wmax_b(i)=Wmax_S
         type_marc(i)=t_serp
       endif 
       if(Tmarc(i).gt.(751.D0+(0.18D0*Pmarc(i)/1.D06)-(3.1D-05*(Pmarc(i)/1.D06)**2)).and.materiale(type_marc(i)).eq.'serpentine')then
         type_marc(i)=t_mantle
       endif
       if(type_marc(i).eq.t_wet.or.type_marc(i).eq.t_mantle)then
         if(Pmarc(i).lt.2.5D09.and.Pmarc(i).ge.(0.042D09*Tmarc(i)-40.54D09).and.(W_b(i)+W_f(i)).gt.H2O_0)then
           Wmax_b(i)=H2O_wet/2.D0
           type_marc(i)=t_wet
         elseif(Pmarc(i).lt.2.5D09.and.Pmarc(i).lt.(0.042D09*Tmarc(i)-40.54D09).and.(W_b(i)+W_f(i)).gt.H2O_0)then
           Wmax_b(i)=H2O_wet/4.D0
           type_marc(i)=t_wet
         elseif(Pmarc(i).ge.2.5D09.and.Pmarc(i).lt.(-0.023D09*Tmarc(i)+25.98D09).and.(W_b(i)+W_f(i)).gt.H2O_0)then
           Wmax_b(i)=H2O_wet
           type_marc(i)=t_wet
         else
           Wmax_b(i)=H2O_0
           type_marc(i)=t_mantle
         endif
       endif
     endif
     if(.not.perple_x.and.(materiale(type_marc(i)).eq.'oceanic crust'.or.materiale(type_marc(i)).eq.'lower oceanic crust'))then
       if(Pmarc(i).lt.2.4D09.and.Pmarc(i).ge.(0.012D09*Tmarc(i)-7.2D09))then
         Wmax_b(i)=H2O_ocean
       elseif(Pmarc(i).lt.2.4D09.and.Pmarc(i).lt.(0.012D09*Tmarc(i)-7.2D09).and.Pmarc(i).ge.(0.024D09*Tmarc(i)-18.55D09))then
         Wmax_b(i)=2.D0
       elseif(Pmarc(i).lt.2.4D09.and.Pmarc(i).lt.(0.024D09*Tmarc(i)-18.55D09))then
         Wmax_b(i)=1.5D0
       elseif(Pmarc(i).ge.2.4D09.and.Pmarc(i).ge.(0.008D09*Tmarc(i)-4.D09).and.Tmarc(i).lt.923.15D0.and.Pmarc(i).lt.7.D09)then
         Wmax_b(i)=1.D0
       elseif(Pmarc(i).ge.7.D09)then
         Wmax_b(i)=0.1D0
       else
         Wmax_b(i)=0.5D0
       endif
     endif
   endif
   if(ic.eq.1)W_b_n(i)=W_b(i)
   if(W_b_n(i).gt.Wmax_b(i))then
     W_e_tmp(elem_marc(n_marc(i)))=W_e_tmp(elem_marc(n_marc(i)))+(W_b_n(i)-Wmax_b(i))
     W_b_n(i)=Wmax_b(i)
   endif
 enddo
 W_e_n=0.D0
 do i=1,marcatori
   W_f(i)=W_e_tmp(elem_marc(n_marc(i)))/C(elem_marc(n_marc(i)))
   if(W_f(i).gt.0.D0.and.W_b_n(i).lt.Wmax_b(i))then
     if(W_f(i).gt.(Wmax_b(i)-W_b_n(i)))then
       W_f(i)=W_f(i)-(Wmax_b(i)-W_b_n(i))
       W_b_n(i)=Wmax_b(i)
     else
       W_b_n(i)=W_b_n(i)+W_f(i)
       W_f(i)=0.D0
     endif
   endif
   W_e_n(elem_marc(n_marc(i)))=W_e_n(elem_marc(n_marc(i)))+W_f(i)
 enddo
 do i=elementi,elem_x+1,-1
   W_e_tmp(i)=W_e_n(i-elem_x)
 enddo
 W_e_n=0.D0
 do i=1,marcatori
   W_f(i)=W_e_tmp(elem_marc(n_marc(i)))/C(elem_marc(n_marc(i)))
   if(W_f(i).gt.0.D0.and.W_b_n(i).lt.Wmax_b(i))then
     if(W_f(i).gt.(Wmax_b(i)-W_b_n(i)))then
       W_f(i)=W_f(i)-(Wmax_b(i)-W_b_n(i))
       W_b_n(i)=Wmax_b(i)
     else
       W_b_n(i)=W_b_n(i)+W_f(i)
       W_f(i)=0.D0
     endif
   endif
   W_e_n(elem_marc(n_marc(i)))=W_e_n(elem_marc(n_marc(i)))+W_f(i)
 enddo

 end subroutine water_markers

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine melt(materiale,W_b_n,W_f,W_e_n,melt_fraction,Pmarc,Tmarc,Rho_melt,Rho_M,alpha_M,Cp_M,Tref_M,&
 num_mat,Cp_melt,alpha_melt,W_m,W_m_n,Wmax_b,kappa,gamma,melting,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol,Hl,&
 t_oc,t_mantleM,t_lmantle,type_marc,beta_press,melt_oc,ocean,melt_ext)

 Implicit none

 integer :: it
 double precision :: water,press,T_sol,T_sl,temp,x0,x1,f0,f1,xc,fc,P_sol,deltaT,deltaT_max,melt_max,dT,dP
 double precision, external :: molten,delta,roots
 logical, intent(in) :: melting,ocean
 character(len=*), intent(in) :: materiale
 integer, intent(in) :: num_mat,t_oc,t_mantleM,t_lmantle
 double precision, intent(in) :: Pmarc,Tmarc,Rho_M(num_mat),alpha_M(num_mat),Tref_M(num_mat),Cp_M(num_mat),W_m,&
 kappa,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol,Hl,beta_press,melt_ext

 double precision, intent(out) :: W_m_n
 logical, intent(inout) :: melt_oc
 integer, intent(inout) :: type_marc
 double precision, intent(inout) :: W_b_n,W_f,W_e_n,melt_fraction,Rho_melt,Cp_melt,alpha_melt,Wmax_b

 melt_oc=.false.
 if(melting)then
   W_m_n=0.D0
   water=W_b_n+W_f+W_m
   if(Pmarc.lt.1.D10)then
     press=Pmarc/1.D09
     temp=Tmarc-273.15D0
     T_sol=A1+(A2*press)+(A3*(press**2))
     T_sl=x_t*log(press)+y_t
     x0=0.D0
     f0=molten(temp,T_sol,T_sl,kappa,water,DH2O,gamma,x0)
     x1=1.D0
     f1=molten(temp,T_sol,T_sl,kappa,water,DH2O,gamma,x1)
     if(f0*f1.le.0.D0.and.(materiale.eq.'upper mantle'.or.materiale.eq.'melt'.or.materiale.eq.'mantle'.or.materiale.eq.'wet mantle'))then
       it=1
       fc=f0
       do while(DABS(fc).gt.tol.and.it.le.100)
         xc=x1-f1*((x1-x0)/(f1-f0))
         fc=molten(temp,T_sol,T_sl,kappa,water,DH2O,gamma,xc)
         if(f0*fc.le.0.D0)then
           x1=xc
         else
           x0=xc
         endif
         f0=molten(temp,T_sol,T_sl,kappa,water,DH2O,gamma,x0)
         f1=molten(temp,T_sol,T_sl,kappa,water,DH2O,gamma,x1)
         it=it+1
       enddo
       deltaT_max=kappa*(chi*(press**lambda)+press)**gamma
       melt_max=(temp-(T_sol-deltaT_max))/T_sl
       melt_fraction=MIN(xc,melt_max)
       if(melt_fraction.lt.0.D0)melt_fraction=0.D0
       if(melt_fraction.gt.1.D0)melt_fraction=1.D0
       deltaT=MIN(kappa*(water/(DH2O*(1.D0-melt_fraction)+melt_fraction))**gamma,deltaT_max)
       dT=DABS(temp-(T_sol-deltaT))
       if(delta(A3,A2,A1-temp-deltaT).gt.0.D0)then
         P_sol=MIN(roots(A3,A2,sqrt(delta(A3,A2,A1-temp-deltaT))),roots(A3,A2,-sqrt(delta(A3,A2,A1-temp-deltaT))))*1.D09
       elseif(delta(A3,A2,A1-temp-deltaT).eq.0.D0)then
         P_sol=roots(A3,A2,0.D0)*1.D09
       endif
       dP=(Tmarc*DABS(Pmarc-P_sol))
       if(dT.lt.1.D0.or.dP.lt.1.D0)melt_fraction=0.D0
       if(melt_fraction.gt.0.D0)then
         if(ocean.and.melt_fraction.gt.melt_ext/100.D0)then
           melt_oc=.true.
           melt_fraction=melt_ext
         endif
         type_marc=t_mantleM
         Rho_melt=Rho_M(type_marc)*(1.D0+Pmarc*beta_press)*(1.D0-alpha_M(type_marc)*(Tmarc-Tref_M(type_marc)))
         Cp_melt=Cp_M(type_marc)+((Hl*melt_fraction)/dT)
         alpha_melt=alpha_M(type_marc)+((Rho_melt*Hl*melt_fraction)/dP)
       elseif(melt_fraction.eq.0.D0.and.type_marc.eq.t_mantleM)then
         type_marc=t_lmantle
       endif
     endif
   endif
   W_m_n=(water*melt_fraction)/(DH2O+(melt_fraction*(1.D0-DH2O)))
   Wmax_b=Wmax_b*(1.D0-melt_fraction)
   W_b_n=MIN(water-W_m_n,Wmax_b)
   W_f=water-W_m_n-W_b_n
   W_e_n=W_e_n+W_f
 else
   Rho_melt=Rho_M(t_mantle)
   Cp_melt=Cp_M(t_mantle)
   alpha_melt=alpha_M(t_mantle)
 endif
    
 end subroutine melt

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine PerpleX(rho_press,Rho_X,Cp_X,alpha_X,Tmarc,Pmarc,materiale,dTX,dPX,minTX,minPX,RhoX,alphaX,CpX,&
 dPXs,RhoXs,CpXs,alphaXs,Rho_XC,celle,celleP,Cpref,alpharef,perple_x,Rho_M,Cp_M,alpha_M,Tref_M,Rho_S,water,beta_press,delta_water)

 Implicit none

 integer :: riga,colonna
 logical, intent(in) :: rho_press
 integer, intent(in) :: celle,celleP
 double precision, intent(in) :: Cpref,alpharef,RhoX(celle,celleP,5),dTX,dPX,minTX,minPX,CpX(celle,celleP,5),&
 alphaX(celle,celleP,5),Pmarc,Tmarc,dPXs,RhoXs(celle,celle),CpXs(celle,celle),alphaXs(celle,celle),Rho_M,Cp_M,&
 alpha_M,Tref_M,Rho_S,water,beta_press,delta_water
 
 logical, intent(in) :: perple_x
 character(len=*), intent(in) :: materiale
 double precision, intent(out) :: Cp_X,alpha_X,Rho_X,Rho_XC

 if(perple_x)then
   colonna=INT(((Tmarc-minTX)/dTX)*celle)+1
   riga=INT(((Pmarc-minPX)/dPX)*celleP)+1
   if(colonna.gt.celle)then
     colonna=celle
   elseif(colonna.lt.1)then
     colonna=1
   endif
   if(riga.gt.celleP)then
     riga=celleP
   elseif(riga.lt.1)then
     riga=1
   endif
   if(materiale.eq.'air')then
     Rho_X=1.D0
     Rho_XC=1.D0
     alpha_X=alpharef
     Cp_X=Cpref
   elseif(materiale.eq.'upper mantle'.or.materiale.eq.'wet mantle'.or.materiale.eq.'mantle'.or.materiale.eq.'melt')then
     Rho_X=RhoX(colonna,riga,1)
     Rho_XC=RhoX(1,riga,1)
     alpha_X=alphaX(colonna,riga,1)
     Cp_X=CpX(colonna,riga,1)
   elseif(materiale.eq.'continental crust')then
     Rho_X=RhoX(colonna,riga,2)
     Rho_XC=RhoX(1,riga,2)
     alpha_X=alphaX(colonna,riga,2)
     Cp_X=CpX(colonna,riga,2)
   elseif(materiale.eq.'lower continental crust')then
     Rho_X=RhoX(colonna,riga,3)
     Rho_XC=RhoX(1,riga,3)
     alpha_X=alphaX(colonna,riga,3)
     Cp_X=CpX(colonna,riga,3)
   elseif(materiale.eq.'oceanic crust'.or.materiale.eq.'lower oceanic crust')then
     Rho_X=RhoX(colonna,riga,4)
     Rho_XC=RhoX(1,riga,4)
     alpha_X=alphaX(colonna,riga,4)
     Cp_X=CpX(colonna,riga,4)
   elseif(materiale.eq.'sediment')then
     Rho_X=RhoX(colonna,riga,5)
     Rho_XC=RhoX(1,riga,5)
     alpha_X=alphaX(colonna,riga,5)
     Cp_X=CpX(colonna,riga,5)
   elseif(materiale.eq.'serpentine')then
     riga=INT(((Pmarc-minPX)/dPXs)*celle)+1
     if(riga.gt.celle)then
       riga=celle
     elseif(riga.lt.1)then
       riga=1
     endif
     Rho_X=RhoXs(colonna,riga)
     Rho_XC=RhoXs(1,riga)
     alpha_X=alphaXs(colonna,riga)
     Cp_X=CpXs(colonna,riga)
   endif
 else
   if(materiale.eq.'melt')then
     Rho_X=Rho_S*(1.D0-alpha_M*(Tmarc-Tref_M))
     Rho_XC=Rho_S
   else
     if(rho_press)then
       Rho_X=Rho_M*(1.D0+Pmarc*beta_press)*(1.D0-alpha_M*(Tmarc-Tref_M)-water*delta_water)
     else
       Rho_X=Rho_M*(1.D0-alpha_M*(Tmarc-Tref_M)-water*delta_water)
     endif
     Rho_XC=Rho_M
   endif
   alpha_X=alpha_M
   Cp_X=Cp_M
 endif

 end subroutine PerpleX

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Elemental_properties(C,Cp_X,Rho_X,alpha_X,Rho_melt,Cp_melt,alpha_melt,Rho_M,Hr_M,Kc_M,Rhoref,Cpref,Kcref,&
 RhoC,Rho,Cp,alpha,Hr,Kc,melt_fraction,alpharef,Href,e2t_p,e2t_v,e2t_p_elem,e2t_v_elem,melt_elem,W_b,w_b_elem)

 Implicit none

 double precision :: Rho_sm,Cp_sm,alpha_sm
 integer, intent(in) :: C
 double precision, intent(in) :: melt_fraction,Cp_X,Rho_X,alpha_X,Rho_melt,Cp_melt,alpha_melt,Rho_M,Hr_M,Kc_M,&
 Rhoref,Cpref,Kcref,alpharef,Href,e2t_p,e2t_v,W_b
 
 double precision, intent(inout) :: RhoC,Rho,Cp,alpha,Hr,Kc,e2t_p_elem,e2t_v_elem,melt_elem,w_b_elem

 Rho_sm=Rho_M*(1.D0-melt_fraction)+Rho_melt*melt_fraction
 RhoC=RhoC+(Rho_sm/DBLE(C))
 Rho_sm=Rho_X*(1.D0-melt_fraction)+Rho_melt*melt_fraction
 Rho=Rho+(Rho_sm/DBLE(C))/Rhoref
 Cp_sm=Cp_X*(1.D0-melt_fraction)+Cp_melt*melt_fraction
 Cp=Cp+(Cp_sm/DBLE(C))/Cpref
 alpha_sm=alpha_X*(1.D0-melt_fraction)+alpha_melt*melt_fraction
 alpha=alpha+(alpha_sm/DBLE(C))/alpharef
 Hr=Hr+((Hr_M/Rho_M)/DBLE(C))/Href
 Kc=Kc+(Kc_M/DBLE(C))/Kcref
 e2t_p_elem=e2t_p_elem+(e2t_p/DBLE(C))
 e2t_v_elem=e2t_v_elem+(e2t_v/DBLE(C))
 melt_elem=melt_elem+(melt_fraction/DBLE(C))
 w_b_elem=w_b_elem+(W_b/DBLE(C))

 end subroutine Elemental_properties

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Viscosita(coeffA,ndisl,Va,Ea,Cohe_0,coeffAdf,Eadf,Vadf,teta_fr_0,num_mat,materiale,Rgas,elem_marc,pi,Visc_M,grain_m,grain_p,&
 Visc_max,Visc_min,C,elementi,media,Viscref,Visc,marcatori,marcatoritot,Pmarc,Tmarc,Visc_eff,n_marc,type_marc,e2t_p,e2t_v,e2t_ist,VV,&
 W_b,melt_fraction,teta_fr_inf,Cohe_inf,Ws,e0,e1,e0v,e1v,alpham_dl,alpham_df,rw,OH_0,Pmarc_lito,diffusion,rheology,t_serp,elem_x,eff_average,f_scale,&
 uniaxial,teta_elem,cohe_elem,idratazione,t_oc,serpentine,peierls,Ea_pe,ndisl_pe,Va_pe,coeffA_pe,melt_oc,colonna_marc)

 Implicit none

 integer :: i,VVm
 double precision :: yield,Visc_pl,Visc_dl,Visc_df,Visc_pe,teta_fr,Cohe,Wv,Visc_cp,Visc_vp
 logical, intent(in) :: diffusion(num_mat),eff_average,uniaxial(num_mat),idratazione,serpentine,peierls(num_mat),melt_oc(elem_x)
 integer, intent(in) :: elementi,num_mat,marcatori,n_marc(marcatori),marcatoritot,elem_marc(marcatoritot),C(elementi),rheology(num_mat),&
 t_serp,elem_x,t_oc,colonna_marc(marcatoritot)

 double precision, intent(in) :: coeffA(num_mat),ndisl(num_mat),Va(num_mat),Ea(num_mat),Cohe_0(num_mat),coeffAdf(num_mat),&
 Eadf(num_mat),Vadf(num_mat),teta_fr_0(num_mat),Visc_M(num_mat),Rgas,pi,Visc_max(num_mat),Visc_min(num_mat),Viscref,&
 W_b(marcatori),Tmarc(marcatori),Pmarc(marcatori),e2t_ist(marcatori),melt_fraction(marcatori),e0,e1,e0v,e1v,alpham_dl,alpham_df,rw,OH_0,&
 teta_fr_inf(num_mat),Cohe_inf(num_mat),Ws(num_mat),f_scale(num_mat),Pmarc_lito(marcatori),Ea_pe(num_mat),ndisl_pe(num_mat),&
 Va_pe(num_mat),coeffA_pe(num_mat),grain_m(num_mat),grain_p(num_mat)

 character(len=*), intent(in) :: media,materiale(num_mat)
 integer, intent(out) :: VV(marcatori)
 double precision, intent(out) :: Visc(elementi),Visc_eff(marcatori),teta_elem(elementi),cohe_elem(elementi)
 integer, intent(inout) :: type_marc(marcatori)
 double precision, intent(inout) :: e2t_p(marcatori),e2t_v(marcatori)

 teta_elem=0.D0
 Cohe_elem=0.D0
 if(media.eq.'aritmetica'.or.media.eq.'armonica')then
   Visc=0.D0
 else
   Visc=1.D0
 endif
 do i=1,marcatori
   if(materiale(type_marc(i)).eq.'air')then
     e2t_p(i)=0.D0
     e2t_v(i)=0.D0
   endif
   if(e2t_v(i).le.e1v.and.e2t_v(i).ge.e0v)then
     Wv=1.D0+(Ws(type_marc(i))-1.D0)*((e2t_v(i)-e0v)/(e1v-e0v))
   elseif(e2t_v(i).gt.e1v)then
     Wv=Ws(type_marc(i))
   else
     Wv=1.D0
   endif
   if(serpentine.and.(materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'mantle'.or.materiale(type_marc(i)).eq.'wet mantle').and.&
       ((VV(i).eq.4.and.e2t_p(i).ge.e1.and.Pmarc_lito(i).lt.3.D08).or.elem_marc(n_marc(i)).gt.elementi-elem_x))type_marc(i)=t_serp
   if(melt_oc(colonna_marc(n_marc(i))).and.(materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'serpentine'.or.materiale(type_marc(i)).eq.&
      'mantle').and.(Pmarc_lito(i).lt.3.D08.or.elem_marc(n_marc(i)).gt.elementi-elem_x))type_marc(i)=t_oc
   Visc_dl=0.5D0*f_scale(type_marc(i))*(1.D0/(Wv*(coeffA(type_marc(i))**(1.D0/ndisl(type_marc(i))))))*&
           (e2t_ist(i)**((1.D0-ndisl(type_marc(i)))/(1.D0*ndisl(type_marc(i)))))*DEXP((Ea(type_marc(i))+&
           Va(type_marc(i))*Pmarc(i))/(ndisl(type_marc(i))*Rgas*Tmarc(i)))*DEXP(alpham_dl*melt_fraction(i))
   Visc_df=0.5D0*((grain_m(type_marc(i))**grain_p(type_marc(i)))/(Wv*coeffAdf(type_marc(i))))*DEXP((Eadf(type_marc(i))+Vadf(type_marc(i))*Pmarc(i))/&
           (Rgas*Tmarc(i)))*DEXP(alpham_df*melt_fraction(i))
   Visc_pe=0.5D0*(1.D0/(Wv*(coeffA_pe(type_marc(i))**(1.D0/ndisl_pe(type_marc(i))))))*&
           (e2t_ist(i)**((1.D0-ndisl_pe(type_marc(i)))/(1.D0*ndisl_pe(type_marc(i)))))*DEXP((Ea_pe(type_marc(i))+&
           Va_pe(type_marc(i))*Pmarc(i))/(ndisl_pe(type_marc(i))*Rgas*Tmarc(i)))
   if(uniaxial(type_marc(i)))then
     Visc_dl=Visc_dl*(2.D0**((1.D0-ndisl(type_marc(i)))/ndisl(type_marc(i))))/(3.D0**((ndisl(type_marc(i))+1.D0)/(2.D0*ndisl(type_marc(i)))))
     Visc_df=Visc_df/3.D0
   endif
   if(idratazione.and.(materiale(type_marc(i)).eq.'upper mantle'.or.materiale(type_marc(i)).eq.'mantle'.or.materiale(type_marc(i)).eq.'wet mantle'.or.&
      materiale(type_marc(i)).eq.'serpentine').and.W_b(i).gt.OH_0)then
     Visc_dl=Visc_dl*(W_b(i)/OH_0)**(-rw)
     Visc_df=Visc_df*(W_b(i)/OH_0)**(-rw)
     Visc_pe=Visc_pe*(W_b(i)/OH_0)**(-rw)
   endif
   if(e2t_p(i).le.e1.and.e2t_p(i).ge.e0)then
     teta_fr=teta_fr_0(type_marc(i))+(teta_fr_inf(type_marc(i))-teta_fr_0(type_marc(i)))*((e2t_p(i)-e0)/(e1-e0))
     Cohe=Cohe_0(type_marc(i))+(Cohe_inf(type_marc(i))-Cohe_0(type_marc(i)))*((e2t_p(i)-e0)/(e1-e0))
   elseif(e2t_p(i).gt.e1)then
     teta_fr=teta_fr_inf(type_marc(i))
     Cohe=Cohe_inf(type_marc(i))
   else
     teta_fr=teta_fr_0(type_marc(i))
     Cohe=Cohe_0(type_marc(i))
   endif
   yield=Cohe*DCOS(teta_fr*pi/180.D0)+Pmarc(i)*DSIN(teta_fr*pi/180.D0)
   Visc_pl=yield/(2.D0*e2t_ist(i))
   if(diffusion(type_marc(i)).and..not.peierls(type_marc(i)))then
     Visc_cp=((1.D0/Visc_dl)+(1.D0/Visc_df))**(-1)
     if(MIN(Visc_dl,Visc_df).eq.Visc_df)then
       VVm=3
     else
       VVm=2
     endif
   elseif(.not.diffusion(type_marc(i)).and.peierls(type_marc(i)))then
     Visc_cp=((1.D0/Visc_dl)+(1.D0/Visc_pe))**(-1)
     if(MIN(Visc_dl,Visc_pe).eq.Visc_pe)then
       VVm=4
     else
       VVm=2
     endif
   elseif(diffusion(type_marc(i)).and.peierls(type_marc(i)))then
     Visc_cp=((1.D0/Visc_dl)+(1.D0/Visc_df)+(1.D0/Visc_pe))**(-1)
     if(MIN(Visc_dl,Visc_df,Visc_pe).eq.Visc_df)then
       VVm=3
     elseif(MIN(Visc_dl,Visc_df,Visc_pe).eq.Visc_dl)then
       VVm=2
     else
       VVm=4
     endif
   else
     Visc_cp=Visc_dl
     VVm=2
   endif
   if(rheology(type_marc(i)).eq.1)then
     Visc_vp=Visc_M(type_marc(i))
     VV(i)=1
   elseif(rheology(type_marc(i)).eq.2)then
     Visc_vp=Visc_cp
     VV(i)=VVm
   elseif((rheology(type_marc(i)).eq.3))then
     Visc_vp=Visc_pl
     VV(i)=5
   elseif(rheology(type_marc(i)).eq.4)then
     if(eff_average)then
       Visc_vp=((1.D0/Visc_cp)+(1.D0/Visc_pl))**(-1)
     else
       Visc_vp=MIN(Visc_pl,Visc_cp)
     endif
     if(Visc_pl.lt.Visc_cp)then
       VV(i)=5
     else
       VV(i)=VVm
     endif
   endif
   if(eff_average)then
     Visc_eff(i)=(Visc_min(type_marc(i))+((1.D0/Visc_max(type_marc(i)))+(1.D0/Visc_vp))**(-1))/Viscref
   else
     Visc_eff(i)=MIN(MAX(Visc_vp,Visc_min(type_marc(i))),Visc_max(type_marc(i)))/Viscref
   endif
   if(media.eq.'aritmetica')then
     Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))+(Visc_eff(i)/DBLE(C(elem_marc(n_marc(i)))))
   elseif(media.eq.'armonica')then
     Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))+(1.D0/Visc_eff(i))
   else
     Visc(elem_marc(n_marc(i)))=Visc(elem_marc(n_marc(i)))*(Visc_eff(i)**(1.D0/DBLE(C(elem_marc(n_marc(i))))))
   endif
   teta_elem(elem_marc(n_marc(i)))=teta_elem(elem_marc(n_marc(i)))+(teta_fr/DBLE(C(elem_marc(n_marc(i)))))
   cohe_elem(elem_marc(n_marc(i)))=cohe_elem(elem_marc(n_marc(i)))+(Cohe/DBLE(C(elem_marc(n_marc(i)))))
 enddo
 if(media.eq.'armonica')Visc=DBLE(C)/Visc

 end subroutine Viscosita

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Output_nl_vtu(OUT_VEL,OUT_DIR,it,ic,elementi,nodi_elem,nodi,connessioni,x,y,Visc,Ptot,dPtot,Press_old,dPress,Vel0,dvel0,res,ndof,e2_nodi,&
 OUT_MODEL,OUT_PHASE)

 Implicit none

 integer :: i,j
 double precision :: Vel(nodi*ndof),dvel(nodi*ndof)
 character(len=100) :: filename
 integer, intent(in) :: it,ic,nodi,elementi,nodi_elem,connessioni(elementi,nodi_elem),ndof
 double precision, intent(in) :: x(nodi),y(nodi),Ptot(nodi),Visc(elementi),dPtot(nodi),res(nodi*ndof),Press_old(elementi),dPress(elementi),&
 e2_nodi(nodi),Vel0(nodi*ndof),dvel0(nodi*ndof)

 logical, intent(in) :: OUT_VEL
 character(len=*), intent(in) :: OUT_DIR,OUT_MODEL,OUT_PHASE

 if(OUT_VEL)then
   Vel=Vel0*1.D02*365.D0*24.D0*3600.D0
   dvel=dvel0*1.D02*365.D0*24.D0*3600.D0
 else
   Vel=Vel0
   dvel=dvel0
 endif

 write(filename,'("/vtu_files/Convergence.",I2.2,".",I3.3,".vtu")')it,ic
 open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown',form='formatted')
 rewind(1)

 write(1,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
 write(1,*) '<UnstructuredGrid>'
 write(1,*) '<Piece NumberOfPoints="',nodi,'" NumberOfCells="',elementi,'">'

 write(1,*) '<PointData Scalars="scalars">'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Name="Velocity" Format="ascii">'
 do i=1,nodi
   write(1,*) Vel(i*2-1),Vel(i*2),0.0
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Name="Diff. Velocity" Format="ascii">'
 do i=1,nodi
   write(1,*) dVel(i*2-1),dVel(i*2),0.0
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Name="Residual" Format="ascii">'
 do i=1,nodi
   write(1,*) res(i*2-1),res(i*2),0.0
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Smoothed Pressure" Format="ascii">'
 do i=1,nodi
   write(1,*) Ptot(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Diff. Smoothed Pressure" Format="ascii">'
 do i=1,nodi
   write(1,*) dPtot(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Strain Rates" Format="ascii">'
 do i=1,nodi
   write(1,*) e2_nodi(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</PointData>'
 write(1,*) '<Points>'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
 do i=1,nodi
   write(1,*) x(i),y(i),0.D0
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</Points>'

 write(1,*) '<CellData Scalars="scalars">'
 write(1,*) '<DataArray type="Float32" Name="Elemental pressure" Format="ascii">'
 do i=1,elementi
   write(1,*) Press_old(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Diff. elem. pressure" Format="ascii">'
 do i=1,elementi
   write(1,*) dPress(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Viscosity" Format="ascii">'
 do i=1,elementi
   write(1,*) Visc(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</CellData>'
 write(1,*) '<Cells>'
 write(1,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
 do i=1,elementi
   write(1,*) (connessioni(i,j)-1,j=1,nodi_elem)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
 write(1,*) (i*4,i=1,elementi)
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="types" Format="ascii">'
 write(1,*) (9,i=1,elementi)
 write(1,*) '</DataArray>'
 write(1,*) '</Cells>'
 write(1,*) '</Piece>'
 write(1,*) '</UnstructuredGrid>'
 write(1,*) '</VTKFile>'
 close(1)

 end subroutine Output_nl_vtu

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Output_nl(OUT_VEL,OUT_DIR,it,ic,it_max,res,vel,press,Pn_min,Pn_max,Pe_min,Pe_max,visc_min,visc_max,Vx0_min,Vx0_max,Vy0_min,Vy0_max,&
 OUT_MODEL,OUT_PHASE)

 implicit none

 double precision :: Vx_min,Vx_max,Vy_min,Vy_max
 integer, intent(in) :: it,ic,it_max
 double precision, intent(in) :: res,press,Pn_min,Pn_max,Pe_min,Pe_max,visc_min,visc_max,vel,Vx0_min,Vx0_max,Vy0_min,Vy0_max
 character(len=*), intent(in) :: OUT_DIR,OUT_MODEL,OUT_PHASE
 logical, intent(in) :: OUT_VEL

 if(OUT_VEL)then
   Vx_min=Vx0_min*1.D02*365.D0*24.D0*3600.D0
   Vy_min=Vy0_min*1.D02*365.D0*24.D0*3600.D0
   Vx_max=Vx0_max*1.D02*365.D0*24.D0*3600.D0
   Vy_max=Vy0_max*1.D02*365.D0*24.D0*3600.D0
 else
   Vx_min=Vx0_min
   Vx_max=Vx0_max
   Vy_min=Vy0_min
   Vy_max=Vy0_max
 endif

 if(it.eq.1.and.ic.eq.1)then
   open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Residuals.txt',status='replace')
   open(unit=2,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/MinMax.txt',status='replace')
 else
   open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Residuals.txt',status='old',position='append')
   open(unit=2,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/MinMax.txt',status='old',position='append')
 endif
 write(1,*)DBLE(it)+DBLE(ic)/(1.5D0*it_max),res,vel,press
 write(2,*)DBLE(it)+DBLE(ic)/(1.5D0*it_max),Pn_min,Pn_max,Pe_min,Pe_max,visc_min,visc_max,Vx_min,Vx_max,Vy_min,Vy_max
 close(1)
 close(2)

 end subroutine Output_nl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine marker_to_element(num_mat,elem_marc,n_marc,marcatori,marcatoritot,elementi,type_marc,type_elem,VV,def_elem)

 Implicit none

 integer :: i,cont(elementi,num_mat),cont_def(elementi,5)
 integer, intent(in) :: num_mat,elem_marc(marcatoritot),n_marc(marcatori),marcatori,marcatoritot,elementi,type_marc(marcatori),VV(marcatori)
 integer, intent(out) :: type_elem(elementi),def_elem(elementi)

 cont=0
 cont_def=0
 do i=1,marcatori
   cont(elem_marc(n_marc(i)),type_marc(i))=cont(elem_marc(n_marc(i)),type_marc(i))+1
   cont_def(elem_marc(n_marc(i)),VV(i))=cont_def(elem_marc(n_marc(i)),VV(i))+1
 enddo
 do i=1,elementi
   type_elem(i)=MAXLOC(cont(i,:),DIM=1)
   def_elem(i)=MAXLOC(cont_def(i,:),DIM=1)
 enddo

 end subroutine marker_to_element

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Output(OUT_VEL,OUT_DIR,elementi,nodi_elem,nodi,connessioni,x,y,x_elem,y_elem,Temp,Vx0,Vy0,v0,Cp,Visc,Rho,Kc,Telem,Ptot,&
 Plito,Ptot_nodi,x_marc,y_marc,n_marc,type_marc,VV,marcatori,e2_nodi,s2,Pmarc,Tmarc,e2t_p,e2t_v,Hr,Hs,Hd,tempo,W_b,W_f,melt_fraction,alpha,W_m,melting,&
 perple_x,e2t_ist,Visc_m,comp,OUT_TYPE,nt,OUT_FILE_nb,it,q_x,q_y,power_law,teta,cohe,e2t_p_elem,e2t_v_elem,markers_vtu_out,markers_txt_out,type_elem,def_elem,&
 shear,adiabatic,solver_T,melt_elem,OUT_MODEL,OUT_PHASE,w_b_elem,W_e)

 Implicit none
 
 character(len=100) :: filename
 integer, parameter :: vertype=1
 integer :: i,j,name
 double precision :: Vx(nodi),Vy(nodi)
 character(len=*), intent(in) :: OUT_DIR,OUT_MODEL,OUT_PHASE
 logical, intent(in) :: OUT_VEL,melting,perple_x,OUT_TYPE,power_law,markers_vtu_out,markers_txt_out,shear,adiabatic,solver_T
 integer, intent(in) :: it,nodi,elementi,nodi_elem,connessioni(elementi,nodi_elem),marcatori,n_marc(marcatori),type_marc(marcatori),&
 VV(marcatori),comp(elementi),OUT_FILE_nb,type_elem(elementi),def_elem(elementi)

 double precision, intent(in) :: x(nodi),y(nodi),x_elem(elementi),y_elem(elementi),Temp(nodi),Ptot(elementi),Vx0(nodi),Vy0(nodi),tempo,&
 v0,Cp(elementi),s2(elementi),Visc(elementi),Rho(elementi),Telem(elementi),Plito(elementi),e2_nodi(nodi),Hr(elementi),Hs(elementi),Hd(elementi),&
 Kc(elementi),x_marc(marcatori),y_marc(marcatori),Pmarc(marcatori),W_m(marcatori),Tmarc(marcatori),e2t_p(marcatori),e2t_v(marcatori),W_b(marcatori),&
 W_f(marcatori),melt_fraction(marcatori),alpha(elementi),Ptot_nodi(nodi),e2t_ist(marcatori),Visc_m(marcatori),q_x(elementi),q_y(elementi),teta(elementi),&
 cohe(elementi),e2t_p_elem(elementi),e2t_v_elem(elementi),melt_elem(elementi),w_b_elem(elementi),W_e(elementi)

 integer, intent(inout) :: nt

 if(OUT_VEL)then
   Vx=Vx0*v0*1.D02*365.D0*24.D0*3600.D0
   Vy=Vy0*v0*1.D02*365.D0*24.D0*3600.D0
 else
   Vx=Vx0
   Vy=Vy0
 endif
 if(OUT_TYPE)name=FLOOR(tempo*1000.D0)
 if(.not.OUT_TYPE)name=it-1

 write(filename,'("/vtu_files/Solution.",I5.5,".vtu")')name
 open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown',form='formatted')
 rewind(1)

 write(1,104) tempo
 write(1,*) '<UnstructuredGrid>'
 write(1,*) '<Piece NumberOfPoints="',nodi,'" NumberOfCells="',elementi,'">'

 write(1,*) '<PointData Scalars="scalars">'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Name="Velocity" Format="ascii">'
 do i=1,nodi
   write(1,*) Vx(i),Vy(i),0.0
 enddo
 if(solver_T)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Temperature" Format="ascii">'
   do i=1,nodi
     write(1,*) Temp(i)
   enddo
 endif
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Strain rate" Format="ascii">'
 do i=1,nodi
   write(1,*) e2_nodi(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Total pressure" Format="ascii">'
 do i=1,nodi
   write(1,*) Ptot_nodi(i)
 enddo
 write(filename,'("/txt_files/Nodi.",i5.5,".txt")')name
 open(unit=2,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown')
 write(2,108) tempo
 write(2,2002)
 write(1,*) '</DataArray>'
 write(1,*) '</PointData>'
 write(1,*) '<Points>'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
 do i=1,nodi
   write(1,*) x(i),y(i),0.D0
   write(2,'(i10,1x,2(f20.9,1x),3(f15.4,1x),2(g15.5,1x))')i,x(i),y(i),Vx(i),Vy(i),Temp(i),Ptot_nodi(i),e2_nodi(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</Points>'
 close(2)

 write(filename,'("/txt_files/Elementi.",i5.5,".txt")')name
 open(unit=2,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown')
 write(2,108) tempo
 write(2,2001)
 write(1,*) '<CellData Scalars="scalars">'
 write(1,*) '<DataArray type="Float32" Name="Density" Format="ascii">'
 do i=1,elementi
   write(1,*) Rho(i)
   write(2,'(i10,1x,2(f20.9,1x),4(f10.4,1x),3(g15.5,1x),2x,f10.4,1x,6(g15.5,1x),2(f15.4,1x),2(i5,1x),f10.4,1x,g15.5,1x,f15.5)')i,x_elem(i),&
        y_elem(i),Hr(i),Kc(i),Rho(i),Cp(i),alpha(i),Plito(i),Ptot(i),Telem(i),Visc(i),Hs(i),Hd(i),s2(i),q_x(i),q_y(i),e2t_p_elem(i),e2t_v_elem(i),&
        type_elem(i),def_elem(i),teta(i),cohe(i),melt_elem(i)
 enddo
 close(2)
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Viscosity" Format="ascii">'
 do i=1,elementi
   write(1,*) Visc(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Stress" Format="ascii">'
 do i=1,elementi
   write(1,*) s2(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Number of Markers" Format="ascii">'
 do i=1,elementi
   write(1,*) comp(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Pressure lith." Format="ascii">'
 do i=1,nodi
   write(1,*) Plito(i)
 enddo
 if(solver_t.and.shear)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heating shear" Format="ascii">'
   do i=1,elementi
     write(1,*) Hs(i)
   enddo
 endif
 if(solver_t.and.adiabatic)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heating adiabatic" Format="ascii">'
   do i=1,elementi
     write(1,*) Hd(i)
   enddo
 endif
 if(solver_T)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heating radiogenic" Format="ascii">'
   do i=1,elementi
     write(1,*) Hr(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heating total" Format="ascii">'
   do i=1,elementi
     write(1,*) Hd(i)+Hs(i)+Hr(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heat flux X" Format="ascii">'
   do i=1,elementi
     write(1,*) q_x(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Heat flux Y" Format="ascii">'
   do i=1,elementi
     write(1,*) q_y(i)
   enddo
 endif
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Type" Format="ascii">'
 do i=1,elementi
   write(1,*) type_elem(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Strain plastic" Format="ascii">'
 do i=1,elementi
   write(1,*) e2t_p_elem(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Strain viscous" Format="ascii">'
 do i=1,elementi
   write(1,*) e2t_v_elem(i)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Float32" Name="Strain total" Format="ascii">'
 do i=1,elementi
   write(1,*) e2t_p_elem(i)+e2t_v_elem(i)
 enddo
 if(perple_x.or.melting)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Thermal expansion" Format="ascii">'
   do i=1,elementi
     write(1,*) alpha(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Specific heat" Format="ascii">'
   do i=1,elementi
     write(1,*) Cp(i)
   enddo
 endif
 if(melting)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Melt fraction" Format="ascii">'
   do i=1,elementi
     write(1,*) melt_elem(i)
   enddo
 endif
 if(power_law)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Cohesion" Format="ascii">'
   do i=1,elementi
     write(1,*) cohe(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Angle of friction" Format="ascii">'
   do i=1,elementi
     write(1,*) teta(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Deformation" Format="ascii">'
   do i=1,elementi
     write(1,*) def_elem(i)
   enddo
 endif
 if(idratazione)then
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Water bound" Format="ascii">'
   do i=1,elementi
     write(1,*) w_b_elem(i)
   enddo
   write(1,*) '</DataArray>'
   write(1,*) '<DataArray type="Float32" Name="Water free" Format="ascii">'
   do i=1,elementi
     write(1,*) W_e(i)
   enddo
 endif
 write(1,*) '</DataArray>'
 write(1,*) '</CellData>'
 write(1,*) '<Cells>'
 write(1,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
 do i=1,elementi
   write(1,*) (connessioni(i,j)-1,j=1,nodi_elem)
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
 write(1,*) (i*4,i=1,elementi)
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="types" Format="ascii">'
 write(1,*) (9,i=1,elementi)
 write(1,*) '</DataArray>'
 write(1,*) '</Cells>'
 write(1,*) '</Piece>'
 write(1,*) '</UnstructuredGrid>'
 write(1,*) '</VTKFile>'
 close(1)

 if(markers_vtu_out)then
   write(filename,'("/vtu_files/Markers.",I5.5,".vtk")')name
   open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown')
   write(1,107) tempo
   write(1,101)
   write(1,109)
   write(1,102)
   write(1,111) marcatori
 endif
 if(markers_txt_out)then
   write(filename,'("/txt_files/Marcatori.",i5.5,".txt")')name
   open(unit=2,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown')
   write(2,108) tempo
   write(2,2003)
 endif
 do i=1,marcatori
   if(markers_vtu_out)write(1,*) x_marc(i),y_marc(i),0.0
   if(markers_txt_out)write(2,'(i10,2(3x,f20.9),1x,i5,1x,g15.5,1x,2(f15.5,1x),g15.5,3x,f8.2,1x,g15.5,7x,4(f6.2,7x))')n_marc(i),x_marc(i),y_marc(i),&
        type_marc(i),e2t_ist(i),e2t_p(i),e2t_v(i),Visc_m(i),Tmarc(i),Pmarc(i),W_b(i),W_f(i),W_m(i),melt_fraction(i)
 enddo
 if(markers_txt_out)close(2)
 if(markers_vtu_out)then
   write(1,103) marcatori,(marcatori)*2
   do i=1,marcatori
     write(1,*) vertype,i
   enddo
   write(1,117) marcatori
   write(1,105)
   write(1,119)
   do i=1,marcatori
     write(1,*) type_marc(i)
   enddo
   write(1,100)
   write(1,119)
   do i=1,marcatori
     write(1,*) VV(i)
   enddo
   if(solver_T)then
     write(1,118)
     write(1,119)
     do i=1,marcatori
       write(1,*) Tmarc(i)
     enddo
   endif
   write(1,121)
   write(1,119)
   do i=1,marcatori
     write(1,*) Pmarc(i)
   enddo
   write(1,128)
   write(1,119)
   do i=1,marcatori
     write(1,*) e2t_p(i)
   enddo
   write(1,129)
   write(1,119)
   do i=1,marcatori
     write(1,*) e2t_v(i)
   enddo
   write(1,130)
   write(1,119)
   do i=1,marcatori
     write(1,*) e2t_p(i)+e2t_v(i)
   enddo
   write(1,135)
   write(1,119)
   do i=1,marcatori
     write(1,*) e2t_ist(i)
   enddo
   write(1,136)
   write(1,119)
   do i=1,marcatori
     write(1,*) Visc_m(i)
   enddo
   if(idratazione)then
     write(1,131)
     write(1,119)
     do i=1,marcatori
       write(1,*) W_b(i)
     enddo
     write(1,132)
     write(1,119)
     do i=1,marcatori
       write(1,*) W_f(i)
     enddo
     write(1,134)
     write(1,119)
     do i=1,marcatori
       write(1,*) W_m(i)
     enddo
   endif
   if(melting)then
     write(1,133)
     write(1,119)
     do i=1,marcatori
       write(1,*) melt_fraction(i)
     enddo
   endif
   close(1)
 endif

 if(OUT_TYPE)nt=nt+1
 if(.not.OUT_TYPE.and.it.eq.nt)nt=nt+OUT_FILE_nb

100 format(26HSCALARS Deformazione float)
101 format(15H2D markers plot)  
102 format(16HDATASET POLYDATA)
103 format(8HVERTICES,1x,i10,1x,i10)
104 format('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">, Total time (Myr):',f11.6)
105 format(26HSCALARS Markers_Type int 1)
107 format(26H# vtk DataFile Version 1.0,' Total time (Myr):',f11.6)
108 format(' Total time (Myr):',f11.6)
109 format(5HASCII)
111 format(6HPOINTS,1x,i10,1x,5Hfloat)
117 format(10HPOINT_DATA,1x,i10)
118 format(25HSCALARS Temperature float)
119 format(20HLOOKUP_TABLE default)
121 format(22HSCALARS Pressure float)
128 format(28HSCALARS Strain_plastic float)
129 format(28HSCALARS Strain_viscous float)
130 format(26HSCALARS Strain_total float)
131 format(25HSCALARS Water_bound float)
132 format(24HSCALARS Water_free float)
133 format(27HSCALARS Melt_fraction float)
134 format(24HSCALARS Water_melt float)
135 format(25HSCALARS Strain_rate float)
136 format(23HSCALARS Viscosity float)
2001 format(2x,'Elemento',9x,'Coordinata X',9x,'Coordinata Y',9x,'Hr',7x,'Cond.',3x,'Densità',9x,'Cp',11x,'alpha',2x,'Pressione lit.',&
     2x,'Pressione tot.',2x,'Temperatura',7x,'Viscosità',3x,'Shear Heating',5x,'Ad. Heating',10x,'Stress',5x,'X Heat flux',5x,'Y Heat flux',&
     2x,'Plastic Strain',2x,'Viscous Strain',2x,'Type',2x,'Def.',6x,'Angle',8x,'Cohesion',3x,'Melt fraction')
2002 format(6x,'Nodo',9x,'Coordinata X',9x,'Coordinata Y',6x,'Velocità X',6x,'Velocità Y',5x,'Temperatura',2x,'Pressione tot.',5x,'Strain Rate')
2003 format(1x,'Marcatore',11x,'Coordinata X',11x,'Coordinata Y',2x,'Tipo',5x,'Strain Rate',2x,'Plastic strain'2x,'Viscous strain',7x,'Viscosity',2x,&
     'Temperatura',5x,'Pressione',2x,'Bound water',3x,'Free water',3x,'Melt water',4x,'Melt fraction')

 end subroutine Output

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Output_lgrid(it,tempo,x_lg,y_lg,OUT_TYPE,nodi,elementi,connessioni,dim_x_lg,dim_y_lg,em_lg)

 Implicit none

 logical :: inside(elementi)
 character(len=100) :: filename
 integer :: i,j,name,elem_g
 double precision :: dist,dist_max
 logical, intent(in) :: OUT_TYPE
 integer, intent(in) :: it,nodi,elementi,connessioni(elementi,4),em_lg(nodi)
 double precision, intent(in) :: x_lg(nodi),y_lg(nodi),tempo,dim_x_lg,dim_y_lg

 inside=.true.

 if(OUT_TYPE)name=FLOOR(tempo*1000.D0)
 if(.not.OUT_TYPE)name=it-1

 do i=1,elementi
   do j=1,nodi_elem
     if(em_lg(connessioni(i,j)).eq.0)then
!       inside(i)=.false.
     endif
     if(j.eq.1)dist=DSQRT((x_lg(connessioni(i,j))-x_lg(connessioni(i,3)))**2+(y_lg(connessioni(i,j))-y_lg(connessioni(i,3)))**2)
     if(j.eq.2)dist=DSQRT((x_lg(connessioni(i,j))-x_lg(connessioni(i,4)))**2+(y_lg(connessioni(i,j))-y_lg(connessioni(i,4)))**2)
     if(j.eq.3)dist=DSQRT((x_lg(connessioni(i,j))-x_lg(connessioni(i,1)))**2+(y_lg(connessioni(i,j))-y_lg(connessioni(i,1)))**2)
     if(j.eq.4)dist=DSQRT((x_lg(connessioni(i,j))-x_lg(connessioni(i,2)))**2+(y_lg(connessioni(i,j))-y_lg(connessioni(i,2)))**2)
     dist_max=10.D0*DSQRT(dim_x_lg**2+dim_y_lg**2)
     if(dist.gt.dist_max)then
       inside(i)=.false.
     endif
   enddo
 enddo
 elem_g=COUNT(inside)

 write(filename,'("/vtu_files/Grid.",I5.5,".vtu")')name
 open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//filename,status='unknown',form='formatted')
 rewind(1)

 write(1,104) tempo
 write(1,*) '<UnstructuredGrid>'
 write(1,*) '<Piece NumberOfPoints="',elem_g*nodi_elem,'" NumberOfCells="',elem_g,'">'

 write(1,*) '<Points>'
 write(1,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
 do i=1,elementi
   do j=1,nodi_elem
     if(inside(i))write(1,*) x_lg(connessioni(i,j)),y_lg(connessioni(i,j)),0.0
   enddo
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</Points>'

 write(1,*) '<Cells>'
 write(1,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
 do i=1,elem_g
   write(1,*) (i-1)*4,(i-1)*4+1,(i-1)*4+2,(i-1)*4+3
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
 do i=1,elem_g
   write(1,*) i*4
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '<DataArray type="Int32" Name="types" Format="ascii">'
 do i=1,elem_g
   write(1,*) 9
 enddo
 write(1,*) '</DataArray>'
 write(1,*) '</Cells>'
 write(1,*) '</Piece>'
 write(1,*) '</UnstructuredGrid>'
 write(1,*) '</VTKFile>'
 close(1)

104 format('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">, Total time (Myr):',f11.6)

 endsubroutine Output_lgrid

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Strain(e2t_ist,dt,VV,e2t_p,e2t_v,marcatori,healing,healing_rate,thermal_act,e1,e1v,T_heal)

 Implicit none

 integer :: i
 double precision :: heal,app_heal
 logical, intent(in) :: healing
 integer, intent(in) :: marcatori,VV(marcatori)
 double precision, intent(in) :: e2t_ist(marcatori),dt,healing_rate,thermal_act,e1,e1v,T_heal(marcatori)
 double precision, intent(inout) :: e2t_p(marcatori),e2t_v(marcatori)

 do i=1,marcatori
   if(healing)then
     heal=healing_rate*DEXP((-thermal_act/2.D0)*((1.D0/(T_heal(i)+1.D0))-0.5D0))
   else
     heal=0.D0
   endif
   if(VV(i).eq.2.or.VV(i).eq.3.or.VV(i).eq.4)then
     app_heal=heal*e2t_v(i)
     e2t_v(i)=e2t_v(i)+(e2t_ist(i)-app_heal)*dt
   elseif(VV(i).eq.5)then
     app_heal=heal*e2t_p(i)
     e2t_p(i)=e2t_p(i)+(e2t_ist(i)-app_heal)*dt
   endif
   if(e2t_p(i).gt.e1)e2t_p(i)=e1
   if(e2t_v(i).gt.e1v)e2t_v(i)=e1v
   if(e2t_p(i).lt.0.D0)e2t_p(i)=0.D0
   if(e2t_v(i).lt.0.D0)e2t_v(i)=0.D0
 enddo

 end subroutine Strain

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Upper_Boundary(x,y,nodi,nodi_x,nodi_y,Vx,Vy,dt,delta_x,delta_y,connessioni,nodi_elem,elementi,chain,x_chain,y_chain,&
 elem_x,elem_y,larghezza,altezza,h0,v0,eps,free_surface,erosione,aria,xmin_elem,ymin_elem,esp_m,esp_n,kappa_f,kappa_d,x_min,x_max,&
 y_min,y_max,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,kappa_fsed,kappa_dsed,&
 layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y,sed_layer)

 Implicit none
 
 integer :: i,j,k,colonna,riga
 double precision :: vx_chain(chain),vy_chain(chain),x_chain_temp(chain),y_chain_temp(chain),xc,yc,N(nodi_elem),x_temp(nodi_x),y_temp(nodi_x),y0(nodi_x),&
 y_up_tmp(nodi_x),y_uplift(nodi_x),dy,dy_sed(nodi_x)

 logical, intent(in) :: erosione,free_surface
 integer, intent(in) :: nodi,nodi_x,nodi_y,elementi,nodi_elem,connessioni(elementi,nodi_elem),chain,elem_x,elem_y,&
 layer_x,layer_y,n_elem_x(layer_x),n_elem_y(layer_y)
 
 double precision, intent(in) :: x(nodi),Vx(nodi),Vy(nodi),dt,larghezza,altezza,delta_x(elementi),eps,h0,v0,aria,xmin_elem(elementi),&
 ymin_elem(elementi),esp_m,esp_n,kappa_f,kappa_d,x_min,x_max,y_min,y_max,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,&
 zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,kappa_fsed,kappa_dsed,thick_x(layer_x),thick_y(layer_y)

 logical, intent(out) :: sed_layer(elem_x)
 double precision, intent(inout) :: x_chain(chain),y_chain(chain),y(nodi),delta_y(elementi)

 sed_layer=.false.
 if(free_surface)then
   do i=nodi_x*(nodi_y-1)+1,nodi
     y_up_tmp(i-nodi_x*(nodi_y-1))=y(i)+Vy(i)*dt
     if(.not.erosione.or.equilib)y_temp(i-nodi_x*(nodi_y-1))=y_up_tmp(i-nodi_x*(nodi_y-1))
     x_temp(i-nodi_x*(nodi_y-1))=x(i)+Vx(i)*dt
     y0(i-nodi_x*(nodi_y-1))=y(i)
   enddo
   if(erosione.and..not.equilib)call Erosion(nodi_x,3,h0,larghezza*h0,2.D03,dt*t0/(365.D0*24.D0*3600.D0),x_temp*h0,y(nodi_x*(nodi_y-1)+1:nodi)*h0,&
                                     Vy(nodi_x*(nodi_y-1)+1:nodi)*v0*365.D0*24.D0*3600.D0,kappa_d,kappa_f,g_coeff,gsed_coeff,p_coeff,sealevel,&
                                     poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,esp_m,esp_n,kappa_fsed,kappa_dsed,y_temp)
   do i=2,nodi_x-1
     if(x_temp(i).ge.x(i+nodi_x*(nodi_y-1)))then
       y(i+nodi_x*(nodi_y-1))=(y_temp(i)-y_temp(i-1))*((x(i+nodi_x*(nodi_y-1))-x_temp(i-1))/(x_temp(i)-x_temp(i-1)))+y_temp(i-1)
       y_uplift(i)=(y_up_tmp(i)-y_up_tmp(i-1))*((x(i+nodi_x*(nodi_y-1))-x_temp(i-1))/(x_temp(i)-x_temp(i-1)))+y_up_tmp(i-1)
     elseif(x_temp(i).lt.x(i+nodi_x*(nodi_y-1)))then
       y(i+nodi_x*(nodi_y-1))=(y_temp(i+1)-y_temp(i))*((x(i+nodi_x*(nodi_y-1))-x_temp(i))/(x_temp(i+1)-x_temp(i)))+y_temp(i)
       y_uplift(i)=(y_up_tmp(i+1)-y_up_tmp(i))*((x(i+nodi_x*(nodi_y-1))-x_temp(i))/(x_temp(i+1)-x_temp(i)))+y_up_tmp(i)
     endif
   enddo
   y(nodi_x*(nodi_y-1)+1)=y_temp(1)
   y_uplift(1)=y_temp(1)
   y(nodi)=y_temp(nodi_x)
   y_uplift(nodi_x)=y_temp(nodi_x)
   do i=1,nodi_x
     dy=y(i+nodi_x*(nodi_y-1))-y0(i)
     dy_sed(i)=y(i+nodi_x*(nodi_y-1))-y_uplift(i)
     do j=2,nodi_y-1
       y(nodi_x*(j-1)+i)=y(nodi_x*(j-1)+i)+dy*((y(nodi_x*(j-1)+i)-y(i))/(y(nodi_x*(nodi_y-1)+i)-y(i)))
     enddo
   enddo
   do j=1,elem_y
     do i=1,elem_x
       if(j.eq.elem_y.and.(dy_sed(i)+dy_sed(i+1))/2.D0.gt.0.D0)sed_layer(i)=.true.
       indice=elem_x*(j-1)+i
       delta_y(indice)=DABS((y(connessioni(indice,4))-y(connessioni(indice,1)))+(y(connessioni(indice,3))-y(connessioni(indice,2))))/2.D0
     enddo
   enddo
 else
   vx_chain=0.D0
   vy_chain=0.D0
   do i=1,chain
     call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_chain(i),y_chain(i),j,x_min,x_max,&
          y_min,y_max,larghezza,altezza,free_surface,colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
     xc=(2.D0/delta_x(j))*(x_chain(i)-xmin_elem(j))-1.D0
     yc=(2.D0/delta_y(j))*(y_chain(i)-ymin_elem(j))-1.D0
     N=(/0.25D0*(1.D0-xc)*(1.D0-yc),0.25D0*(1.D0+xc)*(1.D0-yc),&
         0.25D0*(1.D0+xc)*(1.D0+yc),0.25D0*(1.D0-xc)*(1.D0+yc)/)
     do k=1,nodi_elem
       vx_chain(i)=vx_chain(i)+Vx(connessioni(j,k))*N(k)
       vy_chain(i)=vy_chain(i)+Vy(connessioni(j,k))*N(k)
     enddo
     if(x_chain(i).le.(10.D0*delta_x(1)).or.x_chain(i).ge.larghezza-(10.D0*delta_x(elem_x)))vy_chain(i)=0.D0
     if(.not.erosione)y_chain_temp(i)=y_chain(i)+vy_chain(i)*dt
     x_chain_temp(i)=x_chain(i)+vx_chain(i)*dt
     if(x_chain_temp(i).le.x_chain_temp(i-1))x_chain_temp(i)=x_chain_temp(i-1)+eps
   enddo
   if(erosione)call Erosion(chain,3,h0,larghezza*h0,2.D03,dt*t0/(365.D0*24.D0*3600.D0),x_chain*h0,y_chain,vy_chain*v0*365.D0*24.D0*3600.D0,&
                    kappa_d,kappa_f,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,esp_m,&
                    esp_n,kappa_fsed,kappa_dsed,y_chain_temp)
   do i=1,chain
     if(x_chain(i).lt.x_chain_temp(1))then
       y_chain(i)=(y_chain_temp(2)-y_chain_temp(1))*((x_chain(i)-x_chain_temp(1))/(x_chain_temp(2)-&
                  x_chain_temp(1)))+y_chain_temp(1)
     elseif(x_chain(i).ge.x_chain_temp(chain))then
       y_chain(i)=(y_chain_temp(chain)-y_chain_temp(chain-1))*((x_chain(i)-x_chain_temp(chain-1))/(x_chain_temp(chain)-&
                  x_chain_temp(chain-1)))+y_chain_temp(chain-1)
     else
       do j=1,chain-1
         if(x_chain(i).ge.x_chain_temp(j).and.x_chain(i).lt.x_chain_temp(j+1))then
           y_chain(i)=(y_chain_temp(j+1)-y_chain_temp(j))*((x_chain(i)-x_chain_temp(j))/(x_chain_temp(j+1)-&
                      x_chain_temp(j)))+y_chain_temp(j)
           exit
         endif
       enddo
     endif
   enddo
 endif

 end subroutine Upper_Boundary

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Erosion(nx,ny,h0,xl,yl,dt,xt,ht,ut,kappa_d,kappa_f,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,&
 poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,esp_m,esp_n,kfsed,kdsed,h_out)

 Implicit none

 double precision, parameter :: dtmax=2.5D03
 integer :: istep,nstep
 double precision :: dt0,x(nx*ny),y(nx*ny),h(nx*ny),u(nx*ny),kf(nx*ny),kd(nx*ny)
 integer, intent(in) :: nx,ny
 double precision, intent(in) :: h0,xl,yl,dt,xt(nx),ht(nx),ut(nx),kfsed,esp_m,esp_n,kdsed,g_coeff,gsed_coeff,p_coeff,&
 sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,kappa_d,kappa_f

 double precision, intent(out) :: h_out(nx)

 nstep=INT(dt/dtmax)
 dt0=dt/(nstep+1)

 x(1:nx)=xt
 x(nx+1:nx*(ny-1))=xt
 x(nx*(ny-1)+1:nx*ny)=xt
 y(1:nx)=-1.D03
 y(nx+1:nx*(ny-1))=0.D03
 y(nx*(ny-1)+1:nx*ny)=1.D03
 h(1:nx)=ht
 h(nx+1:nx*(ny-1))=ht
 h(nx*(ny-1)+1:nx*ny)=ht
 u(1:nx)=ut
 u(nx+1:nx*(ny-1))=ut
 u(nx*(ny-1)+1:nx*ny)=ut
 kf=kappa_f
 kd=kappa_d

 call FastScape_Init()

 call FastScape_Set_NX_NY(nx,ny)
 call FastScape_Setup()
 call FastScape_Set_XL_YL(xl,yl)
 call FastScape_Set_DT(dt0)

 call FastScape_Init_H(h)
 call FastScape_Set_Erosional_Parameters(kf,kfsed,esp_m,esp_n,kd,kdsed,g_coeff,gsed_coeff,p_coeff)
 call FastScape_Set_Marine_Parameters(sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh)
 call FastScape_Set_U(u)
 call FastScape_Set_BC(1111)

 call FastScape_Get_Step(istep)

 do while (istep.lt.nstep+1)
   call FastScape_Execute_Step()
   call FastScape_Get_Step(istep)
   call FastScape_Copy_h(h)
 enddo
 h_out=h(nx+1:nx*(ny-1))/h0

! call FastScape_Debug()
 call FastScape_Destroy()

 end subroutine Erosion

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Advezione(elementi,elem_x,elem_y,nodi,nodi_elem,x,y,x_marc,y_marc,n_marc,marcatori,marcatoritot,elem_marc,Vx,Vy,connessioni,dt,C,nmax,larghezza,&
 altezza,x_min,x_max,y_min,y_max,ndime,vel_corr,runge,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y,inside,marc_in,marc_each_elem,marc_agg,e2t_p,e2t_v,&
 type_marc,W_b,W_f,W_m,melt_fraction,free_surface,erosione,chain,x_chain,y_chain,t_air,t_sed,n_marc_tmp,type_marc_tmp,elem_marc_tmp,x_marc_tmp,y_marc_tmp,&
 e2t_p_tmp,e2t_v_tmp,W_b_tmp,W_f_tmp,W_m_tmp,melt_fraction_tmp,colonna_marc,colonna_marc_tmp)

 Implicit none

 integer :: i,colonna,riga,col_chain
 double precision :: ymin
 logical, intent(in) :: vel_corr,runge,free_surface,erosione
 integer, intent(in) :: elementi,nodi_elem,nodi,connessioni(elementi,nodi_elem),elem_x,elem_y,ndime,layer_x,layer_y,n_elem_x(layer_x),&
 n_elem_y(layer_y),marcatori,marcatoritot,n_marc(marcatori),nmax(elementi),chain,t_air,t_sed

 double precision, intent(in) :: Vx(nodi),Vy(nodi),dt,x(nodi),y(nodi),larghezza,altezza,x_min,x_max,y_min,y_max,thick_x(layer_x),thick_y(layer_y),&
 W_b(marcatori),W_f(marcatori),W_m(marcatori),melt_fraction(marcatori),x_chain(chain),y_chain(chain)

 logical, intent(out) :: inside(marcatori)
 integer, intent(out) :: C(elementi),marc_in,marc_agg,marc_each_elem(elementi,MAXVAL(nmax*2)),n_marc_tmp(marcatori),elem_marc_tmp(marcatoritot),&
 type_marc_tmp(marcatori),colonna_marc_tmp(marcatoritot)

 double precision, intent(out) :: x_marc_tmp(marcatori),y_marc_tmp(marcatori),e2t_p_tmp(marcatori),e2t_v_tmp(marcatori),W_b_tmp(marcatori),&
 W_f_tmp(marcatori),W_m_tmp(marcatori),melt_fraction_tmp(marcatori)

 integer, intent(inout) :: elem_marc(marcatoritot),type_marc(marcatori),colonna_marc(marcatoritot)
 double precision, intent(inout) :: x_marc(marcatori),y_marc(marcatori),e2t_p(marcatori),e2t_v(marcatori)

 C=0
 marc_in=0
 marc_agg=0
 marc_each_elem=0
 inside=.false.
 do i=1,marcatori
   if(runge)then
     call Runge_4th(x_marc(i),y_marc(i),elem_marc(n_marc(i)),ndime,elementi,nodi_elem,connessioni,x,y,nodi,vel_corr,Vx,Vy,dt,elem_x,elem_y,&
          x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
   else
     call Runge_2nd(x_marc(i),y_marc(i),elem_marc(n_marc(i)),ndime,elementi,nodi_elem,connessioni,x,y,nodi,vel_corr,Vx,Vy,dt,elem_x,elem_y,&
          x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
   endif
   if(.not.free_surface.and.erosione)then
     col_chain=INT((x_marc(i)/larghezza)*(chain-1))+1
     if(col_chain.gt.chain-1)col_chain=chain-1
     ymin=(y_chain(col_chain+1)-y_chain(col_chain))*((x_marc(i)-x_chain(col_chain))/(x_chain(col_chain+1)-&
          x_chain(col_chain)))+y_chain(col_chain)
     if(type_marc(i).eq.t_air.and.y_marc(i).lt.ymin)type_marc(i)=t_sed
     if(type_marc(i).ne.t_air.and.y_marc(i).gt.ymin)type_marc(i)=t_air
   endif
   call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_marc(i),y_marc(i),elem_marc(n_marc(i)),&
        x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,colonna_marc(n_marc(i)),riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
   if(elem_marc(n_marc(i)).ne.0)then
     C(elem_marc(n_marc(i)))=C(elem_marc(n_marc(i)))+1
     marc_in=marc_in+1
     inside(i)=.true.
     if(check_markers)marc_each_elem(elem_marc(n_marc(i)),C(elem_marc(n_marc(i))))=i
   endif
   n_marc_tmp(i)=n_marc(i)
   type_marc_tmp(i)=type_marc(i)
   elem_marc_tmp(i)=elem_marc(n_marc(i))
   colonna_marc_tmp(i)=colonna_marc(n_marc(i))
   x_marc_tmp(i)=x_marc(i)
   y_marc_tmp(i)=y_marc(i)
   e2t_p_tmp(i)=e2t_p(i)
   e2t_v_tmp(i)=e2t_v(i)
   W_b_tmp(i)=W_b(i)
   W_f_tmp(i)=W_f(i)
   W_m_tmp(i)=W_m(i)
   melt_fraction_tmp(i)=melt_fraction(i)
 enddo

 end subroutine Advezione

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 Subroutine Runge_4th(x_marc,y_marc,elem_marc,ndime,elementi,nodi_elem,connessioni,x,y,nodi,vel_corr,Vx,Vy,dt,elem_x,elem_y,&
 x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)

 Implicit none

 integer :: i,j,e_m,colonna,riga
 double precision :: vx_rg(4),vy_rg(4),x_m,y_m,posgp(ndime),corr_u,corr_v,N(nodi_elem),vx_marc,vy_marc
 logical, intent(in) :: vel_corr,free_surface
 integer, intent(in) :: elem_marc,ndime,elementi,nodi_elem,nodi,connessioni(elementi,nodi_elem),elem_x,elem_y,layer_x,layer_y,&
 n_elem_x(layer_x),n_elem_y(layer_y)
 double precision, intent(in) :: x(nodi),y(nodi),Vx(nodi),Vy(nodi),dt,x_min,x_max,y_min,y_max,larghezza,altezza,thick_x(layer_x),thick_y(layer_y)
 double precision, intent(inout) :: x_marc,y_marc

 vx_rg=0.D0
 vy_rg=0.D0
 x_m=x_marc
 y_m=y_marc
 e_m=elem_marc
 do i=1,4
   posgp(1)=(2.D0/delta_x(e_m))*(x_m-xmin_elem(e_m))-1.D0
   posgp(2)=(4*y_m-((1.D0-posgp(1))*y(connessioni(e_m,1))+(1.D0+posgp(1))*y(connessioni(e_m,2))+&
          (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4))))/&
          (-(1.D0-posgp(1))*y(connessioni(e_m,1))-(1.D0+posgp(1))*y(connessioni(e_m,2))+&
          (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4)))
   if(posgp(1).lt.-1.D0.or.posgp(1).gt.1.D0.or.posgp(2).lt.-1.D0.or.posgp(2).gt.1.D0)then
     call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_m,y_m,e_m,x_min,x_max,y_min,y_max,&
          larghezza,altezza,free_surface,colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
     posgp(1)=(2.D0/delta_x(e_m))*(x_m-xmin_elem(e_m))-1.D0
     posgp(2)=(4*y_m-((1.D0-posgp(1))*y(connessioni(e_m,1))+(1.D0+posgp(1))*y(connessioni(e_m,2))+&
              (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4))))/&
              (-(1.D0-posgp(1))*y(connessioni(e_m,1))-(1.D0+posgp(1))*y(connessioni(e_m,2))+&
              (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4)))
   endif
   N=(/0.25D0*(1.D0-posgp(1))*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(1))*(1.D0-posgp(2)),&
       0.25D0*(1.D0+posgp(1))*(1.D0+posgp(2)),0.25D0*(1.D0-posgp(1))*(1.D0+posgp(2))/)
   corr_u=0.D0
   corr_v=0.D0
   if(vel_corr)call Correzione_Velocita(posgp,nodi_elem,ndime,x(connessioni(e_m,:)),y(connessioni(e_m,:)),&
                    Vx(connessioni(e_m,:)),Vy(connessioni(e_m,:)),corr_u,corr_v)
   do j=1,nodi_elem
     vx_rg(i)=vx_rg(i)+N(j)*Vx(connessioni(e_m,j))
     vy_rg(i)=vy_rg(i)+N(j)*Vy(connessioni(e_m,j))
   enddo
   vx_rg(i)=vx_rg(i)+corr_u
   vy_rg(i)=vy_rg(i)+corr_v
   if(i.eq.1.or.i.eq.2)then
     x_m=x_marc+vx_rg(i)*dt/2.D0
     y_m=y_marc+vy_rg(i)*dt/2.D0
   elseif(i.eq.3)then
     x_m=x_marc+vx_rg(i)*dt
     y_m=y_marc+vy_rg(i)*dt
   endif
 enddo
 vx_marc=1.D0/6.D0*(vx_rg(1)+2.D0*vx_rg(2)+2.D0*vx_rg(3)+vx_rg(4))
 vy_marc=1.D0/6.D0*(vy_rg(1)+2.D0*vy_rg(2)+2.D0*vy_rg(3)+vy_rg(4))
 x_marc=x_marc+vx_marc*dt
 y_marc=y_marc+vy_marc*dt

 end subroutine Runge_4th

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 Subroutine Runge_2nd(x_marc,y_marc,elem_marc,ndime,elementi,nodi_elem,connessioni,x,y,nodi,vel_corr,Vx,Vy,dt,elem_x,elem_y,&
 x_min,x_max,y_min,y_max,larghezza,altezza,free_surface,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)

 Implicit none

 integer :: i,j,e_m,colonna,riga
 double precision :: vx_rg(2),vy_rg(2),x_m,y_m,posgp(ndime),corr_u,corr_v,N(nodi_elem)
 logical, intent(in) :: vel_corr,free_surface
 integer, intent(in) :: elem_marc,ndime,elementi,nodi_elem,nodi,connessioni(elementi,nodi_elem),elem_x,elem_y,layer_x,layer_y,&
 n_elem_x(layer_x),n_elem_y(layer_y)
 double precision, intent(in) :: x(nodi),y(nodi),Vx(nodi),Vy(nodi),dt,x_min,x_max,y_min,y_max,larghezza,altezza,thick_x(layer_x),thick_y(layer_y)
 double precision, intent(inout) :: x_marc,y_marc

 vx_rg=0.D0
 vy_rg=0.D0
 x_m=x_marc
 y_m=y_marc
 e_m=elem_marc
 do i=1,2
   posgp(1)=(2.D0/delta_x(e_m))*(x_m-xmin_elem(e_m))-1.D0
   posgp(2)=(4*y_m-((1.D0-posgp(1))*y(connessioni(e_m,1))+(1.D0+posgp(1))*y(connessioni(e_m,2))+&
          (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4))))/&
          (-(1.D0-posgp(1))*y(connessioni(e_m,1))-(1.D0+posgp(1))*y(connessioni(e_m,2))+&
          (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4)))
   if(posgp(1).lt.-1.D0.or.posgp(1).gt.1.D0.or.posgp(2).lt.-1.D0.or.posgp(2).gt.1.D0)then
     call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_m,y_m,e_m,x_min,x_max,y_min,y_max,&
          larghezza,altezza,free_surface,colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
     posgp(1)=(2.D0/delta_x(e_m))*(x_m-xmin_elem(e_m))-1.D0
     posgp(2)=(4*y_m-((1.D0-posgp(1))*y(connessioni(e_m,1))+(1.D0+posgp(1))*y(connessioni(e_m,2))+&
              (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4))))/&
              (-(1.D0-posgp(1))*y(connessioni(e_m,1))-(1.D0+posgp(1))*y(connessioni(e_m,2))+&
              (1.D0+posgp(1))*y(connessioni(e_m,3))+(1.D0-posgp(1))*y(connessioni(e_m,4)))
   endif
   N=(/0.25D0*(1.D0-posgp(1))*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(1))*(1.D0-posgp(2)),&
       0.25D0*(1.D0+posgp(1))*(1.D0+posgp(2)),0.25D0*(1.D0-posgp(1))*(1.D0+posgp(2))/)
   corr_u=0.D0
   corr_v=0.D0
   if(vel_corr)call Correzione_Velocita(posgp,nodi_elem,ndime,x(connessioni(e_m,:)),y(connessioni(e_m,:)),&
                    Vx(connessioni(e_m,:)),Vy(connessioni(e_m,:)),corr_u,corr_v)
   do j=1,nodi_elem
     vx_rg(i)=vx_rg(i)+N(j)*Vx(connessioni(e_m,j))
     vy_rg(i)=vy_rg(i)+N(j)*Vy(connessioni(e_m,j))
   enddo
   vx_rg(i)=vx_rg(i)+corr_u
   vy_rg(i)=vy_rg(i)+corr_v
   if(i.eq.1)then
     x_m=x_marc+vx_rg(i)*dt/2.D0
     y_m=y_marc+vy_rg(i)*dt/2.D0
   endif
 enddo
 x_marc=x_marc+vx_rg(2)*dt
 y_marc=y_marc+vy_rg(2)*dt

 end subroutine Runge_2nd

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 Subroutine Correzione_Velocita(posgp,nodi_elem,ndime,x_conn,y_conn,u_conn,v_conn,corr_u,corr_v)

 Implicit none

 double precision :: coord(ndime,nodi_elem),dNst(nodi_elem,ndime),jacob(ndime,ndime),u14,u23,u12,u34,v14,v23,v12,v34,corr_b,corr_d
 integer, intent(in) :: nodi_elem,ndime
 double precision, intent(in) :: x_conn(nodi_elem),y_conn(nodi_elem),posgp(ndime),u_conn(nodi_elem),v_conn(nodi_elem)
 double precision, intent(out) :: corr_u,corr_v

 dNst=RESHAPE((/-0.25D0*(1.D0-posgp(2)),0.25D0*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(2)),&
      -0.25D0*(1.D0+posgp(2)),-0.25D0*(1.D0-posgp(1)),-0.25D0*(1.D0+posgp(1)),&
      0.25D0*(1.D0+posgp(1)),0.25D0*(1.D0-posgp(1))/),SHAPE=(/4,2/))
 coord(1,:)=x_conn(:)
 coord(2,:)=y_conn(:)
 jacob=matmul(coord,dNst)
 u14=0.25D0*(u_conn(1)-u_conn(4))
 u23=0.25D0*(u_conn(2)-u_conn(3))
 u12=0.25D0*(u_conn(1)-u_conn(2))
 u34=0.25D0*(u_conn(3)-u_conn(4))
 v14=0.25D0*(v_conn(1)-v_conn(4))
 v23=0.25D0*(v_conn(2)-v_conn(3))
 v12=0.25D0*(v_conn(1)-v_conn(2))
 v34=0.25D0*(v_conn(3)-v_conn(4))
 corr_b=(1.D0/(2.D0*jacob(1,1)))*(jacob(1,2)*(u14-u23)+jacob(2,2)*(v14-v23))
 corr_d=(1.D0/(2.D0*jacob(2,2)))*(jacob(1,1)*(u12+u34)+jacob(2,1)*(v12+v34))
 corr_u=corr_b*(1.D0-posgp(1))*(1.D0+posgp(1))
 corr_v=corr_d*(1.D0-posgp(2))*(1.D0+posgp(2))

 end subroutine Correzione_Velocita

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Add_remove_markers(free_surface,elementi,nmin,nmax,nodi,nodi_elem,connessioni,marc_in_elem,marcatori,marcatoritot,n_marc,type_marc,&
 delta_x,xmin_elem,ymin_elem,ymax_elem,x,y,x_marc,y_marc,larghezza,e2t_p,e2t_v,W_b,W_f,W_m,melt_fraction,C,n_marc_a,type_marc_a,elem_marc_a,&
 x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_f_a,W_m_a,melt_fraction_a,marc_in,marc_agg,erosione,elem_x,elem_y,t_sed,equilib,sed_layer,&
 colonna_marc,colonna_marc_a)

 Implicit none

 logical :: sed_mark
 integer :: i,j,ii,jj,k,kk,pos,in_elem,missed,plus,flag
 double precision :: dist,distmin,xx,yy,ymin,ymax,xmarc,ymarc,y_sed,dy
 logical, intent(in) :: free_surface,erosione,equilib,sed_layer(elem_x)
 integer, intent(in) :: elementi,nmin(elementi),nmax(elementi),nodi,nodi_elem,connessioni(elementi,nodi_elem),marc_in_elem(elementi,MAXVAL(nmax*2)),&
 marcatori,marcatoritot,n_marc(marcatori),type_marc(marcatori),elem_x,elem_y,t_sed,colonna_marc(marcatoritot)

 double precision :: distm(MAXVAL(nmax)),xtmp(MAXVAL(nmax)),ytmp(MAXVAL(nmax))
 double precision, intent(in) :: delta_x(elementi),xmin_elem(elementi),ymin_elem(elementi),ymax_elem(elementi),x(nodi),y(nodi),x_marc(marcatori),&
 y_marc(marcatori),larghezza,e2t_p(marcatori),e2t_v(marcatori),W_b(marcatori),W_f(marcatori),W_m(marcatori),melt_fraction(marcatori)

 integer, intent(out) :: n_marc_a(SUM(nmin)),type_marc_a(SUM(nmin)),elem_marc_a(SUM(nmin)),colonna_marc_a(SUM(nmin))
 double precision, intent(out) :: x_marc_a(SUM(nmin)),y_marc_a(SUM(nmin)),e2t_p_a(SUM(nmin)),e2t_v_a(SUM(nmin)),W_b_a(SUM(nmin)),W_f_a(SUM(nmin)),&
 W_m_a(SUM(nmin)),melt_fraction_a(SUM(nmin))

 integer, intent(inout) :: C(elementi),marc_in,marc_agg

 do jj=1,elem_y
   do ii=1,elem_x 
     i=elem_x*(jj-1)+ii
     if(C(i).lt.nmin(i))then
       in_elem=C(i)
       missed=nmin(i)-C(i)
       y_sed=MIN(y(connessioni(i,4))-y(connessioni(i,1)),y(connessioni(i,3))-y(connessioni(i,2)))
       do k=1,in_elem
         if(type_marc(marc_in_elem(i,k)).ne.t_sed)then
           dy=((y(connessioni(i,3))-y(connessioni(i,4)))*((x_marc(marc_in_elem(i,k))-x(connessioni(i,4)))/(x(connessioni(i,3))-&
                x(connessioni(i,4))))+y(connessioni(i,4)))-y_marc(marc_in_elem(i,k))
         if(dy.le.y_sed)y_sed=dy
         endif
       enddo
       do j=1,missed
         marc_agg=marc_agg+1
         C(i)=C(i)+1
         sed_mark=.false.
         call random_number(xx)
         call random_number(yy)
         xmarc=xx*delta_x(i)+xmin_elem(i)
         if(free_surface)then
           ymax=(y(connessioni(i,3))-y(connessioni(i,4)))*((xmarc-x(connessioni(i,4)))/(x(connessioni(i,3))-&
                 x(connessioni(i,4))))+y(connessioni(i,4))
           ymin=(y(connessioni(i,2))-y(connessioni(i,1)))*((xmarc-x(connessioni(i,1)))/(x(connessioni(i,2))-&
                 x(connessioni(i,1))))+y(connessioni(i,1))
         else
           ymin=ymin_elem(i)
           ymax=ymax_elem(i)
         endif
         if(free_surface.and.erosione.and.sed_layer(ii).and.jj.eq.elem_y.and.ii.ne.1.and.ii.ne.elem_x.and.j.le.&
            CEILING(DBLE(missed)*y_sed/(ymax-ymin)))then
           ymarc=yy*y_sed+(ymax-y_sed)
           sed_mark=.true.
         else
           ymarc=yy*(ymax-ymin)+ymin
         endif
         distmin=larghezza
         do k=1,in_elem
           dist=DSQRT(((xmarc-x_marc(marc_in_elem(i,k)))**2)+((ymarc-y_marc(marc_in_elem(i,k)))**2))
           if(dist.lt.distmin)then
             distmin=dist
             pos=k
           endif
         enddo
         n_marc_a(marc_agg)=marc_agg+marcatoritot
         x_marc_a(marc_agg)=xmarc
         y_marc_a(marc_agg)=ymarc
         elem_marc_a(marc_agg)=i
         colonna_marc_a(marc_agg)=ii
         W_b_a(marc_agg)=W_b(marc_in_elem(i,pos))
         W_m_a(marc_agg)=W_m(marc_in_elem(i,pos))
         W_f_a(marc_agg)=W_f(marc_in_elem(i,pos))
         e2t_p_a(marc_agg)=e2t_p(marc_in_elem(i,pos))
         e2t_v_a(marc_agg)=e2t_v(marc_in_elem(i,pos))
         melt_fraction_a(marc_agg)=melt_fraction(marc_in_elem(i,pos))
         if(sed_mark)then
           type_marc_a(marc_agg)=t_sed
         else
           type_marc_a(marc_agg)=type_marc(marc_in_elem(i,pos))
         endif
       enddo
     elseif(C(i).gt.nmax(i))then
       in_elem=C(i)
       plus=C(i)-nmax(i)
       do j=1,plus
         distmin=larghezza
         do k=1,in_elem
           if(inside(marc_in_elem(i,k)))then
             dist=0.D0
             do kk=1,in_elem
               if(k.ne.kk.and.inside(marc_in_elem(i,kk)))dist=dist+DSQRT(((x_marc(marc_in_elem(i,k))-x_marc(marc_in_elem(i,kk)))**2)+&
                 ((y_marc(marc_in_elem(i,k))-y_marc(marc_in_elem(i,kk)))**2))
             enddo
             if(dist.lt.distmin)then
               distmin=dist
               flag=marc_in_elem(i,k)
             endif
           endif
         enddo
         inside(flag)=.false.
         C(i)=C(i)-1
         marc_in=marc_in-1
       enddo
     endif
   enddo
 enddo

 end subroutine Add_remove_markers

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Refine_markers(x_marc,y_marc,e2t_p,e2t_v,n_marc,marcatori,marcatoritot,elem_marc,type_marc,W_b,W_m,W_f,melt_fraction,inside,&
  marc_in,marc_agg,n_marc_a,type_marc_a,elem_marc_a,x_marc_a,y_marc_a,e2t_p_a,e2t_v_a,W_b_a,W_f_a,W_m_a,melt_fraction_a,n_marc_tmp,type_marc_tmp,&
  elem_marc_tmp,x_marc_tmp,y_marc_tmp,e2t_p_tmp,e2t_v_tmp,W_b_tmp,W_f_tmp,W_m_tmp,melt_fraction_tmp,colonna_marc_a,colonna_marc_tmp)

 Implicit none

 integer :: i,j
 logical, intent(in) :: inside(marcatori)
 integer, intent(in) :: marc_in,marc_agg,n_marc_a(SUM(nmin)),type_marc_a(SUM(nmin)),elem_marc_a(SUM(nmin)),colonna_marc_a(SUM(nmin)),&
 n_marc_tmp(marcatori),type_marc_tmp(marcatori),elem_marc_tmp(marcatoritot),colonna_marc_tmp(marcatoritot)

 double precision, intent(in) :: x_marc_tmp(marcatori),y_marc_tmp(marcatori),e2t_p_tmp(marcatori),e2t_v_tmp(marcatori),W_b_tmp(marcatori),&
 x_marc_a(SUM(nmin)),y_marc_a(SUM(nmin)),e2t_p_a(SUM(nmin)),e2t_v_a(SUM(nmin)),W_b_a(SUM(nmin)),W_m_tmp(marcatori),W_m_a(SUM(nmin)),&
 melt_fraction_tmp(marcatori),melt_fraction_a(SUM(nmin)),W_f_tmp(marcatori),W_f_a(SUM(nmin))

 integer, intent(out) :: n_marc(marc_in+marc_agg),type_marc(marc_in+marc_agg),elem_marc(marcatoritot+marc_agg)
 double precision, intent(out) :: x_marc(marc_in+marc_agg),y_marc(marc_in+marc_agg),e2t_p(marc_in+marc_agg),e2t_v(marc_in+marc_agg),W_b(marc_in+marc_agg),&
 W_m(marc_in+marc_agg),W_f(marc_in+marc_agg),melt_fraction(marc_in+marc_agg)

 integer, intent(inout) :: marcatori,marcatoritot

 j=0
 do i=1,marcatori
   if(inside(i))then
     j=j+1
     n_marc(j)=n_marc_tmp(i)
     x_marc(j)=x_marc_tmp(i)
     y_marc(j)=y_marc_tmp(i)
     type_marc(j)=type_marc_tmp(i)
     e2t_p(j)=e2t_p_tmp(i)
     e2t_v(j)=e2t_v_tmp(i)
     elem_marc(n_marc(j))=elem_marc_tmp(i)
     colonna_marc(n_marc(j))=colonna_marc_tmp(i)
     W_b(j)=W_b_tmp(i)
     W_m(j)=W_m_tmp(i)
     W_f(j)=W_f_tmp(i)
     melt_fraction(j)=melt_fraction_tmp(i)
   endif
 enddo
 do i=1,marc_agg
   n_marc(i+marc_in)=n_marc_a(i)
   x_marc(i+marc_in)=x_marc_a(i)
   y_marc(i+marc_in)=y_marc_a(i)
   type_marc(i+marc_in)=type_marc_a(i)
   e2t_p(i+marc_in)=e2t_p_a(i)
   e2t_v(i+marc_in)=e2t_v_a(i)
   elem_marc(n_marc(i+marc_in))=elem_marc_a(i)
   colonna_marc(n_marc(i+marc_in))=colonna_marc_a(i)
   W_b(i+marc_in)=W_b_a(i)
   W_m(i+marc_in)=W_m_a(i)
   W_f(i+marc_in)=W_f_a(i)
   melt_fraction(i+marc_in)=melt_fraction_a(i)
 enddo 
 marcatori=marc_in+marc_agg
 marcatoritot=marcatoritot+marc_agg

 end subroutine Refine_markers

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Advection_lgrid(free_surface,nodi_lg,elementi,nodi,nodi_elem,connessioni,elem_x,elem_y,layer_x,layer_y,n_elem_x,n_elem_y,x,y,&
 x_min,x_max,y_min,y_max,larghezza,altezza,thick_x,thick_y,Vx,Vy,dt,delta_x,xmin_elem,x_lg,y_lg,e_m,elem_lg,conn_lg,eps,nodi_x_lg,nodi_y_lg,dim_y_lg)

 Implicit none

 integer :: i,j,k,kk,colonna,riga,vx_rg,vy_rg
 double precision :: posgp(2),N(nodi_elem)
 logical, intent(in) :: free_surface
 integer, intent(in) :: nodi_lg,elementi,nodi,nodi_elem,connessioni(elementi,nodi_elem),elem_x,elem_y,layer_x,layer_y,&
 n_elem_x(layer_x),n_elem_y(layer_y),elem_lg,conn_lg(elem_lg,4),nodi_x_lg,nodi_y_lg

 double precision, intent(in) :: x(nodi),y(nodi),x_min,x_max,y_min,y_max,larghezza,altezza,thick_x(layer_x),thick_y(layer_y),&
 Vx(nodi),Vy(nodi),dt,delta_x(elementi),xmin_elem(elementi),eps,dim_y_lg

 integer, intent(out) :: e_m(nodi_lg)
 double precision, intent(inout) :: x_lg(nodi_lg),y_lg(nodi_lg)

 do kk=1,nodi_y_lg
   do k=1,nodi_x_lg
     i=nodi_x_lg*(kk-1)+k
     vx_rg=0.D0
     vy_rg=0.D0
     call Composizione(elementi,nodi,nodi_elem,connessioni,x,y,elem_x,elem_y,x_lg(i),y_lg(i),e_m(i),x_min,x_max,y_min,y_max,larghezza,altezza,&
          free_surface,colonna,riga,layer_x,layer_y,n_elem_x,n_elem_y,thick_x,thick_y)
     if(e_m(i).ne.0)then
       posgp(1)=(2.D0/delta_x(e_m(i)))*(x_lg(i)-xmin_elem(e_m(i)))-1.D0
       posgp(2)=(4*y_lg(i)-((1.D0-posgp(1))*y(connessioni(e_m(i),1))+(1.D0+posgp(1))*y(connessioni(e_m(i),2))+&
                (1.D0+posgp(1))*y(connessioni(e_m(i),3))+(1.D0-posgp(1))*y(connessioni(e_m(i),4))))/&
                (-(1.D0-posgp(1))*y(connessioni(e_m(i),1))-(1.D0+posgp(1))*y(connessioni(e_m(i),2))+&
                (1.D0+posgp(1))*y(connessioni(e_m(i),3))+(1.D0-posgp(1))*y(connessioni(e_m(i),4)))
       N=(/0.25D0*(1.D0-posgp(1))*(1.D0-posgp(2)),0.25D0*(1.D0+posgp(1))*(1.D0-posgp(2)),&
           0.25D0*(1.D0+posgp(1))*(1.D0+posgp(2)),0.25D0*(1.D0-posgp(1))*(1.D0+posgp(2))/)
       do j=1,nodi_elem
         vx_rg=vx_rg+N(j)*Vx(connessioni(e_m(i),j))
         vy_rg=vy_rg+N(j)*Vy(connessioni(e_m(i),j))
       enddo
       x_lg(i)=x_lg(i)+vx_rg*dt
       y_lg(i)=y_lg(i)+vy_rg*dt
     endif
     if(kk.eq.nodi_y_lg.and.e_m(i).eq.0)then
       x_lg(i)=x_lg(i-nodi_x_lg)
       y_lg(i)=y_lg(i-nodi_x_lg)+dim_y_lg
     elseif(kk.eq.nodi_y_lg.and.e_m(i).ne.0)then
       x_lg(i)=x_lg(i-nodi_x_lg)
     endif
   enddo
 enddo
 do i=1,elem_lg
   do j=1,4
     if(j.eq.1.and.x_lg(j).gt.x_lg(j+1))x_lg(j)=x_lg(j+1)-eps
     if(j.eq.1.and.y_lg(j).gt.y_lg(4))y_lg(j)=y_lg(4)-eps
     if(j.eq.2.and.x_lg(j).lt.x_lg(j-1))x_lg(j)=x_lg(j-1)+eps
     if(j.eq.2.and.y_lg(j).gt.y_lg(j+1))y_lg(j)=y_lg(j+1)-eps
     if(j.eq.3.and.x_lg(j).lt.x_lg(j+1))x_lg(j)=x_lg(j+1)+eps
     if(j.eq.3.and.y_lg(j).lt.y_lg(j-1))y_lg(j)=y_lg(j-1)+eps
     if(j.eq.4.and.x_lg(j).gt.x_lg(j-1))x_lg(j)=x_lg(j-1)-eps
     if(j.eq.4.and.y_lg(j).lt.y_lg(1))y_lg(j)=y_lg(1)+eps
   enddo
 enddo

      
 end subroutine Advection_lgrid

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 subroutine Output_final(nodi,elementi,nodi_elem,layer_x,layer_y,thick_x,thick_y,n_elem_x,n_elem_y,connessioni,Temp,x,y,&
 x_marc,y_marc,n_marc,marcatori,type_marc,e2t_p,e2t_v,W_b,W_f,W_m,melt_fraction,OUT_DIR,OUT_MODEL,OUT_PHASE)

 Implicit none
 
 integer :: i,j
 integer, intent(in) :: nodi,elementi,nodi_elem,layer_x,layer_y,n_elem_x(layer_x),n_elem_y(layer_y),connessioni(elementi,nodi_elem),&
 marcatori,type_marc(marcatori),n_marc(marcatori)

 double precision, intent(in) :: x(nodi),y(nodi),thick_x(layer_x),thick_y(layer_y),Temp(nodi),x_marc(marcatori),y_marc(marcatori),&
 e2t_p(marcatori),e2t_v(marcatori),W_b(marcatori),W_f(marcatori),W_m(marcatori),melt_fraction(marcatori)

 character(len=*), intent(in) :: OUT_DIR,OUT_MODEL,OUT_PHASE

 open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Griglia_final',status='unknown') 
 write(1,1000)
 write(1,'(i10,1x,20(f15.5,i10))')layer_x,(thick_x(i),n_elem_x(i),i=1,layer_x)
 write(1,'(i10,1x,20(f15.5,i10))')layer_y,(thick_y(i),n_elem_y(i),i=1,layer_y)
 write(1,*)''
 write(1,2000)
 do i=1,nodi
   write(1,'(i10,1x,2(f20.9,1x),f15.4,1x,i4)')i,x(i),y(i),Temp(i)
 enddo
 write(1,3000)
 do i=1,elementi
     write(1,'(5(i10,1x))')i,(connessioni(i,j),j=1,nodi_elem)
 enddo
 close(1)

 open(unit=1,file=TRIM(OUT_DIR)//'/'//TRIM(OUT_MODEL)//'/'//TRIM(OUT_PHASE)//'/Marcatori_final',status='unknown') 
 rewind(1)
 write(1,*)'Marcatori totali'
 write(1,*)marcatori
 write(1,4000)
 do i=1,marcatori
   write(1,'(i10,2(f20.9,1x),i10,1x,6(f15.4,1x))')n_marc(i),x_marc(i),y_marc(i),type_marc(i),e2t_p(i),e2t_v(i),W_b(i),W_f(i),W_m(i),melt_fraction(i)
 enddo
 close(1)

1000 format(5x,'Layer',6x,'Dimensione',2x,'Elementi')
2000 format(6x,'Nodo',9x,'Coordinata X',9x,'Coordinata Y',5x,'Temperatura')
3000 format(2x,'Elemento',4x,'Conn. 1',4x,'Conn. 2',4x,'Conn. 3',4x,'Conn. 4')
4000 format(1x,'Marcatore',8x,'Coordinata X',9x,'Coordinata Y',7x,'Type',2x,'Plastic strain',2x,'Viscous strain',5x,'Bound water',&
     6x,'Free water',6x,'Melt water',3x,'Melt fraction')

 end subroutine Output_final

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end program FALCON

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     13 - FUNCTIONS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 function molten(a,b,c,d,e,f,g,x)
 double precision :: molten,x,a,b,c,d,e,f,g
 molten=-a+b+c*x-d*(e/(f*(1.D0-x)+x))**g
 end function

 function delta(a,b,c)
 double precision :: delta,a,b,c
 delta=(b**2.D0)-4.D0*a*c
 end function

 function roots(a,b,d)
 double precision :: roots,a,b,d
 roots=(-b+d)/(2.D0*a)
 end function
