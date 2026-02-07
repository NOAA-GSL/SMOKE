!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_anthro_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  private

  public :: mpas_smoke_anthro_emis_driver

contains


  subroutine mpas_smoke_anthro_emis_driver(dt,gmt,julday,kemit,                      &
                           xlat,xlong, chem,num_chem,dz8w,t_phy,rho_phy,             &    
                           e_ant_in, e_ant_out, num_e_ant_in, num_e_ant_out,         &
                           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
                           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday, num_chem, kemit,           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_e_ant_in, num_e_ant_out,       &
           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse

   REAL(RKIND), INTENT(IN    ) :: dt,gmt

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w,rho_phy,t_phy
   REAL(RKIND),DIMENSION(ims:ime,1:kemit,jms:jme,1:num_e_ant_in), INTENT(IN)    :: e_ant_in
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_e_ant_out),INTENT(INOUT) :: e_ant_out
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n
   REAL(RKIND) :: conv_aer, conv_gas, emis

   REAL(RKIND), PARAMETER :: rwc_t_thresh = 283.15 ! [ 50 F]


   do j = jts,jte
   do k = kts, kemit
   do i = its,ite
!  
!     Conversion factor for aerosol emissions (ug/m2/s) --> ug/kg
      conv_aer = dt / (rho_phy(i,k,j) *  dz8w(i,k,j))
!     Conversion factor for gas phase emissions (mol/m2/s) --> ppm/ppm
      conv_gas = 60._RKIND * 1.E6_RKIND * 4.828E-4_RKIND * dt / ( rho_phy(i,k,j) * dz8w(i,k,j) )
!
      if (p_unspc_fine .gt. 0 .and. index_e_ant_in_unspc_fine .gt. 0 ) then
         emis = conv_aer*e_ant_in(i,k,j,index_e_ant_in_unspc_fine)
         chem(i,k,j,p_unspc_fine)   = chem(i,k,j,p_unspc_fine) + emis
         e_ant_out(i,k,j,index_e_ant_out_unspc_fine) = e_ant_out(i,k,j,index_e_ant_out_unspc_fine) + emis
      endif
      if (index_e_ant_in_unspc_coarse .gt. 0 ) then
         emis = conv_aer*e_ant_in(i,k,j,index_e_ant_in_unspc_coarse)
         if ( p_unspc_coarse .gt. 0 ) then
            chem(i,k,j,p_unspc_coarse)   = chem(i,k,j,p_unspc_coarse) + emis
         elseif (p_dust_coarse .gt. 0 ) then
            chem(i,k,j,p_dust_coarse)   = chem(i,k,j,p_dust_coarse) + emis
         endif
         e_ant_out(i,k,j,index_e_ant_out_unspc_coarse) = e_ant_out(i,k,j,index_e_ant_out_unspc_coarse) + emis
      endif
!     
   enddo ! i
   enddo ! k
   enddo ! j

  end subroutine mpas_smoke_anthro_emis_driver

end module module_anthro_emissions
