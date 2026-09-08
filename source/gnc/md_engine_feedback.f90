! Heuristic GL2026 (arXiv:2608.28947) positive feedback:
! Gamma_TDE -> molecular-cloud density via covering + tidal + virial,
! then a cloud-driven boost to two-body relaxation (Coulomb log / DCs).
!
! Disk-fed geometry (paper Sec. IV.4): f_Omega=0.2, beamed UDRs (b=1).
! Covering (eq. 60) + tidal rho_MC=9 Mh/(4 pi r_b^3) (eq. 21) + virial
! tau_c close n_MC(Gamma). Heating (eq. 61) then gives R_MC and M_MC.
! Default: hold M_MC fixed so nm2_cloud = rho_MC * M_fixed ∝ n_MC
! (Sec. III fixed-mass collapse = the actual positive loop). Set
! --cloud mass = 0 to use the solved M_MC instead (nearly Gamma-independent).
module md_engine_feedback
	use com_main_gw
	use constant
	implicit none

	real(8),parameter::mu_mean=1.4d0
		real(8),parameter::Lambda0_co=1.3d-27
	real(8),parameter::E50_norm=1.8d0
	real(8),parameter::RUDR_pc_norm=0.22d0
	real(8),parameter::n_mc_lo=1d0
	real(8),parameter::n_mc_hi=1d9

contains

	subroutine engine_feedback_reset()
		implicit none
		ctl%engine_relax_boost=1d0
		ctl%engine_nm2_cloud=0d0
		ctl%engine_n_mc=0d0
		ctl%engine_r_b=0d0
		ctl%engine_R_mc=0d0
		ctl%engine_M_mc=0d0
		ctl%engine_gamma_tde=0d0
	end subroutine

	subroutine engine_feedback_update()
		use md_dms_saving_data
		implicit none
		real(8) gamma_yr,n_mc,r_b_pc,R_mc_pc,M_mc_sol,M_use,rho_cgs
		real(8) rho_au,nm2_cloud,nm2_rh,boost,cover

		if(ctl%engine_feedback.lt.1)then
			call engine_feedback_reset()
			return
		end if

		gamma_yr=oe_star%se_td%rate/1d6
		ctl%engine_gamma_tde=gamma_yr
		if(gamma_yr.le.0d0)then
			call engine_feedback_reset()
			return
		end if

		call gl2026_cloud_from_gamma(gamma_yr,spp_new%mbh,n_mc,r_b_pc,&
			R_mc_pc,M_mc_sol,cover)
		if(n_mc.le.0d0)then
			if(rid.eq.0)then
				print*, "engine feedback: Gamma_TDE too low to drive belt, Gamma=",gamma_yr
			end if
			call engine_feedback_reset()
			ctl%engine_gamma_tde=gamma_yr
			return
		end if

		ctl%engine_n_mc=n_mc
		ctl%engine_r_b=r_b_pc
		ctl%engine_R_mc=R_mc_pc
		ctl%engine_M_mc=M_mc_sol

		if(ctl%engine_cloud_mass.gt.0d0)then
			M_use=ctl%engine_cloud_mass
		else
			M_use=M_mc_sol
		end if

		rho_cgs=mu_mean*m_proton_GS*n_mc
		rho_au=rho_cgs*AU_GS**3/m_sun_GS
		nm2_cloud=rho_au*M_use

		call engine_nm2_at_rh(nm2_rh)
		boost=1d0
		if(nm2_rh.gt.0d0)then
			boost=1d0+nm2_cloud/nm2_rh
		end if
		if(boost.gt.ctl%engine_max_boost)then
			boost=ctl%engine_max_boost
			if(nm2_rh.gt.0d0)then
				nm2_cloud=(boost-1d0)*nm2_rh
			end if
		end if
		ctl%engine_relax_boost=boost
		ctl%engine_nm2_cloud=nm2_cloud

		if(rid.eq.0)then
			print*, "engine feedback: Gamma(/yr), n_MC(cm-3), r_b(pc), R_MC(pc), M_MC, boost=",&
				gamma_yr,n_mc,r_b_pc,R_mc_pc,M_use,boost
		end if
	end subroutine

	subroutine gl2026_cloud_from_gamma(gamma_yr,Mh,n_mc,r_b_pc,R_mc_pc,M_mc_sol,cover)
		implicit none
		real(8),intent(in)::gamma_yr,Mh
		real(8),intent(out)::n_mc,r_b_pc,R_mc_pc,M_mc_sol,cover
		real(8) lo,hi,mid,c_lo,c_mid
		real(8) Mh6,beta,rho,rb,E_udr,R_mc,f_omega,b_udr
		integer k

		n_mc=0d0
		r_b_pc=0d0
		R_mc_pc=0d0
		M_mc_sol=0d0
		cover=0d0
		if(Mh.le.0d0.or.gamma_yr.le.0d0) return

		f_omega=ctl%engine_f_omega
		if(f_omega.le.0d0) f_omega=0.2d0
		if(ctl%engine_beamed.ge.1)then
			b_udr=1d0
			beta=b_udr/f_omega
		else
			b_udr=f_omega
			beta=1d0
		end if

		c_lo=gl2026_covering_lhs(n_mc_lo,gamma_yr,Mh,beta)
		if(c_lo.lt.1d0) return

		lo=n_mc_lo
		hi=n_mc_hi
		do k=1,80
			mid=sqrt(lo*hi)
			c_mid=gl2026_covering_lhs(mid,gamma_yr,Mh,beta)
			if(c_mid.gt.1d0)then
				lo=mid
			else
				hi=mid
			end if
		end do
		n_mc=sqrt(lo*hi)
		cover=gl2026_covering_lhs(n_mc,gamma_yr,Mh,beta)

		Mh6=Mh/1d6
		rho=mu_mean*m_proton_GS*n_mc
		rb=(9d0*Mh*m_sun_GS/(4d0*pi*rho))**(1d0/3d0)
		r_b_pc=rb/pc_GS

		E_udr=1.8d50*Mh6**(1d0/3d0)
		R_mc=(gamma_yr/one_year)*E_udr*b_udr/ &
			((16d0*pi/3d0)*f_omega*Lambda0_co*n_mc**2*rb**2)
		if(R_mc.le.0d0)then
			R_mc_pc=0d0
			M_mc_sol=0d0
			return
		end if
		R_mc_pc=R_mc/pc_GS
		M_mc_sol=(4d0*pi/3d0)*rho*R_mc**3/m_sun_GS
	end subroutine

	real(8) function gl2026_covering_lhs(n_mc,gamma_yr,Mh,beta)
		implicit none
		real(8),intent(in)::n_mc,gamma_yr,Mh,beta
		real(8) rho,rb,n4,E50,R_udr,f_udr,tau_c,Mh6
		if(n_mc.le.0d0)then
			gl2026_covering_lhs=0d0
			return
		end if
		Mh6=Mh/1d6
		rho=mu_mean*m_proton_GS*n_mc
		rb=(9d0*Mh*m_sun_GS/(4d0*pi*rho))**(1d0/3d0)
		n4=n_mc/1d4
		E50=E50_norm*Mh6**(1d0/3d0)
		R_udr=RUDR_pc_norm*(E50**(5d0/17d0))*(n4**(-7d0/17d0))*pc_GS
		f_udr=R_udr**2/(4d0*rb**2)
		tau_c=2d0/sqrt((4d0*pi/5d0)*G_GS*rho)/one_year
		gl2026_covering_lhs=gamma_yr*f_udr*tau_c*beta
	end function

	subroutine engine_nm2_at_rh(nm2_tot)
		implicit none
		real(8),intent(out)::nm2_tot
		real(8) rh_now,rh_dmless,rho_tmp
		integer k
		nm2_tot=0d0
		call get_rh_now(rh_now)
		rh_dmless=10d0**rh_now
		do k=1,dms%n
			call get_rho_full_range(dms%mb(k)%all%fmden,dms%mb(k)%all%spt_rho_rmin,&
				log10(rh_dmless),rho_tmp,spp_new)
			nm2_tot=nm2_tot+rho_tmp*ctl%n0*dms%mb(k)%mc
		end do
	end subroutine

end module
