! Orbit-averaged star-star collisions.
! Local rate Gamma(r) = n(r) * sigma * v_rel(r), with gravitational focusing.
! A sample is hit with p = 1 - exp(-<Gamma> dt). Outcome vs surface escape:
!   v_rel <  v_esc = sqrt(2G(M1+M2)/(R1+R2))  -> merger: keep the sample,
!     redraw J uniform in (0, Jc], recompute (rp,ra,P). GL2026: collision
!     randomizes the orbit and usually kicks it out of the loss cone.
!   v_rel >= v_esc  -> destroy both (exit_collision), as before.
! Mass is not grown: the operator is one-body vs a field (no second MC
! sample), and discrete mass bins + clone weights cannot host M1+M2.
module md_particle_collision
	use com_main_gw
	use,intrinsic::ieee_arithmetic
	implicit none

	integer,parameter::stellar_collision_method_off=0
	integer,parameter::stellar_collision_method_oa_destroy=1
	integer,parameter::stellar_collision_pair_samebin=0
	integer,parameter::stellar_collision_pair_all=1

	type(s1d_type)::coll_gamma
	type(s1d_type)::coll_n_collid
	type(s1d_type)::coll_vrel
	type(s1d_type)::coll_n_phys
	type(s1d_type)::coll_n_mb(20)
	real(8)::coll_M_mb(20)
	real(8)::coll_R_mb(20)
	integer::coll_nmb=0
	logical::coll_tables_ready=.false.

	real(8),parameter::code_time_to_myr=2d0*pi*1d6
	integer,parameter::nquad_orbit=48

	private::collision_sigma,collision_vesc,get_vrel_at_logr
	private::sample_is_collidable,sample_collision_radius
	private::orbit_average_gamma,gamma_at_r_for_sample
	private::select_fden,apply_collision_outcome,merge_collision_remnant

contains

	real(8) function collision_sigma(r1,r2,m1,m2,vrel)
		implicit none
		real(8),intent(in)::r1,r2,m1,m2,vrel
		real(8) rsum,v2,focus
		real(8),parameter::v2min=1d-30
		rsum=r1+r2
		if(rsum.le.0d0)then
			collision_sigma=0d0
			return
		end if
		v2=max(vrel*vrel,v2min)
		! G=1 AU/Msun (same physical units as get_trlx_* / star_Radius).
		! Cap focusing: a numerically vanishing vrel (1e-11 vs Kepler ~1e3)
		! otherwise makes 2GM/(R v^2) ~ 1e24 and Gamma ~ 1e12 /Myr.
		focus=1d0+2d0*(m1+m2)/(rsum*v2)
		if(focus.gt.1d6) focus=1d6
		collision_sigma=pi*rsum*rsum*focus
	end function

	! G=1 AU/Msun, same units as collision_sigma / star_Radius.
	real(8) function collision_vesc(r1,r2,m1,m2)
		implicit none
		real(8),intent(in)::r1,r2,m1,m2
		real(8) rsum
		rsum=r1+r2
		if(rsum.le.0d0.or.m1+m2.le.0d0)then
			collision_vesc=0d0
			return
		end if
		collision_vesc=sqrt(2d0*(m1+m2)/rsum)
	end function

	subroutine select_fden(so,fden)
		implicit none
		type(dms_stellar_object),intent(in)::so
		type(s1d_type),intent(out)::fden
		if(ctl%source_fden.eq.source_simu)then
			fden=so%fden_simu
		else
			fden=so%fden
		end if
	end subroutine

	logical function sample_is_collidable(sp)
		implicit none
		type(particle_sample_type),intent(in)::sp
		select case(sp%obtype)
		case(star_type_ms,star_type_rg,star_type_nakedHe,star_type_bd)
			sample_is_collidable=.true.
		case default
			sample_is_collidable=.false.
		end select
	end function

	real(8) function sample_collision_radius(sp)
		implicit none
		type(particle_sample_type),intent(in)::sp
		real(8),external::star_Radius
		if(sp%byot%ms%radius.gt.0d0)then
			sample_collision_radius=sp%byot%ms%radius
		else
			sample_collision_radius=star_Radius(sp%m)
		end if
	end function

	subroutine get_vrel_at_logr(logr,vrel)
		implicit none
		real(8),intent(in)::logr
		real(8),intent(out)::vrel
		vrel=0d0
		if(.not.coll_tables_ready) return
		if(logr.lt.coll_vrel%xmin.or.logr.gt.coll_vrel%xmax) return
		call coll_vrel%get_value_l(logr,vrel)
		if(vrel.lt.0d0) vrel=0d0
	end subroutine

	real(8) function gamma_at_r_for_sample(r_xy,sp,idx)
		implicit none
		real(8),intent(in)::r_xy
		type(particle_sample_type),intent(in)::sp
		integer,intent(in)::idx
		real(8) logr,n_j,vrel,sig,Rs,gamma_code
		integer j
		gamma_at_r_for_sample=0d0
		if(r_xy.le.0d0) return
		logr=log10(r_xy)
		call get_vrel_at_logr(logr,vrel)
		if(vrel.le.0d0) return
		Rs=sample_collision_radius(sp)
		if(Rs.le.0d0) return

		if(ctl%stellar_collision_pair_mode.eq.stellar_collision_pair_samebin)then
			if(idx.ge.1.and.idx.le.coll_nmb)then
				if(logr.ge.coll_n_mb(idx)%xmin.and.logr.le.coll_n_mb(idx)%xmax)then
					call coll_n_mb(idx)%get_value_l(logr,n_j)
					if(n_j.gt.0d0)then
						sig=collision_sigma(Rs,coll_R_mb(idx),sp%m,coll_M_mb(idx),vrel)
						gamma_code=n_j*sig*vrel
						gamma_at_r_for_sample=gamma_code*code_time_to_myr
					end if
				end if
			end if
		else
			do j=1,coll_nmb
				if(logr.ge.coll_n_mb(j)%xmin.and.logr.le.coll_n_mb(j)%xmax)then
					call coll_n_mb(j)%get_value_l(logr,n_j)
					if(n_j.gt.0d0)then
						sig=collision_sigma(Rs,coll_R_mb(j),sp%m,coll_M_mb(j),vrel)
						gamma_code=n_j*sig*vrel
						gamma_at_r_for_sample=gamma_at_r_for_sample+gamma_code*code_time_to_myr
					end if
				end if
			end do
		end if
	end function

	subroutine orbit_average_gamma(sp,idx,gamma_avg)
		implicit none
		type(particle_sample_type),intent(in)::sp
		integer,intent(in)::idx
		real(8),intent(out)::gamma_avg
		real(8) rp,ra,jc_xy,theta,r,vr,vt,vm,dth,sint,cost
		real(8) gsum,g_r,pd_xy
		integer i
		gamma_avg=0d0
		rp=sp%rp/r0_cl
		ra=sp%ra/r0_cl
		if(rp.le.0d0.or.ra.le.rp)then
			if(rp.gt.0d0)then
				gamma_avg=gamma_at_r_for_sample(rp,sp,idx)
			end if
			return
		end if
		if(sp%jm.ge.0.99999d0)then
			gamma_avg=gamma_at_r_for_sample(0.5d0*(rp+ra),sp,idx)
			return
		end if
		jc_xy=sp%jc/(r0_cl*ctl%v0)
		pd_xy=sp%period*ctl%v0/r0_cl
		if(pd_xy.le.0d0.or..not.ieee_is_finite(pd_xy))then
			gamma_avg=gamma_at_r_for_sample(0.5d0*(rp+ra),sp,idx)
			return
		end if

		dth=0.5d0*pi/dble(nquad_orbit)
		gsum=0d0
		do i=1,nquad_orbit
			theta=(dble(i)-0.5d0)*dth
			sint=sin(theta)
			cost=cos(theta)
			r=rp+(ra-rp)*sint*sint
			if(r.le.0d0) cycle
			call get_v_from_r(spp_new,r,sp%x,sp%jm,jc_xy,vm,vr,vt)
			if(vr.le.1d-20) cycle
			! dt ~ dr/vr ; dr = 2(ra-rp) sin theta cos theta dtheta
			g_r=gamma_at_r_for_sample(r,sp,idx)
			gsum=gsum+g_r*2d0*(ra-rp)*sint*cost/vr
		end do
		gsum=gsum*dth
		! <Gamma> = (2/P) * int Gamma dr/vr ; the 2 is inbound+outbound
		gamma_avg=2d0*gsum/pd_xy
		if(gamma_avg.lt.0d0.or..not.ieee_is_finite(gamma_avg)) gamma_avg=0d0
	end subroutine

	subroutine update_stellar_collision_rates()
		use md_star_pot
		implicit none
		type(s1d_type)::fden
		real(8) rho_tmp,rho_tot,phi_tmp,vh,n_phys,logr,r_xy
		real(8) gamma_max,gamma_pop,gamma_r,sig,vrel,n_j,v_char,wmax,wbin
		real(8),external::star_Radius
		integer i,k,j,igmax

		if(ctl%stellar_collision_method.lt.stellar_collision_method_oa_destroy)then
			ctl%tmax_collision=1d99
			coll_tables_ready=.false.
			return
		end if
		if(dms%n.le.0)then
			ctl%tmax_collision=1d99
			coll_tables_ready=.false.
			return
		end if

		call select_fden(dms%all%all,fden)
		if(fden%nbin.le.0)then
			ctl%tmax_collision=1d99
			coll_tables_ready=.false.
			return
		end if

		coll_nmb=dms%n
		do k=1,dms%n
			coll_M_mb(k)=dms%mb(k)%mc
			coll_R_mb(k)=star_Radius(dms%mb(k)%mc)
			call select_fden(dms%mb(k)%all,coll_n_mb(k))
			! convert dimensionless fden to physical number density (AU^-3)
			coll_n_mb(k)%fx=coll_n_mb(k)%fx*ctl%n0
		end do

		coll_vrel=fden
		coll_n_phys=fden
		coll_gamma=fden
		coll_n_collid=fden
		coll_vrel%fx=0d0
		coll_n_phys%fx=0d0
		coll_gamma%fx=0d0
		coll_n_collid%fx=0d0

		do i=1,fden%nbin
			logr=fden%xb(i)
			r_xy=10**logr
			rho_tot=0d0
			n_phys=0d0
			do k=1,dms%n
				call get_rho_full_range(dms%mb(k)%all%fmden,dms%mb(k)%all%spt_rho_rmin,&
					logr,rho_tmp)
				if(rho_tmp.gt.0d0) rho_tot=rho_tot+rho_tmp
				if(logr.ge.coll_n_mb(k)%xmin.and.logr.le.coll_n_mb(k)%xmax)then
					call coll_n_mb(k)%get_value_l(logr,n_j)
					if(n_j.gt.0d0) n_phys=n_phys+n_j
				end if
			end do
			coll_n_phys%fx(i)=n_phys

			call get_phi_star_full_range(spp_new,logr,phi_tmp)
			phi_tmp=10**phi_tmp+spp_new%mbh_dmless/r_xy
			vh=0d0
			if(rho_tot.gt.0d0)then
				call get_v_dispersion_one_ir(dms%all%all%barge_ir,phi_tmp,rho_tot,vh)
			end if
			! get_v_dispersion can return ~1e-11 when barge_ir has no support
			! but fmden is extrapolated (C0c snap1 i=3). Relaxation already
			! skips vh<1d-5; collisions cannot — Gamma_focus ~ 1/v.
			! Require vh to be a non-tiny fraction of the local circular speed.
			v_char=0d0
			if(phi_tmp.gt.0d0) v_char=sqrt(phi_tmp)*ctl%v0
			if(vh.gt.0d0.and.ieee_is_finite(vh).and.vh.gt.1d-4*v_char)then
				! mean relative speed of two Maxwellians: 4/sqrt(pi) * sigma_1D
				vrel=4d0/sqrt(pi)*vh
			else
				vrel=0d0
			end if
			coll_vrel%fx(i)=vrel

			gamma_r=0d0
			if(vrel.gt.0d0.and.n_phys.gt.0d0)then
				if(ctl%stellar_collision_pair_mode.eq.stellar_collision_pair_samebin)then
					do k=1,coll_nmb
						n_j=0d0
						if(logr.ge.coll_n_mb(k)%xmin.and.logr.le.coll_n_mb(k)%xmax)then
							call coll_n_mb(k)%get_value_l(logr,n_j)
						end if
						if(n_j.gt.0d0.and.coll_R_mb(k).gt.0d0)then
							sig=collision_sigma(coll_R_mb(k),coll_R_mb(k),&
								coll_M_mb(k),coll_M_mb(k),vrel)
							gamma_r=max(gamma_r,n_j*sig*vrel*code_time_to_myr)
						end if
					end do
				else
					! fiducial all-pairs rate for the most massive MS-like bin
					do k=1,coll_nmb
						n_j=0d0
						if(logr.ge.coll_n_mb(k)%xmin.and.logr.le.coll_n_mb(k)%xmax)then
							call coll_n_mb(k)%get_value_l(logr,n_j)
						end if
						if(n_j.gt.0d0.and.coll_R_mb(k).gt.0d0)then
							sig=collision_sigma(coll_R_mb(k),coll_R_mb(1),&
								coll_M_mb(k),coll_M_mb(1),vrel)
							gamma_r=gamma_r+n_j*sig*vrel*code_time_to_myr
						end if
					end do
				end if
			end if
			if(.not.ieee_is_finite(gamma_r).or.gamma_r.lt.0d0) gamma_r=0d0
			coll_gamma%fx(i)=gamma_r
			! expected collision-rate density (n * Gamma) on the r-grid
			coll_n_collid%fx(i)=n_phys*gamma_r
		end do

		! Local Gamma in the empty inner cusp can be ~1/Myr even when
		! almost no mass lives there (C0c: M(<10 AU)~0.02 Msun). The
		! destruction operator already uses orbit-averaged <Gamma>, so
		! timestep only on populated bins: n r^3 above 1e-4 of its max.
		gamma_max=0d0
		do i=1,coll_gamma%nbin
			if(coll_gamma%fx(i).gt.gamma_max) gamma_max=coll_gamma%fx(i)
		end do
		wmax=0d0
		do i=1,coll_n_phys%nbin
			r_xy=10**coll_n_phys%xb(i)
			wbin=coll_n_phys%fx(i)*r_xy*r_xy*r_xy
			if(wbin.gt.wmax) wmax=wbin
		end do
		gamma_pop=0d0
		igmax=1
		do i=1,coll_gamma%nbin
			r_xy=10**coll_gamma%xb(i)
			wbin=coll_n_phys%fx(i)*r_xy*r_xy*r_xy
			if(wmax.gt.0d0.and.wbin.ge.1d-4*wmax)then
				if(coll_gamma%fx(i).gt.gamma_pop)then
					gamma_pop=coll_gamma%fx(i)
					igmax=i
				end if
			end if
		end do
		if(gamma_pop.gt.0d0.and.ieee_is_finite(gamma_pop))then
			ctl%tmax_collision=ctl%tfractor_collision/gamma_pop
		else
			ctl%tmax_collision=1d99
		end if
		coll_tables_ready=.true.

		if(rid.eq.0)then
			print*, "stellar collision: gamma_max_all, gamma_pop(1/Myr), tmax_collision(Myr)=", &
				gamma_max, gamma_pop, ctl%tmax_collision
			if(igmax.ge.1.and.igmax.le.coll_gamma%nbin)then
				print*, "stellar collision: logr, r/r0, vrel, n_phys at gamma_pop=", &
					coll_gamma%xb(igmax), 10**coll_gamma%xb(igmax), &
					coll_vrel%fx(igmax), coll_n_phys%fx(igmax)
			end if
		end if
	end subroutine

	subroutine merge_collision_remnant(sp)
		implicit none
		type(particle_sample_type),intent(inout)::sp
		real(8) jm_new
		real(8),external::rnd
		! Uniform J in (0, Jc]: jm = J/Jc. Does not change M, weight, or bin.
		jm_new=rnd(0d0,1d0)
		if(jm_new.le.0d0) jm_new=jmin_value
		call set_jm_bound(jm_new)
		sp%jm=jm_new
		sp%jph=sp%jm*sp%jc
		call update_sample_para(sp,spp_new)
	end subroutine

	subroutine apply_collision_outcome(sp,idx,n_merge,n_destroy)
		implicit none
		type(particle_sample_type),intent(inout)::sp
		integer,intent(in)::idx
		integer,intent(inout)::n_merge,n_destroy
		real(8) r_xy,vrel,vesc,r1,r2,m1,m2
		r1=sample_collision_radius(sp)
		m1=sp%m
		if(idx.ge.1.and.idx.le.coll_nmb)then
			m2=coll_M_mb(idx)
			r2=coll_R_mb(idx)
		else
			m2=m1
			r2=r1
		end if
		vesc=collision_vesc(r1,r2,m1,m2)
		r_xy=0.5d0*(sp%rp+sp%ra)/r0_cl
		if(r_xy.le.0d0) r_xy=sp%rp/r0_cl
		vrel=0d0
		if(r_xy.gt.0d0) call get_vrel_at_logr(log10(r_xy),vrel)
		! Slow / focused: keep the star, randomize J out of the loss cone.
		! Fast / head-on, or unknown vrel: shred (exit_collision).
		if(vrel.gt.0d0.and.vesc.gt.0d0.and.vrel.lt.vesc)then
			call merge_collision_remnant(sp)
			n_merge=n_merge+1
		else
			sp%exit_flag=exit_collision
			sp%exit_time=ctl%run_snap_time_f
			n_destroy=n_destroy+1
		end if
	end subroutine

	subroutine apply_stellar_collision_operator(dt)
		implicit none
		real(8),intent(in)::dt
		type(chain_pointer_type),pointer::ps
		real(8) gamma_avg,p_hit,surv
		real(8),external::rnd
		integer idx,n_merge,n_destroy
		if(ctl%stellar_collision_method.lt.stellar_collision_method_oa_destroy) return
		if(.not.coll_tables_ready) return
		if(dt.le.0d0) return

		n_merge=0
		n_destroy=0
		ps=>bksams%head
		do while(associated(ps))
			select type(ca=>ps%ob)
			type is (particle_sample_type)
				if(ca%exit_flag.eq.exit_normal.and.sample_is_collidable(ca))then
					call get_mass_idx(ca%m,idx)
					if(idx.lt.1) idx=1
					call orbit_average_gamma(ca,idx,gamma_avg)
					if(gamma_avg.gt.0d0)then
						p_hit=1d0-exp(-gamma_avg*dt)
						if(p_hit.gt.1d0) p_hit=1d0
						if(ctl%collision_consider_weight.ge.1)then
							surv=exp(-gamma_avg*dt)
							ca%weight_n=ca%weight_n*surv
							ca%weight_real=ca%weight_clone*ca%weight_n*ctl%n_basic
							if(ca%weight_real.lt.1d-12)then
								ca%exit_flag=exit_collision
								ca%exit_time=ctl%run_snap_time_f
								n_destroy=n_destroy+1
							elseif(rnd(0d0,1d0).lt.p_hit)then
								call apply_collision_outcome(ca,idx,n_merge,n_destroy)
							end if
						else
							if(rnd(0d0,1d0).lt.p_hit)then
								call apply_collision_outcome(ca,idx,n_merge,n_destroy)
							end if
						end if
					end if
				end if
			end select
			ps=>ps%next
		end do
		if(rid.eq.0)then
			print*, "stellar collision operator: n_merge, n_destroy (rank0)=", &
				n_merge, n_destroy
		end if
	end subroutine

end module
