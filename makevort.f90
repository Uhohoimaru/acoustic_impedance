
!===============================================================
! Reading signinstTcos0/cos1/sin1 ... light.bin
! to compile : ifort -r8 -O3 -convert big_endian -mcmodel=medium -shared-intel -o read_signalTcos0light read_signalTcos0light.f90
!===============================================================

!!read_signalT90light.f90
subroutine makevort

	 use variables

    implicit none


	open(25,file="filename.txt")
	read(25,*)
	read(25,*)flowinst
	read(25,*)
	read(25,*)flowave
	read(25,*)
	read(25,*)gridname
	read(25,*)
	read(25,*)outfile


	uj = a0 * 0.9d0

	pi = 3.14159265358979323846d0
	
    ! eg for Tcos0 (change for Tcos1 ...)
    ! Instantaneous fields
    open(804,file=flowinst,form='unformatted',status='old')
	write(*,*)flowinst
    rewind(804)
    read(804) nbinst                ! nb of sampings
    read(804) nr                    ! nb of points in the radial direction
    read(804) nz                    ! nb of points in the axial direction
    read(804) thetagT0              ! azimuth
    
    print *, nbinst,nr,nz
    allocate(rg(nr))
    allocate(zg(nz))
    read(804) (rg(i),i=1,nr)        ! grid in the radial direction
    read(804) (zg(n),n=1,nz)        ! grid in the axial direction
    read(804) deltat
    write(*,*) deltat                ! time step
    read(804) r0                    ! jet radius
    write(*,*) r0
    deltat_nondim = deltat/(r0/a0)
    write(*,*) 'deltat_nonidim=',deltat_nondim 
    

    nb =11
 !	nbinst = 10 !for test 
    write(*,*)'flow inport start'
    
    !------
	jmax = nz
	kmax = 1
	lmax = nr
	!------

	!----------allocate matrix--------------
	!--------------------------------------------
	!data for cylindrical coordinate (instantaneous)
	 allocate(rhoT0(nz,nr))
    allocate(urT0(nz,nr))
    allocate(utT0(nz,nr))
    allocate(uzT0(nz,nr))
    allocate(pT0(nz,nr))
    allocate(buf(nz,nr,nbinst,5))
	 allocate(buf_ave(nz,nr,5))	
	!-----------------------------------------------
	!data for cylindrical cordinate (time-averaged)
	  allocate(rhoT0moy(nz,nr))
     allocate(urT0moy(nz,nr))
     allocate(utT0moy(nz,nr))
     allocate(uzT0moy(nz,nr))
     allocate(pT0moy(nz,nr))
	!-------------------------------------------------
	!data for cartecian cordinate (instantaneous)
	allocate(rho(jmax,lmax))
	allocate(u(jmax,lmax))
	allocate(v(jmax,lmax))
	allocate(w(jmax,lmax))
	!allocate(e(jmax,lmax))
	allocate(p(jmax,lmax))
allocate(p_total(jmax,lmax))
allocate(kinetic_e(jmax,lmax))
	!-------------------------------------------------
	!data for cartecian cordinate (time-averaged)
	allocate(rho_ave(jmax,lmax))
	allocate(u_ave(jmax,lmax))
	allocate(v_ave(jmax,lmax))
	allocate(w_ave(jmax,lmax))
	!allocate(e_ave(jmax,lmax))
	allocate(p_ave(jmax,lmax))
	allocate(vdvd_ave(jmax,lmax))
	!--------------------------------------------------
	!fluctuation value u'
	allocate(rhod(jmax,lmax))
	allocate(ud(jmax,lmax))
	allocate(vd(jmax,lmax))
	allocate(wd(jmax,lmax))
	!allocate(ed(jmax,lmax))
	allocate(pd(jmax,lmax))
	!-----------------------------------------------------
	!fluctuation value u'u'
	allocate(vdvd(jmax,lmax))

	!-----------------------------------------------------
	!gradient of p' 
	allocate(dpdx(jmax,lmax))
	allocate(d2pdx2(jmax,lmax))
	!-----------------------------------------------------
	!vorticity of xj' 
	allocate(vortxd(jmax,lmax))
	allocate(vortyd(jmax,lmax))
	allocate(vortzd(jmax,lmax))	
        allocate(aimp(jmax,lmax))
        allocate(tke(jmax,lmax))
        allocate(ptotal_all(jmax,lmax))
        allocate(ptotal_lami(jmax,lmax))
        allocate(ptotal_turb(jmax,lmax))
        allocate(aimp_all(jmax,lmax))
        allocate(aimp_lami(jmax,lmax))
        allocate(aimp_turb(jmax,lmax))

	!----------------------------------------------------
	!for out data
	allocate(q(jmax,kmax,lmax,20))
	!--------------------------------------------
	vdvd_ave = 0.0d0
	q=0.0d0


   allocate(x(nz,nr,nb))
    allocate(y(nz,nr,nb))
    allocate(z(nz,nr,nb))
    
    ds = 1.0d0
    do ii = 1,nb      ! samplings
       do j = 1,nz    ! axial direction
          do l = 1,nr ! radial direction
             x(j,l,ii) = zg(j)/r0
             y(j,l,ii) = rg(l) * dsin(0.d0)/r0 + ds*dble(ii)
             z(j,l,ii) = rg(l) * dcos(0.d0)/r0 
          enddo
       enddo
    enddo
    rg(:) = rg(:) / r0
    zg(:) = zg(:) / r0  
	
	dx = dabs( x(2,1,1) - x(1,1,1) ) 
	dy = dabs( y(1,1,2) - y(1,1,1) )
	dz = dabs( z(1,2,1) - z(1,1,1) )

	write(*,*)dx,dy,dz
	write(*,*)"read grid"


!-----------------------------------------------average flow 

   write(*,*)'average flow start'
    
    ! Time-averaged fields
    

    
     open(805,file=flowave,form='unformatted',status='old')
     rewind(805)
     read(805) ((rhoT0moy(j,l),l=1,nr),j=1,nz)  ! time-averaged density, velocities and pressure, vorticity norm and pressure
     read(805) ((urT0moy(j,l),l=1,nr),j=1,nz)!
     read(805) ((utT0moy(j,l),l=1,nr),j=1,nz)
     read(805) ((uzT0moy(j,l),l=1,nr),j=1,nz)
     read(805) ((pT0moy(j,l),l=1,nr),j=1,nz)

     close(805)
    
    !!! IMPLEMENT your own post-treatment
     !do ii=1,nb
       do j=1,nz
       do l=1,nr
    
!    !!--------------Nondimensionalization------------
     
          rhoT0moy(j,l) = rhoT0moy(j,l) / rho0
          urT0moy(j,l) = urT0moy(j,l) / uj
          utT0moy(j,l) = utT0moy(j,l) / uj
          uzT0moy(j,l) = uzT0moy(j,l) / uj
          pT0moy(j,l) = pT0moy(j,l) / (rho0 * uj * uj)
    !------------------------------------------------
    
    !!    cartesian verocity u v w when theta=0
          uu = (uzT0moy(j,l))
          vv = ( rg(l)*utT0moy(j,l)*dcos(0.0d0) + urT0moy(j,l)*dsin(0.0d0) )!
          ww = ( urT0moy(j,l)*dcos(0.0d0) -rg(l)*utT0moy(j,l)*dsin(0.0d0) )
	   !vv = 0.0e0
    !!!------------------------
    
                        rho_ave(j,l) = rhoT0moy(j,l)
    	         	u_ave(j,l) = uu!
			v_ave(j,l) = vv
			w_ave(j,l) = ww
			p_ave(j,l) = pT0moy(j,l)
	
       end do
   	end do
   

    write(6,*)'average flow inport ok'
!    

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------



!-------------read instanaous parameter-------------------------------------

!----initialization--------------------------------
 call system("mkdir ./gnuplot/data_rho")
 call system("mkdir ./gnuplot/data_rho/data")
 call system("cp -r ./plot ./gnuplot/data_rho/plot") 
 

 call system("mkdir ./gnuplot/data_u")
 call system("mkdir ./gnuplot/data_u/data")
 call system("cp -r ./plot ./gnuplot/data_u/plot")


 call system("mkdir ./gnuplot/data_v")
 call system("mkdir ./gnuplot/data_v/data")
 call system("cp -r ./plot ./gnuplot/data_v/plot")


 call system("mkdir ./gnuplot/data_w")
 call system("mkdir ./gnuplot/data_w/data")
 call system("cp -r ./plot ./gnuplot/data_w/plot")
 

 call system("mkdir ./gnuplot/data_p")
 call system("mkdir ./gnuplot/data_p/data")
 call system("cp -r ./plot ./gnuplot/data_p/plot")


 call system("mkdir ./gnuplot/data_ptotal_all")
 call system("mkdir ./gnuplot/data_ptotal_all/data")
 call system("cp -r ./plot ./gnuplot/data_ptotal_all/plot")
 

 call system("mkdir ./gnuplot/data_ptotal_turb")
 call system("mkdir ./gnuplot/data_ptotal_turb/data")
 call system("cp -r ./plot ./gnuplot/data_ptotal_turb/plot")


 call system("mkdir ./gnuplot/data_ptotal_lami")
 call system("mkdir ./gnuplot/data_ptotal_lami/data")
 call system("cp -r ./plot ./gnuplot/data_ptotal_lami/plot")
 

 call system("mkdir ./gnuplot/data_kinetic_e")
 call system("mkdir ./gnuplot/data_kinetic_e/data")
 call system("cp -r ./plot ./gnuplot/data_kinetic_e/plot")
 

 call system("mkdir ./gnuplot/data_aimp_all")
 call system("mkdir ./gnuplot/data_aimp_all/data")
 call system("cp -r ./plot ./gnuplot/data_aimp_all/plot")
 

 call system("mkdir ./gnuplot/data_aimp_turb")
 call system("mkdir ./gnuplot/data_aimp_turb/data")
 call system("cp -r ./plot ./gnuplot/data_aimp_turb/plot")


 call system("mkdir ./gnuplot/data_tke")
 call system("mkdir ./gnuplot/data_tke/data")
 call system("cp -r ./plot ./gnuplot/data_tke/plot")


 call system("mkdir ./gnuplot/data_aimp_lami")
 call system("mkdir ./gnuplot/data_aimp_lami/data")
 call system("cp -r ./plot ./gnuplot/data_aimp_lami/plot")


 call system("mkdir ./gnuplot/data_ptotal_tke")
 call system("mkdir ./gnuplot/data_ptotal_tke/data")
 call system("cp -r ./plot ./gnuplot/data_ptotal_tke/plot")


 call system("mkdir ./gnuplot/data_average")

!nbinst = 2
!---------------------------------------------------
    do ii=1,nbinst !nbinst
	 nbg = ii*interval
    write(*,*)nbg,'/',nbinst,'ok'
    ! density, velocity components and pressure, vorticity norm and pressure

       read(804) ((rhoT0(j,l),l=1,nr),j=1,nz)
       read(804) ((urT0(j,l),l=1,nr),j=1,nz)
       read(804) ((utT0(j,l),l=1,nr),j=1,nz)
       read(804) ((uzT0(j,l),l=1,nr),j=1,nz)
       read(804) ((pT0(j,l),l=1,nr),j=1,nz)

     !!! here your may need to store the fields in 3D tables (nr,nz,nbtinst)
      
       do j=1,nz
       do l=1,nr
    
    !--------------Nondimensionalization------
          rhoT0(j,l) = rhoT0(j,l) / rho0
          urT0(j,l) = urT0(j,l) / uj
          utT0(j,l) = utT0(j,l) / uj
          uzT0(j,l) = uzT0(j,l) / uj
          pT0(j,l) = pT0(j,l) / (rho0 * uj * uj)   
    !-----------------------------------------------

    !------cartesian verocity u v w 
          uu = uzT0(j,l)
          vv = utT0(j,l)
          ww = urT0(j,l)
		   !vv = 0.0e0
    !!----------------------------------
   
          rho(j,l) = rhoT0(j,l) 
          u(j,l) = uu  
          v(j,l) = vv
          w(j,l) = ww
          p(j,l) = pT0(j,l)
          kinetic_e(j,l) = 0.5d0*rho(j,l)*( uu*uu + vv*vv + ww*ww )	       
          
          !  aimp(j,l) = sqrt(rho(j,l)*p_total(j,l)*1.4d0)                       
          aimp(j,l) = sqrt(rho(j,l)*p(j,l)*1.4d0) 
 
 
       end do
    	end do


     do j=1,jmax
        do l=1,lmax
           
           !--------variable component------------
           
           rhod(j,l) = rho(j,l) - rho_ave(j,l)
           ud(j,l) = u(j,l) - u_ave(j,l)
           vd(j,l) = v(j,l) - v_ave(j,l)
           wd(j,l) = w(j,l) - w_ave(j,l)
           !ed(j,l) = e(j,l) - e_ave(j,l)
           pd(j,l) = p(j,l) - p_ave(j,l)
           
           tke(j,l) = 0.5d0*rho(j,l)*( ud(j,l)*ud(j,l) + vd(j,l)*vd(j,l) + wd(j,l)*wd(j,l) )
           ptotal_all(j,l) = p(j,l) + kinetic_e(j,l)
           ptotal_turb(j,l) =  tke(j,l)
           ptotal_lami(j,l) = kinetic_e(j,l)
 !+ 0.5d0*rho(j,l)*( u_ave(j,l)*u_ave(j,l) + v_ave(j,l)*v_ave(j,l) + w_ave(j,l)*w_ave(j,l) )
           
           aimp_all(j,l) = sqrt(rho(j,l)*ptotal_all(j,l)*1.4d0) 
           aimp_turb(j,l) = sqrt(rho(j,l)*ptotal_turb(j,l)*1.4d0)
           aimp_lami(j,l) = sqrt(rho(j,l)*ptotal_lami(j,l)*1.4d0)
    !----------test! vdvd means v'*v'
           vdvd(j,l) = vd(j,l) * vd(j,l)
           vdvd_ave(j,l) = vdvd(j,l) + vdvd_ave(j,l)
	!-------------------------------------------

        end do
     end do

		do j=2,jmax-1
		do l=2,lmax-1


	!--------vorticity fluctuation value-----------------
	
			dudx = -( ud(j+1,l) - ud(j-1,l) ) / (2.0d0*dx)
			dudy = -( ud(j,l) - ud(j,l) ) / (2.0d0*dy)
			dudz = -( ud(j,l+1) - ud(j,l-1) ) / (2.0d0*dz)

			dvdx = -( vd(j+1,l) - vd(j-1,l) ) / (2.0d0*dx)
			dvdy = -( vd(j,l) - vd(j,l) ) / (2.0d0*dy)
			dvdz = -( vd(j,l+1) - vd(j,l-1) ) / (2.0d0*dz)

			dwdx = -( wd(j+1,l) - wd(j-1,l) ) / (2.0d0*dx)
			dwdy = -( wd(j,l) - wd(j,l) ) / (2.0d0*dy)
			dwdz = -( wd(j,l+1) - wd(j,l-1) ) / (2.0d0*dz)

			vortxd(j,l) = dwdy - dvdz
			vortyd(j,l) = dudz - dwdx
			vortzd(j,l) = dvdx - dudy

	!---------------------------------------------


                end do
                end do

!		do l=1,lmax
!		do j=1+1,jmax-1


	!-------pressure fluctuation gradient

	!	dpdx(j,l) = -( pd(j+1,l) - 2.0d0*pd(j,l) + pd(j-1,l) ) / (dx*dx)
	!	d2pdx2(j,l) = -( pd(j+1,l) - pd(j-1,l) ) / (2.0d0*dx)

	!-----------------------------------

!		end do
!		end do


	!----outpostkun
  
!    write (filename,'(A,"flow_z00001_",i8.8)')trim(adjustl(outfile)),ii
!	filename = ''//trim(adjustl(filename))//''
!    open(55,file=''//trim(adjustl(filename))//'',form='unformatted',status='unknown')
  

!	open(55,file=filename,form="unformatted")

                do j=1,jmax
                   do l=1,lmax
                      
                      q(j,1,l,1) = vortxd(j,l) !fn1001 omegax'
                      q(j,1,l,2) = vortyd(j,l) !fn1002 omegay'
                      q(j,1,l,3) = vortzd(j,l) !fn1003 omegaz'
                      q(j,1,l,4) = dpdx(j,l) !fn1004 dp'dx
                      q(j,1,l,5) = d2pdx2(j,l) !fn1005 d2p'dx2!!
                      q(j,1,l,6) = pd(j,l) !fn1006 p'
                      q(j,1,l,7) = rhod(j,l) !fn1007 rho'
                      q(j,1,l,8) = rho(j,l)
                      q(j,1,l,9) = u(j,l)
                      q(j,1,l,10) = v(j,l)
                      q(j,1,l,11) = w(j,l)
                      q(j,1,l,12) = p(j,l)
		
                      q(j,1,l,13) = ptotal_all(j,l) + q(j,1,l,13)	
                      q(j,1,l,14) = ptotal_turb(j,l) + q(j,1,l,14)
                      q(j,1,l,15) = ptotal_lami(j,l) + q(j,1,l,15)
                      
                      q(j,1,l,16) = aimp_all(j,l) + q(j,1,l,16)	
                      q(j,1,l,17) = aimp_turb(j,l) + q(j,1,l,17)
                      q(j,1,l,18) = aimp_lami(j,l) + q(j,1,l,18)
	
                      q(j,1,l,19) = kinetic_e(j,l) + q(j,1,l,19)
                      !    q(j,1,l,17) = aimp(j,l) + q(j,1,l,17)
                      q(j,1,l,20) = tke(j,l) + q(j,1,l,20)

                   end do
                end do
	
        
	
	write(55) jmax,kmax,lmax
!	write(55) 7,0
!	do n=1,7
!	write(55)  (((real(q(j,1,l,n)),j=1,jmax),k=1,kmax),l=1,lmax)
!	end do 


!-----------------------------------------------------------------
!------------for gnuplot----------------------------------
!-----------------------------------------------------------------
 
 if(ii.gt.500) then
    cycle
    
 else 
    

    !!make grid!!
    write(filename,'("./gnuplot/data_rho/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,8)
          
       end do
       write(75,*)
       
    end do

    close(75)


    write(filename,'("./gnuplot/data_u/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,9)

       end do
       write(75,*)
       
    end do

    close(75)

    
    write(filename,'("./gnuplot/data_v/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,10)
          
       end do
       write(75,*)
       
    end do
    
    close(75)

    
    
    write(filename,'("./gnuplot/data_w/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,11)
          
       end do
       write(75,*)
       
    end do

    close(75)


    write(filename,'("./gnuplot/data_p/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,12)

       end do
       write(75,*)
       
    end do

    close(75)


    
    write(filename,'("./gnuplot/data_ptotal_all/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)

    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),ptotal_all(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)


  
    write(filename,'("./gnuplot/data_ptotal_turb/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),ptotal_turb(j,l)
          
       end do
       write(75,*)
       
    end do
 
    close(75)


  
    write(filename,'("./gnuplot/data_ptotal_lami/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),ptotal_lami(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)
    

  
    write(filename,'("./gnuplot/data_kinetic_e/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax

          write(75,*) x(j,l,1),z(j,l,1),kinetic_e(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)
  

  
    
    write(filename,'("./gnuplot/data_aimp_all/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),aimp_all(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)

    
  
    write(filename,'("./gnuplot/data_aimp_turb/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
  
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),aimp_turb(j,l)
        
       end do
       write(75,*)
       
    end do
    
    close(75)


    
    write(filename,'("./gnuplot/data_aimp_lami/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),aimp_lami(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)

    
    
    write(filename,'("./gnuplot/data_tke/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),tke(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)
 

    
    write(filename,'("./gnuplot/data_ptotal_tke/data/data_pm3d_",i8.8,"")')ii
    open(75,file=filename)
    
    do j=1,jmax
       do l=1,lmax
          
          write(75,*) x(j,l,1),z(j,l,1),p_total(j,l)+tke(j,l)
          
       end do
       write(75,*)
       
    end do
    
    close(75)
 


	end if
    enddo

     close(804)





!--------------------------------------------------------------------------
!-------------------------------------------------------------------------    
    
     q(:,1,:,13) = q(:,1,:,13) / dble(nbinst)
     q(:,1,:,14) = q(:,1,:,14) / dble(nbinst)
     q(:,1,:,15) = q(:,1,:,15) / dble(nbinst)
     q(:,1,:,16) = q(:,1,:,16) / dble(nbinst)    
     q(:,1,:,17) = q(:,1,:,17) / dble(nbinst)
     q(:,1,:,18) = q(:,1,:,18) / dble(nbinst)
     q(:,1,:,19) = q(:,1,:,19) / dble(nbinst)
     q(:,1,:,20) = q(:,1,:,20) / dble(nbinst)    
     

 	write(filename,'("./gnuplot/data_average/data_pm3d_rho")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),rho_ave(j,l)


		end do
                write(75,*)

	end do

	close(75)

 	write(filename,'("./gnuplot/data_average/data_pm3d_u")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),u_ave(j,l)


		end do
                write(75,*)

	end do

	close(75)


 	write(filename,'("./gnuplot/data_average/data_pm3d_v")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),v_ave(j,l)


		end do
                write(75,*)

	end do

	close(75)


 	write(filename,'("./gnuplot/data_average/data_pm3d_w")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),w_ave(j,l)


		end do
                write(75,*)

	end do

	close(75)


 	write(filename,'("./gnuplot/data_average/data_pm3d_p")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),p_ave(j,l)


		end do
                write(75,*)

	end do

	close(75)

	write(filename,'("./gnuplot/data_average/data_pm3d_p_total_all")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,13)


		end do
                write(75,*)

	end do

	close(75)

	write(filename,'("./gnuplot/data_average/data_pm3d_p_total_turb")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,14)


		end do
                write(75,*)

	end do

	close(75)


	write(filename,'("./gnuplot/data_average/data_pm3d_p_total_lami")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,15)


		end do
                write(75,*)

	end do

	close(75)


	write(filename,'("./gnuplot/data_average/data_pm3d_aimp_all")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,16)


		end do
                write(75,*)

	end do

	close(75)

	write(filename,'("./gnuplot/data_average/data_pm3d_aimp_turb")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,17)


		end do
                write(75,*)

	end do

	close(75)
	write(filename,'("./gnuplot/data_average/data_pm3d_aimp_lami")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,18)


		end do
                write(75,*)

	end do

	close(75)


	write(filename,'("./gnuplot/data_average/data_pm3d_kinetic_e")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,19)


		end do
                write(75,*)

	end do

	close(75)

	write(filename,'("./gnuplot/data_average/data_pm3d_tke")')
 	open(75,file=filename)

	do j=1,jmax
		do l=1,lmax

		write(75,*) x(j,l,1),z(j,l,1),q(j,1,l,20)


		end do
                write(75,*)

	end do

	close(75)



     deallocate(rg)
     deallocate(zg)
    
     deallocate(rhoT0)
     deallocate(urT0)
     deallocate(utT0)
     deallocate(uzT0)
     deallocate(pT0)

     deallocate(rhoT0moy)
     deallocate(urT0moy)
     deallocate(utT0moy)
     deallocate(uzT0moy)
     deallocate(pT0moy)
!     deallocate(vortT0moy)
!     deallocate(dilatT0moy)
    
    write(6,*) 'read data done'
    write(6,*) 'start gnuplot'

!-------------read instanaous parameter-------------------------------------

!----write data--------------------------------
call system("pwd") 
! call system("cd ./gnuplot/data_rho/plot/") 
 call system("cd ./gnuplot/data_rho/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_u/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_v/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_w/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_p/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_ptotal_all/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_ptotal_turb/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_ptotal_lami/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_kinetic_e/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_aimp_all/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_aimp_turb/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_tke/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_aimp_lami/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_ptotal_tke/plot/ ; gnuplot ./p_gif.plt") 

 call system("cd ./gnuplot/data_average/ ; gnuplot ./p_gif.plt") 


 call system("gnuplot ./gnuplot/data_average/p_gif.plt") 

    
    
    end subroutine
