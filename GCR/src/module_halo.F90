module module_halo
        use h1
        use mpi
        interface glob_updatehalo_inc
            module procedure glob_updatehalo_inc1
            module procedure glob_updatehalo_inc2
        end interface
        interface grapes_sum_ew
            module procedure grapes_allsum_ew_real_1d
            module procedure grapes_allsum_ew_double_1d
        end interface
    
        INTEGER, PARAMETER :: UPDATE_ALL = 0
        INTEGER, PARAMETER :: EW_ONLY = 1
        INTEGER, PARAMETER :: SN_ONLY = 2

!#define tmp


contains

      subroutine glob_updatehalo_real_3d(f,m,iter_max)  
          implicit none
          logical ,external  ::    On_Monitor
          !type( type_model_partition ) :: grid 
          real(kind = 4) , dimension(its-1:ite+1,kts-1:kte+1,jts-1:jte+1,0:iter_max-1),intent(inout)     :: f
          integer  :: halo=1, iter_max,m
          !logical,optional,intent(in)  :: vpolar, isu, diagv
          !integer,optional,intent(in)  :: update_dir
          !local var:

          !include 'glob_updatehalo_3d.inc'
          call  glob_updatehalo_inc(f(:,:,:,m),halo)
        end subroutine glob_updatehalo_real_3d

        subroutine glob_updatehalo_double_3d(f)

            !    !, vpolar,isu,diagv,update_dir)  
            implicit none
            logical ,external  ::    On_Monitor
            !type( type_model_partition ) :: grid 
            !real(kind = 8) , dimension(grid%ims:grid%ime,grid%kms:grid%kme,grid%jms:grid%jme),intent(inout)     :: f
            real(kind = 8) , dimension(its:ime,kts-1:kte+1,jms:jme),intent(inout)     :: f
            integer  :: halo=1
            !local var:

            ! write(*,*)'mytid=',myprcid,',jte=',jte
            !include 'glob_updatehalo_3d.inc'
            call glob_updatehalo_inc(f,halo)
          end subroutine glob_updatehalo_double_3d


          subroutine glob_updatehalo_inc1(f,halo,vpolar,isu,diagv,update_dir)
            use h1
            logical,optional,intent(in)  :: vpolar, isu, diagv
            integer,optional,intent(in)  :: update_dir
            integer,optional,intent(in)  :: halo
            real(kind = 4) , dimension(its-1:ite+1,kts-1:kte+1,jts-1:jte+1),intent(inout)     :: f
            !local var:
            real(kind = 4) , dimension(:), allocatable     :: ws, wr, es, er, ss, sr, ns, nr
            real(kind = 4) , dimension(:,:,:), allocatable :: sbuf, rbuf
            real(kind = 4) , dimension(:,:,:), allocatable :: buf
            real(kind = 4) , dimension(kts-1:kte+1)  :: u0, v0
            integer :: len_s, len_r, len_sn, bdy_h, bdy, bdy_s, bdy_r
            integer :: ierr, status(MPI_STATUS_SIZE),tag
            integer :: i,j,k,idx,nx,ii
            integer :: mpidatatype
            logical :: update_ew_dir
            logical :: update_sn_dir 
            integer:: iims,iime,jjms,jjme
            iims=ims
            iime=ime
            jjms=jms
            jjme=jme
            ims=its-1
            ime=ite+1
            jms=jts-1
            jme=jte+1
            mpidatatype = mpi_real
            myprcid=mytid(1)
            Sid=mytid(2)
            Nid=mytid(3)
            Wid=mytid(4)
            Eid=mytid(5)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)

        include '../inc/glob_updatehalo.inc'
            ims=iims
            ime=iime
            jms=jjms
            jme=jjme
        
        end SUBROUTINE

        subroutine glob_updatehalo_inc2(f,halo,vpolar,isu,diagv,update_dir)
          use h1
          logical,optional,intent(in)  :: vpolar, isu, diagv
          integer,optional,intent(in)  :: update_dir
          integer,optional,intent(in)  :: halo
          real(kind = 8) , dimension(its-1:ite+1,kts-1:kte+1,jts-1:jte+1),intent(inout)     :: f
          !local var:
          real(kind = 8) , dimension(:), allocatable     :: ws, wr, es, er, ss, sr, ns, nr
          real(kind = 8) , dimension(:,:,:), allocatable :: sbuf, rbuf
          real(kind = 8) , dimension(:,:,:), allocatable :: buf
          real(kind = 8) , dimension(kts-1:kte+1)  :: u0, v0
          integer :: len_s, len_r, len_sn, bdy_h, bdy, bdy_s, bdy_r
          integer :: ierr, status(MPI_STATUS_SIZE),tag
          integer :: i,j,k,idx,nx,ii
          integer :: mpidatatype
          logical :: update_ew_dir 
          logical :: update_sn_dir 
          mpidatatype = mpi_real8
          myprcid=mytid(1)
          Sid=mytid(2)
          Nid=mytid(3)
          Wid=mytid(4)
          Eid=mytid(5)
!#ifdef tmp
          !include 'tmp.inc'
          include '../inc/glob_updatehalo.inc'
!#else
!          include 'inc/glob_updatehalo_3d.inc'
!#endif
      end SUBROUTINE
      SUBROUTINE grapes_abort
          USE h1, ONLY : local_communicator
          IMPLICIT NONE
          INTEGER ierr
          CALL mpi_abort(local_communicator,1,ierr)
          stop
      END SUBROUTINE grapes_abort
      SUBROUTINE grapes_allsum_ew_real_1d(Data,n)
          IMPLICIT NONE
          integer ,intent(in) :: n
          real(kind=4),intent(inout):: data(n)
          real(kind=4):: data_tmp(n)
          INTEGER ierr
         
          CALL MPI_ALLREDUCE(data,data_tmp,n,MPI_REAL, MPI_SUM,EWcomm,ierr)
          data=data_tmp
 
      END SUBROUTINE grapes_allsum_ew_real_1d

      SUBROUTINE grapes_allsum_ew_double_1d(Data,n)
          IMPLICIT NONE
          integer ,intent(in) :: n
          real(kind=8) ,intent(inout):: data(n)
          real(kind=8) :: data_tmp(n)
          INTEGER ierr
       
          CALL MPI_ALLREDUCE(data,data_tmp,n,MPI_REAL8, MPI_SUM,EWcomm,ierr)
          data=data_tmp

      END SUBROUTINE grapes_allsum_ew_double_1d

!      SUBROUTINE  svrasr( yl, m,iter_max,b,  NI, NK, NJ)
!!
!
!      IMPLICIT NONE
!      INTEGER    NI,NK,NJ,NG1,ILUK,NFWPP
!      integer :: i,j,k,halo,m,iter_max
!      real :: yl(ni,nk,nj,iter_max),b(ni,nk,nj,7)
!   
!        !!$OMP PARALLEL
!        !!$OMP DO PRIVATE(j,i,k)
!        do j=1,nj
!          do i=2,ni
!            yl(i,1,j,m)=yl(i,1,j,m)-b(i,1,j,2)*yl(i-1,1,j,m)
!          enddo
!           !k>1
!          do k=2,nk
!            do i=1,ni
!              yl(i,k,j,m)=yl(i,k,j,m)-b(i,k,j,4)*yl(i,k-1,j,m)
!            enddo
!
!            do i=2,ni
!              yl(i,k,j,m)=yl(i,k,j,m)-b(i,k,j,2)*yl(i-1,k,j,m)
!            enddo
!          enddo
!        enddo
!        !!$OMP END DO
!
!          !j=1
!           !k=1
!           !!$OMP MASTER
!           yl(ni,nk,nj,m)=yl(ni,nk,nj,m)*b(ni,nk,nj,1)
!           !!$OMP END MASTER
!           !!$OMP BARRIER
!           !!$OMP DO PRIVATE(j,i,k)
!         do j=nj,1,-1
!           do i=ni-1,1,-1
!             yl(i,nk,j,m)=(yl(i,nk,j,m)-b(i,nk,j,3)*yl(i+1,nk,j,m))*b(i,nk,j,1)
!           enddo
!           !k>1
!           do k=nk-1,1,-1
!             do i=1,ni
!               yl(i,k,j,m)=yl(i,k,j,m)-b(i,k,j,5)*yl(i,k+1,j,m)
!             enddo
!             yl(ni,k,j,m)=yl(ni,k,j,m)*b(ni,k,j,1)
!             do i=ni-1,1,-1
!                yl(i,k,j,m)=(yl(i,k,j,m)-b(i,k,j,3)*yl(i+1,k,j,m))*b(i,k,j,1)
!             enddo
!           enddo
!         enddo
!         !!$OMP END DO NOWAIT
!         !!$OMP END PARALLEL
!
!      END SUBROUTINE  SVRASR
      subroutine getijk(ijk_index)
        use h1
        integer::ijk_index(20)
        
        idep=ijk_index(1)
	jdep=ijk_index(2)
	ids=ijk_index(3)
	ide=ijk_index(4)
	jds=ijk_index(5)
	jde=ijk_index(6)
	kds=ijk_index(7)
	kde=ijk_index(8)
	
        ims=ijk_index(9)
	ime=ijk_index(10)
	jms=ijk_index(11)
	jme=ijk_index(12)
	kms=ijk_index(13)
	kme=ijk_index(14)
	
        its=ijk_index(15)
	ite=ijk_index(16)
	jts=ijk_index(17)
	jte=ijk_index(18)
	kts=ijk_index(19)
	kte=ijk_index(20)
	
      end subroutine
end module
