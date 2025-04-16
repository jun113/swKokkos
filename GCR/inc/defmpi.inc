        CALL MPI_INIT(ierr)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_COMM_RANK (MPI_COMM_WORLD, myprcid, ierr)
        CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nproc, ierr)
         isperiodic(1) = .false.
         isperiodic(2) = .true.
         dims(1)=prcnumb_a
         dims(2)=prcnumb_b
        CALL mpi_cart_create(MPI_COMM_WORLD, 2, dims, isperiodic, .false., new_comm, ierr )
        CALL mpi_comm_rank(new_comm,myprcid, ierr )
        CALL mpi_cart_coords(new_comm,myprcid, 2, coords, ierr )
        myprcid_a = coords(1)   ! row task (y)
        myprcid_b = coords(2)   ! col task (x)

    allocate ( procmap(0:prcnumb_a-1,0:prcnumb_b-1))
        DO ii=0,nproc-1
            CALL mpi_cart_coords(new_comm, ii, 2, coords, ierr )
            procmap(coords(1),coords(2)) = ii
        ENDDO
        CALL mpi_cart_shift(new_comm,0,1,Sid,Nid,ierr)
        CALL mpi_cart_shift(new_comm,1,1,Wid,Eid,ierr)
        if ( global_opt ) then
            if (  myprcid_a  == prcnumb_a-1 ) THEN
               jj = mod( (prcnumb_b/2+myprcid_b), prcnumb_b)
               Nid = procmap(  myprcid_a , jj ) 
            end if
            if (  myprcid_a == 0 ) THEN
               jj = mod( (prcnumb_b/2+myprcid_b), prcnumb_b)
               Sid = procmap(  myprcid_a , jj ) 
            end if
        else
            if (  myprcid_a  == prcnumb_a-1 ) THEN
               Nid = MPI_PROC_NULL 
            end if
            if (  myprcid_a == 0 ) THEN
               Sid = MPI_PROC_NULL
            end if
        end if
         
!    allocate(NSmember(0:prcnumb_a-1))
!    allocate(EWmember(0:prcnumb_b-1))
!    allocate(NScomm_lat(0:prcnumb_b-1))
!    allocate(EWcomm_lon(0:prcnumb_a-1))
!    allocate(NSgroup_lat(0:prcnumb_b-1))
!    allocate(EWgroup_lon(0:prcnumb_a-1))
!
!call MPI_Comm_group(new_comm, world_group,ierr)
!
!    do j=0,prcnumb_a-1
!      do i=0,prcnumb_b-1
!        EWmember(i)=procmap(j,i)
!      end do
!      call MPI_group_incl(world_group, prcnumb_b, EWmember , EWgroup_lon(j),ierr)
!      call MPI_Comm_create(new_comm,EWgroup_lon(j),EWcomm_lon(j),ierr)
!    end do
!    EWgroup = EWgroup_lon(myprcid_a)
!    EWcomm = EWcomm_lon(myprcid_a)
!
!    do j=0,prcnumb_b-1
!      do i=0,prcnumb_a-1
!        NSmember(i)=procmap(i,j)
!      end do
!      call MPI_group_incl(world_group, prcnumb_a, NSmember , NSgroup_lat(j),ierr)
!      call MPI_Comm_create(new_comm,NSgroup_lat(j),NScomm_lat(j),ierr)
!    end do
!    NSgroup = NSgroup_lat(myprcid_b)
!    NScomm = NScomm_lat(myprcid_b)
!
!
!    local_communicator = new_comm
!    call MPI_Comm_rank(NScomm,i,ierr)
!    call MPI_Comm_rank(EWcomm,i,ierr)
!    call mpi_COMM_SIZE(EWCOMM,ewcommsize,ierr)
!    call mpi_COMM_SIZE(NSCOMM,nscommsize,ierr)


        mytid(1)=myprcid
        mytid(2)=Sid
        mytid(3)=Nid
        mytid(4)=Wid
        mytid(5)=Eid
        mytid(6)=myprcid_a
        mytid(7)= myprcid_b
        mytid(8)=prcnumb_a
        mytid(9)= prcnumb_b


