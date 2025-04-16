module h1
    integer :: local_communicator
        integer:: myprcid,mytid(9), nproc,  new_comm, EWcomm, World_group, EWgroup
        logical, dimension(2) :: isperiodic
        integer, dimension(2) :: dims, coords
        integer,parameter:: prcnumb_a=3,prcnumb_b=2!a=n_procy,b=n_procx=4
        !integer:: prcnumb_a,prcnumb_b!a=n_procy,b=n_procx=4
        integer:: myprcid_a,myprcid_b,Sid,Nid,Wid,Eid
        logical::global_opt=.true.
        integer,dimension(:),allocatable:: NSmember,EWmember,NScomm_lat,EWcomm_lon,NSgroup_lat,EWgroup_lon
        integer,dimension(:,:),allocatable::procmap,NSgroup,NScomm,nscommsize,ewcommsize
        integer:: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme, its, ite, kts, kte, jts, jte
        integer, parameter,public ::r8 = selected_real_kind(13),r4 =selected_real_kind(6),r3=selected_real_kind(5)
end module
