      bdy = its-ims
      if ((bdy .ne. ime-ite).or.(bdy .ne. jme-jte).or.(bdy .ne. jts-jms)) then
              !write(*,*) ' something error about the halo width, set the halo being minimal width'
              write(*,*)ime-ite,jme-jte,jts-jms,its-ims
              bdy = min(bdy,ime-ite)
              bdy = min(bdy,jme-jte)
              bdy = min(bdy,jts-jms)
      endif
      bdy_h = bdy
      if(present(halo)) then
              bdy_h = halo
              if (bdy_h > bdy) then
                bdy_h = bdy
              endif
      endif
      !bdy_h=1

            update_ew_dir = .true.
            update_sn_dir = .true.

            !*******************************************************************
            ! for west-east halo exchange and period boundary:
            !-----------------------------------------------------------------
!write(*,*)'begin ew'
if(update_ew_dir) then
      if( myprcid == Eid .and. myprcid == Wid ) then
            bdy = its-ims
            do j=jts,jte
            do k=kms,kme
            do i=bdy,1,-1
            f(ids-i,k,j) = f(ide-i+1,k,j)
            enddo
            enddo
            enddo
            do j=jts,jte
            do k=kms,kme
            do i=1,bdy
            f(ide+i,k,j) = f(ids+i-1,k,j)
            enddo
            enddo
            enddo

      else
            ! send data from west to east:
            bdy_s = bdy_h
            len_s = bdy_s*(jte-jts+1)*(kme-kms+1)
            bdy_r = bdy_h
            len_r = len_s
            if( myprcid_b == 0 ) then
              bdy_r = its-ims
              len_r = bdy_r*(jte-jts+1)*(kme-kms+1)
            endif
            if( myprcid_b == prcnumb_b-1 ) then 
              bdy_s = its-ims
              len_s = bdy_s*(jte-jts+1)*(kme-kms+1)
            endif
            allocate(sbuf(bdy_s,kms:kme,jts:jte))
            allocate(rbuf(bdy_r,kms:kme,jts:jte))

            sbuf(1:bdy_s,kms:kme,jts:jte) = f(ite-bdy_s+1:ite,kms:kme,jts:jte)

            tag = 880
            call mpi_sendrecv(sbuf,len_s,mpidatatype,eid,tag,       &
              rbuf,len_r,mpidatatype,wid,tag,MPI_COMM_WORLD,status,ierr)

            f(its-bdy_r:its-1,kms:kme,jts:jte) = rbuf(1:bdy_r,kms:kme,jts:jte)

            deallocate(sbuf)
            deallocate(rbuf)

            ! send data from east to west:
            bdy_s = bdy_h
            len_s = bdy_s*(jte-jts+1)*(kme-kms+1)
            bdy_r = bdy_h
            len_r = len_s
            if( myprcid_b == 0 ) then
              bdy_s = its-ims
              len_s = bdy_s*(jte-jts+1)*(kme-kms+1)
            endif
            if( myprcid_b == prcnumb_b-1 ) then 
              bdy_r = its-ims
              len_r = bdy_r*(jte-jts+1)*(kme-kms+1)
            endif
            allocate(sbuf(bdy_s,kms:kme,jts:jte))
            allocate(rbuf(bdy_r,kms:kme,jts:jte))

            sbuf(1:bdy_s,kms:kme,jts:jte) = f(its:its+bdy_s-1,kms:kme,jts:jte)

            tag = 881
            call mpi_sendrecv(sbuf,len_s,mpidatatype,wid,tag,       &
              rbuf,len_r,mpidatatype,eid,tag,MPI_COMM_WORLD,status,ierr)

            f(ite+1:ite+bdy_r,kms:kme,jts:jte) = rbuf(1:bdy_r,kms:kme,jts:jte)

            deallocate(sbuf)
            deallocate(rbuf)
      endif
endif !update ew halo

          !*******************************************************************
          ! for north-south halo exchange and symmetric boundary:
          !-----------------------------------------------------------------
if(update_sn_dir) then
      nx = ide-ids+1
      !write(*,*)'ids=',ids,', ide=',ide,', nx=',nx
      if(mod(nx,2) /= 0) then
            stop 'grapes_abort'
      endif
     
            !-----------------------------------------------------------------
            ! south-north only one node  *************************************
      if( prcnumb_a ==1 ) then
            bdy = jts-jms
     
            ! east-west only one node -------------------*********************
            if ( prcnumb_b == 1) then
     
              !south polar
     
              do j=1,bdy 
              do i=ims,ime
              ii = mod(i+nx/2-1,nx) + 1
              f(i,:,jds-j) = f(ii,:,jds+j-1)
              enddo   
              enddo   
     
              !nouth polar
     
              do j=1,bdy
              do i=ims,ime
              ii = mod(i+nx/2-1,nx) + 1
              f(i,:,jde-1+j) = f(ii,:,jde-j)
              enddo
              enddo
     
              ! east-west more than one node -------------------*********************
            else
     
              len_sn = bdy*(ime-ims+1)*(kme-kms+1)
              allocate(sbuf(ims:ime,kms:kme,1:bdy))
              allocate(rbuf(ims:ime,kms:kme,1:bdy))
     
              ! north polar
              sbuf = 0.
              do j=1,bdy
              sbuf(:,:,j) = f(:,:,jde-j)
              enddo
              tag = 980
              call mpi_sendrecv(sbuf,len_sn,mpidatatype,nid,tag,       &
                rbuf,len_sn,mpidatatype,nid,tag,MPI_COMM_WORLD,status,ierr)
     
              do j=1,bdy
              f(:,:,jde-1+j) = rbuf(:,:,j)
              enddo
     
              ! south polar
              sbuf = 0.
              ! sbuf(:,:,1:bdy) = f(:,:,jds:jds+bdy-1)
     
              do j=1,bdy
              sbuf(:,:,j) = f(:,:,jds+j-1)
              enddo
     
     
              tag = 981
              call mpi_sendrecv(sbuf,len_sn,mpidatatype,sid,tag,       &
                rbuf,len_sn,mpidatatype,sid,tag,MPI_COMM_WORLD,status,ierr)
     
              !         if(present(isu).or.(present(vpolar).and.(vpolar))) then
     
              do j=1,bdy 
              f(:,:,jds-j) = rbuf(:,:,j)
              enddo   
              deallocate(sbuf)
              deallocate(rbuf)
     
            endif
     
            !---------------------------------------------------------------------
            ! south-north more then one node  ************************************
      else
            if( myprcid_a == prcnumb_a-1 ) then                  ! north node
              !exchange halo:
              len_sn = bdy_h*(ime-ims+1)*(kme-kms+1)
              allocate(sbuf(ims:ime,kms:kme,1:bdy_h))
              allocate(rbuf(ims:ime,kms:kme,1:bdy_h))
     
              sbuf(:,:,1:bdy_h) = f(:,:,jts:jts+bdy_h-1)
     
              tag = 891
              call mpi_send(sbuf,len_sn,mpidatatype,sid,tag,      &           ! send data to south node
                MPI_COMM_WORLD,ierr)
     
              tag = 892
              call mpi_recv(rbuf,len_sn,mpidatatype,sid,tag,      &           ! recv data from south node
                MPI_COMM_WORLD,status,ierr)
              f(:,:,jts-bdy_h:jts-1) = rbuf(:,:,1:bdy_h)
     
              deallocate(sbuf)
              deallocate(rbuf)
     
              !-----------------
              !for polar symmetric boundary:
              bdy = jts - jms
              if( myprcid == nid ) then                     ! east-west only one node
                do j=1,bdy
                do i=ims,ime
                ii = mod(i+nx/2-1,nx) + 1
                f(i,:,jde-1+j) = f(ii,:,jde-j)
                enddo
                enddo
     
              else                                          ! east-west more than one node
                len_sn = bdy*(ime-ims+1)*(kme-kms+1)
                allocate(sbuf(ims:ime,kms:kme,1:bdy))
                allocate(rbuf(ims:ime,kms:kme,1:bdy))
     
                do j=1,bdy
                sbuf(:,:,j) = f(:,:,jte-j)
                enddo
     
                tag = 893
                call mpi_sendrecv(sbuf,len_sn,mpidatatype,nid,tag,       &
                  rbuf,len_sn,mpidatatype,nid,tag,MPI_COMM_WORLD,status,ierr)
     
                do j=1,bdy
                f(:,:,jde-1+j) = rbuf(:,:,j)
                enddo
                deallocate(sbuf)
                deallocate(rbuf)
              endif
     
              !-----------------------------------------------------------------
            else if( myprcid_a /= 0 ) then                       ! middle node
              len_sn = bdy_h*(ime-ims+1)*(kme-kms+1)
              allocate(sbuf(ims:ime,kms:kme,1:bdy_h))
              allocate(rbuf(ims:ime,kms:kme,1:bdy_h))
     
              sbuf(:,:,1:bdy_h) = f(:,:,jts:jts+bdy_h-1)
     
              tag = 891
              call mpi_sendrecv(sbuf,len_sn,mpidatatype,sid,tag,       &      ! sendrecv data with south node
                rbuf,len_sn,mpidatatype,nid,tag,MPI_COMM_WORLD,status,ierr)
     
              f(:,:,jte+1:jte+bdy_h) = rbuf(:,:,1:bdy_h)
     
              !-------------
              sbuf(:,:,1:bdy_h) = f(:,:,jte-bdy_h+1:jte)
     
              tag = 892
              call mpi_sendrecv(sbuf,len_sn,mpidatatype,nid,tag,       &      ! sendrecv data with north node
                rbuf,len_sn,mpidatatype,sid,tag,MPI_COMM_WORLD,status,ierr)
     
              f(:,:,jts-bdy_h:jts-1) = rbuf(:,:,1:bdy_h)
     
              deallocate(sbuf)
              deallocate(rbuf)
     
              !-----------------------------------------------------------------
            else                                                  ! south node
              !exchange halo:
              len_sn = bdy_h*(ime-ims+1)*(kme-kms+1)
              allocate(sbuf(ims:ime,kms:kme,1:bdy_h))
              allocate(rbuf(ims:ime,kms:kme,1:bdy_h))
     
              tag = 891
              call mpi_recv(rbuf,len_sn,mpidatatype,nid,tag,      &          ! recv data from north node
                MPI_COMM_WORLD,status,ierr)
     
              sbuf(:,:,1:bdy_h) = f(:,:,jte-bdy_h+1:jte)
     
              tag = 892
              call mpi_send(sbuf,len_sn,mpidatatype,nid,tag,       &
                MPI_COMM_WORLD,ierr)                             ! send data to north node
     
              f(:,:,jte+1:jte+bdy_h) = rbuf(:,:,1:bdy_h)
     
              deallocate(sbuf)
              deallocate(rbuf)
     
              !for polar symmetric boundary:
              bdy = jts - jms
              if( myprcid == sid ) then                      ! east-west only one node
                !            if(present(isu).or.(present(vpolar).and.(vpolar))) then
                do j=1,bdy 
                do i=ims,ime
                ii = mod(i+nx/2-1,nx) + 1
                f(i,:,jds-j) = f(ii,:,jds+j-1)
                enddo   
                enddo   
     
              else                                           ! east-west more than one node
                len_sn = bdy*(ime-ims+1)*(kme-kms+1)
                allocate(sbuf(ims:ime,kms:kme,1:bdy))
                allocate(rbuf(ims:ime,kms:kme,1:bdy))
     
                sbuf = 0.
                !  sbuf(:,:,1:bdy) = f(:,:,jts:jts+bdy-1)
                do j=1,bdy
                sbuf(:,:,j) = f(:,:,jts+j-1)
                enddo
     
     
                tag = 890
                call mpi_sendrecv(sbuf,len_sn,mpidatatype,sid,tag,       &
                  rbuf,len_sn,mpidatatype,sid,tag,MPI_COMM_WORLD,status,ierr)
     
                !            if(present(isu).or.(present(vpolar).and.(vpolar))) then
                do j=1,bdy 
                f(:,:,jts-j) = rbuf(:,:,j)
                enddo   
                deallocate(sbuf)
                deallocate(rbuf)
              endif
            endif
     
            !-----------------------------------------------------------------
      endif
endif  !update sn halo

