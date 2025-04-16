        PROGRAM main
        use h1
        implicit none
        INCLUDE 'mpif.h'
        !character*14 ::postname
        integer::ii,jj,i,j,ierr
        
        include '../inc/defmpi.inc'
        !write(postname, '(a,i2.2,a)') "postvar",myprcid,".txt"
        CALL mygetvar_c(mytid)
        CALL MPI_FINALIZE ( ierr )
        
        END PROGRAM
