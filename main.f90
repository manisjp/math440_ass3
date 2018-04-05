program main 

	use mpi

	! DECLARES VARIABLES
	implicit none
	integer :: ierr, rank, i, k, P, J, L, i_max = 0, i_min = 0, k_max = 0, k_min = 0
	character(len=32) :: arg
	double precision :: stime, etime, C_max = -huge(0.0d0), C_min = huge(0.0d0), core_max, core_min, S, rand_num, r_sum
	double precision, dimension(:), allocatable :: rk, w, xi, xi_all
	double precision, dimension(:,:), allocatable :: X, C, X_all, C_all
	double precision, dimension(2) :: maxvals, minvals
	integer, dimension(mpi_status_size) :: mystatus

	! INITIALIZES MPI AND GETS RANK NUMBERS AND NUMBER OF CORES
	call mpi_init(ierr)
	stime = mpi_wtime()
	call mpi_comm_rank(mpi_comm_world,rank,ierr)
	call mpi_comm_size(mpi_comm_world,P,ierr)

	! BEGIN MAIN PROGRAM

	! GETS THE COMMAND LINE ARGUMENT FOR J AND COMPUTES L
	call getarg(1,arg)
	read(arg,*) J
	L = J*P

	! CREATES THE RANDOM WEIGHT VECTOR rk
	allocate(rk(1:J))
	do i = 1, J
		call random_number(rand_num)
		rk(i) = real(rand_num*(real(i*(rank+1)+1)-real(rank+1)/real(i)))+real(rank+1)/real(i)
	enddo

	if (rank == 0) then
		! GATHERS THE rk AND CREATES UNNORMALIZED w
		allocate(w(0))
		w = (/ w, rk /)
		do i = 1, P-1
			call mpi_recv(rk,J,mpi_double_precision,i,0,mpi_comm_world,mystatus,ierr)
			w = (/ w, rk /)
		enddo
		! CREATES NORMALIZED w
		w = real(w)/real(sum(w))
	else
		! SEND rk TO MASTER
		call mpi_send(rk,J,mpi_double_precision,0,0,mpi_comm_world,ierr)
		allocate(w(L))
	endif
	! BROADCASTS w AND ABORTS IF ERROR IS TOO LARGE
	call mpi_bcast(w,L,mpi_double_precision,0,mpi_comm_world,ierr)
	if (abs(real(sum(w)) - 1) > EPSILON(1.0d0)) then
		if (rank == 0) then
			write(*,*) "Weighted vector sum was not 1."
		endif
		call mpi_abort(mpi_comm_world,1,ierr)
	endif

	! CREATES THE ROWS OF X SEPERATELY
	allocate(X(J*rank+1:J*(rank+1),L),C(J*rank+1:J*(rank+1),L),xi(J*rank+1:J*(rank+1)))
	do i = J*rank+1, J*(rank+1)
		do k = 1, L
			call random_number(rand_num)
			X(i,k) = (rand_num*real(real(i)/real(k)+real(i*k)))-real(i*k)
		enddo
		xi(i) = sum(X(i,:))/real(L)
	enddo
	! STITCHES TOGETHER THE INDIVIDUAL X MATRICES AND xi VECTORS IN ORDER TO COMPUTE C
	allocate(X_all(L,L))
	if (rank == 0) then
		allocate(xi_all(0))
		X_all(1:J,1:L) = X
		xi_all = (/ xi_all, xi /)
		do i = 1, P-1
			call mpi_recv(X,J*L,mpi_double_precision,i,0,mpi_comm_world,mystatus,ierr)
			call mpi_recv(xi,J,mpi_double_precision,i,0,mpi_comm_world,mystatus,ierr)
			X_all(J*i+1:J*(i+1),1:L) = X
			xi_all = (/ xi_all, xi /)
		enddo
	else
		allocate(xi_all(L))
		call mpi_send(X,J*L,mpi_double_precision,0,0,mpi_comm_world,ierr)
		call mpi_send(xi,J,mpi_double_precision,0,0,mpi_comm_world,ierr)
	endif
	! BROADCASTS THE FULL X MATRIX AND xi VECTOR
	call mpi_bcast(X_all,L*L,mpi_double_precision,0,mpi_comm_world,ierr)
	call mpi_bcast(xi_all,L,mpi_double_precision,0,mpi_comm_world,ierr)
	! CREATES THE ROWS OF C SEPARATELY
	do i = J*rank+1, J*(rank+1)
		do k = 1, L
			C(i,k) = sum(w*(X_all(i,:)-xi_all(i))*(X_all(k,:)-xi_all(k)))/(1.0d0-sum(w*w))
			if (C(i,k) > C_max) then
				C_max = C(i,k)
				i_max = i
				k_max = k
			elseif (C(i,k) < C_min) then
				C_min = C(i,k)
				i_min = i
				k_min = k
			endif
		enddo
	enddo
	write(*,*) C_max, C_min, sum(X_all(J*rank+1,:))
	! SEARCHES FOR THE LARGEST C_max AND C_min AND GETS rank OF CORRESPONDING CORES	
	call mpi_reduce((/ C_max, real(rank,kind=kind(0.0d0)) /),maxvals,2,mpi_2double_precision,mpi_maxloc,0,mpi_comm_world,ierr)
	call mpi_reduce((/ C_min, real(rank,kind=kind(0.0d0)) /),minvals,2,mpi_2double_precision,mpi_minloc,0,mpi_comm_world,ierr)
	call mpi_bcast(maxvals,2,mpi_double_precision,0,mpi_comm_world,ierr)
	call mpi_bcast(minvals,2,mpi_double_precision,0,mpi_comm_world,ierr)	
	
	! SENDS INDEX INFORMATION
	if (rank == int(maxvals(2))) then
		call mpi_send(i_max,1,mpi_integer,0,0,mpi_comm_world,ierr)
		call mpi_send(k_max,1,mpi_integer,0,0,mpi_comm_world,ierr)
	endif
	if (rank == int(minvals(2))) then
		call mpi_send(i_min,1,mpi_integer,0,0,mpi_comm_world,ierr)
		call mpi_send(k_min,1,mpi_integer,0,0,mpi_comm_world,ierr)
	endif
	
	! GATHERS INDEX INFORMATION AND PRINTS OUTPUT
	if (rank == 0) then
		call mpi_recv(i_max,1,mpi_integer,int(maxvals(2)),0,mpi_comm_world,mystatus,ierr)
		call mpi_recv(k_max,1,mpi_integer,int(maxvals(2)),0,mpi_comm_world,mystatus,ierr)
		call mpi_recv(i_min,1,mpi_integer,int(minvals(2)),0,mpi_comm_world,mystatus,ierr)
		call mpi_recv(k_min,1,mpi_integer,int(minvals(2)),0,mpi_comm_world,mystatus,ierr)
		write(*,"(A14,I5,A1,I5,A1,I5,A1)")   "(P, J, L) = (", P, ",", J, ",", L, ")"
		write(*,"(A22,ES23.15E3,A1,I22,A1)") "(C_max, core_max) = (", maxvals(1), ",", int(maxvals(2)), ")"
		write(*,"(A22,ES23.15E3,A1,I22,A1)") "(C_min, core_min) = (", minvals(1), ",", int(minvals(2)), ")"
		write(*,"(A22,I23,A1,I22,A1)") "(i_max, j_max)    = (", i_max, ",", k_max, ")"
		write(*,"(A22,I23,A1,I22,A1)") "(i_min, j_min)    = (", i_min, ",", k_min, ")"
	endif 
	! GATHERS AND CREATES MATRIX C THEN PRINTS IT (SPECIAL CASE)
	allocate(C_all(L,L))
	if (J == 2 .and. P == 4) then
		if (rank == 0) then
			C_all(1:J,1:L) = C
			do i = 1, P-1
				call mpi_recv(C,J*L,mpi_double_precision,i,0,mpi_comm_world,mystatus,ierr)
				C_all(J*i+1:J*(i+1),1:L) = C
			enddo
			write(*,*) "C matrix:"
			do i = 1, L
				write(*,*) C_all(i,:)
			enddo
		else
			call mpi_send(C,J*L,mpi_double_precision,0,0,mpi_comm_world,ierr)
		endif
	endif
	
	deallocate(rk,w,xi,xi_all,X,C,X_all,C_all)

	! END MAIN PROGRAM

	! FINALIZES MPI AND PRINTS TIME TO COMPUTE
	etime = mpi_wtime()
	if (rank == 0) then
		write(*,*) "Time taken: ", etime-stime
	endif
	call mpi_finalize(ierr)
	
endprogram main
