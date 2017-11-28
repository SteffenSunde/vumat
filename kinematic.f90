! Abaqus Explicit material model.
! Von Mises yield criterion, linear kinematic hardening
! Written in fortran 90/95: Add /free compiler flag in abaqus environment file.
! Draft model for fretting fatigue crack initiation
! (Compiles, converges and overlaps with abaqus intrinsic)
! Todo:
! - Check energies
! - Add damage parameter for elastic and plastic damage
! - Clean up code
! - Add more test cases

include 'modules.f90'

SUBROUTINE vumat ( &
  nBlock, nDir, nShr, nStateV, nFieldV, nProps, lAnneal, stepTime, &
  totalTime, dt, cmName, coordMp, charLength, props, density, &
  strainInc, relSpinInc, tempOld, stretchOld, defGradOld, fieldOld, &
  stressOld, stateOld, enerInternOld, enerInelasOld, tempNew, &
  stretchNew, defGradNew, fieldNew, &
  stressNew, stateNew, enerInternNew, enerInelasNew)

  use isohard
  use utilities, only: dp
  use ifport, only: sleepqq
  INCLUDE 'VABA_PARAM.INC'

  DIMENSION props(nProps), density(nBlock), coordMp(nBlock,*), &
  charLength(*), strainInc(nBlock, nDir+nShr), relSpinInc(*), &
  tempOld(*), stretchOld(*), defGradOld(*), fieldOld(*), &
  stressOld(nBlock, nDir+nShr), stateOld(nBlock, nStateV), &
  enerInternOld(nBlock), enerInelasOld(nBlock), tempNew(*), &
  stretchNew(*), defGradNew(*), fieldNew(*), &
  stressNew(nBlock, nDir+nShr), stateNew(nBlock, nStateV), &
  enerInternNew(nblock), enerInelasNew(nBlock)

  CHARACTER(80) CMNAME
  real(dp) :: C(3)
  integer :: i,iter, maxIter
  real(dp) :: eq_stress(6)
  real(dp) :: back_stress(6)
  real(dp) :: dLambda, ddLambda, phi, E, v, sigma_y, H, tol, mu, dfdl
  real(dp) :: properties(nProps)
  real(dp) :: eldmg, pldmg ! Not implemented yet
  real(dp) :: de(nDir+nShr)
  real(dp), dimension(nDir+nShr) :: stress_old, strain_incr, stress, stress_tr, dfds
  real(dp), dimension(nDir+nShr) :: state_old, n_hat, eta
  real(dp) :: eq_stress_norm

  ! For debugging purposes...
  character(*), parameter :: path = "set absolute path"

  ! Workaround for double precision
  properties = real(props, kind=dp)

  ! Read parameters from ABAQUS material card
  E = properties(1)
  v = properties(2)
  sigma_y = properties(3)
  H = properties(4)
  tol = properties(5)
  maxIter = properties(6)


! Compute elasticity matrix
  C = isotropic_stiffness(E, v)
  mu = E / (1-v) / 2.0

  ! Loop over integration points
  do i=1,nBlock
    ! Work-around due to abaqus using both single and double precision
    de = real(strainInc(i, :), kind=dp)
    state_old = real(stateOld(i, :), kind=dp)
    stress_old = real(stressOld(i, :), kind=dp)

    if (totalTime == 0.0) then ! Data check
      ! Elastic operator for abaqus data check (wave speeds etc.)
      stressNew(i,:) = stress_old + elastic_increment(C, de)

      !call sleepqq(500) ! For debugging attachment: sleep for 500 milliseconds
    else
      ! Elastic predictor
      stress_tr = stress_old + elastic_increment(C, de)
      !back_stress = deviatoric(state_old(1:6))
      !back_stress_norm = sqrt(1.d5) * l2norm(back_stress)

      ! Calculate elastic damage scalar
      eldmg = state_old(7)
      pldmg = state_old(8)

      ! Debugging in file. Slows down analysis considerably
      open(unit=19, file=path//"debug.txt", action="WRITE")

      ! "Equivalent" stress
      eq_stress = deviatoric(stress_tr) - back_stress

      ! Calculate Von Mises equivalent stress
      phi = sqrt(1.5)*l2norm(eq_stress)

      ! Yield function
      F = phi - sigma_y

      ! Plastic corrector
      if (F > 0.0) then ! Yielding
        ! Initialize plastic multiplier
        dLambda = F / (3.0*mu + H)

        ! Iteration variables
        iter = 0
        error = abs(F/sigma_y)

        stress = stress_tr

        ! Newton iterations on plastic multiplier
        do while (error > tol)
          if (iter < maxIter) then ! Perform iteration on plastic multiplier
            ! Yield surface unit normal
            n_hat =  sqrt(1.5) * eq_stress / l2norm(eq_stress)
            !n_hat = sqrt(1.5) * eq_stress / sigma_y

            dfdl = -3.0*mu - H
            ddLambda = F/dfdl
            dLambda = dLambda - ddLambda

            ! Update back stress
            back_stress = state_old(1:6) + h * dLambda * n_hat

            ! Update stresses (Isotropic stiffness)
            stress = stress_tr - 2.0*mu*dLambda*n_hat

            ! New yield function
            eq_stress = deviatoric(stress) - back_stress
            phi = sqrt(1.5)*l2norm(eq_stress)
            F = phi - sigma_y

            ! Update iteration variables
            error = abs(F/sigma_y)
            iter = iter + 1

            ! For debugging purposes. Slows things down considerably
            write(19, *) iter, ",", F
            !write(*,*) iter, F, dLambda
            ! if (iter == 20) then
            !   write(*,*) "WARNING: 20 iterations reached!"
            !   write(*,*) n_hat
            ! end if
          else ! Unable to reach convergence. Print info and stop analysis
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(*,*) "Error: Did not converge in ", maxIter, " iterations!"
            write(*,*) "Integration point: ", i
            write(*,*) "Error: ", error
            write(*,*) "Number of iterations: ", iter
            write(*,*) "Plastic multiplier: ", dLambda
            write(*,*) "Step time: ", stepTime
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

            stop
          end if ! End of check if max iterations reached
        end do ! End loop until equilibrium

        close(19) ! Close file for debugging

        stressNew(i, :) = stress
        stateNew(i, 1:6) = back_stress
        stateNew(i, 7) = eldmg
        stateNew(i, 8) = pldmg
        stateNew(i, 9) = iter
        !write(*,*) "Integration point, iterations"
        !write(*,*) i, iter

        ! TODO TODO TODO
        ! Update energy
      else ! Not yielding
        ! Update stresses
        stressNew(i, :) = stress_tr

        ! Update state dependent variables
        stateNew(i, 1:6) = stateOld(i, 1:6) ! No movement of yield surface
        stateNew(i, 7) = stateOld(i, 7) + calc_eldmg(stress_tr) !
        stateNew(i, 8) = 0.d0 ! No plastic damage
        stateNew(i, 9) = 0 ! number of iterations
        ! No need for p ? Yield surf. pos. is determined by stateNew(1:6)
      end if ! End check for yield
    end if ! End check for time being 0.0
end do ! End loop over integration points

return
end subroutine vumat
