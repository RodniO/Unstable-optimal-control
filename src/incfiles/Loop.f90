
!Constructs a chain with 1900 segments from (0;1) to (1;0) of shortest length with integral of x*y under the curve equal to 0
!Integral under the curve is denoted by z
subroutine Loop(verbose)
  USE ModVec
  Integer(4), intent(in) :: verbose !3 - print everything, 2 - print only results, 1 - do not print anything from UpdateCurve, 0 - print nothing
  Integer(4) ch !curve direction
  Integer(4) mult !multiplier for z in loss function (used for interpolation)
  Integer(4) i, j !for loop indices
  Double precision A !curve length lower bound (should be close to solution)
  Double precision curx, cury, tmp !temporaries
  Double precision dsecnd, dstmp !time calculation
  Type(Vector) sf !start and final values of x, y and z
  Type(Vector) mults !multipliers from 1 to 10^7
  Type(Vector) vz !interpolation data for z
  
  !DATA FILE NOTATION
  !First .15 and Second.15 contain segment angles of the first and second part of the curve
  !First .25 and Second.25 contain chain (curve) vertices
  !First .5 and Second.5 contain final verices and integrals under both parts of the curve
  !First  part contains 1000 segments; second part contains 900 segments (both exclude the segment between the two parts)
  !Old data ends in 4 instead of 5 and can be used to draw the loop close to the one in the figure

  !Setting curve length lower bound
  A = 1.986d0
  
  !Setting curve direction
  ch = 1
  
  !Setting initial and final values of x, y and z (first 3 - initial; last 3 - final)
  call sf%init(6)
  sf%d(1) = 0.0d0
  sf%d(2) = 1.0d0
  sf%d(3) = 0.0d0
  sf%d(4) = 0.0d0
  sf%d(5) = 0.0d0
  sf%d(6) = 0.0d0
  
  !Initializing part of the curve with loop
  !All file names end in letter 5. Rewrite the code and use different number as "fnamenum" in "UpdateCurve" to save different files.
  open(2, file='First .15', action='write')
  do i = 1, 40
    write(2, *) 0.5d0*pi
  end do
  do i = 1, 40
    write(2, *) pi
  end do
  do i = 1, 40
    write(2, *) -0.5d0*pi
  end do
  do i = 1, 880
    write(2, *) 0.0d0
  end do
  close(2)
  !Calculating and saving the cost of the initial curve (performing 0 iterations)
  !We save information about first part of the curve in "First " and about second in "Second"
  call UpdateCurve(verbose-1, A*0.1d0, 1000, ch, 0, 1, sf, 'First ', 1, 10, 5)

  !Initializing part of the curve without loop
  open(2, file='Second.15', action='write')
  do i = 1, 500
    write(2, *) -0.5d0*pi
  end do
  do i = 1, 400
    write(2, *) -pi
  end do
  close(2)
  if (verbose > 0) then
    print *, 'PREPARATION COMPLETE'
  end if

  !Finding the solution. 1 iteration should be enough
  do j = 1, 1
  
    !SECOND PART OF THE CURVE (INTERPOLATION DATA)
    !Setting start to (1;0)
    sf%d(1) = 1.0d0
    sf%d(2) = 0.0d0
    sf%d(3) = 0.0d0
    !Recreating InterpolationData file
    open(2, file='InterpolationData', action='write', status='replace')
    close(2)
    !Solving for 8 different multipliers in the loss function
    mult = 100000000
    do i = 1, 8
      open(1, file='First .5', status='old')
      read(1, *) sf%d(4)
      read(1, *) sf%d(5)
      sf%d(6) = 0.0d0
      close(1)
      mult = mult/10
      if (verbose > 0) then
        print *, 'CURRENT MULTIPLIER', mult
      end if
      !Converging to z = 0 with the current multiplier
      call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 0, sf, 'Second', 0, mult, 5)
      !Reading final point and integral and saving as interpolation data
      open(1, file='Second.5', status='old')
      read(1, *) sf%d(4)
      read(1, *) sf%d(5)
      read(1, *) sf%d(6)
      close(1)
      open(2, file='InterpolationData', position='append')
      write(2, *) mult, sf%d(4), sf%d(5), sf%d(6)
      close(2)
    end do

    !SECOND PART OF THE CURVE (CHOOSING BEST INTERPOLATION)
    !Reading values for the end of the first part (with loop)
    open(1, file='First .5', status='old')
    read(1, *) sf%d(4)
    read(1, *) sf%d(5)
    read(1, *) sf%d(6)
    close(1)
    !Reading interpolation data
    open(1, file='InterpolationData', status='old')
    call mults%init(8)
    call vz%init(8)
    do i = 1, 8
      read(1, *) mults%d(i), curx, cury, vz%d(i)
    end do
    close(1)
    !Finding the best multiplier
    if (-vz%d(1) > sf%d(6)) then
      mult = int(mults%d(1))
    else if (-vz%d(8) < sf%d(6)) then
      mult = int(mults%d(8))
    else
      do i = 1, 8
        mults%d(i) = i
      end do
      mult = int(10.0d0**(8.0d0-quadratic(vz,mults,-sf%d(6))))
    end if
    if (verbose > 0) then
      print *, 'FINAL MULTIPLIER', mult
      print *, 'Setting Z to zero'
    end if
    tmp = sf%d(6)
    !Converging to z = 0
    sf%d(6) = 0.0d0
    call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 0, sf, 'Second', 0, mult, 5)
    !Converging to z in the end of the first part of the curve
    sf%d(6) = tmp
    call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 1, sf, 'Second', 0, mult, 5)
    call mults%deinit()
    call vz%deinit()

    !FIRST PART OF THE CURVE
    !Set start point to (0;1)
    sf%d(1) = 0.0d0
    sf%d(2) = 1.0d0
    sf%d(3) = 0.0d0
    !Reading values in the end of the second part of the curve
    open(1, file='Second.5', status='old')
    read(1, *) sf%d(4)
    read(1, *) sf%d(5)
    read(1, *) sf%d(6)
    close(1)
    !Creating the first part of the curve (with loop)
    if (verbose > 0) then
      print *, 'LOOP IS CREATING'
    end if
    dstmp = dsecnd()
    call UpdateCurve(verbose-1, A*0.1d0, 1000, ch, 500, 2, sf, 'First ', 0, 10, 5)
    if (verbose > 0) then
      print *, 'Time (current iteration)', dsecnd() - dstmp
    end if
  end do

  !Final iteration of the algorithm
  !SECOND PART OF THE CURVE
  !Setting start to (1;0)
  sf%d(1) = 1.0d0
  sf%d(2) = 0.0d0
  sf%d(3) = 0.0d0
  !Reading values for the end of the first part (with loop)
  open(1, file='First .5', status='old')
  read(1, *) sf%d(4)
  read(1, *) sf%d(5)
  read(1, *) sf%d(6)
  close(1)
  !Reading interpolation data
  open(1, file='InterpolationData', status='old')
  call mults%init(8)
  call vz%init(8)
  do i = 1, 8
    read(1, *) mults%d(i), curx, cury, vz%d(i)
  end do
  close(1)
  !Finding the best multiplier
  if (-vz%d(1) > sf%d(6)) then
    mult = int(mults%d(1))
  else if (-vz%d(8) < sf%d(6)) then
    mult = int(mults%d(8))
  else
    do i = 1, 8
      mults%d(i) = i
    end do
    mult = int(10.0d0**(8.0d0-quadratic(vz,mults,-sf%d(6))))
  end if
  if (verbose > 0) then
    print *, 'FINAL MULTIPLIER', mult
    print *, 'Setting Z to zero'
  end if
  tmp = sf%d(6)
  !Converging to z = 0
  sf%d(6) = 0.0d0
  call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 0, sf, 'Second', 0, mult, 5)
  !Converging to z in the end of the first part of the curve
  sf%d(6) = tmp
  call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 1, sf, 'Second', 0, mult, 5)
  call mults%deinit()
  call vz%deinit()

  !FIRST PART OF THE CURVE
  !Setting start to (0;1)
  sf%d(1) = 0.0d0
  sf%d(2) = 1.0d0
  sf%d(3) = 0.0d0
  !Reading values for the end of the second part
  open(1, file='Second.5', status='old')
  read(1, *) sf%d(4)
  read(1, *) sf%d(5)
  read(1, *) sf%d(6)
  close(1)
  if (verbose > 0) then
    print *, 'FINAL ITERATION'
  end if
  !Creating the first part of the curve (with loop)
  dstmp = dsecnd()
  call UpdateCurve(verbose-1, A*0.1d0, 1000, ch, 500, 1, sf, 'First ', 0, 10, 5)
  if (verbose > 0) then
    print *, 'Time (final iteration)', dsecnd() - dstmp
  end if
end

!Improves existing curve, obtained from Loop
subroutine LoopImprove(verbose)
  USE ModVec
  Integer(4), intent(in) :: verbose !3 - print everything, 2 - print only results, 1 - do not print anything from UpdateCurve, 0 - print nothing
  
  Integer(4) ch !curve direction
  Integer(4) mult !multiplier for z in loss function (used for interpolation)
  Integer(4) i !for loop index
  Double precision A !curve length lower bound (should be close to solution)
  Double precision curx, cury, tmp !temporaries
  Type(Vector) sf !start and final values of x, y and z
  Type(Vector) mults !multipliers from 1 to 10^7
  Type(Vector) vz !interpolation data for z

  !Setting curve length lower bound
  A = 1.986d0
  
  !Setting curve direction
  ch = 1
  
  !Setting initial and final values of x, y and z (first 3 - initial; last 3 - final)
  call sf%init(6)
  sf%d(1) = 1.0d0
  sf%d(2) = 0.0d0
  sf%d(3) = 0.0d0
  open(1, file='First .5', status='old')
  read(1, *) sf%d(4)
  read(1, *) sf%d(5)
  read(1, *) sf%d(6)
  close(1)
  
  !Initializing part of the curve without loop
  open(2, file='Second.15', action='write')
  do i = 1, 500
    write(2, *) -0.5d0*pi
  end do
  do i = 1, 400
    write(2, *) -pi
  end do
  close(2)
  if (verbose > 0) then
    print *, 'PREPARATION COMPLETE'
  end if

  !SECOND PART OF THE CURVE (INTERPOLATION DATA)
  !Setting start to (1;0)
  sf%d(1) = 1.0d0
  sf%d(2) = 0.0d0
  sf%d(3) = 0.0d0
  mult = 100000000
  !Recreating InterpolationData file
  open(2, file='InterpolationData', action='write', status='replace')
  close(2)
  !Solving for 8 different multipliers in the loss function
  do i = 1, 8
    open(1, file='First .5', status='old')
    read(1, *) sf%d(4)
    read(1, *) sf%d(5)
    sf%d(6) = 0.0d0
    close(1)
    mult = mult/10
    if (verbose > 0) then
      print *, 'MULTIPLIER', mult
    end if
    !Converging to z = 0 with the current multiplier
    call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 0, sf, 'Second', 0, mult, 5)
    !Reading final point and integral and saving as interpolation data
    open(1, file='Second.5', status='old')
    read(1, *) sf%d(4)
    read(1, *) sf%d(5)
    read(1, *) sf%d(6)
    close(1)
    open(2, file='InterpolationData', position='append')
    write(2, *) mult, sf%d(4), sf%d(5), sf%d(6)
    close(2)
  end do

  !SECOND PART OF THE CURVE (CHOOSING BEST INTERPOLATION)
  !Reading values for the end of the first part (with loop)
  open(1, file='First .5', status='old')
  read(1, *) sf%d(4)
  read(1, *) sf%d(5)
  read(1, *) sf%d(6)
  close(1)
  !Reading interpolation data
  open(1, file='InterpolationData', status='old')
  call mults%init(8)
  call vz%init(8)
  do i = 1, 8
    read(1, *) mults%d(i), curx, cury, vz%d(i)
  end do
  close(1)
  !Finding the best multiplier
  if (-vz%d(1) > sf%d(6)) then
    mult = int(mults%d(1))
  else if (-vz%d(8) < sf%d(6)) then
    mult = int(mults%d(8))
  else
    do i = 1, 8
      mults%d(i) = i
    end do
    mult = int(10.0d0**(8.0d0-quadratic(vz,mults,-sf%d(6))))
  end if
  if (verbose > 0) then
    print *, 'FINAL MULTIPLIER', mult
    print *, 'Setting Z to zero'
  end if
  tmp = sf%d(6)
  !Converging to z = 0
  sf%d(6) = 0.0d0
  call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 0, sf, 'Second', 0, mult, 5)
  !Converging to z in the end of the first part of the curve
  sf%d(6) = tmp
  call UpdateCurve(verbose-1, A*0.9d0, 900, -ch, 500, 1, sf, 'Second', 0, mult, 5)
  call mults%deinit()
  call vz%deinit()

  !FIRST PART OF THE CURVE
  !Set start point to (0;1)
  sf%d(1) = 0.0d0
  sf%d(2) = 1.0d0
  sf%d(3) = 0.0d0
  !Reading values for the end of the second part
  open(1, file='Second.5', status='old')
  read(1, *) sf%d(4)
  read(1, *) sf%d(5)
  read(1, *) sf%d(6)
  close(1)
  !Creating the first part of the curve (with loop)
  if (verbose > 0) then
    print *, 'FINAL LOOP IS CREATING'
  end if
  call UpdateCurve(verbose-1, A*0.1d0, 1000, ch, 500, 2, sf, 'First ', 1, 10, 5)
  if (verbose > 0) then
    print *, 'DONE'
  end if
end

!Updating the corresponding part of the curve using coordinate descent
subroutine UpdateCurve(verbose, A, N, cht, maxit, ddz, sf, fname, squares, multin, fnamenum)
  USE ModVec
  Integer(4), intent(in) :: verbose !2 - print everything, 1 - only results, 0 - nothing
  Double precision, intent(in) :: A !Lower bound for the curve length
  Integer(4), intent(in) :: N !Number of chain segments
  Integer(4), intent(in) :: cht !Size and direction of integral under the curve (includes curve direction)
  Integer(4), intent(in) :: maxit !Limit on the number of iterations
  Integer(4), intent(in) :: ddz !Loss function number (type)
  Type(Vector), intent(in) :: sf !start and final values of x, y and z
  Character(6), intent(in) :: fname !Name of the data file
  Integer(4), intent(in) :: squares !2 - allow to generate squares on curve of arbbitrary size, 1 - allow to generate squares of size "dt", 0 - forbid squares
  Integer(4), intent(in), optional :: multin !multiplier for integral z (only as input, used name is "mult")
  Integer(4), intent(in), optional :: fnamenum !Number, added to file name
  
  Character(1) filenum !file number as character
  Double precision dt !minimum rotation angle
  Double precision tmp1, tmp2 !temporaries
  Double precision mult !multiplier for integral z
  Double precision itnum !current number of iterations
  Double precision dz !Possible increase of z from creating dt by dt square
  Type(Vector) u !Segment angles
  Type(Vector) uu !Temporary segment angles
  Type(Vector) ux, uy !x and y segment directions, calculated from angles u
  Type(Vector) tv !temporary values of loss function for different curve updates
  Double precision dist !modified distance till the end of curve (including weighted z)
  Type(Vector) x, y !Chain vertices
  Type(Vector) z !Values of the integral under the curve
  Type(Vector) uv, uv2 !Temporary values (usually for angles, but rarely for z)
  Type(Vector) vx, vy, vz !Interpolation data
  Integer(4) changes !Total number of changes to the chain
  Integer(4) changed !1 - best choice is to change the rotation of the segment; 0 - better to keep the segment rotation
  Integer(4) i, j, k !indices
  Integer(4) curi !remembers the appropriate index i in for loops
  Integer(4) forcedbsize !Force creating big square in the beginning of the curve
  Integer(4) vxn !Number of interpolation points
  
  !Перемещать из "до петли" в "после петли" и наоборот.
  !В точку перед u, меньшего по модулю (напр. 0.11 переместить рядом к -0.11)
  !mult table:
  !100: <4.5*10^-9
  !1000: 10^-8
  !10000: 10^-6
  !100000: 10^-5
  !1000000: 6*10^-5
  !10000000: 10^-4
  !100000000: 1.3*10^-4
  
  !Setting values of optional inputs
  if (present(fnamenum)) then
    write(filenum,'(i0)') fnamenum
  else
    filenum = '5'
  end if
  mult = 700*cht
  if (present(multin)) then
    mult = multin*cht
  end if
  forcedbsize = 0
  if (squares .eq. 2) then
    forcedbsize = 1
  end if
  
  !Setting angle change size
  dt = 1.0d0/N
  
  !Initializing vector arrays
  call u%init(N)
  call uu%init(N)
  call ux%init(N+1)
  call uy%init(N+1)
  call tv%init(3)
  call uv%init(3)
  call uv2%init(3)
  call x%init(N+2)
  call y%init(N+2)
  call z%init(N+2)

  !Reading angle data
  open(1, file=fname//'.1'//filenum, status='old')
  read(1, *) u%d
  close(1)
  
  !Reading interpolation data
  vxn = 8
  call vx%init(vxn)
  call vy%init(vxn)
  call vz%init(vxn)
  if (ddz .eq. 2) then
    open(1, file='InterpolationData', status='old')
    do i = 1, vxn
      read(1, *) j, vx%d(i), vy%d(i), vz%d(i)
    end do
    close(1)
  end if
  x%d(1) = sf%d(1)
  y%d(1) = sf%d(2)
  z%d(1) = sf%d(3)*abs(mult)
  x%d(N+2) = sf%d(4)
  y%d(N+2) = sf%d(5)
  z%d(N+2) = sf%d(6)*abs(mult)
  do i = 1, N
    ux%d(i) = A * sin(u%d(i))
    uy%d(i) = - A * cos(u%d(i))
  end do
  
  !Update angles to be from -pi to pi
  do i = 1, N
    if (u%d(i) > pi) then
      u%d(i) = u%d(i) - 2*pi
    end if
    if (u%d(i) < -pi) then
      u%d(i) = u%d(i) + 2*pi
    end if
  end do
  
  !Calculate segment directions from angles
  do i = 1, N
    ux%d(i) = A * sin(u%d(i))
    uy%d(i) = - A * cos(u%d(i))
  end do
  
  if (verbose > 0) then
    print *, 'Initial cost', recount(x, y, z, ux, uy, mult, dt, A, tmp1)
  end if
  
  !Output minimum values of x and y
  tmp1 = 0
  tmp2 = 0
  do i = 1, N+1
    if (x%d(i) < tmp1) then 
      tmp1 = x%d(i)
    end if
    if (y%d(i) < tmp2) then 
      tmp2 = y%d(i)
    end if
  end do
  if (verbose > 1) then
    print *, 'min(x), min(y)', tmp1, tmp2
  end if
  
  !Finding dz
  !Setting dz to the gain from "dt by dt" square at (0;1)
  dz = 12*abs(mult)*dt**3
  !Finding best position for the square
  changes = 1
  if (changes > 0) then
    do k = 1, 3
      if (k == 2) then
        tmp1 = 0
        curi = 1
        do i = 1, N/2
          tmp2 = (x%d(i)+dt)*(y%d(i)+dt)
          if (tmp2 > tmp1) then
            tmp1 = tmp2
            curi = i
          end if
        end do
        uu = u%shift(4, curi)
        uu%d(curi) = pi/2.0d0
        uu%d(curi+1) = pi
        uu%d(curi+2) = -pi/2.0d0
        uu%d(curi+3) = 0.0d0
        do i = 1, N
          ux%d(i) = A * sin(uu%d(i))
          uy%d(i) = - A * cos(uu%d(i))
        end do
      else
        do i = 1, N
          ux%d(i) = A * sin(u%d(i))
          uy%d(i) = - A * cos(u%d(i))
        end do
      end if
      tv%d(k) = recount(x, y, z, ux, uy, mult, dt, A, tmp1)
      uv%d(k) = tmp1
      if (k == 2) then
        dz = abs(uv%d(1) - uv%d(2))
        if (verbose > 1) then
          print *, 'New Distance (with square)', recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, 1)
        end if
      end if
    end do
    if (verbose > 1) then
      print *, 'Old Distance', recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, 1)
    end if
  end if
  if (verbose > 1) then
    print *, 'Current Cost', tv%d(3)
  end if
  
  !Creating large square at (0;1)
  if (forcedbsize == 1) then
    tmp1 = 0
    curi = 1
    do i = 1, N/2
      tmp2 = (x%d(i)+dt)*(y%d(i)+dt)
      if (tmp2 > tmp1) then
        tmp1 = tmp2
        curi = i
      end if
    end do
    uu = u%shift(4, curi)
    uu%d(curi) = pi/2.0d0
    uu%d(curi+1) = pi
    uu%d(curi+2) = -pi/2.0d0
    uu%d(curi+3) = 0.0d0
    do i = 1, N
      ux%d(i) = A * sin(uu%d(i))
      uy%d(i) = - A * cos(uu%d(i))
    end do
    u = 1.0d0*uu
  end if
  
  !Iterative curve update (coordinate descent)
  itnum = 0
  do while ((changes > 0) .and. (itnum < maxit))
    itnum = itnum + 1
    changes = 0
    
    !In case we deal with second part of the curve, all angles should be at least -pi/2
    if (cht < 0) then
      do k = 1, N-1
        u%d(k) = min(u%d(k),-pi/2)
      end do
      k = 1
    else
      k = 0
    end if
    !Try rotating each segment by +dt and -dt
    do while (k < N)
      k = k + 1
      uv%d(1) = u%d(k) - dt
      uv%d(3) = u%d(k)
      uv%d(2) = u%d(k) + dt
      if ((mult < 0) .and. (k < N)) then
        if (abs(uv%d(1) - u%d(k+1)) < pi) then
          uv%d(1) = max(uv%d(1), u%d(k+1))
        end if
      end if
      if ((mult < 0) .and. (k > 1)) then
        if (abs(uv%d(2) - u%d(k-1)) < pi) then
          uv%d(2) = min(uv%d(2), u%d(k-1))
        end if
      end if
      do j = 1, 3
        ux%d(k) = A * sin(uv%d(j))
        uy%d(k) = - A * cos(uv%d(j))
        tv%d(j) = recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, ddz, vx, vy, vz)
      end do
      changed = 0
      if ((tv%d(1) <= tv%d(2)) .and. (tv%d(1) < tv%d(3))) then
        u%d(k) = uv%d(1)
        changed = 1
      elseif ((tv%d(2) <= tv%d(1)) .and. (tv%d(2) < tv%d(3))) then
        u%d(k) = uv%d(2)
        changed = 1
      endif
      ux%d(k) = A * sin(u%d(k))
      uy%d(k) = - A * cos(u%d(k))
      if (changed == 1) then
        changes = changes + 1
        if (changes < 10000) then
          k = k - 1
        end if
      end if
    end do
    
    !Simultaneously rotate current and next segment
    k = 0
    do while (k < N-1)
      k = k + 1
      uv%d(1) = u%d(k) - dt
      uv%d(3) = u%d(k)
      uv%d(2) = u%d(k) + dt
      uv2%d(1) = u%d(k+1) + dt
      uv2%d(3) = u%d(k+1)
      uv2%d(2) = u%d(k+1) - dt
      if ((mult < 0) .and. (k > 1)) then
        if (abs(uv2%d(1) - u%d(k-1)) < pi) then
          uv2%d(1) = min(uv2%d(1), u%d(k-1))
        end if
        if (abs(uv%d(2) - u%d(k-1)) < pi) then
          uv%d(2) = min(uv%d(2), u%d(k-1))
        end if
      end if
      if (mult < 0) then
        if (abs(uv%d(1) - uv2%d(1)) < pi) then
          uv%d(1) = max(uv%d(1), uv2%d(1))
        end if
      end if
      if ((mult < 0) .and. (k < N-1)) then
        if (abs(uv2%d(2) - u%d(k+2)) < pi) then
          uv2%d(2) = max(uv2%d(2), u%d(k+2))
        end if
      end if
      do j = 1, 3
        ux%d(k) = A * sin(uv%d(j))
        uy%d(k) = - A * cos(uv%d(j))
        ux%d(k+1) = A * sin(uv2%d(j))
        uy%d(k+1) = - A * cos(uv2%d(j))
        tv%d(j) = recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, ddz, vx, vy, vz)
      end do
      changed = 0
      if ((tv%d(1) <= tv%d(2)) .and. (tv%d(1) < tv%d(3))) then
        u%d(k) = uv%d(1)
        u%d(k+1) = uv2%d(1)
        changed = 1
      elseif ((tv%d(2) <= tv%d(1)) .and. (tv%d(2) < tv%d(3))) then
        u%d(k) = uv%d(2)
        u%d(k+1) = uv2%d(2)
        changed = 1
      endif
      ux%d(k) = A * sin(u%d(k))
      uy%d(k) = - A * cos(u%d(k))
      ux%d(k+1) = A * sin(u%d(k+1))
      uy%d(k+1) = - A * cos(u%d(k+1))
      if (changed == 1) then
        changes = changes + 1
        if (changes < 10000) then
          k = k - 1
        end if
      end if
    end do
    
    !Switch current and next segment (TURNED OFF, NEVER USED)
    k = 0
    do while (k > N)!(k < N-1)
      k = k + 1
      uv%d(1) = u%d(k+1)
      uv%d(2) = u%d(k)
      uv2%d(1) = u%d(k)
      uv2%d(2) = u%d(k+1)
      if (mult < 0) then
        if (abs(uv%d(1) - uv2%d(1)) < pi) then
          uv%d(1) = max(uv%d(1), uv2%d(1))
        end if
      end if
      do j = 1, 2
        ux%d(k) = A * sin(uv%d(j))
        uy%d(k) = - A * cos(uv%d(j))
        ux%d(k+1) = A * sin(uv2%d(j))
        uy%d(k+1) = - A * cos(uv2%d(j))
        tv%d(j) = recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, ddz, vx, vy, vz)
      end do
      changed = 0
      if ((tv%d(1) < tv%d(2)) .or. ((mult > 0) .and. (abs(u%d(k) - u%d(k+1)) < pi) .and. (u%d(k) > u%d(k+1)))) then
        u%d(k) = uv%d(1)
        u%d(k+1) = uv2%d(1)
        changed = 1
      endif
      ux%d(k) = A * sin(u%d(k))
      uy%d(k) = - A * cos(u%d(k))
      ux%d(k+1) = A * sin(u%d(k+1))
      uy%d(k+1) = - A * cos(u%d(k+1))
      if (changed == 1) then
        changes = changes + 1
      end if
    end do
    
    !Rotate segment by 0.1 of the current angle (TURNED OFF, NEVER USED)
    k = 0
    do while (k > N)!(k < N-1)
      k = k + 1
      tmp1 = 0.1d0 * max(dt/10,abs(u%d(k)-u%d(k+1)))
      uv%d(1) = u%d(k) - tmp1
      uv%d(3) = u%d(k)
      uv%d(2) = u%d(k) + tmp1
      if ((mult < 0) .and. (k < N)) then
        if (abs(uv%d(1) - u%d(k+1)) < pi) then
          uv%d(1) = max(uv%d(1), u%d(k+1))
        end if
      end if
      if ((mult < 0) .and. (k > 1)) then
        if (abs(uv%d(2) - u%d(k-1)) < pi) then
          uv%d(2) = min(uv%d(2), u%d(k-1))
        end if
      end if
      do j = 1, 3
        ux%d(k) = A * sin(uv%d(j))
        uy%d(k) = - A * cos(uv%d(j))
        tv%d(j) = recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, ddz, vx, vy, vz)
      end do
      changed = 0
      if ((tv%d(1) <= tv%d(2)) .and. (tv%d(1) < tv%d(3))) then
        u%d(k) = uv%d(1)
        changed = 1
      elseif ((tv%d(2) <= tv%d(1)) .and. (tv%d(2) < tv%d(3))) then
        u%d(k) = uv%d(2)
        changed = 1
      endif
      ux%d(k) = A * sin(u%d(k))
      uy%d(k) = - A * cos(u%d(k))
      if (changed == 1) then
        changes = changes + 1
        if (changes < 10000) then
          k = k - 1
        end if
      end if
    end do
    
    !Try creating 1 by 1 square at best position, taking 4 segments from the end of the curve
    if (squares == 1) then
      do k = 1, 3
        if (k == 2) then
          tmp1 = 0
          curi = 1
          do i = 1, N/2
            tmp2 = (x%d(i)+dt)*(y%d(i)+dt)
            if (tmp2 > tmp1) then
              tmp1 = tmp2
              curi = i
            end if
          end do
          uu = u%shift(4, curi)
          uu%d(curi) = pi/2.0d0
          uu%d(curi+1) = pi
          uu%d(curi+2) = -pi/2.0d0
          uu%d(curi+3) = 0.0d0
          do i = 1, N
            ux%d(i) = A * sin(uu%d(i))
            uy%d(i) = - A * cos(uu%d(i))
          end do
        else
          do i = 1, N
            ux%d(i) = A * sin(u%d(i))
            uy%d(i) = - A * cos(u%d(i))
          end do
        end if
        tv%d(k) = recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, ddz, vx, vy, vz)
        uv%d(k) = tmp1
        if ((k == 2) .and. (tv%d(2) < tv%d(1))) then
          u = 1.0d0*uu
          changes = changes + 1
        end if
      end do
    end if
    
    !Update dz in case better square is found
    dz = max(dz, uv%d(1) - uv%d(2))
    !Loss function (modified distance)
    dist = 4*A*dt*(z%d(N+1)-z%d(N+2))/dz + sqrt((x%d(N+1) - x%d(N+2))**2 + (y%d(N+1) - y%d(N+2))**2)
    if (verbose > 1) then
      print *, 'Itnum, dz, dist', itnum, dz, dist
      print *, 'Changes, cost  ', changes, tv%d(3)
    end if
  end do
  
  if (verbose > 0) then
    print *, 'Final cost', recount(x, y, z, ux, uy, mult, dt, A, tmp1)
  end if
  
  !Save chain vertices
  open(2, file=fname//'.2'//filenum)
  do i = 1, N+1
    write(2, *) x%d(i), y%d(i)
  end do
  write(2, *) x%d(N+2), y%d(N+2)
  close(2)
  
  !Save final point and integral value
  open(2, file=fname//'.'//filenum)
  write(2, *) x%d(N+1)
  write(2, *) y%d(N+1)
  write(2, *) -z%d(N+1)/abs(mult)
  close(2)
  if (verbose > 0) then
    print *, 'Final data', x%d(N+1), y%d(N+1), z%d(N+1)/mult
  end if
  
  !Save segment directions (angles)
  open(1, file=fname//'.1'//filenum)
  do i = 1, N
    write(1, *) u%d(i)
  end do
  close(1)
  
  !Output total curve length
  if (verbose > 0) then
    if (A > 1.0d0) then
      print *, 'Total length', recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, max(1, ddz), vx, vy, vz)+A/0.9d0
    else
      print *, 'Total length', recount(x, y, z, ux, uy, mult, dt, A/dz, tmp1, max(1, ddz), vx, vy, vz)+A*10.0d0
    end if
  end if
end

!Calculating quadratic interpolation of the function "vy" in the point "point"
function quadratic(vx, vy, point) Result(res)
  USE ModVec
  Type(Vector), intent(in) :: vx !x coordinates
  Type(Vector), intent(in) :: vy !y coordinates (function values)
  Type(Vector) :: vx1, vy1 !coordinates for quadratic interpolation
  Double precision, intent(in) :: point !current point
  Double precision :: res !result
  
  Double precision p1, p2 !temporary points
  Integer(4) i, j !point indices
  
  j = vx%n-2
  
  !Find the segment j, containing the point
  do i = 1, vx%n-2
    p1 = vx%d(i)
    p2 = vx%d(i+1)
    if ((point <= p1) .and. (point >= p2)) then
      j = i
      exit
    end if
    if ((point >= p1) .and. (point <= p2)) then
      j = i
      exit
    end if
  end do
  
  !Calculating quadratic interpolation as a special case of lagrange interpolation
  call vx1%init(3)
  call vy1%init(3)
  do i = 1, 3
    vx1%d(i) = vx%d(j+i-1)
    vy1%d(i) = vy%d(j+i-1)
  end do
  res = lagrange(vx1, vy1, point)
  
  !If the result is smaller or larger than in the ends, return the result from one of the end points
  if (res > max(vy%d(j), vy%d(j+1))) then
    res = max(vy%d(j), vy%d(j+1))
  end if
  if (res < min(vy%d(j), vy%d(j+1))) then
    res = min(vy%d(j), vy%d(j+1))
  end if
end

!Calculating lagrange interpolation polynomial
function lagrange(vx, vy, point) Result(res)
  USE ModVec
  Type(Vector), intent(in) :: vx, vy !interpolation points
  Double precision, intent(in) :: point !current point
  Double precision :: res !result
  
  Double precision pres !result for the products (basis polynomials)
  Integer(4) i, j !loop indices
  
  res = 0.0d0
  do i = 1, vx%n
    pres = vy%d(i)
    do j = 1, vx%n
      if (j .ne. i) then
        pres = pres * (point - vx%d(j)) / (vx%d(i) - vx%d(j))
      end if
    end do
    res = res + pres
  end do
end

!Calculates the appropriate loss function under the chain
function recount(x, y, z, u1, u2, mult, dt, A, zn, ddz, vx, vy, vz) Result(res)
  USE ModVec
  Type(Vector) :: x, y, z !chain vertices and the integral under the curve
  Type(Vector) :: u1, u2 !directions of chain segments
  Double precision, intent(in) :: mult !multiplier for integral z
  Double precision, intent(in) :: dt !segment length
  Double precision, intent(in) :: A !multiplier in the cost function
  Double precision, intent(out), optional :: zn !total integral (last value of z)
  Integer(4), intent(in), optional :: ddz !loss function number (type)
  Type(Vector), intent(in), optional :: vx, vy, vz !interpolation points and functions
  Double precision :: res !result
  
  Double precision lagrx, lagry !interpolation results
  Integer(4) N !number of chain segments, excluding last
  Double precision tmp1, tmp2 !temporaries
  Integer(4) i !loop index
  
  N = x%n-2
  
  !Calculating position x, y and integral z under the chain (exact formula)
  do i = 2, N+1
    x%d(i) = x%d(i-1) + u1%d(i-1) * dt
    y%d(i) = y%d(i-1) + u2%d(i-1) * dt
    tmp1 = 8 * (x%d(i)**3 - x%d(i-1)**3) * u2%d(i-1)*(y%d(i-1) - x%d(i-1)*u2%d(i-1)*dt)
    tmp1 = tmp1 + 6 * (x%d(i) + x%d(i-1)) * u1%d(i-1) * (y%d(i-1) - x%d(i-1)*u2%d(i-1)*dt)**2
    tmp1 = tmp1 + 3 * (x%d(i)**4 - x%d(i-1)**4) * u2%d(i-1)**2 * dt
    z%d(i) = z%d(i-1) + tmp1*(mult*dt)
  end do
  
  !Adding the integral under the last segment to the final point outside the chain
  u1%d(N+1) = (x%d(N+2)-x%d(N+1))/dt
  u2%d(N+1) = (y%d(N+2)-y%d(N+1))/dt
  tmp1 = 8 * (x%d(N+2)**3 - x%d(N+1)**3) * u2%d(N+1)*(y%d(N+1) - x%d(N+1)*u2%d(N+1)*dt)
  tmp1 = tmp1 + 6 * (x%d(N+2) + x%d(N+1)) * u1%d(N+1) * (y%d(N+1) - x%d(N+1)*u2%d(N+1)*dt)**2
  tmp1 = tmp1 + 3 * (x%d(N+2)**4 - x%d(N+1)**4) * u2%d(N+1)**2 * dt
  tmp2 = z%d(N+1) + tmp1*(mult*dt)
  if (present(zn)) then
    zn = z%d(N+1)
  end if
  if (tmp2 < z%d(N+2)) then
    tmp2 = z%d(N+2)
  end if
  
  !Calculating loss function
  if (present(ddz)) then
    if (ddz == 1) then
      !Case 1: distance to the final point + weighted integral under the curve
      res = 4*A*dt*(tmp2-z%d(N+2)) + sqrt((x%d(N+1) - x%d(N+2))**2 + (y%d(N+1)-y%d(N+2))**2)
    else if (ddz == 2) then
      !Case 2: distance to the second part of the curve using interpolation for it's end depending on the integral under it
      !Ошибка в шестом знаке. z%d(N+2) или N+1? Был n+2.
      !Посмотреть, насколько быстро сойдется.
      lagrx = quadratic(vz, vx, z%d(N+2)/abs(mult))
      lagry = quadratic(vz, vy, z%d(N+2)/abs(mult))
      if (z%d(N+2)/abs(mult) > vz%d(1)) then
        lagrx = vx%d(1)
        lagry = vy%d(1)
      end if
      if (z%d(N+2)/abs(mult) < vz%d(vz%n)) then
        lagrx = vx%d(vz%n)
        lagry = vy%d(vz%n)
      end if
      u1%d(N+1) = (lagrx-x%d(N+1))/dt
      u2%d(N+1) = (lagry-y%d(N+1))/dt
      tmp1 = 8 * (lagrx**3 - x%d(N+1)**3) * u2%d(N+1)*(y%d(N+1) - x%d(N+1)*u2%d(N+1)*dt)
      tmp1 = tmp1 + 6 * (lagrx + x%d(N+1)) * u1%d(N+1) * (y%d(N+1) - x%d(N+1)*u2%d(N+1)*dt)**2
      tmp1 = tmp1 + 3 * (lagrx**4 - x%d(N+1)**4) * u2%d(N+1)**2 * dt
      tmp2 = z%d(N+1) + tmp1*(mult*dt)
      res = sqrt((x%d(N+1) - lagrx)**2 + (y%d(N+1) - lagry)**2)
      lagrx = quadratic(vz, vx, tmp2/abs(mult))
      lagry = quadratic(vz, vy, tmp2/abs(mult))
      if (tmp2/abs(mult) > vz%d(1)) then
        lagrx = vx%d(1)
        lagry = vy%d(1)
      end if
      if (tmp2/abs(mult) < vz%d(vz%n)) then
        lagrx = vx%d(vz%n)
        lagry = vy%d(vz%n)
      end if
      res = sqrt((x%d(N+1) - lagrx)**2 + (y%d(N+1) - lagry)**2)
    else
      res = (x%d(N+1) - x%d(N+2))**2 + (y%d(N+1)-y%d(N+2))**2 + (tmp2-z%d(N+2))**2
    end if
  else
    !Case 3: simple quadratic loss function
    res = (x%d(N+1) - x%d(N+2))**2 + (y%d(N+1)-y%d(N+2))**2 + (tmp2-z%d(N+2))**2
  end if
end
