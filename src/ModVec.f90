Module ModVec
  DOUBLE PRECISION, parameter :: eps = 1.0d-15
  DOUBLE PRECISION, parameter :: pi = acos(-1.0d0)
  Type Vector
    Integer(4) n
    DOUBLE PRECISION,Allocatable :: d(:)
    Contains
      Procedure :: init => Vector_constructor
      Procedure :: deinit => Vector_destructor
      Procedure :: set => Vector_set
      Procedure :: norm => vec_norm
      Procedure :: cnorm => vec_cnorm
      Procedure :: random => Vec_random
      Procedure :: badrandom => Vec_badrandom
      Procedure :: pertapp => Vec_pertapp
      Procedure :: pert => Vec_pert
      Procedure :: swap => Vec_swap
      Procedure :: normalize => Vec_normalize
      Procedure :: maxelement => Vec_maxelement
      Procedure :: house => Vec_house
      Procedure :: intarray => Vec_toint
      Procedure :: reverse => vec_reverse
      Procedure :: shift => vec_shift
      Procedure :: copy => vec_copy
      Procedure :: update1 => vec_update1
      Procedure :: perminv => vec_perminv
  End type
  
  interface operator (.dot.)
    module procedure vec_dotmul
  end interface
  
  interface operator (.dd.)
    module procedure vec_dotdiv
  end interface
  
  interface assignment(=)
    module procedure IntArray_transform, VecArray_transform
  end interface
  
  interface operator(+)
    module procedure vec_sum
  end interface
  
  interface operator(-)
    module procedure vec_sub
  end interface
  
  interface operator(*)
    module procedure mul_vec, mul_tnum, mul_tnumr, mul_numt, mul_numtr
  end interface
  
  interface operator(/)
    module procedure div_num, div_numr
  end interface
  
  Contains
    
    subroutine vec_pertapp(this, pert, c)
      Class(Vector) :: this
      Type(Vector), intent(in) :: pert
      Integer(4), intent(in) :: c
      Type(Vector) :: res
      Integer(4) i
      call res%init(this%n)
      if (c == 1) then
        do i = 1, pert%n
          res%d(i) = this%d(floor(pert%d(i)+0.5d0))
        end do
      else
        do i = 1, pert%n
          res%d(floor(pert%d(i)+0.5d0)) = this%d(i)
        end do
      end if
      call this%copy(res)
    end
    
    function vec_perminv(this) Result(res)
      Class(Vector) :: this
      Type(Vector) :: res
      Integer(4) i
      res%n = this%n
      Allocate(res%d(this%n))
      do i = 1, this%n
        res%d(floor(this%d(i)+0.5d0)) = i
      end do
    end
    
    subroutine vec_swap(this, a, b)
      Class(Vector) :: this
      Integer(4), intent(in) :: a, b
      DOUBLE PRECISION tmp
      if ((a > this%n) .or. (b > this%n)) then
        print *, "error in swap_vec"
        !call backtrace()
        return
      end if
      tmp = this%d(a)
      this%d(a) = this%d(b)
      this%d(b) = tmp
    end
    
    function vec_shift(this, sh, start) Result(res)
      Class(Vector) :: this
      Integer(4), intent(in) :: sh
      Integer(4), optional :: start
      Type(Vector) :: res
      Integer(4) i
      if (.not. present(start)) start = 1
      call res%init(this%n)
      if (sh >= 0) then
        call dcopy(start-1, this%d, 1, res%d, 1)
        do i = this%n, start+sh, -1
          res%d(i) = this%d(i-sh)
        end do
      else
        do i = this%n+sh-start+2, this%n
          res%d(i) = this%d(i)
        end do
        do i = 1, this%n+sh-start+1
          res%d(i) = this%d(i-sh)
        end do
      end if
    end
    
    function vec_reverse(this) Result(res)
      Class(Vector) :: this
      Type(Vector) :: res
      Integer(4) i
      res%n = this%n
      Allocate(res%d(res%n))
      do i = 1, res%n
        res%d(i) = this%d(this%n-i+1)
      end do
    end
    
    subroutine IntArray_transform(this, array)
      Integer(4), dimension(:), intent(in) :: array
      Class(Vector), intent(out) :: this
      Integer(4) i
      this%n = size(array)
      Allocate(this%d(this%n))
      do i = 1, this%n
        this%d(i) = array(i)
      end do
    end
    
    subroutine VecArray_transform(this, array)
      DOUBLE PRECISION, dimension(:), intent(in) :: array
      Class(Vector), intent(out) :: this
      this%n = size(array)
      this%d = array
    end
    
    subroutine vec_copy(this, v)
      Class(Vector) :: this
      Type(Vector), intent(in) :: v
      if (.not. allocated(this%d)) allocate(this%d(v%n))
      if (this%n < v%n) then
        Deallocate(this%d)
        Allocate(this%d(v%n))
      end if
      this%n = v%n
      call dcopy(v%n, v%d, 1, this%d, 1)
    end
    
    subroutine vec_update1(this, alpha, x)
      Class(Vector) :: this
      DOUBLE PRECISION, intent(in) :: alpha
      Type(Vector), intent(in) :: x
      call daxpy(this%n, alpha, x%d, 1, this%d, 1)
    end
    
    function Vec_toint(this) Result(array)
      Class(Vector), intent(in) :: this
      Integer(4) :: array(this%n)
      Integer(4) i
      do i = 1, this%n
        array(i) = floor(this%d(i) + 0.5)
      end do
    end
  
    subroutine Vec_normalize(this)
      Class(Vector) :: this
      call dscal(this%n, 1.0d0/this%norm(), this%d, 1)
    end
    
    subroutine Vec_maxelement(this, k)
      Class(Vector) :: this
      Integer(4) k, idamax
      k = idamax(this%n,this%d,1)
    end
    
    subroutine Vec_house(this, v, beta)
      Class(Vector) :: this
      Type(Vector), intent(out) :: v
      DOUBLE PRECISION, intent(out) :: beta
      DOUBLE PRECISION sigma
      sigma = this%norm()
      sigma = sigma**2 - this%d(1)**2
      call v%copy(this)
      v%d(1) = 1.0d0
      if ((sigma == 0) .and. (this%d(1) >= 0)) then
        beta = 0.0d0
      else if ((sigma == 0) .and. (this%d(1) < 0)) then
        beta = -2.0d0
      else
        v%d(1) = this%d(1) - this%norm()
        beta = 2.0d0*v%d(1)**2/(sigma + v%d(1)**2)
        v = v/v%d(1)
      end if
    end
    
    subroutine Vec_pert(this, n)
      Class(Vector) :: this
      Integer(4), intent(in) :: n
      Integer(4) i
      this%n = n
      Allocate(this%d(n))
      do i = 1, n
        this%d(i) = i
      end do
    end
    
    subroutine Vec_random(this, n, nz)
      Class(Vector) :: this
      Integer(4) n, i
      Integer(4), Optional, intent(in) :: nz
      DOUBLE PRECISION r1, r2, r
      this%n = n
      Allocate(this%d(n))
      do i = 1, n
        r = 1
        do while (r >= 1)
          call random_number(r1)
          call random_number(r2)
          r1 = r1 * 2 - 1
          r2 = r2 * 2 - 1
          r = r1 * r1 + r2 * r2
        end do
        r2 = r1 * sqrt((-2) * log(r) / r)
        this%d(i) = r2
      end do
      if (.not. (present(nz))) then
        call this%normalize()
      end if
    end
    
    subroutine Vec_badrandom(this, n)
      Class(Vector) :: this
      Integer(4) n, i
      DOUBLE PRECISION r
      this%n = n
      Allocate(this%d(n))
      do i = 1, n
        call random_number(r)
        this%d(i) = 2*r-1
      end do
    end
  
    subroutine Vector_constructor(this, n)
      Class(Vector) :: this
      Integer(4) n
      this%n = n
      Allocate(this%d(n))
      this%d = 0
    end
    
    subroutine Vector_set(this, d)
      Class(Vector) :: this
      DOUBLE PRECISION :: d(:)
      this%d = d
      this%n = size(d)
    end
    
    subroutine Vector_destructor(this)
      Class(Vector) :: this
      this%n = 0
      Deallocate(this%d)
    end
    
    function mul_vec(this, v2) Result(res)
      Type(Vector), intent(in) :: this
      DOUBLE PRECISION res, ddot
      Type(Vector), intent(in) :: v2
      if (this%n == v2%n) then
        res = ddot(this%n, this%d, 1, v2%d, 1)
      else
        print *, "error mult_vec"
      endif
    end
      
    function mul_tnum(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    function mul_tnumr(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    function mul_numt(num, this) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    function mul_numtr(num, this) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    function div_num(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    function div_numr(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    function vec_norm(this) Result(res)
      Class(Vector) :: this
      DOUBLE PRECISION res, dnrm2
      res = dnrm2(this%n, this%d, 1)
    end
    
    function vec_cnorm(this) Result(res)
      Class(Vector) :: this
      DOUBLE PRECISION res, idamax
      res = idamax(this%n, this%d, 1)
    end
    
    function vec_sum(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      if (v1%n == v2%n) then
        res%d = v1%d + v2%d
      else
        print *, "error sum_vec", v1%n, v2%n
      endif
    end
    
    function vec_sub(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      if (v1%n == v2%n) then
        res%d = v1%d - v2%d
      else
        print *, "error sub_vec", v1%n, v2%n
      endif
    end
    
    function vec_dotmul(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      if (v1%n == v2%n) then
        res%d = v1%d * v2%d
      else
        print *, "error dotmul_vec", v1%n, v2%n
      endif
    end
    
    function vec_dotdiv(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      if (v1%n == v2%n) then
        res%d = v1%d / v2%d
      else
        print *, "error dotdiv_vec", v1%n, v2%n
      endif
    end
    
    function evec(n, m1) Result(res)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: m1
      Type(Vector) :: res
      Integer(4) i
      call res%init(n)
      if (present(m1)) then
        res%d(m1) = 1.0d0
      else
        do i = 1, n
          res%d(i) = 1.0d0
        end do
      end if
    end
    
end
