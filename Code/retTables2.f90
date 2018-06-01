subroutine bisection(a,n,x, ifind)
  integer :: n, i, ifind
  real    :: a(n), x
  integer :: i1, i2, imid
  i1=1
  i2=n
  if(x>=a(n)) then
     ifind=n
     return
  endif
  if(x<=a(1)) then
     ifind=1
     return
  endif
  do while (i2-i1>1)
     imid=(i1+i2)/2
     if(a(imid)>x) then
        i2=imid
     else 
        if (a(imid)<x) then
           i1=imid
        else
           ifind=imid
           return
        endif
     endif
  end do
  
  ifind=i1
  
end subroutine bisection


subroutine bisectiond0(a,d0, ad0, n,x, ifind)
  integer :: n, i, ifind
  real    :: a(n), d0(n), ad0, x
  integer :: i1, i2, imid
  i1=1
  i2=n
  
  if(x>=a(n)+ad0*d0(n)) then
     ifind=n
     return
  endif
  if(x<=a(1)+ad0*d0(1)) then
     ifind=1
     return
  endif
  do while (i2-i1>1)
     imid=(i1+i2)/2
     if(a(imid)+ad0*d0(imid)>x) then
        i2=imid
     else 
        if (a(imid)+ad0*d0(imid)<x) then
           i1=imid
        else
           ifind=imid
           return
        endif
     endif
  end do
  
  ifind=i1
  
end subroutine bisectiond0
