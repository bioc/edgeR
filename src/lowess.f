      subroutine lowess(x, y, n, f, nsteps, delta, ys, rw, res)
      integer n
      integer nsteps
      double precision x(n), y(n), f, delta, ys(n), rw(n)
      double precision res(n)
      double precision sc
      integer nright, min0, max0, i, j, ifix
      integer iter, last, m1, m2, ns, nleft
      double precision abs, cut, cmad, r, d1, d2
      double precision c1, c9, alpha, denom, float
      logical ok
      if (n .ge. 2) goto 1
         ys(1) = y(1)
         return
c at least two, at most n points
   1  ns = max0(min0(int(f*dble(n) + 1.0d-7), n), 2)
      iter = 1
         goto  3
   2     iter = iter+1
   3     if (iter .gt. nsteps+1) goto  22
c robustness iterations
         nleft = 1
         nright = ns
c index of prev estimated point
         last = 0
c index of current point
         i = 1
   4        if (nright .ge. n) goto  5
c move nleft, nright to right if radius decreases
               d1 = x(i)-x(nleft)
c if d1<=d2 with x(nright+1)==x(nright), lowest fixes
               d2 = x(nright+1)-x(i)
               if (d1 .le. d2) goto  5
c radius will not decrease by move right
               nleft = nleft+1
               nright = nright+1
               goto  4
c fitted value at x(i)
   5        call lowest(x, y, n, x(i), ys(i), nleft, nright, res, iter
     +     .gt. 1, rw, ok)
            if (.not. ok) ys(i) = y(i)
c all weights zero - copy over value (all rw==0)
            if (last .ge. i-1) goto 9
               denom = x(i)-x(last)
c skipped points -- interpolate
c non-zero - proof?
               j = last+1
                  goto  7
   6              j = j+1
   7              if (j .ge. i) goto  8
                  alpha = (x(j)-x(last))/denom
                  ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
                  goto  6
   8           continue
c last point actually estimated
   9        last = i
c x coord of close points
            cut = x(last)+delta
            i = last+1
               goto  11
  10           i = i+1
  11           if (i .gt. n) goto  13
c find close points
               if (x(i) .gt. cut) goto  13
c i one beyond last pt within cut
               if (x(i) .ne. x(last)) goto 12
                  ys(i) = ys(last)
c exact match in x
                  last = i
  12           continue
               goto  10
c back 1 point so interpolation within delta, but always go forward
  13        i = max0(last+1, i-1)
  14        if (last .lt. n) goto  4
c residuals
         do  15 i = 1, n
            res(i) = y(i)-ys(i)
  15        continue
c overall scale estimate (from R clowess.c)
         sc = 0.0d0
         do  151 i = 1, n
            sc = sc + abs(res(i))
  151    continue
         sc = sc/dble(n)
         if (iter .gt. nsteps) goto  22
c compute robustness weights except last time
         do  16 i = 1, n
            rw(i) = abs(res(i))
  16        continue
         m1 = n/2+1
         m2 = n-m1+1
c partial selection of the two median order statistics (matches
c the rPsort partial sort in R clowess.c; O(n) vs the old O(n**2) sort)
         call qselect(rw, n, m1)
         if (m2 .ne. m1) call qselect(rw, n, m2)
c 6 median abs resid
         cmad = 3.0*(rw(m1)+rw(m2))
c effectively zero (from R clowess.c)
         if (cmad .lt. 1.0d-7*sc) goto 22
         c9 = 0.999d0*cmad
         c1 = 0.001d0*cmad
         do  21 i = 1, n
            r = abs(res(i))
            if (r .gt. c1) goto 17
               rw(i) = 1.
c near 0, avoid underflow
               goto  20
  17           if (r .le. c9) goto 18
                  rw(i) = 0.
c near 1, avoid underflow
                  goto  19
  18              rw(i) = (1.0-(r/cmad)**2)**2
  19        continue
  20        continue
  21        continue
         goto  2
  22  return
      end
      subroutine lowest(x, y, n, xs, ys, nleft, nright, w, userw
     +, rw, ok)
      integer n
      integer nleft, nright
      double precision x(n), y(n), xs, ys, w(n), rw(n)
      logical userw, ok
      integer nrt, j
      double precision abs, a, b, c, h, r
      double precision h1, sqrt, h9, amax1, range
      range = x(n)-x(1)
      h = max(xs-x(nleft), x(nright)-xs)
      h9 = 0.999d0*h
      h1 = 0.001d0*h
c sum of weights
      a = 0.0
      j = nleft
         goto  2
   1     j = j+1
   2     if (j .gt. n) goto  7
c compute weights (pick up all ties on right)
         w(j) = 0.
         r = abs(x(j)-xs)
         if (r .gt. h9) goto 5
            if (r .le. h1) goto 3
               w(j) = (1.0-(r/h)**3)**3
c small enough for non-zero weight
               goto  4
   3           w(j) = 1.
   4        if (userw) w(j) = rw(j)*w(j)
            a = a+w(j)
            goto  6
   5        if (x(j) .gt. xs) goto  7
c get out at first zero wt on right
   6     continue
         goto  1
c rightmost pt (may be greater than nright because of ties)
   7  nrt = j-1
      if (a .gt. 0.0) goto 8
         ok = .false.
         goto  16
   8     ok = .true.
c weighted least squares
         do  9 j = nleft, nrt
c make sum of w(j) == 1
            w(j) = w(j)/a
   9        continue
         if (h .le. 0.) goto 14
            a = 0.0
c use linear fit
            do  10 j = nleft, nrt
c weighted center of x values
               a = a+w(j)*x(j)
  10           continue
            b = xs-a
            c = 0.0
            do  11 j = nleft, nrt
               c = c+w(j)*(x(j)-a)**2
  11           continue
            if (sqrt(c) .le. 0.001d0*range) goto 13
               b = b/c
c points are spread out enough to compute slope
               do  12 j = nleft, nrt
                  w(j) = w(j)*(b*(x(j)-a)+1.0)
  12              continue
  13        continue
  14     ys = 0.0
         do  15 j = nleft, nrt
            ys = ys+w(j)*y(j)
  15        continue
  16  return
      end
      subroutine qselect(a, n, k)
c in-place quickselect: put the k-th smallest (1-indexed) at a(k), with
c a(lo:k-1) <= a(k) <= a(k+1:hi). average O(n). deterministic median-of-
c three pivot (no RNG) guards the O(n**2) worst case on near-sorted input.
c only the k-th order statistic is needed (median of |residuals|), matching
c the partial select used by R clowess.c (rPsort).
      integer n, k, lo, hi, i, j, mid
      double precision a(n), piv, t
      lo = 1
      hi = n
    1 if (lo .ge. hi) return
c median-of-three: order a(lo) <= a(mid) <= a(hi); pivot = a(mid)
      mid = (lo + hi) / 2
      if (a(mid) .lt. a(lo)) then
         t = a(lo)
         a(lo) = a(mid)
         a(mid) = t
      endif
      if (a(hi) .lt. a(lo)) then
         t = a(lo)
         a(lo) = a(hi)
         a(hi) = t
      endif
      if (a(hi) .lt. a(mid)) then
         t = a(mid)
         a(mid) = a(hi)
         a(hi) = t
      endif
      piv = a(mid)
c Hoare partition of a(lo:hi) about piv
      i = lo - 1
      j = hi + 1
    2 continue
    3 i = i + 1
      if (a(i) .lt. piv) goto 3
    4 j = j - 1
      if (a(j) .gt. piv) goto 4
      if (i .lt. j) then
         t = a(i)
         a(i) = a(j)
         a(j) = t
         goto 2
      endif
c a(lo:j) <= a(j+1:hi); keep the side that holds index k
      if (k .le. j) then
         hi = j
      else
         lo = j + 1
      endif
      goto 1
      end
