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
   5        call lowest(x, y, n, x(i), ys(i), nleft, nright, iter
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
      subroutine lowest(x, y, n, xs, ys, nleft, nright, userw, rw, ok)
c local linear fit at xs with tricube weights (times rw if userw).
c one sweep accumulates sum w, w*d, w*d**2, w*y and w*d*y with
c d = x(j)-xs, and the weighted least-squares line is evaluated at xs in
c closed form. this equals the former five-sweep version (weights,
c normalise, mean, variance, weighted sum) up to rounding (~1e-13).
      integer n
      integer nleft, nright
      double precision x(n), y(n), xs, ys, rw(n)
      logical userw, ok
      integer j
      double precision r, h, h1, h9, hinv, q, t, wj, d, range
      double precision s0, s1, s2, t0, t1, md, my, c, b
      range = x(n)-x(1)
      h = max(xs-x(nleft), x(nright)-xs)
      h9 = 0.999d0*h
      h1 = 0.001d0*h
      hinv = 0.0d0
      if (h .gt. 0.0d0) hinv = 1.0d0/h
      s0 = 0.0d0
      s1 = 0.0d0
      s2 = 0.0d0
      t0 = 0.0d0
      t1 = 0.0d0
c weights from nleft on (pick up all ties on right)
      do 6 j = nleft, n
         d = x(j)-xs
         r = abs(d)
         if (r .gt. h9) then
c get out at first zero wt on right
            if (x(j) .gt. xs) goto 7
            goto 6
         endif
         if (r .le. h1) then
c near 0, avoid underflow
            wj = 1.0d0
         else
            q = r*hinv
            t = 1.0d0-q*q*q
            wj = t*t*t
         endif
         if (userw) wj = rw(j)*wj
         s0 = s0+wj
         s1 = s1+wj*d
         s2 = s2+wj*d*d
         t0 = t0+wj*y(j)
         t1 = t1+wj*d*y(j)
   6  continue
c all weights zero
   7  if (s0 .gt. 0.0d0) goto 8
         ok = .false.
         return
   8  ok = .true.
c weighted mean of y; used as is unless the x values are spread out
      my = t0/s0
      ys = my
      if (h .le. 0.0d0) return
      md = s1/s0
      c = s2/s0-md*md
c c <= 0 can only come from rounding when the x values coincide
      if (c .le. 0.0d0) return
      if (sqrt(c) .le. 0.001d0*range) return
c linear fit: slope cov(d,y)/c, evaluated at d = 0
      b = -md/c
      ys = my+b*(t1/s0-md*my)
      return
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
