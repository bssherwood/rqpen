C Output from Public domain Ratfor, version 1.0
      subroutine rqbr(m,nn,m5,n3,n4,a,b,t,toler,ift,x,e,s,wa,wb,nsol,nds
     *ol,sol,dsol,lsol,h,qn,cutoff,ci,tnmat,big,lci1)
      integer i,j,k,kl,kount,kr,l,lsol,m,m1,m2,m3,m4,m5,ift
      integer n,n1,n2,n3,nsol,ndsol,out,s(m),h(nn,nsol)
      integer nn,n4,idxcf
      logical stage,test,init,iend,lup
      logical lci1,lci2,skip
      double precision a1,aux,b1,big,d,dif,pivot,smax,t,t0,t1,tnt
      double precision min,max,toler,zero,half,one,two
      double precision b(m),sol(n3,nsol),a(m,nn),x(nn),wa(m5,n4),wb(m)
      double precision sum,e(m),dsol(m,ndsol)
      double precision qn(nn),cutoff,ci(4,nn),tnmat(4,nn),tnew,told,tn
      parameter( zero = 0.d0)
      parameter( one = 1.d0)
      parameter( two = 2.d0)
      n=nn
      ift = 0
      wa(m+2,nn+1) = one
      if(m5.ne.m+5)then
      ift = 3
      endif
      if(n3.ne.n+3)then
      ift = 4
      endif
      if(n4.ne.n+4)then
      ift = 5
      endif
      if(m.le.zero.or.n.le.zero)then
      ift = 6
      endif
      if(ift.le.two)then
      half = one/two
      iend = .true.
      lci2 = .false.
      lup = .true.
      skip = .false.
      idxcf = 0
      tnew = zero
      tn = zero
      m1 = m+1
      n1 = n+1
      n2 = n+2
      m2 = m+2
      m3 = m+3
      m4 = m+4
      do23010 j = 1,n
      x(j) = zero
23010 continue
23011 continue
      do23012 i = 1,m
      e(i) = zero
23012 continue
23013 continue
      if(t.lt.zero.or.t.gt.one)then
      t0 = one/(two*float(m))
      t1 = one-t0
      t = t0
      iend = .false.
      lci1 = .false.
      endif
23016 continue
      do23019 i = 1,m 
      k = 1
      do23021 j = 1,nn
      if(k.le.nn)then
      if(j.eq.idxcf)then
      skip = .true.
      else
      skip = .false.
      endif
      if(.not.skip)then
      wa(i,k) = a(i,j)
      k = k+1
      endif
      endif
23021 continue
23022 continue
      wa(i,n4) = n+i
      wa(i,n2) = b(i)
      if(idxcf .ne. 0)then
      wa(i,n3) = tnew*a(i,idxcf)
      else
      wa(i,n3) = zero
      endif
      wa(i,n1) = wa(i,n2)-wa(i,n3)
      if(wa(i,n1).lt.zero)then
      do23033 j = 1,n4
      wa(i,j) = -wa(i,j)
23033 continue
23034 continue
      endif
23019 continue
23020 continue
      do23035 j = 1,n 
      wa(m4,j) = j
      wa(m2,j) = zero
      wa(m3,j) = zero
      do23037 i = 1,m 
      aux = sign(one,wa(m4,j))*wa(i,j)
      wa(m2,j) = wa(m2,j)+aux*(one-sign(one,wa(i,n4)))
      wa(m3,j) = wa(m3,j)+aux*sign(one,wa(i,n4))
23037 continue
23038 continue
      wa(m3,j) = two*wa(m3,j)
23035 continue
23036 continue
      dif = zero
      init = .false.
      if(.not.lci2)then
      do23041 k = 1,n 
      wa(m5,k) = zero
      do23043 i = 1,m
      wa(m5,k) = wa(m5,k)+a(i,k)
23043 continue
23044 continue
      wa(m5,k) = wa(m5,k)/float(m)
23041 continue
23042 continue
      endif
      lsol = 1
      kount = 0
23045 continue
      do23048 j = 1,n
      wa(m1,j) = wa(m2,j)+wa(m3,j)*t
23048 continue
23049 continue
      if(.not.init)then
      stage = .true.
      kr = 1
      kl = 1
      go to 30
      endif
23052 continue
      stage = .false.
23055 continue
      max = -big
      do23058 j = kr,n 
      d = wa(m1,j)
      if(d.lt.zero)then
      if(d.gt.(-two))then
      goto 23058
      endif
      d = -d-two
      endif
      if(d.gt.max)then
      max = d
      in = j
      endif
23058 continue
23059 continue
      if(max.le.toler)then
      goto 23054
      endif
      if(wa(m1,in).le.zero)then
      do23070 i = 1,m4
      wa(i,in) = -wa(i,in)
23070 continue
23071 continue
      wa(m1,in) = wa(m1,in)-two
      wa(m2,in) = wa(m2,in)-two
      endif
23072 continue
      k = 0
      do23075 i = kl,m 
      d = wa(i,in)
      if(d.gt.toler)then
      k = k+1
      wb(k) = wa(i,n1)/d
      s(k) = i
      test = .true.
      endif
23075 continue
23076 continue
23079 continue
      if(k.le.0)then
      test = .false.
      else
      min = big
      do23084 i = 1,k
      if(wb(i).lt.min)then
      j = i
      min = wb(i)
      out = s(i)
      endif
23084 continue
23085 continue
      wb(j) = wb(k)
      s(j) = s(k)
      k = k-1
      endif
      if(.not.test.and.stage)then
      goto 23081
      endif
      if(.not.test)then
      goto 23047
      endif
      pivot = wa(out,in)
      if(wa(m1,in)-pivot-pivot.le.toler)then
      go to 10
      endif
      do23094 j = kr,n3 
      d = wa(out,j)
      wa(m1,j) = wa(m1,j)-d-d
      wa(m2,j) = wa(m2,j)-d-d
      wa(out,j) = -d
23094 continue
23095 continue
      wa(out,n4) = -wa(out,n4)
23080 goto 23079
23081 continue
      do23096 i = 1,m4 
      d = wa(i,kr)
      wa(i,kr) = wa(i,in)
      wa(i,in) = d
23096 continue
23097 continue
      kr = kr+1
      go to 20
10    do23098 j = kr,n3
      if(j.ne.in)then
      wa(out,j) = wa(out,j)/pivot
      endif
23098 continue
23099 continue
      do23102 i = 1,m3
      if(i.ne.out)then
      d = wa(i,in)
      do23106 j = kr,n3
      if(j.ne.in)then
      wa(i,j) = wa(i,j)-d*wa(out,j)
      endif
23106 continue
23107 continue
      endif
23102 continue
23103 continue
      do23110 i = 1,m3
      if(i.ne.out)then
      wa(i,in) = -wa(i,in)/pivot
      endif
23110 continue
23111 continue
      wa(out,in) = one/pivot
      d = wa(out,n4)
      wa(out,n4) = wa(m4,in)
      wa(m4,in) = d
      kount = kount+1
      if(.not.stage)then
      goto 23074
      endif
      kl = kl+1
      do23116 j = kr,n4 
      d = wa(out,j)
      wa(out,j) = wa(kount,j)
      wa(kount,j) = d
23116 continue
23117 continue
20    if(kount+kr.eq.n1)then
      goto 23057
      endif
30    max = -one
      do23120 j = kr,n
      if(abs(wa(m4,j)).le.n)then
      d = abs(wa(m1,j))
      if(d.gt.max)then
      max = d
      in = j
      endif
      endif
23120 continue
23121 continue
      if(wa(m1,in).lt.zero)then
      do23128 i = 1,m4
      wa(i,in) = -wa(i,in)
23128 continue
23129 continue
      endif
23073 goto 23072
23074 continue
23056 goto 23055
23057 continue
23053 goto 23052
23054 continue
      if(kr.eq.1)then
      do23132 j = 1,n 
      d = abs(wa(m1,j))
      if(d.le.toler.or.two-d.le.toler)then
      ift = 1
      wa(m2,nn+1) = zero
      go to 80
      endif
23132 continue
23133 continue
      endif
80    kount = 0
      sum = zero
      if(.not.lci2)then
      do23138 i = 1,kl-1 
      k = wa(i,n4)*sign(one,wa(i,n4))
      x(k) = wa(i,n1)*sign(one,wa(i,n4))
23138 continue
23139 continue
      endif
      do23140 i = 1,n 
      kd = abs(wa(m4,i))-n
      dsol(kd,lsol) = one+wa(m1,i)/two
      if(wa(m4,i).lt.zero)then
      dsol(kd,lsol) = one-dsol(kd,lsol)
      endif
      if(.not.lci2)then
      sum = sum + x(i)*wa(m5,i)
      sol(i+3,lsol) = x(i)
      h(i,lsol) = kd
      endif
23140 continue
23141 continue
      do23146 i = kl,m 
      kd = abs(wa(i,n4))-n
      if(wa(i,n4).lt.zero)then
      dsol(kd,lsol) = zero
      endif
      if(wa(i,n4).gt.zero)then
      dsol(kd,lsol) = one
      endif
23146 continue
23147 continue
      if(.not.lci2)then
      sol(1,lsol) = smax
      sol(2,lsol) = sum
      sum = zero
      do23154 j=kl,m
      d = wa(j,n1)*sign(one,wa(j,n4))
      sum = sum + d*(smax + half*(sign(one,d) - one))
23154 continue
23155 continue
      sol(3,lsol) = sum
      do23156 i=1,m
      dsol(i,lsol+1) = dsol(i,lsol)
23156 continue
23157 continue
      endif
      if(lci2)then
      a1 = zero
      do23160 i = 1,m 
      a1 = a1+a(i,idxcf)*(dsol(i,lsol)+t-one)
23160 continue
23161 continue
      tn = a1/sqrt(qn(idxcf)*t*(one-t))
      if(abs(tn).lt.cutoff)then
      if(lup)then
      smax = big
      else
      smax = -big
      endif
      do23166 i =1,kl-1 
      k = wa(i,n4)*sign(one,wa(i,n4))
      sol(k,1) = wa(i,n2)*sign(one,wa(i,n4))
      sol(k,2) = wa(i,n3)*sign(one,wa(i,n4))/tnew
23166 continue
23167 continue
      do23168 i = kl,m 
      a1 = zero
      b1 = zero
      k = wa(i,n4)*sign(one,wa(i,n4))-n
      l = 1
      do23170 j = 1,n
      if(j.eq.idxcf)then
      l = l+1
      endif
      a1 = a1 + a(k,l)*sol(j,1)
      b1 = b1 + a(k,l)*sol(j,2)
      l = l+1
23170 continue
23171 continue
      tnt = (b(k)-a1)/(a(k,idxcf)-b1)
      if(lup)then
      if(tnt.gt.tnew)then
      if(tnt.lt.smax)then
      smax = tnt
      out = i
      endif
      endif
      else
      if(tnt.lt.tnew)then
      if(tnt.gt.smax)then
      smax = tnt
      out = i
      endif
      endif
      endif
23168 continue
23169 continue
      if(lup)then
      told = tnew
      tnew = smax + toler
      ci(3,idxcf) = told - toler
      tnmat(3,idxcf) = tn
      if(.not.(tnew .lt. big-toler))then
      ci(3,idxcf) = big
      ci(4,idxcf) = big
      tnmat(3,idxcf) = tn
      tnmat(4,idxcf) = tn
      lup = .false.
      go to 70
      endif
      else
      told = tnew
      tnew = smax - toler
      ci(2,idxcf) = told + toler
      tnmat(2,idxcf) = tn
      if(.not.(tnew .gt. -big+toler))then
      ci(2,idxcf) = -big
      ci(1,idxcf) = -big
      tnmat(2,idxcf) = tn
      tnmat(1,idxcf) = tn
      lup = .true.
      go to 60
      endif
      endif
      do23190 i = 1,m
      wa(i,n3) = wa(i,n3)/told*tnew
      wa(i,n1) = wa(i,n2) - wa(i,n3)
23190 continue
23191 continue
      do23192 j = kr,n3
      d = wa(out,j)
      wa(m1,j) = wa(m1,j) -d -d
      wa(m2,j) = wa(m2,j) -d -d
      wa(out,j) = -d
23192 continue
23193 continue
      wa(out,n4) = -wa(out,n4)
      init = .true.
      else
      if(lup)then
      ci(4,idxcf) = tnew - toler
      tnmat(4,idxcf) = tn
      lup = .false.
      go to 70
      else
      ci(1,idxcf) = tnew + toler
      tnmat(1,idxcf) = tn
      lup = .true.
      go to 60
      endif
      endif
      endif
      if((iend).and.(.not.lci2))then
      go to 40
      endif
      if(.not.lci2)then
      init = .true.
      lsol = lsol+1
      do23200 i = 1,m
      s(i) = zero
23200 continue
23201 continue
      do23202 j = 1,n
      x(j) = zero
23202 continue
23203 continue
      smax = two
      do23204 j = 1,n 
      b1 = wa(m3,j)
      a1 = (-two-wa(m2,j))/b1
      b1 = -wa(m2,j)/b1
      if(a1.ge.t)then
      if(a1.lt.smax)then
      smax = a1
      dif = (b1-a1)/two
      endif
      endif
      if(b1.gt.t)then
      if(b1.lt.smax)then
      smax = b1
      dif = (b1-a1)/two
      endif
      endif
23204 continue
23205 continue
      tnt = smax+toler*(one+abs(dif))
      if(tnt.ge.t1+toler)then
      iend = .true.
      endif
      t = tnt
      if(iend)then
      t = t1
      endif
      endif
23046 goto 23045
23047 continue
      wa(m2,nn+1) = two
      ift = 2
      go to 50
40    if(lsol.gt.2)then
      sol(1,1) = zero
      sol(3,1) = zero
      sol(1,lsol) = one
      sol(3,lsol) = zero
      do23220 i = 1,m 
      dsol(i,1) = one
      dsol(i,lsol) = zero
      dsol(i,lsol+1) = zero
23220 continue
23221 continue
      endif
      l = kl-1
      do23222 i = 1,l
      if(wa(i,n1).lt.zero)then
      do23226 j = kr,n4
      wa(i,j) = -wa(i,j)
23226 continue
23227 continue
      endif
23222 continue
23223 continue
50    sum = zero
      if(.not.lci2)then
      do23230 i = kl,m 
      k = wa(i,n4)*sign(one,wa(i,n4))
      d = wa(i,n1)*sign(one,wa(i,n4))
      sum = sum+d*sign(one,d)*(half+sign(one,d)*(t-half))
      k = k-n
      e(k) = d
23230 continue
23231 continue
      wa(m2,n2) = kount
      wa(m1,n2) = n1-kr
      wa(m1,n1) = sum
      endif
      if(wa(m2,nn+1).eq.two)then
      goto 23018
      endif
      if(.not.lci1)then
      goto 23018
      endif
      if(.not.(.not.lci2))goto 23234
      lci2 = .true.
      n = nn-1
      n1 = n+1
      n2 = n+2
      n3 = n+3
      n4 = n+4
60    idxcf = idxcf+1
      if(.not.(idxcf.gt.nn))goto 23236
      goto 23018
23236 continue 
70    if(.not.(lup))goto 23238
      tnew = x(idxcf)+toler
      told = tnew
      ci(3,idxcf) = x(idxcf)
      tnmat(3,idxcf) = zero
      goto 23239
23238 continue
      tnew = x(idxcf)-toler
      told = tnew
      ci(2,idxcf) = x(idxcf)
      tnmat(2,idxcf) = zero
23239 continue
23234 continue
23017 goto 23016
23018 continue
      do23242 i=1,m
      dsol(i,lsol) = dsol(i,lsol+1)
23242 continue
23243 continue
      endif
      return
      end
