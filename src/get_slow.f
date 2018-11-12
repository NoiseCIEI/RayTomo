c============================================================
c Estimate slovness for phase veloocity using model file
c============================================================
      subroutine GET_SLOW(type,type1,period,slowness,ierr)
      implicit none
      include "tomo.h"
      character*4 type,type1
      real*8    period,slowness,per(100),r_ph(100),l_ph(100)
      real*8    r_gr(100),l_gr(100)
      integer*4 n_t,i,ierr,lnblnk
      logical   tf
      ierr=0
c
c Check that file fname5 exists
c
      inquire(file=fname5,exist=tf)
      if(.not.tf) then
        i=lnblnk(fname5)
        write(*,*) '(R0010) File ',fname5(1:i), ' does not exists'
        STOP
      endif
c
c Read input data
c
      open(26,file=fname5,status='old')
      read(26,*) n_t
      do i=1,n_t
        read(26,*) per(i),r_gr(i),r_ph(i),l_gr(i),l_ph(i)
      enddo
      close(26)
      if(type.eq.'R'.or.type.eq.'r') then
         do i=1,n_t-1
           if(period.ge.per(i).and.period.le.per(i+1)) then
             if(type1.eq.'G'.or.type1.eq.'g') then
             slowness=1.0/((r_gr(i+1)-r_gr(i))/(per(i+1)-per(i))*
     *       (period-per(i))+r_gr(i))
             else
             slowness=1.0/((r_ph(i+1)-r_ph(i))/(per(i+1)-per(i))*
     *       (period-per(i))+r_ph(i))
            endif
           goto 2
           endif
         enddo
         ierr=1
    2    continue
      else
         do i=1,n_t-1
           if(period.ge.per(i).and.period.le.per(i+1)) then
             if(type1.eq.'G'.or.type1.eq.'g') then
             slowness=1.0/((l_gr(i+1)-l_gr(i))/(per(i+1)-per(i))*
     *       (period-per(i))+l_gr(i))
             else
             slowness=1.0/((l_ph(i+1)-l_ph(i))/(per(i+1)-per(i))*
     *       (period-per(i))+l_ph(i))
             endif
           goto 3
           endif
         enddo
         ierr=1
    3    continue
      endif
      return
      end
