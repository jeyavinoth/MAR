      integer function zext(logvar)
      logical logvar
C +
      if                   (logvar) then
                       zext=-1
      else
                       zext= 0
      end if
C +
      return
      end
