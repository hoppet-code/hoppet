      subroutine printer2(na,nb) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
