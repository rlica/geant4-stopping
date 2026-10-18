c      Subroutine condenses a character string and gives its length
c       word is input, cword is condenced word
c--------------------------------------------------------
      subroutine condens_s(word,cword,length)
      character*1 letter
      character*(*) word
      character*(*) cword
      l = len(word)
      cword = word(1:1)
      length = 1
      do i=2,l
       letter = word(i:i)
       if(letter.ne.' '.and.letter.ne.char(26))then
        length = length + 1       
        cword(length:length)=letter
       endif
      enddo
      return
      end
