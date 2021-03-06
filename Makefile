## Copyright (C) 2010 Riad S Wahby <rsw@jfet.org>
## 
## This file is part of resgrad
##
##  resgrad is free software.  It comes without any warranty, to
##  to the extent permitted by applicable law.  You can redistribute it
##  and/or modify it under the terms of the Do What The Fuck You Want To
##  Public License, Version 2, as published by Sam Hocevar.  See
##  the COPYING file or http://sam.zoy.org/wtfpl/COPYING for more details
##

all: resgrad

resgrad: resgrad.hs
	ghc -O6 --make -threaded resgrad.hs

clean:
	rm -rf *.o *.hi resgrad *.svg

sclean:
	rm -rf *.svg
